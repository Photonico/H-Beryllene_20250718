#!/usr/bin/env python3
"""Adaptive fixed-density SOC zoom into explicitly chosen direct-gap basins.

Example: python -B refine_gap_targeted.py STRUCTURE_DIR --vasp /path/to/vasp_ncl \
             --center 0.32338,0.36190,0.000595 --center 0.3283,0.3135,0.00476 --levels 4
Each --center is k1,k2,halfwidth (fractional). Level l writes gap_refine_3<a,b,...>
with a (2*steps+1)^2 patch per basin, re-centred on that basin's previous minimum;
the halfwidth shrinks by --shrink only when the previous minimum was interior.
Every level is validated like the earlier stages; gap_refine_3_summary.json
records each basin's trajectory, a local slope bound and the aggregate extrema.
Existing stages are only read. Local sampling, not a proof of a continuous-BZ gap.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import string
import subprocess
import sys

import h5py
import numpy as np

sys.dont_write_bytecode = True
from check_topology_outputs import fingerprint, manifest_entry
from refine_gap import aggregate, dump, read_stage, set_incar


def patch(center, halfwidth, steps):
    step = halfwidth / steps
    return [((center[0] + ix * step) % 1, (center[1] + iy * step) % 1, ix, iy)
            for ix in range(-steps, steps + 1) for iy in range(-steps, steps + 1)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('structure_dir', type=Path)
    parser.add_argument('--vasp', type=Path, required=True)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--ranks', type=int, default=42)
    parser.add_argument('--center', action='append', required=True)
    parser.add_argument('--levels', type=int, default=4)
    parser.add_argument('--steps', type=int, default=4)
    parser.add_argument('--shrink', type=float, default=4.0)
    args = parser.parse_args()
    root = args.structure_dir.resolve(strict=True)
    vasp = args.vasp.resolve(strict=True)
    if not (0 < args.levels <= 8 and args.steps > 0 and args.shrink > 1):
        raise ValueError('Need 0 < levels <= 8, steps > 0 and shrink > 1')
    manifest = root / 'manifest.json'
    entry, nelect = manifest_entry(manifest, None)
    candidates = [nelect] if nelect % 2 == 0 else [nelect - 1, nelect + 1]
    scf, trim = root / 'scf', root / 'trim'
    summary = json.loads((root / 'gap_refinement_summary.json').read_text())
    if fingerprint(scf / 'CHGCAR') != summary['fixed_density']:
        raise ValueError('SCF charge density differs from the one used by earlier refinements')
    names = ['scf', 'trim', 'bands'] + [level['stage'] for level in summary['levels']]
    stages = [read_stage(root / name, nelect, candidates) for name in names]
    initial = aggregate(stages, candidates)
    basins = []
    for text in args.center:
        k1, k2, halfwidth = (float(v) for v in text.split(','))
        if not 0 < halfwidth <= 0.05:
            raise ValueError('Basin halfwidth must be in (0, 0.05]')
        basins.append({'start': [k1 % 1, k2 % 1], 'center': [k1 % 1, k2 % 1],
                       'halfwidth': halfwidth, 'levels': []})
    reciprocal = 2 * np.pi * np.linalg.inv(np.array(json.loads(
        (trim / 'topology_validation.json').read_text())['structure']['lattice'])).T[:2, :2]
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'structure_id': entry['id'],
              'candidate_spinor_band_counts': candidates, 'fixed_density': summary['fixed_density'],
              'previous_stages': names, 'steps': args.steps, 'shrink': args.shrink,
              'basins': basins, 'topology_certified': False}
    for level in range(args.levels):
        stage = root / ('gap_refine_3' + string.ascii_lowercase[level])
        points = [(b, p) for b, basin in enumerate(basins)
                  for p in patch(basin['center'], basin['halfwidth'], args.steps)]
        if len(points) > 1000:
            raise ValueError('Refinement exceeded the 1000-k-point budget')
        stage.mkdir()
        for name in ('POSCAR', 'POTCAR'):
            shutil.copyfile(trim / name, stage / name)
        shutil.copyfile(scf / 'CHGCAR', stage / 'CHGCAR')
        incar = set_incar((trim / 'INCAR').read_text(),
                          {'ICHARG': 11, 'ISTART': 0, 'ISYM': -1, 'KPAR': 1,
                           'LWAVE': '.FALSE.', 'LCHARG': '.FALSE.', 'LORBIT': 11})
        (stage / 'INCAR').write_text(incar)
        (stage / 'KPOINTS').write_text('Adaptive SOC direct-gap zoom\n' + str(len(points)) +
                                       '\nReciprocal\n' + ''.join('%.12f %.12f 0.0 1\n' % p[:2]
                                                                  for _, p in points))
        dump(stage / 'refinement_input.json',
             {'stage': stage.name, 'basins': [{'center': b['center'], 'halfwidth': b['halfwidth'],
                                               'spacing': b['halfwidth'] / args.steps} for b in basins]})
        print('Running', entry['id'], stage.name, len(points), 'k points', flush=True)
        with (stage / 'vasp.log').open('x') as log:
            subprocess.run(['mpirun', '-np', str(args.ranks), str(vasp)], cwd=stage,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([args.python, '-B', str(Path(__file__).with_name('check_topology_outputs.py')),
                        str(stage), '--manifest', str(manifest), '--structure-id', entry['id']], check=True)
        stages.append(read_stage(stage, nelect, candidates))
        with h5py.File(stage / 'vaspout.h5', 'r') as stream:
            energies = stream['/results/electron_eigenvalues/eigenvalues'][()][0]
            kpoints = stream['/results/electron_eigenvalues/kpoint_coords'][()]
        written = np.array([p[:2] for _, p in points])
        if len(kpoints) != len(points) or np.abs((kpoints[:, :2] - written + 0.5) % 1 - 0.5).max() > 1e-8:
            raise ValueError('vaspout.h5 k points do not match the written KPOINTS order')
        gaps = np.min([energies[:, n] - energies[:, n - 1] for n in candidates], axis=0)
        for b, basin in enumerate(basins):
            index = [i for i, (owner, _) in enumerate(points) if owner == b]
            local = gaps[index]
            best = index[int(np.argmin(local))]
            ix, iy = points[best][1][2:]
            spacing = basin['halfwidth'] / args.steps
            grid = local.reshape(2 * args.steps + 1, 2 * args.steps + 1)
            dx = np.linalg.norm(spacing * reciprocal[0])
            dy = np.linalg.norm(spacing * reciprocal[1])
            slope = max(np.abs(np.diff(grid, axis=0)).max() / dx, np.abs(np.diff(grid, axis=1)).max() / dy)
            radius = 0.5 * max(np.linalg.norm(spacing * (reciprocal[0] + reciprocal[1])),
                               np.linalg.norm(spacing * (reciprocal[0] - reciprocal[1])))
            interior = abs(ix) < args.steps and abs(iy) < args.steps
            basin['levels'].append({'stage': stage.name, 'center': list(basin['center']),
                                    'halfwidth': basin['halfwidth'], 'spacing': spacing,
                                    'min_direct_gap_ev': float(local.min()),
                                    'min_k': [float(v) for v in points[best][1][:2]],
                                    'min_interior': bool(interior),
                                    'max_slope_ev_angstrom': float(slope),
                                    'sampling_radius_inverse_angstrom': float(radius),
                                    'slope_bound_ev': float(local.min() - slope * radius)})
            basin['center'] = [float(v) for v in points[best][1][:2]]
            if interior:
                basin['halfwidth'] /= args.shrink
    if fingerprint(scf / 'CHGCAR') != summary['fixed_density']:
        raise ValueError('SCF charge density changed during refinement')
    report['stages'] = stages[len(names):]
    report['sampled_extrema_before'] = initial
    report['sampled_extrema_after'] = aggregate(stages, candidates)
    report['note'] = ('slope_bound_ev = min gap - (max neighbour slope) x (sampling radius) is a '
                      'heuristic local bound; a positive value supports, but does not prove, a gap.')
    dump(root / 'gap_refine_3_summary.json', report)
    print(root / 'gap_refine_3_summary.json', flush=True)


if __name__ == '__main__':
    main()
