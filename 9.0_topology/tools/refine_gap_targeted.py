#!/usr/bin/env python3
"""One further fixed-density SOC gap refinement around explicitly chosen basins.

Example: python -B refine_gap_targeted.py STRUCTURE_DIR --vasp /path/to/vasp_ncl \
             --center 0.3232,0.3619 --halfwidth 0.0006 --steps 8
Creates STRUCTURE_DIR/gap_refine_3 (or --stage) with a (2*steps+1)^2 square patch
per center, validates it and writes refinement_summary.json there, aggregating
all earlier stages. Existing stages and gap_refinement_summary.json are only read.
Local sampling, not a proof of a continuous-BZ gap.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
from check_topology_outputs import fingerprint, manifest_entry
from refine_gap import aggregate, dump, read_stage, set_incar


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('structure_dir', type=Path)
    parser.add_argument('--vasp', type=Path, required=True)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--ranks', type=int, default=42)
    parser.add_argument('--stage', default='gap_refine_3')
    parser.add_argument('--center', action='append', required=True,
                        help='fractional k1,k2 of a basin; repeat for several basins')
    parser.add_argument('--halfwidth', type=float, required=True)
    parser.add_argument('--steps', type=int, default=8)
    args = parser.parse_args()
    root = args.structure_dir.resolve(strict=True)
    vasp = args.vasp.resolve(strict=True)
    if not (0 < args.halfwidth <= 0.05 and args.steps > 0):
        raise ValueError('Need 0 < halfwidth <= 0.05 and steps > 0')
    manifest = root / 'manifest.json'
    entry, nelect = manifest_entry(manifest, None)
    candidates = [nelect] if nelect % 2 == 0 else [nelect - 1, nelect + 1]
    scf, trim = root / 'scf', root / 'trim'
    summary = json.loads((root / 'gap_refinement_summary.json').read_text())
    if fingerprint(scf / 'CHGCAR') != summary['fixed_density']:
        raise ValueError('SCF charge density differs from the one used by earlier refinements')
    names = ['scf', 'trim', 'bands'] + [level['stage'] for level in summary['levels']]
    stages = [read_stage(root / name, nelect, candidates) for name in names]
    previous = aggregate(stages, candidates)
    centers = [tuple(float(v) for v in text.split(',')) + (0.0,) for text in args.center]
    if any(len(center) != 3 for center in centers):
        raise ValueError('--center takes two fractional coordinates')
    step = args.halfwidth / args.steps
    points = set()
    for center in centers:
        for ix in range(-args.steps, args.steps + 1):
            for iy in range(-args.steps, args.steps + 1):
                points.add((round((center[0] + ix * step) % 1, 12),
                            round((center[1] + iy * step) % 1, 12), 0.0))
    if len(points) > 1000:
        raise ValueError('Refinement exceeded the 1000-k-point budget')
    stage = root / args.stage
    stage.mkdir()
    for name in ('POSCAR', 'POTCAR'):
        shutil.copyfile(trim / name, stage / name)
    shutil.copyfile(scf / 'CHGCAR', stage / 'CHGCAR')
    incar = set_incar((trim / 'INCAR').read_text(),
                      {'ICHARG': 11, 'ISTART': 0, 'ISYM': -1, 'KPAR': 1,
                       'LWAVE': '.FALSE.', 'LCHARG': '.FALSE.', 'LORBIT': 11})
    (stage / 'INCAR').write_text(incar)
    (stage / 'KPOINTS').write_text('Targeted SOC direct-gap refinement\n' + str(len(points)) +
                                   '\nReciprocal\n' + ''.join('%.12f %.12f %.12f 1\n' % point
                                                              for point in sorted(points)))
    record = {'stage': stage.name, 'centers_fractional': centers,
              'halfwidth_fractional': [args.halfwidth] * 2, 'points': len(points),
              'spacing_fractional': [step] * 2, 'previous_stages': names}
    dump(stage / 'refinement_input.json', record)
    print('Running', entry['id'], stage.name, len(points), 'k points', flush=True)
    with (stage / 'vasp.log').open('x') as log:
        subprocess.run(['mpirun', '-np', str(args.ranks), str(vasp)], cwd=stage,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run([args.python, '-B', str(Path(__file__).with_name('check_topology_outputs.py')),
                    str(stage), '--manifest', str(manifest), '--structure-id', entry['id']], check=True)
    stages.append(read_stage(stage, nelect, candidates))
    combined = aggregate(stages, candidates)
    if fingerprint(scf / 'CHGCAR') != summary['fixed_density']:
        raise ValueError('SCF charge density changed during refinement')
    record['stage_manifolds'] = stages[-1]['manifolds']
    record['sampled_extrema_after_level'] = combined
    record['direct_gap_change_ev'] = {str(after['N']): after['direct_gap_ev'] - before['direct_gap_ev']
                                      for before, after in zip(previous, combined)}
    record['topology_certified'] = False
    dump(stage / 'refinement_summary.json', record)
    print(stage / 'refinement_summary.json', flush=True)


if __name__ == '__main__':
    main()
