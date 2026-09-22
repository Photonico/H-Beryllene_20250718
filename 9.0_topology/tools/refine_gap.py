#!/usr/bin/env python3
"""Two bounded fixed-density SOC gap refinements around sampled extrema.

Example: python -B refine_gap.py STRUCTURE_DIR --vasp /path/to/vasp_ncl
Creates gap_refine_1 and gap_refine_2, then gap_refinement_summary.json.
This is local k-point refinement, not a proof of a continuous-BZ gap.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
from check_topology_outputs import eigenval_scan, fingerprint, manifest_entry


def dump(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def read_stage(stage, nelect, candidates):
    validation_path = stage / 'topology_validation.json'
    validation = json.loads(validation_path.read_text())
    if not validation.get('validation_passed'):
        raise ValueError('Stage did not pass validation: ' + str(stage))
    for name in ('POSCAR', 'OUTCAR', 'EIGENVAL'):
        if validation['provenance'][name] != fingerprint(stage / name):
            raise ValueError('Validated file changed: ' + str(stage / name))
    eigen = eigenval_scan(stage / 'EIGENVAL', nelect, candidates, 1e-6)
    return {'stage': stage.name,
            'eigenval': fingerprint(stage / 'EIGENVAL'),
            'manifolds': [item for item in eigen['manifolds'] if item['N'] in candidates]}


def aggregate(stages, candidates):
    result = []
    for n in candidates:
        choices = [(stage['stage'], item) for stage in stages
                   for item in stage['manifolds'] if item['N'] == n]
        d_stage, direct = min(choices, key=lambda pair: pair[1]['direct_gap_ev'])
        v_stage, valence = max(choices, key=lambda pair: pair[1]['valence_max_ev'])
        c_stage, conduction = min(choices, key=lambda pair: pair[1]['conduction_min_ev'])
        result.append({'N': n, 'bands_one_based': [1, n],
                       'direct_gap_ev': direct['direct_gap_ev'],
                       'direct_gap_k': direct['direct_gap_k'], 'direct_gap_stage': d_stage,
                       'valence_max_ev': valence['valence_max_ev'],
                       'valence_max_k': valence['valence_max_k'], 'valence_max_stage': v_stage,
                       'conduction_min_ev': conduction['conduction_min_ev'],
                       'conduction_min_k': conduction['conduction_min_k'],
                       'conduction_min_stage': c_stage,
                       'indirect_gap_ev': conduction['conduction_min_ev'] - valence['valence_max_ev']})
    return result


def centers_for(stages, candidates):
    # Preserve one direct-gap basin from each original/new sampled set and
    # refine the combined valence/conduction extrema. At most 12 centers.
    points = [item['direct_gap_k'] for stage in stages for item in stage['manifolds']]
    for item in aggregate(stages, candidates):
        points.extend([item['valence_max_k'], item['conduction_min_k']])
    return sorted({tuple(round(float(v) % 1.0, 12) for v in point[:2]) + (0.0,)
                   for point in points})


def base_spacing(scf, override):
    if override is not None:
        if override <= 0:
            raise ValueError('--base-grid must be positive')
        return [1.0 / override] * 2
    lines = (scf / 'KPOINTS').read_text().splitlines()
    if int(lines[1].split()[0]) != 0 or lines[2].strip().lower()[0] not in ('g', 'm'):
        raise ValueError('Supply --base-grid for a nonautomatic SCF KPOINTS')
    mesh = [int(v) for v in lines[3].split()[:3]]
    if len(mesh) != 3 or min(mesh[:2]) <= 0 or mesh[2] != 1:
        raise ValueError('Expected a two-dimensional SCF mesh')
    return [1.0 / value for value in mesh[:2]]


def set_incar(text, updates):
    # Prepared campaign INCARs have one setting per line.
    for key, value in updates.items():
        text = re.sub(r'^\s*' + re.escape(key) + r'\s*=.*\n?', '', text,
                      flags=re.MULTILINE | re.IGNORECASE)
        text += '\n' + key + ' = ' + str(value) + '\n'
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('structure_dir', type=Path)
    parser.add_argument('--vasp', type=Path, required=True)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--ranks', type=int, default=42)
    parser.add_argument('--base-grid', type=int)
    args = parser.parse_args()
    root = args.structure_dir.resolve(strict=True)
    vasp = args.vasp.resolve(strict=True)
    if args.ranks <= 0:
        raise ValueError('--ranks must be positive')
    manifest = root / 'manifest.json'
    entry, nelect = manifest_entry(manifest, None)
    sid = entry['id']
    candidates = [nelect] if nelect % 2 == 0 else [nelect - 1, nelect + 1]
    scf, trim = root / 'scf', root / 'trim'
    steps = base_spacing(scf, args.base_grid)
    stages = [read_stage(root / name, nelect, candidates) for name in ('scf', 'trim', 'bands')]
    initial = aggregate(stages, candidates)
    fixed_density = fingerprint(scf / 'CHGCAR')
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'structure_id': sid,
              'expected_nelect': nelect, 'candidate_spinor_band_counts': candidates,
              'neutral_odd_filling': bool(nelect % 2),
              'fixed_density_source': str(scf / 'CHGCAR'), 'fixed_density': fixed_density,
              'initial_sampled_extrema': initial, 'levels': [], 'topology_certified': False,
              'scope': 'Sampled 2D BZ plus two local refinements; no continuous-BZ proof, '
                       'ENCUT/charge-density convergence claim, or final Z2 certification.'}
    for level in (1, 2):
        stage = root / ('gap_refine_' + str(level))
        centers = centers_for(stages, candidates)
        halfwidth = [value / (4 ** (level - 1)) for value in steps]
        points = set()
        for center in centers:
            for ix in range(-4, 5):
                for iy in range(-4, 5):
                    points.add((round((center[0] + ix * halfwidth[0] / 4) % 1, 12),
                                round((center[1] + iy * halfwidth[1] / 4) % 1, 12), 0.0))
        points.update((x, y, 0.0) for x in (0.0, 0.5) for y in (0.0, 0.5))
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
        (stage / 'KPOINTS').write_text('Local SOC direct/indirect-gap refinement\n' +
                                     str(len(points)) + '\nReciprocal\n' +
                                     ''.join('%.12f %.12f %.12f 1\n' % point for point in sorted(points)))
        level_record = {'level': level, 'stage': stage.name, 'centers_fractional': centers,
                        'halfwidth_fractional': halfwidth, 'points': len(points),
                        'spacing_fractional': [v / 4 for v in halfwidth]}
        dump(stage / 'refinement_input.json', level_record)
        print('Running', sid, stage.name, len(points), 'k points', flush=True)
        with (stage / 'vasp.log').open('x') as log:
            subprocess.run(['mpirun', '-np', str(args.ranks), str(vasp)], cwd=stage,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([args.python, '-B', str(Path(__file__).with_name('check_topology_outputs.py')),
                        str(stage), '--manifest', str(manifest), '--structure-id', sid], check=True)
        refined = read_stage(stage, nelect, candidates)
        previous = aggregate(stages, candidates)
        stages.append(refined)
        combined = aggregate(stages, candidates)
        level_record['sampled_extrema_after_level'] = combined
        level_record['direct_gap_change_ev'] = {
            str(after['N']): after['direct_gap_ev'] - before['direct_gap_ev']
            for before, after in zip(previous, combined)}
        report['levels'].append(level_record)
        dump(stage / 'refinement_summary.json', level_record)
    if fingerprint(scf / 'CHGCAR') != fixed_density:
        raise ValueError('SCF charge density changed during refinement')
    report['stages'] = stages
    report['final_sampled_extrema'] = aggregate(stages, candidates)
    report['further_refinement_needed_if'] = [
        'The direct gap is near EIGENVAL printing precision (approximately 1 microelectronvolt).',
        'The minimum is comparable to its change between refinement levels.',
        'A new lower gap or other candidate minimum is found elsewhere in the 2D BZ.']
    dump(root / 'gap_refinement_summary.json', report)
    print(root / 'gap_refinement_summary.json', flush=True)


if __name__ == '__main__':
    main()
