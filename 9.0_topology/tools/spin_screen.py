#!/usr/bin/env python3
"""Collinear spin-polarised stability screen of the nonmagnetic SOC reference.

Example: python -B spin_screen.py STRUCTURE_DIR --vasp /path/to/vasp_std
Creates STRUCTURE_DIR/spin_screen/{fm,afm}/ and spin_screen/spin_screen_summary.json.
Every SOC run started from MAGMOM=0, so it cannot reveal a lower-energy
magnetic solution. Here ISPIN=2 SCFs start from finite ferromagnetic and,
with two or more Be atoms, layer-alternating moments on the same geometry,
mesh, cutoff, smearing and dipole setup (no SOC: Be/H SOC is meV-scale).
Collapse of every seed supports, but does not prove, a nonmagnetic ground
state; supercell magnetic orders are not sampled.
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
from check_topology_outputs import fingerprint, manifest_entry, read_outcar
from refine_gap import dump, set_incar

MOMENT_TOLERANCE_MUB = 1e-3


def seeds_for(poscar):
    lines = poscar.read_text().splitlines()
    species, counts = lines[5].split(), [int(v) for v in lines[6].split()]
    labels = [name for name, count in zip(species, counts) for _ in range(count)]
    heights = [float(line.split()[2]) for line in lines[8:8 + len(labels)]]
    beryllium = sorted((z, i) for i, (z, name) in enumerate(zip(heights, labels)) if name == 'Be')
    fm = [1.0 if name == 'Be' else 0.5 for name in labels]
    seeds = {'fm': fm}
    if len(beryllium) >= 2:
        # Layer-alternating Be moments; for two Be this is also inversion-odd.
        afm = [0.0] * len(labels)
        for rank, (_, index) in enumerate(beryllium):
            afm[index] = 1.0 if rank % 2 == 0 else -1.0
        seeds['afm'] = afm
    if len(beryllium) >= 3 and len(beryllium) % 2:
        # Odd layer count: alternating seed is inversion-even, so add the
        # inversion-odd (PT-symmetric) pattern +, 0, ..., 0, -.
        odd = [0.0] * len(labels)
        odd[beryllium[0][1]], odd[beryllium[-1][1]] = 1.0, -1.0
        seeds['p_odd'] = odd
    return labels, seeds


def site_moments(outcar):
    """Last LORBIT=11 'magnetization (x)' table: per-ion and total 'tot' column."""
    lines = outcar.read_text(errors='replace').splitlines()
    starts = [i for i, line in enumerate(lines) if line.strip() == 'magnetization (x)']
    if not starts:
        return None, None
    rows = []
    for line in lines[starts[-1] + 1:]:
        fields = line.split()
        if fields and fields[0] == 'tot':
            return rows, float(fields[-1])
        if fields and fields[0].isdigit():
            rows.append(float(fields[-1]))
    return None, None


def oszicar_last(oszicar):
    line = [row for row in oszicar.read_text().splitlines() if ' F= ' in row][-1]
    values = dict(re.findall(r'(\w+)=\s*(\S+)', line))
    return {'F_ev': float(values['F']), 'E0_ev': float(values['E0']), 'mag_muB': float(values['mag'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('structure_dir', type=Path)
    parser.add_argument('--vasp', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=42)
    args = parser.parse_args()
    root = args.structure_dir.resolve(strict=True)
    vasp = args.vasp.resolve(strict=True)
    entry, nelect = manifest_entry(root / 'manifest.json', None)
    scf = root / 'scf'
    labels, seeds = seeds_for(scf / 'POSCAR')
    screen = root / 'spin_screen'
    screen.mkdir()
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'structure_id': entry['id'],
              'expected_nelect': nelect, 'vasp': str(vasp), 'atoms': labels,
              'inputs': {name: fingerprint(scf / name) for name in ('POSCAR', 'POTCAR', 'KPOINTS', 'INCAR')},
              'moment_tolerance_muB': MOMENT_TOLERANCE_MUB, 'seeds': {}}
    for name, magmom in seeds.items():
        run = screen / name
        run.mkdir()
        for source in ('POSCAR', 'POTCAR', 'KPOINTS'):
            shutil.copyfile(scf / source, run / source)
        text = (scf / 'INCAR').read_text()
        for tag in ('LSORBIT', 'LNONCOLLINEAR', 'SAXIS', 'ISYM'):
            text = re.sub(r'^\s*' + tag + r'\s*=.*\n?', '', text, flags=re.MULTILINE | re.IGNORECASE)
        text = text.replace('SOC topology', 'collinear spin screen (' + name + ')')
        text = set_incar(text, {'ISPIN': 2, 'MAGMOM': ' '.join('%g' % m for m in magmom),
                                'NBANDS': 14, 'ISTART': 0, 'ICHARG': 2, 'KPAR': 1, 'NELM': 300,
                                'EDIFF': '1E-7', 'LWAVE': '.FALSE.', 'LCHARG': '.FALSE.', 'LORBIT': 11})
        (run / 'INCAR').write_text(text)
        print('Running', entry['id'], name, 'MAGMOM', magmom, flush=True)
        with (run / 'vasp.log').open('x') as log:
            subprocess.run(['mpirun', '-np', str(args.ranks), str(vasp)], cwd=run,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        metadata, errors = read_outcar(run / 'OUTCAR')
        errors = [e for e in errors if not e.startswith(('LSORBIT', 'LNONCOLLINEAR'))]
        sites, site_total = site_moments(run / 'OUTCAR')
        final = oszicar_last(run / 'OSZICAR')
        collapsed = (not errors and sites is not None and abs(final['mag_muB']) < MOMENT_TOLERANCE_MUB
                     and max(abs(v) for v in sites) < MOMENT_TOLERANCE_MUB)
        report['seeds'][name] = {'initial_magmom_muB': magmom, 'errors': errors,
                                 'electronic_iterations': metadata.get('final_electronic_iteration'),
                                 **final, 'site_moments_muB': sites, 'site_moment_total_muB': site_total,
                                 'moment_collapsed': collapsed}
    report['all_seeds_collapsed'] = all(s['moment_collapsed'] for s in report['seeds'].values())
    report['interpretation'] = ('All finite-moment seeds relaxed to a nonmagnetic state: no ferromagnetic or '
                                'in-cell layer-alternating instability at this level of theory.'
                                if report['all_seeds_collapsed'] else
                                'At least one seed kept a moment or failed: compare its energy with an ISPIN=1 '
                                'reference before any TR-based Z2 statement.')
    dump(screen / 'spin_screen_summary.json', report)
    print(screen / 'spin_screen_summary.json', flush=True)


if __name__ == '__main__':
    main()
