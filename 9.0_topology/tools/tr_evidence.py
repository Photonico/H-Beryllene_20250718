#!/usr/bin/env python3
"""Electronic time-reversal evidence of the converged SOC state (read-only).

Example: python -B tr_evidence.py CAMPAIGN_ROOT
Writes <structure>/tr_evidence.json for every manifest structure from
scf/CHGCAR (magnetisation density and PAW one-centre occupancies), the
float64 eigenvalues in trim/vaspout.h5 (Kramers pairs at the four TRIM) and,
for structures without inversion, E_n(k) vs E_n(-k) on the scf grid.
This characterises the converged state only: every SCF started from m=0 and
ISYM=2 removed inversion-odd m, so stability against magnetic order is a
separate question (spin_screen.py).
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import h5py
import numpy as np

sys.dont_write_bytecode = True
from check_topology_outputs import fingerprint
from refine_gap import dump


def chgcar_magnetisation(path):
    """Grids are rho*V (x fastest); each is followed by PAW occupancy blocks."""
    lines = path.read_text().split('\n')
    lattice = np.array([[float(v) for v in lines[i].split()] for i in (2, 3, 4)]) * float(lines[1].split()[0])
    volume = abs(np.linalg.det(lattice))
    header = 8 + sum(int(v) for v in lines[6].split())
    while not lines[header].strip():
        header += 1
    key = lines[header].split()
    size = int(np.prod([int(v) for v in key]))
    starts = [i for i, line in enumerate(lines) if line.split() == key]
    if len(starts) != 4:
        raise ValueError('Expected four noncollinear CHGCAR grids: ' + str(path))
    grids, real, imaginary = [], [], []
    for index, start in enumerate(starts):
        rows = -(-size // 5)
        grids.append(np.array(' '.join(lines[start + 1:start + 1 + rows]).split(), float))
        end = starts[index + 1] if index + 1 < len(starts) else len(lines)
        re_max = im_max = 0.0
        row = start + 1 + rows
        while row < end:
            if lines[row].startswith('augmentation occupancies'):
                count, imag = int(lines[row].split()[-1]), 'imaginary' in lines[row]
                values, row = [], row + 1
                while len(values) < count:
                    values += [float(v) for v in lines[row].split()]
                    row += 1
                if imag:
                    im_max = max(im_max, float(np.abs(values).max()))
                else:
                    re_max = max(re_max, float(np.abs(values).max()))
                continue
            row += 1
        real.append(re_max)
        imaginary.append(im_max)
    rho, m = grids[0], np.vstack(grids[1:])
    norm = np.sqrt((m ** 2).sum(0))
    return {'electrons_from_grid': float(rho.sum() / size),
            'max_abs_m_muB_per_A3': float(norm.max() / volume),
            'max_abs_m_over_max_rho': float(norm.max() / rho.max()),
            'integral_abs_m_muB': float(norm.sum() / size),
            'integral_m_muB': (m.sum(1) / size).tolist(),
            'paw_occupancy_max_abs_re_m': max(real[1:]),
            'paw_occupancy_max_abs_im_m': max(imaginary[1:]),
            'note': 'With real PAW projectors TR requires Re(one-centre m)=0; Im(one-centre m) '
                    'is the TR-even spin-orbit correlation and sets the SOC scale.'}


def eigenvalues(path):
    with h5py.File(path, 'r') as stream:
        return (stream['/results/electron_eigenvalues/eigenvalues'][()][0],
                stream['/results/electron_eigenvalues/kpoint_coords'][()])


def minus_k_deviation(energies, kpoints, bands):
    """max |E_n(k) - E_n(-k)| over +/-k pairs of the sampled set (in-plane, mod G)."""
    scale = 10 ** 7
    keys = {tuple(np.mod(np.rint(k[:2] * scale).astype(np.int64), scale)): i
            for i, k in enumerate(kpoints)}
    worst, pairs = 0.0, 0
    for i, k in enumerate(kpoints):
        j = keys.get(tuple(np.mod(np.rint(-k[:2] * scale).astype(np.int64), scale)))
        if j is not None and j > i:
            pairs += 1
            worst = max(worst, float(np.abs(energies[i, :bands] - energies[j, :bands]).max()))
    return pairs, worst


def structure_evidence(folder, entry):
    nelect = entry['expected_nelect']
    subspaces = [nelect] if nelect % 2 == 0 else [nelect - 1, nelect + 1]
    bands = max(subspaces) + 2
    energies, _ = eigenvalues(folder / 'trim' / 'vaspout.h5')
    kramers = float(np.abs(energies[:, 1:bands:2] - energies[:, 0:bands:2]).max())
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'structure_id': entry['id'],
              'inputs': {name: fingerprint(folder / name) for name in
                         ('scf/CHGCAR', 'trim/vaspout.h5', 'scf/vaspout.h5')},
              'magnetisation_density': chgcar_magnetisation(folder / 'scf' / 'CHGCAR'),
              'trim_kramers_max_splitting_ev': kramers,
              'trim_kramers_bands_one_based': [1, bands]}
    if not entry['inversion_expected']:
        energies, kpoints = eigenvalues(folder / 'scf' / 'vaspout.h5')
        pairs, worst = minus_k_deviation(energies, kpoints, bands)
        report['scf_minus_k'] = {'kpoints': int(len(kpoints)), 'pairs': pairs,
                                 'bands_one_based': [1, bands], 'max_abs_dE_ev': worst}
    report['scope'] = ('Converged SOC state only. SCFs started from m=0 and ISYM=2 symmetrisation '
                       'removed inversion-odd m; stability against magnetic order is tested separately.')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    args = parser.parse_args()
    campaign = args.campaign.resolve(strict=True)
    for entry in json.loads((campaign / 'manifest.json').read_text())['structures']:
        folder = campaign / entry['directory']
        report = structure_evidence(folder, entry)
        target = folder / 'tr_evidence.json'
        if target.exists():
            target.unlink()
        dump(target, report)
        print(target, flush=True)


if __name__ == '__main__':
    main()
