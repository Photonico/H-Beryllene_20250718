#!/usr/bin/env python3
"""Recompute the 2H-β-beryllene PBE-D3(BJ) phonon bands along the thesis path.

The force sets in ../2H_beryllene_phonons are used unchanged. Their band.yaml follows
S–X–Γ–S–X1–Γ–Y–X1, whose labels differ from the Brillouin-zone convention of the PBE
dispersions (3.7_phonon_dispersion_with_hydrogen_pbe/b-Beryllene_bb/band.conf): there,
X = b1/2, Y = b2/2, M = (b2 - b1)/2, and K, K′ = -K + b2 are the zone corners.
This script writes band.conf and band.yaml here for the PBE path Γ–X–K–Γ–Y–K′–M–Γ,
with K recomputed for the D3(BJ) cell. The other settings follow the original band step:
the unit cell as the primitive cell, no force-constant symmetrization, the default
symmetry tolerance, 201 points per segment, and band connection.

    python make_band_yaml.py
"""

from pathlib import Path

import numpy as np
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "2H_beryllene_phonons"
DIM = (12, 12, 1)
NPOINTS = 201
LABELS = ["$\\Gamma$", "X", "K", "$\\Gamma$", "Y", "K′", "M", "$\\Gamma$"]


def zone_corner(lattice):
    """Corner K = k1 b1 + k2 b2 where the Bragg lines of b1 and b2 meet."""
    reciprocal = np.linalg.inv(np.asarray(lattice, dtype=float)).T[:2, :2]
    metric = reciprocal @ reciprocal.T
    return np.linalg.solve(metric, np.diag(metric) / 2)


def main():
    phonon = phonopy.load(
        unitcell_filename=SOURCE / "POSCAR",
        supercell_matrix=DIM,
        primitive_matrix=np.eye(3),
        force_sets_filename=SOURCE / "FORCE_SETS",
        symmetrize_fc=False,
        log_level=0,
    )
    k1, k2 = zone_corner(phonon.unitcell.cell)
    corners = [[0, 0, 0], [0.5, 0, 0], [k1, k2, 0], [0, 0, 0], [0, 0.5, 0], [-k1, 1 - k2, 0], [-0.5, 0.5, 0], [0, 0, 0]]
    qpoints, connections = get_band_qpoints_and_path_connections([corners], npoints=NPOINTS)
    phonon.run_band_structure(qpoints, path_connections=connections, labels=LABELS, is_band_connection=True)
    phonon.write_yaml_band_structure(filename=HERE / "band.yaml")
    band = "  ".join(" ".join(f"{x:.8f}" for x in point) for point in corners)
    (HERE / "band.conf").write_text(
        f"DIM = {' '.join(map(str, DIM))}\n"
        f"BAND = {band}\n"
        f"BAND_POINTS = {NPOINTS}\n"
        f"BAND_LABELS = {' '.join(LABELS)}\n"
        "BAND_CONNECTION = .TRUE.\n",
        encoding="utf-8",
    )
    print(f"K = ({k1:.8f}, {k2:.8f}); wrote {HERE / 'band.yaml'} and {HERE / 'band.conf'}")


if __name__ == "__main__":
    main()
