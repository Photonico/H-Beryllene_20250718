#!/usr/bin/env python3
"""Colour versions of the nine published phonon dispersion figures (article and thesis).

Every figure is drawn exactly as in its source - 3.0_phonon_dispersion.ipynb,
3.5_phonon_dispersion_with_hydrogen_selected.ipynb and
3.9_phonon_dispersion_d3bj/plot_phonon_dispersion.py - except that branch i takes
colour i of vmatplot.phonon.PHONON_BRANCH_COLORS. The files are saved next to the
monochrome ones as figures/phonon/<name>_color.pdf and figures/phonon/<name>_thesis_color.pdf;
the monochrome PDFs are not touched.

    python scripts/plot_phonon_color.py              # write the 18 colour PDFs
    python scripts/plot_phonon_color.py --mono DIR   # redraw the monochrome PDFs into DIR (reproducibility check)
"""

import argparse
from pathlib import Path
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vmatplot.output_settings import figure_version, save_figure
from vmatplot.phonon import PHONON_BRANCH_COLORS, plot_phonons

# (file name, title, data directory, frequency range in THz, figure size set after plotting or None)
FIGURES = (
    # 3.0_phonon_dispersion.ipynb
    ("phonon_hcp-beryllium", "bulk hcp beryllium", "3.8_phonon_dispersion_pbe/a-Beryllium", [0, 30], None),
    ("phonon_alpha-beryllene", "α-beryllene", "3.8_phonon_dispersion_pbe/a-Beryllene", [0, 30], None),
    ("phonon_beta-beryllene", "β-beryllene", "3.8_phonon_dispersion_pbe/b-Beryllene_p3m1", [0, 30], None),
    ("phonon_bcc-beryllium", "bulk bcc beryllium", "3.8_phonon_dispersion_pbe/c0-Beryllium", [-5, 30], None),
    ("phonon_cubic-beryllene-trilayer", "trilayer cubic beryllene", "3.8_phonon_dispersion_pbe/c3-Beryllene_trilayer", [0, 30], None),
    # 3.5_phonon_dispersion_with_hydrogen_selected.ipynb
    ("phonon_2H-alpha-beryllene", "2H-α-beryllene", "3.7_phonon_dispersion_with_hydrogen_pbe/a-Beryllene_hh", [-10, 60], None),
    ("phonon_1H-beta-beryllene", "1H-β-beryllene", "3.7_phonon_dispersion_with_hydrogen_pbe/b-Beryllene_b", [-10, 60], (9, 6)),
    ("phonon_2H-beta-beryllene", "2H-β-beryllene", "3.7_phonon_dispersion_with_hydrogen_pbe/b-Beryllene_bb", [-10, 60], (9, 6)),
    # 3.9_phonon_dispersion_d3bj/plot_phonon_dispersion.py
    ("phonon_2H-beta-beryllene_d3bj", "2H-β-beryllene, PBE+D3", "3.9_phonon_dispersion_d3bj/2H_beryllene_phonons_thesis_path", [-10, 60], (9, 6)),
)

def draw(title, directory, eigen_range, size, colored):
    matters = [[title, directory, "blue", "solid", 1.5, 1.0, 0]]
    plot_phonons(title, matters, eigen_range, None, branch_colors=PHONON_BRANCH_COLORS if colored else None)
    figure = plt.gcf()
    if size is not None:
        figure.set_size_inches(*size)
    return figure

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mono", type=Path, metavar="DIR", help="redraw the monochrome figures into DIR instead")
    args = parser.parse_args()
    os.chdir(ROOT)  # the data directories are relative to the repository, as in the notebooks
    directory = args.mono.resolve() if args.mono else ROOT / "figures" / "phonon"
    for version in ("article", "thesis"):
        figure_version(version)
        for name, title, data, eigen_range, size in FIGURES:
            figure = draw(title, data, eigen_range, size, colored=args.mono is None)
            print(save_figure(name, directory=directory, variant="" if args.mono else "_color"))
            plt.close(figure)
    figure_version("article")

if __name__ == "__main__":
    main()
