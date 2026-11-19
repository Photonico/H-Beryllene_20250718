#!/usr/bin/env python3
"""Plot the supplied 2H-β-beryllene PBE-D3(BJ) phonon bands."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from vmatplot.output_settings import figure_version, save_figure
from vmatplot.phonon import plot_phonons


def main():
    data_dir = HERE / "2H_beryllene_phonons"
    title = "2H-β-beryllene"
    # Match the 2H-β panel in 3.5_phonon_dispersion_with_hydrogen_selected.ipynb.
    matters = [[title, data_dir, "blue", "solid", 1.5, 1.0, 0]]
    for version in ("article", "thesis"):
        figure_version(version)
        with plt.rc_context():
            plot_phonons(title, matters, [-10, 60], legend_loc=None)
            figure = plt.gcf()
            figure.set_size_inches(9, 6)
            print(save_figure("phonon_2H-beta-beryllene_d3bj", directory=ROOT / "figures"))
            plt.close(figure)


if __name__ == "__main__":
    main()
