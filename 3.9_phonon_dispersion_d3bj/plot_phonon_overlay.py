#!/usr/bin/env python3
"""Overlay the PBE and PBE-D3(BJ) phonon bands of 2H-β-beryllene.

PBE: 3.7_phonon_dispersion_with_hydrogen_pbe/b-Beryllene_bb (blue, alpha 0.75, below).
PBE+D3: 2H_beryllene_phonons_thesis_path (purple, drawn on top).
Both follow Γ–X–K–Γ–Y–K′–M–Γ. The D3 cell is 2.4 % smaller, so each D3 segment is mapped
linearly onto the corresponding PBE segment; the high-symmetry points then coincide.
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from vmatplot.output_settings import figure_version, save_figure
from vmatplot.phonon import plot_phonons

PBE_DIR = ROOT / "3.7_phonon_dispersion_with_hydrogen_pbe" / "b-Beryllene_bb"
D3_DIR = HERE / "2H_beryllene_phonons_thesis_path"
PBE_COLOR = "#1e78e1"
D3_COLOR = "#aa3cb9"
PBE_ALPHA = 0.75
LINE_WIDTH = 1.5


def read_segments(directory):
    """Return per-segment (distances, frequencies[nq, nband]) and the segment labels from band.yaml."""
    loader = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
    with open(directory / "band.yaml", encoding="utf-8") as handle:
        data = yaml.load(handle, Loader=loader)
    points = data["phonon"]
    segments, start = [], 0
    for count in data["segment_nqpoint"]:
        block = points[start:start + count]
        start += count
        distance = np.array([p["distance"] for p in block])
        frequency = np.array([[b["frequency"] for b in p["band"]] for p in block])
        segments.append((distance, frequency))
    return segments, data["labels"]


def main():
    d3_segments, d3_labels = read_segments(D3_DIR)
    _, pbe_labels = read_segments(PBE_DIR)
    if d3_labels != pbe_labels:
        raise ValueError(f"paths differ: {pbe_labels} vs {d3_labels}")
    for version in ("article", "thesis"):
        figure_version(version)
        with plt.rc_context():
            # PBE through the standard routine (axes, ticks, zero line), every branch in PBE_COLOR at zorder 4.
            plot_phonons("2H-β-beryllene", [["PBE", PBE_DIR, "blue", "solid", LINE_WIDTH, PBE_ALPHA, 0]], [-10, 60],
                         legend_loc=None, branch_colors=[PBE_COLOR] * 12)
            figure, axis = plt.gcf(), plt.gca()
            ticks = list(axis.get_xticks())
            if len(ticks) != len(d3_segments) + 1:
                raise ValueError(f"{len(ticks)} ticks for {len(d3_segments)} segments")
            # PBE+D3 on top: map each segment onto the PBE tick interval.
            pbe_line = next(line for line in axis.get_lines() if line.get_label() == "PBE")
            for (distance, frequency), x0, x1 in zip(d3_segments, ticks[:-1], ticks[1:]):
                x = x0 + (distance - distance[0]) / (distance[-1] - distance[0]) * (x1 - x0)
                d3_lines = axis.plot(x, frequency, color=D3_COLOR, lw=pbe_line.get_linewidth(), zorder=5)
            d3_lines[0].set_label("PBE+D3")
            handles = [pbe_line, d3_lines[0]]
            axis.legend(handles=handles, loc="upper right", ncol=2, frameon=True,
                        fontsize=plt.rcParams["xtick.labelsize"])
            figure.set_size_inches(9, 6)
            print(save_figure("phonon_2H-beta-beryllene_pbed3", directory=ROOT / "figures" / "phonon"))
            plt.close(figure)


if __name__ == "__main__":
    main()
