#! python
"""Draw the Eliashberg spectral functions and cumulative coupling strengths.

Run from any directory: python /path/to/10_superconductivity/supercond_a2f.py
Dependencies: matplotlib, numpy.

The curves are read directly from the EPW ph.a2f files of the 40x40x1 and
80x80x1 fine meshes, using the first phonon smearing (0.5 meV) as in the
reported coupling strengths. The cumulative coupling is the integrated column
written by EPW. ST-beryllene is trilayer cubic beryllene; beta-H-2Be and
beta-2H-2Be are 1H-beta-beryllene and 2H-beta-beryllene.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np

HERE = Path(__file__).resolve().parent
EPW = HERE / "5_electron_phonon_EPW"
BLUE = "#1464cd"
CYAN = "#0a8c8c"
WINE = "#C82364"
GREY = "#787878"
LINE_WIDTH = 1.5
THRESHOLD = 40 * 0.12398419843320026  # 40 cm^-1 in meV

# Same presentation settings as supercond.py; standalone by design.
STYLE = {
    "font.family": "serif", "mathtext.fontset": "cm", "text.usetex": False,
    "axes.titlesize": 18, "axes.labelsize": 16, "axes.titlepad": 6,
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 12,
    "lines.linewidth": LINE_WIDTH, "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "round", "lines.solid_joinstyle": "round",
    "lines.dash_joinstyle": "round", "figure.facecolor": "white",
    "axes.facecolor": "white", "pdf.fonttype": 42,
    "path.simplify": False,
}

# (folder, title, colour, frequency limit, frequency tick, a2F limit, lambda limit)
STRUCTURES = (
    ("ST-beryllene", "Trilayer cubic beryllene", CYAN, 100, 20, 0.5, 1.2),
    ("beta-H-2Be", "1H-β-beryllene", BLUE, 220, 50, 0.7, 0.8),
    ("beta-2H-2Be", "2H-β-beryllene", WINE, 220, 50, 1.8, 3.5),
)
MESHES = (("2_80x80", "solid"), ("1_coarse_to_40x40", "dashed"))

def read_a2f(path):
    """Return frequency (meV), a2F, and cumulative lambda for the first smearing."""
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 21:
            continue
        try:
            rows.append([float(value) for value in parts])
        except ValueError:
            continue
    data = np.asarray(rows)
    if data.ndim != 2 or not np.isfinite(data).all():
        raise ValueError(f"Missing or non-finite data in {path}")
    return data[:, 0], data[:, 1], data[:, 11]

def legend(ax, handles, loc):
    return ax.legend(handles=handles, loc=loc, frameon=True, fancybox=True,
                     facecolor="white", edgecolor="#CCCCCC", framealpha=.75,
                     handlelength=1.6, handletextpad=.6, borderpad=.4)

def make_figure():
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(3, 2, figsize=(12, 14), dpi=196)
        fig.subplots_adjust(left=.085, right=.98, bottom=.055, top=.965,
                            wspace=.25, hspace=.36)
        letters = iter("abcdef")
        for row, (folder, title, color, wmax, wtick, fmax, lmax) in enumerate(STRUCTURES):
            ax_f, ax_l = axes[row]
            handles = []
            for mesh, style in MESHES:
                w, a2f, lam = read_a2f(EPW / folder / mesh / "ph.a2f")
                line_color = color if style == "solid" else GREY
                ax_f.plot(w, a2f, color=line_color, ls=style, zorder=3)
                ax_l.plot(w, lam, color=line_color, ls=style, zorder=3)
                label = r"$80\times80\times1$" if mesh.endswith("80x80") else r"$40\times40\times1$"
                handles.append(Line2D([], [], color=line_color, ls=style, label=label))
            for ax, ymax, ylabel in ((ax_f, fmax, r"$\alpha^2F(\omega)$"), (ax_l, lmax, r"$\lambda(\omega)$")):
                # Highlight the low-frequency window below 40 cm^-1.
                ax.axvspan(0, THRESHOLD, color="#698cd7", alpha=.20, lw=0, zorder=0)
                ax.set(xlim=(0, wmax), ylim=(0, ymax), ylabel=ylabel)
                ax.set_title(f"({next(letters)}) {title}")
                ax.xaxis.set_major_locator(MultipleLocator(wtick))
                ax.xaxis.set_minor_locator(AutoMinorLocator(5))
                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
                ax.tick_params(which="both", direction="in", top=True, right=True)
                ax.tick_params(which="major", length=5)
                ax.tick_params(which="minor", length=2.5)
                if row == len(STRUCTURES) - 1:
                    ax.set_xlabel(r"$\hbar\omega$ (meV)")
            legend(ax_l, handles, "lower right")
        fig.align_ylabels(axes[:, 0])
        fig.align_ylabels(axes[:, 1])
        return fig

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE.parent / "figures/superconductivity_a2f.pdf")
    args = parser.parse_args()
    fig = make_figure()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        fig.savefig(args.output, metadata={"CreationDate": None})
    plt.close(fig)
    print(args.output.resolve())

if __name__ == "__main__":
    main()
