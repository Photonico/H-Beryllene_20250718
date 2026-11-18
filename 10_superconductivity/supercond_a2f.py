#! python
"""Draw the Eliashberg spectral functions and cumulative coupling strengths.

Run from any directory: python /path/to/10_superconductivity/supercond_a2f.py
It writes figures/superconductivity_a2f.pdf, the one version used by both the paper and the
thesis. The columns are titled above the top row with a clear gap to the panels, the subplot
titles sit in rounded boxes inside the upper-left corners of the panels, and an upper y limit
is raised only as far as needed to keep the curves clear of the box. The two meshes are solid
lines, the 40x40x1 reference in grey beneath the 80x80x1 curve.
Dependencies: matplotlib, numpy.

The curves are read directly from the EPW ph.a2f files of the 40x40x1 and
80x80x1 fine meshes, using the first phonon smearing (0.5 meV) as in the
reported coupling strengths. The cumulative coupling is the integrated column
written by EPW. ST-beryllene is trilayer cubic beryllene; beta-H-2Be and
beta-2H-2Be are 1H-beta-beryllene and 2H-beta-beryllene.
"""
import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np

HERE = Path(__file__).resolve().parent
EPW = HERE / "5_electron_phonon_EPW"
# One colour per material, shared with the other superconductivity figure and the topology gap-zoom figure
TRILAYER = "#fa8c00"   # trilayer cubic beryllene
BETA_1H = "#1478e1"    # 1H-β-beryllene
BETA_2H = "#19a0a0"    # 2H-β-beryllene
GREY = "#787878"
LINE_WIDTH = 2.0  # thicker than supercond.py: the dense a2F spectra read better at print size
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
    ("ST-beryllene", "Trilayer cubic beryllene", TRILAYER, 100, 20, 0.5, 1.2),
    ("beta-H-2Be", "1H-β-beryllene", BETA_1H, 220, 50, 0.7, 0.8),
    ("beta-2H-2Be", "2H-β-beryllene", BETA_2H, 220, 50, 1.8, 3.5),
)
# (folder, legend label, drawn in the material colour); the 40x40x1 reference is grey
MESHES = (("2_80x80", r"$80\times80\times1$", True), ("1_coarse_to_40x40", r"$40\times40\times1$", False))

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

def panel_label(ax, text):
    """Subplot title inside the upper-left corner, styled like the chapter 6 thesis figures."""
    return ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=(8, -8),
                       textcoords="offset points", ha="left", va="top", zorder=6,
                       fontsize=STYLE["axes.titlesize"],
                       bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=.75,
                                 edgecolor="#B4B4B4", linewidth=1.5))

def legend(ax, handles, loc, **kwargs):
    return ax.legend(handles=handles, loc=loc, frameon=True, fancybox=True,
                     facecolor="white", edgecolor="#CCCCCC", framealpha=.75,
                     handlelength=1.6, handletextpad=.6, borderpad=.4, **kwargs)

def clear_panel_label(fig, ax, label, gap=4):
    """Raise the upper y limit, by whole minor-tick steps, just enough for the curves to clear the label."""
    box = label.get_bbox_patch().get_window_extent(fig.canvas.get_renderer())
    pad = gap * fig.dpi / 72
    lower, upper = ax.get_ylim()
    left, right = ax.transData.inverted().transform([(box.x0 - pad, 0), (box.x1 + pad, 0)])[:, 0]
    peak = max(max(np.interp([left, right], line.get_xdata(), line.get_ydata()).max(),
                   np.max(line.get_ydata()[(line.get_xdata() >= left) & (line.get_xdata() <= right)], initial=-np.inf))
               for line in ax.lines)
    fraction = (box.y0 - pad - ax.bbox.y0) / ax.bbox.height
    needed = lower + (peak - lower) / fraction
    if needed > upper:
        ticks = ax.yaxis.get_majorticklocs()
        step = (ticks[1] - ticks[0]) / 5
        ax.set_ylim(lower, math.ceil(needed / step - 1e-9) * step)

COLUMN_TITLES = ("Eliashberg spectral function", "Cumulative coupling strength")
COLUMN_TITLE_PAD = 16  # points of clear space between the column titles and the top row

def make_figure():
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(3, 2, figsize=(12, 14), dpi=196)
        fig.subplots_adjust(left=.085, right=.98, bottom=.055, top=.945,
                            wspace=.25, hspace=.16)
        letters = iter("abcdef")
        labels = []
        for row, (folder, title, color, wmax, wtick, fmax, lmax) in enumerate(STRUCTURES):
            ax_f, ax_l = axes[row]
            handles = []
            for mesh, label, fine in MESHES:
                w, a2f, lam = read_a2f(EPW / folder / mesh / "ph.a2f")
                # Solid lines; the grey 40x40x1 reference lies beneath the 80x80x1 curve
                line_color, order = (color, 3) if fine else (GREY, 2)
                ax_f.plot(w, a2f, color=line_color, zorder=order)
                ax_l.plot(w, lam, color=line_color, zorder=order)
                handles.append(Line2D([], [], color=line_color, label=label))
            for ax, ymax, ylabel in ((ax_f, fmax, r"$\alpha^2F(\omega)$"), (ax_l, lmax, r"$\lambda(\omega)$")):
                # Highlight the low-frequency window below 40 cm^-1.
                ax.axvspan(0, THRESHOLD, color="#698cd7", alpha=.25, lw=0, zorder=0)
                ax.set(xlim=(0, wmax), ylim=(0, ymax), ylabel=ylabel)
                labels.append((ax, panel_label(ax, f"({next(letters)}) {title}")))
                ax.xaxis.set_major_locator(MultipleLocator(wtick))
                ax.xaxis.set_minor_locator(AutoMinorLocator(5))
                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
                ax.tick_params(which="both", direction="in", top=True, right=True)
                ax.tick_params(which="major", length=5)
                ax.tick_params(which="minor", length=2.5)
                if row == len(STRUCTURES) - 1:
                    ax.set_xlabel(r"$\hbar\omega$ (meV)")
            legend(ax_l, handles, "lower right", fontsize=STYLE["xtick.labelsize"])
        for ax, column_title in zip(axes[0], COLUMN_TITLES):
            ax.set_title(column_title, pad=COLUMN_TITLE_PAD)
        fig.align_ylabels(axes[:, 0])
        fig.align_ylabels(axes[:, 1])
        fig.canvas.draw()
        for ax, label in labels:
            clear_panel_label(fig, ax, label)
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
