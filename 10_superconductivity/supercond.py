#! python
"""Redraw the four superconductivity panels from supercond_data.xlsx.

Run from any directory: python /path/to/10_superconductivity/supercond.py
Dependencies: matplotlib, numpy, openpyxl.

The draft's ST-beryllene is trilayer cubic beryllene; beta-H-2Be and
beta-2H-2Be are called 1H-beta-beryllene and 2H-beta-beryllene in this repo.
The supplied fit samples are drawn directly, without refitting or extending
them to the anisotropic critical temperatures. Diamonds are extrapolated
anisotropic limits, not the final values of the iteration traces.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, LogLocator, MultipleLocator, NullFormatter
import numpy as np
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
BLUE = "#1478E1"
CYAN = "#19A0A0"
WINE = "#C82364"
GREY = "#787878"
LINE_WIDTH = 1.5

# Local copies of the repository's article presentation settings; standalone
# by design, with no imports from or changes to vmatplot.
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

def read_pair(sheet, x_column, y_column, first_row=5):
    """Read one independent series; zero is data and blank tails are absent."""
    values = []
    for row in sheet.iter_rows(min_row=first_row, values_only=True):
        x, y = row[x_column - 1], row[y_column - 1]
        if x is None and y is None:
            continue
        if y is None and first_row == 4:
            continue  # Iteration columns have different numbers of samples.
        if x is None or y is None:
            raise ValueError(f"Incomplete pair in {sheet.title}, columns {x_column}/{y_column}")
        values.append((float(x), float(y)))
    data = np.asarray(values)
    if data.ndim != 2 or data.shape[1] != 2 or not np.isfinite(data).all():
        raise ValueError(f"Missing or non-finite data in {sheet.title}")
    return data.T

def isotropic(ax, sheet, raw_columns, fit_columns, color, marker, label):
    temperature, gap = read_pair(sheet, *raw_columns)
    fit_temperature, fit_gap = read_pair(sheet, *fit_columns)
    ax.plot(fit_temperature, fit_gap, color=color, zorder=2)
    ax.plot(temperature, gap, ls="none", marker=marker, ms=6.5, color=color,
            mec="white", mew=.8, zorder=3)
    return Line2D([], [], color=color, marker=marker, ms=6.5, mec="white", mew=.8, label=label)

def anisotropic(ax, sheet, columns, color, label):
    temperature, gap = read_pair(sheet, *columns)
    ax.plot(temperature, gap, ls="none", marker="D", ms=8, mfc="none",
            mec=color, mew=LINE_WIDTH, zorder=4)
    return Line2D([], [], ls="none", marker="D", ms=8, mfc="none", mec=color,
                  mew=LINE_WIDTH, label=label)

def critical_window(ax, sheet, columns, color, estimate):
    lower, upper = (float(sheet.cell(5, column).value) for column in columns)
    ax.axvspan(lower, upper, color=color, alpha=.12, lw=0, zorder=0)
    ax.text(upper + .35, .57, rf"anisotropic $T_\mathrm{{c}} \approx {estimate}$ K",
            transform=ax.get_xaxis_transform(), rotation=90, ha="center", va="center",
            fontsize=12, color="#3C3C3C")

def legend(ax, handles=None, loc="upper right"):
    return ax.legend(handles=handles, loc=loc, frameon=True, fancybox=True,
                     facecolor="white", edgecolor="#CCCCCC", framealpha=.75,
                     handlelength=1.6, handletextpad=.6, borderpad=.4)

def iteration_panel(ax, sheet, temperatures, colors, ymax):
    for column, (temperature, color) in enumerate(zip(temperatures, colors), start=2):
        iteration, gap = read_pair(sheet, 1, column, first_row=4)
        if np.any(gap <= 0):
            raise ValueError(f"Non-positive gap on logarithmic axis in {sheet.title}")
        ax.semilogy(iteration, gap, color=color, label=f"{temperature} K")
    ax.set_xlabel("Migdal–Eliashberg iteration")
    ax.set_ylim(1e-2 if ymax == 3 else 1e-3, ymax)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    legend(ax, loc="lower right")


def make_figure(workbook_path):
    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    try:
        with plt.rc_context(STYLE):
            fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=196)
            fig.suptitle("Superconductivity in pristine and hydrogen-functionalized beryllenes",
                         fontsize=20, y=.98)
            fig.subplots_adjust(left=.085, right=.98, bottom=.075, top=.905,
                                wspace=.25, hspace=.28)
            for ax in axes.flat:
                ax.set_ylabel(r"$\Delta(i\omega_1)$ (meV)")
                ax.tick_params(which="both", direction="in", top=True, right=True)
                ax.tick_params(which="major", length=5)
                ax.tick_params(which="minor", length=2.5)

            ax = axes[0, 0]
            sheet = workbook["Fig. 1a"]
            handles = [isotropic(ax, sheet, (7, 8), (10, 11), BLUE, "o", "1H-β-beryllene, isotropic"),
                       isotropic(ax, sheet, (1, 2), (4, 5), CYAN, "s", "Trilayer cubic beryllene, isotropic"),
                       anisotropic(ax, sheet, (13, 14), BLUE, "1H-β-beryllene, anisotropic")]
            critical_window(ax, sheet, (17, 18), BLUE, "14.2")
            ax.set(xlim=(0, 16), ylim=(0, 3), xlabel=r"$T$ (K)")
            ax.set_title("(a) Trilayer cubic and 1H-β-beryllene")
            ax.text(1.1, .8, r"$T_\mathrm{c}^{\mathrm{iso}}=3.7$ K", color=CYAN, fontsize=12)
            ax.text(8.4, .25, r"$T_\mathrm{c}^{\mathrm{iso}}=11.6$ K", color=BLUE, fontsize=12)
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            legend(ax, handles, loc="upper left")

            ax = axes[0, 1]
            iteration_panel(ax, workbook["Fig. 1b"], ("14.0", "14.1", "14.2", "14.5"),
                            ("#145AAA", BLUE, GREY, "#FA8C00"), 3)
            ax.set_xlim(0, 5200)
            ax.set_title("(b) 1H-β-beryllene: anisotropic")

            ax = axes[1, 0]
            sheet = workbook["Fig. 1c"]
            handles = [isotropic(ax, sheet, (1, 2), (4, 5), WINE, "o", "2H-β-beryllene, isotropic"),
                       anisotropic(ax, sheet, (7, 8), WINE, "2H-β-beryllene, anisotropic")]
            critical_window(ax, sheet, (11, 12), WINE, "31.7")
            ax.set(xlim=(0, 35), ylim=(0, 6), xlabel=r"$T$ (K)")
            ax.set_title("(c) 2H-β-beryllene")
            ax.text(20, .45, r"$T_\mathrm{c}^{\mathrm{iso}}=28.3$ K", color=WINE, fontsize=12)
            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            legend(ax, handles, loc="lower left")

            ax = axes[1, 1]
            iteration_panel(ax, workbook["Fig. 1d"], ("30", "31", "32", "34"),
                            ("#145AAA", BLUE, "#FA8C00", "#C81423"), 10)
            ax.set_xlim(0, 3100)
            ax.set_title("(d) 2H-β-beryllene: anisotropic")
            fig.align_ylabels(axes[:, 0])
            fig.align_ylabels(axes[:, 1])
            return fig
    finally:
        workbook.close()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "supercond_data.xlsx")
    parser.add_argument("--output", type=Path, default=HERE.parent / "figures/superconductivity.pdf")
    args = parser.parse_args()
    fig = make_figure(args.data)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        fig.savefig(args.output, metadata={"CreationDate": None})
    plt.close(fig)
    print(args.output.resolve())

if __name__ == "__main__":
    main()
