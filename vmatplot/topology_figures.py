"""Shared presentation settings for the topology notebook."""
from functools import lru_cache
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from vmatplot.output_settings import LINE_WIDTH

STYLES = {
    "article": dict(label=16, tick=14, title=20, panel=18, note=12, inside=12, outside=14, scale=1.),
    "thesis": dict(label=13, tick=11, title=13, panel=11, note=11, inside=11, outside=11, scale=1.),
}
SIZES = {
    "article": {"bands": (18, 9), "gap_maps": (18, 6), "wcc": (12, 5), "gap_zoom": (10, 6)},
    "thesis": {"bands": (8, 9), "gap_maps": (8, 7), "wcc": (8, 4), "gap_zoom": (8, 5)},
}


def topology_canvas(kind, version):
    s = STYLES[version]
    params = {"font.family": "serif", "mathtext.fontset": "cm", "text.usetex": False,
              "axes.labelsize": s["label"], "axes.titlesize": s["panel"],
              "xtick.labelsize": s["tick"], "ytick.labelsize": s["tick"],
              "legend.fontsize": s["inside"], "lines.linewidth": LINE_WIDTH, "axes.titlepad": 6., "pdf.fonttype": 42,
              "figure.facecolor": "white", "axes.facecolor": "white", "path.simplify": False}
    return SIZES[version][kind], 196, params, (s["title"], s["panel"]), "upper right"


def cache_topology_readers(module):
    names = ("extract_json", "extract_band_path", "extract_fermi_topology", "extract_gap_extrema",
             "extract_parity", "extract_wcc", "extract_reciprocal_2d", "extract_point_group_2d", "extract_direct_gap_grid")
    readers = []
    for name in names:
        reader = getattr(module, name)
        if not hasattr(reader, "cache_clear"):
            reader = lru_cache(maxsize=None)(reader)
            setattr(module, name, reader)
        readers.append(reader)
    return readers


def all_axes(fig):
    axes = list(fig.axes)
    for ax in axes:
        axes.extend(child for child in ax.child_axes if child not in axes)
    return axes


def style_topology_figure(fig, version):
    s = STYLES[version]
    axes = all_axes(fig)
    panels = [ax for ax in axes if ax.axison and not hasattr(ax, "_colorbar") and ax.get_label() != "<colorbar>"]
    legends = list(fig.legends)
    for ax in axes:
        ax.xaxis.label.set_fontsize(s["label"])
        ax.yaxis.label.set_fontsize(s["label"])
        ax.tick_params(which="both", direction="in", labelsize=s["tick"])
        for text in (ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()):
            text.set_fontsize(s["tick"])
        if ax in panels:
            ax.tick_params(which="both", top=True, right=True)
            title = ax.get_title()
            if title:
                panel = len(panels) > 1 or fig._suptitle is not None
                box = version == "thesis" and panel
                inside = box and (ax.get_ylabel().startswith("WCC") or getattr(ax, "_topology_band_panel", False))
                ax.set_title(title, fontsize=s["panel"] if panel else s["title"],
                             y=.97 if inside else 1., pad=0 if inside else 6.)
                if box:
                    ax.title.set(x=.025, ha="left", va="top" if inside else "bottom", zorder=10,
                                 bbox=dict(boxstyle="round", facecolor="white", alpha=.75,
                                           edgecolor=plt.rcParams["legend.edgecolor"]))
            for text in ax.texts:
                text.set_fontsize(s["note"])
            for line in ax.lines:
                line.set_linewidth(line.get_linewidth()*s["scale"])
                line.set_solid_capstyle("round"); line.set_dash_capstyle("round")
                line.set_solid_joinstyle("round"); line.set_dash_joinstyle("round")
            for collection in ax.collections:
                if isinstance(collection, LineCollection):
                    collection.set_linewidths(collection.get_linewidths()*s["scale"])
        legend = ax.get_legend()
        if legend:
            legends.append(legend)
    for legend in legends:
        ax = legend.axes
        outside = ax is None or not ax.axison
        if ax is not None and ax.axison:
            b = legend.get_bbox_to_anchor().transformed(ax.transAxes.inverted())
            outside = b.x0 < 0 or b.x1 > 1 or b.y0 < 0 or b.y1 > 1
        for text in [*legend.get_texts(), legend.get_title()]:
            text.set_fontsize(s["outside"] if outside else s["inside"])
        legend.set_frame_on(True)
        legend.get_frame().set(facecolor="white", alpha=.75, edgecolor=plt.rcParams["legend.edgecolor"])
    if fig._suptitle:
        fig._suptitle.set_fontsize(s["title"])
    return fig


def draw_topology_figure(plotter, filename, *args, figure_version, readers=(), preview=True, **kwargs):
    """Read fresh inputs once per pair, redraw both layouts, and save both PDFs."""
    for reader in readers:
        reader.cache_clear()
    path = Path(filename).with_suffix(".pdf")
    paths = (path, path.with_name(path.stem + "_thesis.pdf"))
    previous = figure_version()
    try:
        for version, output in zip(("article", "thesis"), paths):
            figure_version(version)
            with plt.rc_context():
                plotter(*args, **kwargs)
                fig = style_topology_figure(plt.gcf(), version)
                output.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(output, dpi=196, metadata={"CreationDate": None})
                if preview:
                    plt.show()
                plt.close(fig)
    finally:
        figure_version(previous)
    return paths
