#### Output settings
# pylint: disable = C0103, C0114, C0116, C0301, C0321, R0913, R0914

import os
import math
import matplotlib.pyplot as plt

import matplotlib as mpl

mpl.rcParams["lines.solid_capstyle"] = "round"
mpl.rcParams["lines.dash_capstyle"]  = "round"
mpl.rcParams["lines.solid_joinstyle"] = "round"
mpl.rcParams["lines.dash_joinstyle"]  = "round"

def vasprun_directory(directory="."):
    """Find folders with complete vasprun.xml and print incomplete ones."""
    complete_folders = []
    incomplete_detected = False  # Flag to track if any incomplete files are detected

    for dirpath, _, filenames in os.walk(directory):
        if "vasprun.xml" in filenames:
            file_name_xml = os.path.join(dirpath, "vasprun.xml")

            # Check if vasprun.xml is complete
            try:
                with open(file_name_xml, "r", encoding="utf-8") as f:
                    # Check the last few lines for the closing tag
                    last_lines = f.readlines()[-10:]  # read the last 10 lines
                    for line in last_lines:
                        if "</modeling>" in line or "</vasp>" in line:
                            complete_folders.append(dirpath)
                            break
                    else:
                        print(f"vasprun.xml in {dirpath} is incomplete.")
                        incomplete_detected = True
            except IOError as e:
                print(f"Error reading {file_name_xml}: {e}")

    if not incomplete_detected:
        print("All works are finished.")

    return complete_folders

def _canvas_setting_article(*args):
    help_info = "Usage: canvas_setting(length, width, dpi, font)\n" + \
                "The default setting is length: 10, width: 6, dpi: 196, font: 'Serif'" + \
                "The return values are :" + \
                "\t[0]: size_setting"+\
                "\t[1]: dpi"+\
                "\t[2]: figure params (font)"+\
                "\t[3]: Suptitle and subtitle fontsize"+\
                "\t[4]: legend location"
    default_style = {"length": 10,
                     "width": 6,
                     "dpi": 256,
                     "font_style": "serif"}
    default_params = {"text.usetex": False,
                      "font.family": "serif",
                      "mathtext.fontset": "cm",
                      "axes.titlesize": 20,
                      "axes.labelsize": 16,
                      "xtick.labelsize": 14,
                      "ytick.labelsize": 14,
                      "legend.fontsize": 12,
                      "figure.facecolor": "w"}
    default_suptitle = 20
    default_subtitle = 18
    default_legend = "upper right"
    if len(args) == 0:
        return (default_style["length"],default_style["width"]), default_style["dpi"], default_params, (default_suptitle, default_subtitle), default_legend
    if len(args) == 1:
        if args[0] == "help":
            print(help_info)
            return
        else:
            return (args[0], 6), 196, default_params, (default_suptitle, default_subtitle), default_legend
    if len(args) == 2:
        return (args[0], args[1]), 196, default_params, (default_suptitle, default_subtitle), default_legend
    if len(args) == 3:
        return (args[0], args[1]), args[2], default_params, (default_suptitle, default_subtitle), default_legend
    if len(args) == 4:
        customized_params = {"text.usetex": False,
                             "font.family": args[3],
                             "mathtext.fontset": "cm",
                             "axes.titlesize": 20,
                             "axes.labelsize": 16,
                             "xtick.labelsize": 14,
                             "ytick.labelsize": 14,
                             "legend.fontsize": 12,
                             "figure.facecolor": "w"}
        return (args[0], args[1]), args[2], customized_params, (default_suptitle, default_subtitle), default_legend
    if len(args) == 6:
        customized_params = {"text.usetex": False,
                             "font.family": args[3],
                             "mathtext.fontset": "cm",
                             "axes.titlesize": 20,
                             "axes.labelsize": 16,
                             "xtick.labelsize": 14,
                             "ytick.labelsize": 14,
                             "legend.fontsize": 12,
                             "figure.facecolor": "w"}
        return (args[0],args[1]), args[2], customized_params, (args[4],args[5]), default_legend
    if len(args) == 7:
        customized_params = {"text.usetex": False,
                             "font.family": args[3],
                             "mathtext.fontset": "cm",
                             "axes.titlesize": 20,
                             "axes.labelsize": 16,
                             "xtick.labelsize": 14,
                             "ytick.labelsize": 14,
                             "legend.fontsize": 12,
                             "figure.facecolor": "w"}
        return (args[0],args[1]), args[2], customized_params, (args[4],args[5]), args[6]

## Figure versions: article (default) and thesis
# figure_version("thesis") switches every later figure to the thesis typography below (the thesis settings of
# o-B14_20241024), and save_figure() then writes <name>_thesis.pdf instead of <name>.pdf. Some plotting functions
# also change their layout for the thesis, e.g. a 1x3 row becomes a 2x2 grid with the legend in the fourth panel.
figure_settings = {"version": "article"}
thesis_params = {"axes.titlesize": 24, "axes.labelsize": 18, "xtick.labelsize": 18, "ytick.labelsize": 18,
                 "legend.fontsize": 12, "axes.titlepad": 10, "pdf.fonttype": 42}
article_params = {"axes.titlepad": 6.0, "pdf.fonttype": 3}      # matplotlib defaults, restored when switching back
thesis_titles = (24, 22)                                        # suptitle and subtitle fontsize

def figure_version(version=None):
    help_info = "Usage: figure_version(version)\n" + \
                "Without argument it returns the current figure version; \"article\" (default) or \"thesis\" sets it.\n" + \
                "The thesis version uses the thesis typography and saves files as <name>_thesis.pdf."
    if version is None:
        return figure_settings["version"]
    if version == "help":
        print(help_info)
        return None
    if version not in ("article", "thesis"):
        raise ValueError('figure_version() takes "article" or "thesis"')
    figure_settings["version"] = version
    return version

def canvas_setting(*args):
    # Article canvas and typography, replaced by the thesis typography when figure_version() is "thesis"
    settings = _canvas_setting_article(*args)
    if settings is None:
        return settings
    size, dpi, params, titles, legend = settings
    if figure_settings["version"] == "thesis":
        params = {**params, **thesis_params}
        if len(args) < 6:
            titles = thesis_titles
    else:
        params = {**params, **article_params}
    return size, dpi, params, titles, legend

def fit_thesis_text(fig=None, margin=0.01):
    # The thesis typography is applied to the article canvases, so a long title can leave the canvas and dense rotated
    # tick labels can overlap. Those texts are scaled down until they fit; all other text keeps the thesis sizes.
    fig = plt.gcf() if fig is None else fig
    try:
        renderer = fig.canvas.get_renderer()
    except AttributeError:
        return fig
    left, right = fig.bbox.x0 + margin*fig.bbox.width, fig.bbox.x1 - margin*fig.bbox.width
    titles = [fig._suptitle] + [title for ax in fig.axes for title in (ax.title, ax._left_title, ax._right_title)]
    for text in titles:
        if text is None or not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        align = text.get_horizontalalignment()
        if align == "left":
            room = right - box.x0
        elif align == "right":
            room = box.x1 - left
        else:
            room = 2*min((box.x0 + box.x1)/2 - left, right - (box.x0 + box.x1)/2)
        if 0 < room < box.width:
            text.set_fontsize(text.get_fontsize()*room/box.width)
    for ax in fig.axes:
        labels = [label for label in ax.get_xticklabels() if label.get_visible() and label.get_text()]
        if len(labels) < 2 or abs(math.sin(math.radians(labels[0].get_rotation()))) < 0.1:
            continue
        lower, upper = sorted(ax.get_xlim())
        positions = sorted(ax.transData.transform((x, ax.get_ylim()[0]))[0] for x in ax.get_xticks() if lower <= x <= upper)
        spacing = min(b - a for a, b in zip(positions, positions[1:])) if len(positions) > 1 else None
        size = labels[0].get_fontsize()
        needed = 1.1*size*fig.dpi/72/abs(math.sin(math.radians(labels[0].get_rotation())))
        if spacing and spacing < needed:
            ax.tick_params(axis="x", which="both", labelsize=size*spacing/needed)
    return fig

def save_figure(name, directory="figures", **kwargs):
    # Save the current figure as <directory>/<name>.pdf (article) or <directory>/<name>_thesis.pdf (thesis);
    # the PDF carries no creation date, so rerunning a notebook reproduces the file byte for byte
    suffix = "_thesis" if figure_settings["version"] == "thesis" else ""
    path = os.path.join(directory, f"{name}{suffix}.pdf")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if suffix:
        fit_thesis_text(plt.gcf())
    plt.savefig(path, metadata={"CreationDate": None}, **kwargs)
    return path

def color_sampling(color_family):
    help_info = "Usage: color_family(color_family)\n" + \
                "Input the name of color family will return a series colors." + \
                "Color families: Grey, Red, Orange, Yellow, Green, Cyan, Blue, Violet, Purple, Wine, Brown, Orbit\n" + \
                "Return values:\n" + \
                "color[0]: deep color\n" + \
                "color[1]: major color\n" + \
                "color[2]: shallow color\n" + \
                "color[3]: comparison color 1\n" + \
                "color[4]: comparison color 2\n" + \
                "color[5]: comparison color 3\n" + \
                "color[6]: comparison color 4\n" + \
                "color[7]: comparison color 5\n" + \
                "color[8]: comparison color 6\n" + \
                "color[9]: comparison color 7\n" + \
                "color[10]: comparison color 8\n"
    # Check if the user asked for help
    if color_family == "help":
        print(help_info)
        return

    color_set = []
    if color_family in ("Default", "default", "Normal", "normal", "Orbital", "orbital", "Orbitals", "orbitals"):
        color_set.append("#145AAA") # color[0]: Base
        color_set.append("#1478E1") # color[1]: total
        color_set.append("#14A0FF") # color[2]: Integral

        color_set.append("#8C64F0") # color[3]: s-orbital
        color_set.append("#D25ADC") # color[4]: px-orbital
        color_set.append("#F03C64") # color[5]: py-orbital
        color_set.append("#FA8C00") # color[6]: pz-orbital
        color_set.append("#FAC828") # color[7]
        color_set.append("#96C800") # color[8]: d-orbital
        color_set.append("#14AFAF") # color[9]: f-orbital
        color_set.append("#1E8CA0") # color[10]
        return color_set

    if color_family in ("Grey", "grey", "Gray", "grey"):
        color_set.append("#3C3C3C")
        color_set.append("#787878")
        color_set.append("#B4B4B4")

        color_set.append("#E1322D")
        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        return color_set

    if color_family in ("Silver", "silver"):
        color_set.append("#787D8C")
        color_set.append("#AAAFBE")
        color_set.append("#C8CDD7")

        color_set.append("#E1322D")
        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        return color_set

    if color_family in ("Red", "red"):
        color_set.append("#C81423")
        color_set.append("#E1322D")
        color_set.append("#FF644B")

        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        return color_set

    if color_family in ("Orange", "orange"):
        color_set.append("#EB731E")
        color_set.append("#FA8C00")
        color_set.append("#FFA03C")

        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        color_set.append("#F03C64")
        return color_set

    if color_family in ("Yellow", "yellow", "Gold", "gold"):
        color_set.append("#EBC31E")
        color_set.append("#FAC828")
        color_set.append("#FFD732")

        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        color_set.append("#F03C64")
        color_set.append("#FA8C00")
        return color_set
    
    if color_family in ("Lime", "lime"):
        color_set.append("#8CB423")
        color_set.append("#96BE2D")
        color_set.append("#A0C837")

        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        color_set.append("#F03C64")
        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        return color_set

    if color_family in ("Green", "green"):
        color_set.append("#238C4B")
        color_set.append("#28AF3C")
        color_set.append("#73C81E")

        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        color_set.append("#F03C64")
        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        return color_set

    if color_family in ("Cyan","cyan"):
        color_set.append("#1E7878")
        color_set.append("#19A0A0")
        color_set.append("#14AFAF")

        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        color_set.append("#F03C64")
        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        return color_set

    if color_family in ("Blue", "blue", "Azure", "azure"):
        color_set.append("#145AAA") # color[0]
        color_set.append("#1478E1") # color[1]
        color_set.append("#14A0FF") # color[2]

        color_set.append("#8C64F0") # color[3]
        color_set.append("#D25ADC") # color[4]
        color_set.append("#F03C64") # color[5]
        color_set.append("#FA8C00") # color[6]
        color_set.append("#FAC828") # color[7]
        color_set.append("#96BE2D") # color[8]
        color_set.append("#28AF3C") # color[9]
        color_set.append("#19A0A0") # color[10]
        return color_set

    if color_family in ("Violet", "violet"):
        color_set.append("#643CC3")
        color_set.append("#8C64E1")
        color_set.append("#AF96FF")

        color_set.append("#D25ADC")
        color_set.append("#F03C64")
        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        return color_set

    if color_family in ("Purple", "purple"):
        color_set.append("#AA3CB9")
        color_set.append("#D25ADC")
        color_set.append("#F078FF")

        color_set.append("#F03C64")
        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        return color_set

    if color_family in ("Wine", "wine"):
        color_set.append("#AA1E64")
        color_set.append("#C82364")
        color_set.append("#F03C64")

        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        return color_set

    if color_family in ("Brown", "brown"):
        color_set.append("#966450")
        color_set.append("#B47D50")
        color_set.append("#D29650")

        color_set.append("#FA8C00")
        color_set.append("#FAC828")
        color_set.append("#96BE2D")
        color_set.append("#28AF3C")
        color_set.append("#19A0A0")
        color_set.append("#1478E1")
        color_set.append("#8C64F0")
        color_set.append("#D25ADC")
        return color_set

    if color_family == "all_families":
        return ["Silver", "Grey", "Red", "Orange", "Yellow", "Lime", 
                "Green", "Cyan", "Blue", "Violet", "Purple", "Wine", "Brown", "Default"]
    
    else:
        color_set.extend([color_family] * 10)
        return color_set

def plot_color_families():
    color_families = color_sampling("all_families")
    all_colors = [color_sampling(family) for family in color_families]

    # Figure Settings
    fig_setting = canvas_setting(10,8)
    params = fig_setting[2]; plt.rcParams.update(params)
    plt.rcParams.update(params)

    plt.figure(figsize=fig_setting[0], dpi = fig_setting[1])
    plt.title("Color families")

    for row, color_row in enumerate(all_colors):
        for col, color in enumerate(color_row):
            plt.gca().add_patch(plt.Rectangle((col, row), 1, 1, color=color))
            plt.text(col + 0.5, row + 0.5, color, ha="center", va="center", fontsize=8, color="white")

    plt.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)
    max_length = max([len(colors) for colors in all_colors])
    plt.xlim(0, max_length)
    plt.ylim(0, len(all_colors))
    plt.xticks([])
    plt.yticks([])

    yaxis_offset = 0.5
    for i, label in enumerate(color_families):
        plt.text(-max_length*0.01, i + yaxis_offset, label, ha="right", va="center")

    plt.show()
    plt.tight_layout()
