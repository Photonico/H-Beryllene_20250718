#### Band topology: SOC Z2 screening results and figures
# pylint: disable = C0103, C0114, C0116, C0301, C0321, R0913, R0914

"""Read-only tables and figures for the 9.0_topology SOC campaign.

Reads the campaign manifest and the per-structure reports written by
9.0_topology/tools (parity, gap refinement, Wilson-loop WCC, time-reversal
evidence, spin screen) plus the VASP outputs they validated. Every Z2 value
is a conditional screening result for a fixed lowest-N spinor-band subspace;
nothing here certifies a global topological invariant.
"""

import json
import math
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.gridspec import GridSpec

from vmatplot.output_settings import canvas_setting

LABELS = {"alpha": "α", "beta": "β", "st": "ST", "alpha_2h": "2H-α",
          "beta_1h": "1H-β", "beta_2h": "2H-β"}
TRIMS = {"Gamma": (0.0, 0.0), "X": (0.5, 0.0), "Y": (0.0, 0.5), "M": (0.5, 0.5)}
TRIM_TEX = {"Gamma": "Γ", "X": "X", "Y": "Y", "M": "M"}

# House blue/orange plus aqua: validated all-pairs for CVD (worst ΔE 8.9);
# orange and aqua are below 3:1 on white, so they always carry direct labels.
SUBSPACE_COLORS = ("#1478E1", "#EB731E")
STRUCTURE_COLORS = {"alpha": "#1478E1", "beta": "#EB731E", "st": "#1BAF7A"}
CONTEXT_COLOR = "#B4B4B4"
FERMI_COLOR = "#643CC3"
INK = "#3C3C3C"
MUTED = "#787878"
# One-hue sequential ramp (blue 100 -> 700), reversed so small gaps are dark.
GAP_CMAP = LinearSegmentedColormap.from_list(
    "direct_gap", ["#0d366b", "#104281", "#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5",
                   "#5598e7", "#6da7ec", "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"])
GAP_NORM = LogNorm(vmin=1e-3, vmax=10.0)


# ---------------------------------------------------------------- data access

def campaign(topology=None):
    """Return (campaign directory, manifest entries with a 'folder' Path)."""
    if topology is None:
        cwd = Path.cwd().resolve()
        repo = next((p for p in (cwd, *cwd.parents) if (p / "9.0_topology/manifest.json").is_file()), None)
        if repo is None:
            raise FileNotFoundError("Run from inside the project, or pass the 9.0_topology directory.")
        topology = repo / "9.0_topology"
    topology = Path(topology).resolve()
    entries = json.loads((topology / "manifest.json").read_text())["structures"]
    for entry in entries:
        entry["folder"] = topology / entry["directory"]
    return topology, entries


def read_json(path):
    path = Path(path) if path is not None else None
    return json.loads(path.read_text()) if path is not None and path.is_file() else None


def latest(folder, pattern):
    found = sorted(Path(folder).glob(pattern))
    return found[-1] if found else None


def subspace_sizes(entry):
    n = entry["expected_nelect"]
    return [n] if n % 2 == 0 else [n - 1, n + 1]


def fermi_level(folder):
    return float((Path(folder) / "scf" / "DOSCAR").read_text().splitlines()[5].split()[3])


def gap_extrema(folder):
    """Final sampled extrema over all stages; True if the adaptive zoom ran."""
    zoom = read_json(Path(folder) / "gap_refine_3_summary.json")
    if zoom:
        return zoom["sampled_extrema_after"], True
    return read_json(Path(folder) / "gap_refinement_summary.json")["final_sampled_extrema"], False


def parity_report(folder):
    analysis = latest(Path(folder) / "trim", "parity-analysis-*")
    return read_json(analysis / "parity_summary.json") if analysis else None


def wcc_status(folder):
    return read_json(latest(Path(folder) / "wcc_direct_v3", "status_*.json"))


def lattice(folder):
    return np.array(read_json(Path(folder) / "trim" / "topology_validation.json")["structure"]["lattice"])


def reciprocal_2d(folder):
    return (2 * np.pi * np.linalg.inv(lattice(folder)).T)[:2, :2]


def trim_name(k):
    for name, target in TRIMS.items():
        if all(abs((k[i] - target[i] + 0.5) % 1 - 0.5) < 1e-6 for i in (0, 1)):
            return name
    return None


def sign(value):
    return "+" if value > 0 else "−"


def gap_text(ev):
    return "%.3g eV" % ev if ev >= 0.1 else "%.3g meV" % (1000 * ev)


# ---------------------------------------------------------------- tables

def show_table(header, rows):
    """Render a Markdown table in Jupyter (plain text elsewhere)."""
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in rows]
    text = "\n".join(lines)
    try:
        from IPython.display import Markdown, display
        display(Markdown(text))
    except ImportError:
        print(text)


def summary_table(topology=None):
    _, entries = campaign(topology)
    rows = []
    for entry in entries:
        folder = entry["folder"]
        parity, wcc = parity_report(folder), wcc_status(folder)
        if parity:
            invariant = "ν = %d (parity)" % parity["conditional_fu_kane_nu"]
        elif wcc:
            invariant = ", ".join("N=%s: Z₂ = %d (WCC)" % (n, m["z2"]) for n, m in wcc["manifolds"].items())
        else:
            invariant = "not validated"
        fermi = fermi_level(folder)
        gaps, _ = gap_extrema(folder)
        insulating = all(g["valence_max_ev"] < fermi < g["conduction_min_ev"] for g in gaps)
        spin = read_json(folder / "spin_screen" / "spin_screen_summary.json")
        tr = read_json(folder / "tr_evidence.json")
        rows.append([LABELS[entry["id"]], entry["directory"], entry["expected_nelect"], invariant,
                     ", ".join(gap_text(g["direct_gap_ev"]) for g in gaps),
                     "insulator" if insulating else "metal",
                     "%.1e" % tr["magnetisation_density"]["max_abs_m_muB_per_A3"] if tr else "pending",
                     ("collapsed" if spin["all_seeds_collapsed"] else "MOMENT KEPT") if spin else "pending"])
    return (["Structure", "Directory", "NELECT", "Conditional invariant", "Min. direct gap E(N+1)−E(N)",
             "Neutral E_F", "max \\|m(r)\\| (μB/Å³)", "Spin seeds"], rows)


def parity_tables(topology=None):
    """(parity table, float64 wavefunction-validation table) for centrosymmetric structures."""
    _, entries = campaign(topology)
    rows, checks = [], []
    for entry in entries:
        if not entry["inversion_expected"]:
            continue
        report = parity_report(entry["folder"])
        if report is None or "wavefunction_validation" not in report:
            rows.append([LABELS[entry["id"]], "not validated"] + [""] * 6)
            continue
        cells = []
        for name in TRIMS:
            point, check = report["trims"][name], report["wavefunction_validation"][name]
            pairs = " ".join(sign(p) for block in point["blocks"]
                             for p in block["kramers_pair_parities_unordered_within_block"])
            above = check["next_block_above_N_parity"]
            cells.append("%s \\| %s (%s)" % (pairs, sign(above) if above else "?",
                                             gap_text(check["next_block_above_N_gap_ev"])))
        rows.append([LABELS[entry["id"]], report["occupied_spinor_bands"], *cells,
                     " ".join(sign(report["trims"][n]["delta"]) for n in TRIMS),
                     report["conditional_fu_kane_nu"]])
        values = report["wavefunction_validation"].values()
        checks.append([LABELS[entry["id"]], "%.1e" % max(v["lowdin_inversion_singular_value_max_error"] for v in values),
                       "%.1e" % max(v["plain_metric_opposite_parity_max_overlap"] for v in values),
                       "%.1e" % max(v["inversion_residual_max"] for v in values),
                       "%.1e" % max(v["plain_metric_same_parity_max_overlap"] for v in values),
                       ", ".join("%.3g" % v for v in report["irrep_plain_metric_orthogonality_messages"]) or "none"])
    return ((["Structure", "N", "Γ", "X", "Y", "M", "δ (Γ X Y M)", "ν"], rows),
            (["Structure", "Löwdin unitarity error", "Opposite-parity overlap", "Inversion residual",
              "Same-parity overlap (PAW metric)", "IrRep plain-metric message"], checks))


def gap_table(topology=None):
    _, entries = campaign(topology)
    rows = []
    for entry in entries:
        fermi = fermi_level(entry["folder"])
        gaps, _ = gap_extrema(entry["folder"])
        for g in gaps:
            rows.append([LABELS[entry["id"]], g["N"], gap_text(g["direct_gap_ev"]),
                         "(%.5f, %.5f)" % tuple(g["direct_gap_k"][:2]), g["direct_gap_stage"],
                         "%+.3f" % (g["valence_max_ev"] - fermi), "%+.3f" % (g["conduction_min_ev"] - fermi),
                         "%+.3f" % g["indirect_gap_ev"]])
    return (["Structure", "N", "Min. direct gap", "k (fractional)", "Stage", "max E_N − E_F (eV)",
             "min E_N+1 − E_F (eV)", "Indirect gap (eV)"], rows)


def zoom_table(topology=None):
    _, entries = campaign(topology)
    rows = []
    for entry in entries:
        zoom = read_json(entry["folder"] / "gap_refine_3_summary.json")
        if zoom is None:
            continue
        for index, basin in enumerate(zoom["basins"]):
            final = basin["levels"][-1]
            rows.append([LABELS[entry["id"]], "%d: (%.5f, %.5f)" % (index + 1, *basin["start"]),
                         len(basin["levels"]), " → ".join("%.3f" % (1000 * l["min_direct_gap_ev"]) for l in basin["levels"]),
                         "(%.7f, %.7f)" % tuple(final["min_k"]), "%.1e" % final["spacing"],
                         "yes" if final["min_interior"] else "no", "%+.3f" % (1000 * final["slope_bound_ev"])])
    return (["Structure", "Basin (start k)", "Levels", "Min. gap per level (meV)", "Final k", "Final spacing",
             "Interior", "Slope bound (meV)"], rows)


def wcc_table(topology=None):
    _, entries = campaign(topology)
    rows = []
    for entry in entries:
        status = wcc_status(entry["folder"])
        if status is None:
            continue
        for n, m in status["manifolds"].items():
            gap = min(stage["minimum_direct_gap_ev"] for stage in m["sampled_gaps"].values())
            rows.append([LABELS[entry["id"]], n, m["status"], m["z2"], gap_text(gap), len(m["line_positions"]),
                         m["topology_certified"]])
    return (["Structure", "Bands", "Status", "Conditional Z₂", "Min. sampled direct gap", "Wilson loops",
             "Certified"], rows)


def time_reversal_tables(topology=None):
    """(converged-state TR evidence table, spin-screen table)."""
    _, entries = campaign(topology)
    rows, spins = [], []
    for entry in entries:
        folder, label = entry["folder"], LABELS[entry["id"]]
        tr = read_json(folder / "tr_evidence.json")
        if tr:
            m, minus_k = tr["magnetisation_density"], tr.get("scf_minus_k")
            rows.append([label, "%.1e" % m["max_abs_m_muB_per_A3"], "%.1e" % m["integral_abs_m_muB"],
                         "%.1e / %.1e" % (m["paw_occupancy_max_abs_re_m"], m["paw_occupancy_max_abs_im_m"]),
                         "%.1e" % tr["trim_kramers_max_splitting_ev"],
                         "%.1e (%d pairs)" % (minus_k["max_abs_dE_ev"], minus_k["pairs"]) if minus_k else "—"])
        spin = read_json(folder / "spin_screen" / "spin_screen_summary.json")
        if spin is None:
            spins.append([label, "pending", "", "", "", ""])
            continue
        for seed, result in spin["seeds"].items():
            sites = result["site_moments_muB"] or [float("nan")]
            spins.append([label, seed, " ".join("%g" % v for v in result["initial_magmom_muB"]),
                          "%.4f" % result["mag_muB"], "%.4f" % max(abs(v) for v in sites),
                          "yes" if result["moment_collapsed"] else "NO"])
    return ((["Structure", "max \\|m(r)\\| (μB/Å³)", "∫\\|m\\| (μB)", "One-centre m: max \\|Re\\| / max \\|Im\\|",
              "TRIM Kramers splitting (eV)", "max \\|E(k)−E(−k)\\| (eV)"], rows),
            (["Structure", "Seed", "Initial moments (μB)", "Final total (μB)", "Final max \\|site\\| (μB)",
              "Collapsed"], spins))


# ---------------------------------------------------------------- band path and BZ grids

def band_path(folder):
    """Fixed-density SOC bands on the line-mode path: x (1/Å), E (eV), fractional k, ticks."""
    stage = Path(folder) / "bands"
    with h5py.File(stage / "vaspout.h5", "r") as stream:
        energies = stream["/results/electron_eigenvalues/eigenvalues"][()][0]
        kpoints = stream["/results/electron_eigenvalues/kpoint_coords"][()][:, :2]
    lines = stage.joinpath("KPOINTS").read_text().splitlines()
    per_segment = int(lines[1].split()[0])
    ends = [line.split() for line in lines[4:] if len(line.split()) >= 4]
    labels = [(float(v[0]), float(v[1]), v[3]) for v in ends]
    cart = kpoints @ reciprocal_2d(folder)
    step = np.linalg.norm(np.diff(cart, axis=0), axis=1)
    step[per_segment - 1::per_segment] = 0.0  # segment ends are repeated or jump
    x = np.concatenate([[0.0], np.cumsum(step)])
    ticks, tick_labels = [], []
    for segment in range(len(kpoints) // per_segment):
        for position, label in ((segment * per_segment, labels[2 * segment][2]),
                                ((segment + 1) * per_segment - 1, labels[2 * segment + 1][2])):
            if ticks and abs(x[position] - ticks[-1]) < 1e-9:
                if label != tick_labels[-1]:
                    tick_labels[-1] += "|" + label
                continue
            ticks.append(x[position])
            tick_labels.append(label)
    return x, energies, kpoints, ticks, tick_labels


def _point_group(folder, symprec=1e-5):
    import spglib
    lines = (Path(folder) / "scf" / "POSCAR").read_text().splitlines()
    scale = float(lines[1].split()[0])
    cell = np.array([[float(v) for v in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    counts = [int(v) for v in lines[6].split()]
    start = 8 if lines[7].strip()[0].upper() in "DC" else 9
    positions = np.array([[float(v) for v in lines[start + i].split()[:3]] for i in range(sum(counts))])
    numbers = [i for i, c in enumerate(counts) for _ in range(c)]
    rotations = spglib.get_symmetry((cell, positions, numbers), symprec=symprec)["rotations"]
    # In-plane action of every operation that does not mix z with the plane.
    blocks = {tuple(r[:2, :2].ravel()) for r in rotations if not r[2, :2].any() and not r[:2, 2].any()}
    return [np.array(b).reshape(2, 2) for b in sorted(blocks)]


def scf_direct_gap_grid(folder, n):
    """E_(n+1) - E_n on the full Gamma-centred SCF grid, unfolded with the crystal point group and k -> -k."""
    folder = Path(folder)
    mesh = [int(v) for v in (folder / "scf" / "KPOINTS").read_text().splitlines()[3].split()[:2]]
    with h5py.File(folder / "scf" / "vaspout.h5", "r") as stream:
        energies = stream["/results/electron_eigenvalues/eigenvalues"][()][0]
        kpoints = stream["/results/electron_eigenvalues/kpoint_coords"][()][:, :2]
    gap = energies[:, n] - energies[:, n - 1]
    grid = np.full(mesh, np.nan)
    worst = 0.0
    for rotation in _point_group(folder):
        for parity in (1, -1):
            image = parity * kpoints @ rotation
            index = np.rint(image * mesh).astype(int) % mesh
            if np.abs(image * mesh - np.rint(image * mesh)).max() > 1e-4:
                raise ValueError("Symmetry image is off the SCF grid: " + str(folder))
            known = ~np.isnan(grid[index[:, 0], index[:, 1]])
            if known.any():
                worst = max(worst, float(np.abs(grid[index[known, 0], index[known, 1]] - gap[known]).max()))
            grid[index[:, 0], index[:, 1]] = gap
    if np.isnan(grid).any() or worst > 1e-4:
        raise ValueError("Unfolding incomplete or inconsistent (%.2e eV): %s" % (worst, folder))
    return grid, mesh


def _images(folder, k):
    points = {tuple(np.round(((s * np.asarray(k[:2]) @ r) + 0.5) % 1 - 0.5, 9))
              for r in _point_group(folder) for s in (1, -1)}
    return np.array(sorted(points))


def _style(fontsize=10):
    params = dict(canvas_setting()[2])
    params.update({"axes.titlesize": fontsize + 2, "axes.labelsize": fontsize + 1,
                   "xtick.labelsize": fontsize, "ytick.labelsize": fontsize, "legend.fontsize": fontsize - 1})
    plt.rcParams.update(params)


def _finish(fig, save):
    if save:
        Path(save).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------- figures

def plot_topology_bands(topology=None, ids=None, window=(-11.0, 5.0), save=None):
    """SOC bands (E − E_F) with the lowest-N subspace coloured, and E_(N+1) − E_N along the path."""
    _, entries = campaign(topology)
    entries = [e for e in entries if ids is None or e["id"] in ids]
    _style()
    columns = 3
    rows = math.ceil(len(entries) / columns)
    fig = plt.figure(figsize=(4.4 * columns, 5.4 * rows))
    grid = GridSpec(2 * rows, columns, figure=fig, height_ratios=[3.2, 1.3] * rows, hspace=0.08, wspace=0.28)
    for index, entry in enumerate(entries):
        row, column = divmod(index, columns)
        ax = fig.add_subplot(grid[2 * row, column])
        axg = fig.add_subplot(grid[2 * row + 1, column], sharex=ax)
        folder, sizes = entry["folder"], subspace_sizes(entry)
        x, energies, kpoints, ticks, labels = band_path(folder)
        energies = energies - fermi_level(folder)
        lower = 0
        for color, n in zip(SUBSPACE_COLORS, sizes):
            for band in range(lower, n):
                ax.plot(x, energies[:, band], color=color, lw=1.6, zorder=3,
                        label="bands %d–%d" % (lower + 1, n) if band == lower else None)
            lower = n
        for band in range(lower, energies.shape[1]):
            ax.plot(x, energies[:, band], color=CONTEXT_COLOR, lw=1.0, zorder=2)
        ax.axhline(0.0, color=FERMI_COLOR, lw=1.0, ls="--", zorder=1)
        parity = parity_report(folder)
        for position in ticks:
            for a in (ax, axg):
                a.axvline(position, color=MUTED, lw=0.6, ls="--", zorder=0)
        if parity:
            _annotate_parities(ax, x, energies, kpoints, parity, sizes[0], x[-1])
        for color, n in zip(SUBSPACE_COLORS, sizes):
            axg.semilogy(x, np.maximum(energies[:, n] - energies[:, n - 1], 1e-6), color=color, lw=1.4)
        gaps, _ = gap_extrema(folder)
        for color, g in zip(SUBSPACE_COLORS, gaps):
            axg.axhline(g["direct_gap_ev"], color=color, lw=0.9, ls=":")
            axg.text(x[-1] * 0.99, g["direct_gap_ev"] * 1.6, "BZ min " + gap_text(g["direct_gap_ev"]),
                     ha="right", va="bottom", fontsize=8, color=INK)
        if parity:
            title = "%s: ν = %d (conditional)" % (LABELS[entry["id"]], parity["conditional_fu_kane_nu"])
        else:
            status = wcc_status(folder)
            title = "%s: Z₂ = %s (WCC, conditional)" % (
                LABELS[entry["id"]], " / ".join(str(m["z2"]) for m in status["manifolds"].values())) if status else LABELS[entry["id"]]
        ax.set_title(title)
        ax.set_ylim(*window)
        ax.set_xlim(x[0], x[-1])
        axg.set_ylim(3e-4, 30)
        ax.tick_params(direction="in", which="both", top=True, right=True, labelbottom=False)
        axg.tick_params(direction="in", which="both", top=True, right=True)
        axg.set_xticks(ticks)
        axg.set_xticklabels(labels)
        if column == 0:
            ax.set_ylabel(r"$E-E_\mathrm{F}$ (eV)")
            axg.set_ylabel(r"$E_{N+1}-E_N$ (eV)")
        ax.legend(loc="lower right", frameon=False, fontsize=8)
    return _finish(fig, save)


def _annotate_parities(ax, x, energies, kpoints, parity, n, right_end):
    """Kramers-pair parities of the lowest-N subspace at the TRIM on the path."""
    for index in range(len(kpoints)):
        name = trim_name(kpoints[index])
        if name is None or (index and abs(x[index] - x[index - 1]) < 1e-12):
            continue
        pairs = [p for block in parity["trims"][name]["blocks"]
                 for p in block["kramers_pair_parities_unordered_within_block"]]
        groups = []
        for m, value in enumerate(pairs):
            energy = energies[index, 2 * m:2 * m + 2].mean()
            if groups and energy - groups[-1][0] < 0.45:
                groups[-1][1].append(value)
            else:
                groups.append([energy, [value]])
        left = x[index] > 0.97 * right_end
        for energy, values in groups:
            ax.text(x[index] + (-0.03 if left else 0.03) * right_end, energy, ",".join(sign(v) for v in values),
                    ha="right" if left else "left", va="center", fontsize=9, color=INK, zorder=5,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))


def plot_direct_gap_maps(topology=None, ids=None, save=None):
    """E_(N+1) − E_N over the sampled reciprocal cell (SCF grid) with the refined minima marked."""
    _, entries = campaign(topology)
    panels = [(e, n) for e in entries if ids is None or e["id"] in ids for n in subspace_sizes(e)]
    _style()
    columns = 4
    rows = math.ceil(len(panels) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(3.6 * columns, 3.5 * rows), squeeze=False)
    image = None
    for ax, (entry, n) in zip(axes.flat, panels):
        folder = entry["folder"]
        grid, mesh = scf_direct_gap_grid(folder, n)
        shift = [m // 2 for m in mesh]
        values = np.roll(grid, shift, axis=(0, 1))
        edges = [(np.arange(m + 1) - shift[i] - 0.5) / m for i, m in enumerate(mesh)]
        f1, f2 = np.meshgrid(edges[0], edges[1], indexing="ij")
        b = reciprocal_2d(folder)
        X, Y = f1 * b[0, 0] + f2 * b[1, 0], f1 * b[0, 1] + f2 * b[1, 1]
        image = ax.pcolormesh(X, Y, np.clip(values, GAP_NORM.vmin, None), cmap=GAP_CMAP, norm=GAP_NORM,
                              shading="flat", rasterized=True)
        for name, k in TRIMS.items():
            p = np.array(k) @ b
            ax.plot(*p, marker="o", ms=5, mfc="white", mec=INK, mew=0.8, zorder=4)
            if name == "Gamma":
                ax.annotate("Γ", p, xytext=(4, 4), textcoords="offset points", fontsize=9, color=INK)
        zoom = read_json(folder / "gap_refine_3_summary.json")
        minima = ([(b_["levels"][-1]["min_k"], b_["levels"][-1]["min_direct_gap_ev"]) for b_ in zoom["basins"]]
                  if zoom else [])
        gaps, _ = gap_extrema(folder)
        best = [g for g in gaps if g["N"] == n][0]
        if not zoom:
            minima = [(best["direct_gap_k"][:2], best["direct_gap_ev"])]
        for k, value in minima:
            points = _images(folder, k) @ b
            ax.plot(points[:, 0], points[:, 1], ls="none", marker="x", ms=8, mew=1.6, color="#F03C64", zorder=5)
        ax.set_title("%s, N = %d: min %s" % (LABELS[entry["id"]], n, gap_text(best["direct_gap_ev"])))
        ax.set_aspect("equal")
        ax.tick_params(direction="in", which="both", top=True, right=True)
        ax.set_xlabel(r"$k_x$ (Å$^{-1}$)")
        ax.set_ylabel(r"$k_y$ (Å$^{-1}$)")
    for ax in list(axes.flat)[len(panels):]:
        ax.axis("off")
    bar = fig.colorbar(image, ax=axes, shrink=0.8, pad=0.02)
    bar.set_label(r"$E_{N+1}-E_N$ (eV), 105×105 SCF grid")
    return _finish(fig, save)


def plot_wcc(topology=None, save=None):
    """Wilson-loop Wannier charge centres of the 1H-β subspaces (surface [t, s/2, 0])."""
    _, entries = campaign(topology)
    panels = [(e, n, m) for e in entries for status in [wcc_status(e["folder"])] if status
              for n, m in status["manifolds"].items()]
    _style()
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.8), squeeze=False)
    for ax, (entry, n, m), color in zip(axes.flat, panels, SUBSPACE_COLORS):
        for s, centres in zip(m["line_positions"], m["wcc"]):
            ax.plot([s / 2] * len(centres), centres, ls="none", marker="o", ms=4.5, color=color)
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$k_2$ (units of $\mathbf{b}_2$)")
        ax.set_ylabel(r"WCC along $\mathbf{a}_1$")
        ax.set_title("%s lowest %s bands: Z₂ = %d" % (LABELS[entry["id"]], n, m["z2"]))
        ax.tick_params(direction="in", which="both", top=True, right=True)
    return _finish(fig, save)


def plot_gap_zoom(topology=None, save=None):
    """Minimum direct gap per adaptive-zoom level versus sampling radius, with the local slope bound."""
    _, entries = campaign(topology)
    _style()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    markers = ["o", "s", "^", "D"]
    for entry in entries:
        zoom = read_json(entry["folder"] / "gap_refine_3_summary.json")
        if zoom is None:
            continue
        color = STRUCTURE_COLORS.get(entry["id"], INK)
        for index, basin in enumerate(zoom["basins"]):
            radius = np.array([l["sampling_radius_inverse_angstrom"] for l in basin["levels"]])
            gap = np.array([1000 * l["min_direct_gap_ev"] for l in basin["levels"]])
            bound = np.array([1000 * l["slope_bound_ev"] for l in basin["levels"]])
            ax.plot(radius, gap, color=color, lw=1.6, marker=markers[index], ms=6, zorder=3)
            positive = bound > 0
            ax.vlines(radius[positive], bound[positive], gap[positive], color=color, lw=1.0, zorder=2)
            ax.plot(radius[positive], bound[positive], ls="none", marker="_", ms=10, color=color, zorder=2)
            name = LABELS[entry["id"]] + (" basin %d" % (index + 1) if len(zoom["basins"]) > 1 else "")
            ax.annotate("%s: %.3f meV" % (name, gap[-1]), (radius[-1], gap[-1]), xytext=(6, 0),
                        textcoords="offset points", va="center", fontsize=8, color=INK)
    ax.set_xscale("log")
    radii = [l["sampling_radius_inverse_angstrom"] for e in entries
             for z in [read_json(e["folder"] / "gap_refine_3_summary.json")] if z
             for b_ in z["basins"] for l in b_["levels"]]
    ax.set_xlim(max(radii) * 1.6, min(radii) / 60)  # inverted: finer sampling to the right, room for labels
    ax.set_ylim(0, None)
    ax.set_xlabel(r"Sampling radius (Å$^{-1}$), finer →")
    ax.set_ylabel(r"Min. $E_{N+1}-E_N$ in patch (meV)")
    ax.set_title("Adaptive zoom: minimum and local slope bound")
    ax.tick_params(direction="in", which="both", top=True, right=True)
    return _finish(fig, save)
