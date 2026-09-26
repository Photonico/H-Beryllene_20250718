#### Band topology
# This module provides tables and plots for the SOC band-topology screening in 9.0_topology.
# It reads the campaign reports (parity, gap refinement, Wilson-loop WCC, time-reversal evidence, spin screen)
# and the validated VASP outputs; it never writes into the calculation directories.
# Every Z2 value is a conditional screening result for a fixed lowest-N spinor-band subspace, not a certified invariant.
# pylint: disable = C0103, C0114, C0116, C0301, C0302, C0321, R0913, R0914, R0915, W0612

# Necessary packages invoking
import os
import glob
import json
import math
import h5py
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D

from vmatplot.output_settings import color_sampling, canvas_setting

import matplotlib as mpl

mpl.rcParams["lines.solid_capstyle"] = "round"
mpl.rcParams["lines.dash_capstyle"]  = "round"
mpl.rcParams["lines.solid_joinstyle"] = "round"
mpl.rcParams["lines.dash_joinstyle"]  = "round"

# Four 2D TRIM in the reciprocal basis of the DFT cell
trim_points = {"Gamma": (0.0, 0.0), "X": (0.5, 0.0), "Y": (0.0, 0.5), "M": (0.5, 0.5)}
structure_labels = {"alpha": "α", "beta": "β", "st": "ST", "alpha_2h": "2H-α", "beta_1h": "1H-β", "beta_2h": "2H-β"}

# Blue[1], Orange[1] and Cyan[1] pass the colour-blind all-pairs check (worst ΔE 15.4, tritan 6.1);
# Orange is below 3:1 on white, so it is always direct-labelled.
# One-hue sequential blue ramp for the direct-gap maps, dark = small gap.
gap_colormap = LinearSegmentedColormap.from_list(
    "direct_gap", ["#0d366b", "#104281", "#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5",
                   "#5598e7", "#6da7ec", "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"])

## Data extraction

def extract_json(file_path):
    if file_path is None or not os.path.isfile(file_path):
        return None
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)

def extract_latest(directory, pattern):
    matches = sorted(glob.glob(os.path.join(directory, pattern)))
    return matches[-1] if matches else None

def extract_topology_manifest(topology_dir):
    # Campaign manifest entries, each with its structure directory
    manifest = extract_json(os.path.join(topology_dir, "manifest.json"))
    entries = manifest["structures"]
    for entry in entries:
        entry["path"] = os.path.join(topology_dir, entry["directory"])
    return entries

def extract_topology_entry(directory):
    # Manifest entry of one structure directory, e.g. "9.0_topology/a-Beryllene"
    directory = os.path.normpath(directory)
    for entry in extract_topology_manifest(os.path.dirname(directory)):
        if entry["directory"] == os.path.basename(directory):
            return entry
    raise ValueError(f"{directory} is not in the campaign manifest")

def extract_subspaces(entry):
    # Even filling: the lowest NELECT spinor bands; odd filling: NELECT-1 and NELECT+1
    nelect = entry["expected_nelect"]
    return [nelect] if nelect % 2 == 0 else [nelect - 1, nelect + 1]

def extract_fermi_topology(directory):
    # SCF Fermi energy from DOSCAR (line 6, fourth column)
    with open(os.path.join(directory, "scf", "DOSCAR"), "r", encoding="utf-8") as f:
        lines = f.readlines()
    return float(lines[5].split()[3])

def extract_gap_extrema(directory):
    # Final sampled extrema over scf, trim, bands and every refinement; True if the adaptive zoom ran
    zoom = extract_json(os.path.join(directory, "gap_refine_3_summary.json"))
    if zoom:
        return zoom["sampled_extrema_after"], True
    return extract_json(os.path.join(directory, "gap_refinement_summary.json"))["final_sampled_extrema"], False

def extract_parity(directory):
    analysis = extract_latest(os.path.join(directory, "trim"), "parity-analysis-*")
    return extract_json(os.path.join(analysis, "parity_summary.json")) if analysis else None

def extract_wcc(directory):
    return extract_json(extract_latest(os.path.join(directory, "wcc_direct_v3"), "status_*.json"))

def extract_reciprocal_2d(directory):
    # In-plane reciprocal vectors (1/Å, rows b1 and b2) of the validated cell
    lattice = np.array(extract_json(os.path.join(directory, "trim", "topology_validation.json"))["structure"]["lattice"])
    return (2 * np.pi * np.linalg.inv(lattice).T)[:2, :2]

def identify_trim(kpoint):
    for name, target in trim_points.items():
        if all(abs((kpoint[i] - target[i] + 0.5) % 1 - 0.5) < 1e-6 for i in (0, 1)):
            return name
    return None

def extract_band_path(directory):
    # Fixed-density SOC bands on the line-mode path: kpath (1/Å), eigenvalues (eV), fractional k, ticks
    bands_dir = os.path.join(directory, "bands")
    with h5py.File(os.path.join(bands_dir, "vaspout.h5"), "r") as f:
        eigenvalues = f["/results/electron_eigenvalues/eigenvalues"][()][0]
        kpoints = f["/results/electron_eigenvalues/kpoint_coords"][()][:, :2]
    with open(os.path.join(bands_dir, "KPOINTS"), "r", encoding="utf-8") as f:
        lines = f.readlines()
    per_segment = int(lines[1].split()[0])
    ends = [line.split()[3] for line in lines[4:] if len(line.split()) >= 4]
    cartesian = kpoints @ extract_reciprocal_2d(directory)
    steps = np.linalg.norm(np.diff(cartesian, axis=0), axis=1)
    steps[per_segment - 1::per_segment] = 0.0  # segment ends are repeated (or jump)
    kpath = np.concatenate([[0.0], np.cumsum(steps)])
    positions, labels = [], []
    for segment in range(len(kpoints) // per_segment):
        for index, label in ((segment * per_segment, ends[2 * segment]), ((segment + 1) * per_segment - 1, ends[2 * segment + 1])):
            if positions and abs(kpath[index] - positions[-1]) < 1e-9:
                if label != labels[-1]:
                    labels[-1] += "|" + label
                continue
            positions.append(kpath[index])
            labels.append(label)
    return kpath, eigenvalues, kpoints, positions, labels

def extract_point_group_2d(directory, symprec=1e-5):
    # In-plane action of every space-group operation that does not mix z with the plane
    import spglib
    with open(os.path.join(directory, "scf", "POSCAR"), "r", encoding="utf-8") as f:
        lines = f.readlines()
    scale = float(lines[1].split()[0])
    cell = np.array([[float(v) for v in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    counts = [int(v) for v in lines[6].split()]
    start = 8 if lines[7].strip()[0].upper() in "DC" else 9
    positions = np.array([[float(v) for v in lines[start + i].split()[:3]] for i in range(sum(counts))])
    numbers = [i for i, count in enumerate(counts) for _ in range(count)]
    rotations = spglib.get_symmetry((cell, positions, numbers), symprec=symprec)["rotations"]
    blocks = {tuple(r[:2, :2].ravel()) for r in rotations if not r[2, :2].any() and not r[:2, 2].any()}
    return [np.array(block).reshape(2, 2) for block in sorted(blocks)]

def extract_direct_gap_grid(directory, subspace):
    # E_(N+1) - E_N on the full Gamma-centred SCF grid, unfolded with the point group and k -> -k
    with open(os.path.join(directory, "scf", "KPOINTS"), "r", encoding="utf-8") as f:
        mesh = [int(v) for v in f.readlines()[3].split()[:2]]
    with h5py.File(os.path.join(directory, "scf", "vaspout.h5"), "r") as f:
        eigenvalues = f["/results/electron_eigenvalues/eigenvalues"][()][0]
        kpoints = f["/results/electron_eigenvalues/kpoint_coords"][()][:, :2]
    gaps = eigenvalues[:, subspace] - eigenvalues[:, subspace - 1]
    grid = np.full(mesh, np.nan)
    disagreement = 0.0
    for rotation in extract_point_group_2d(directory):
        for sign in (1, -1):
            images = sign * kpoints @ rotation * mesh
            if np.abs(images - np.rint(images)).max() > 1e-4:
                raise ValueError(f"Symmetry image off the SCF grid in {directory}")
            index = np.rint(images).astype(int) % mesh
            known = ~np.isnan(grid[index[:, 0], index[:, 1]])
            if known.any():
                disagreement = max(disagreement, float(np.abs(grid[index[known, 0], index[known, 1]] - gaps[known]).max()))
            grid[index[:, 0], index[:, 1]] = gaps
    if np.isnan(grid).any() or disagreement > 1e-4:
        raise ValueError(f"Unfolding incomplete or inconsistent ({disagreement:.2e} eV) in {directory}")
    return grid, mesh

def extract_symmetry_images(directory, kpoint):
    images = {tuple(np.round(((sign * np.asarray(kpoint[:2]) @ rotation) + 0.5) % 1 - 0.5, 9))
              for rotation in extract_point_group_2d(directory) for sign in (1, -1)}
    return np.array(sorted(images))

def format_gap(gap_ev):
    return f"{gap_ev:.3g} eV" if gap_ev >= 0.1 else f"{1000*gap_ev:.3g} meV"

def format_sign(value):
    return "+" if value > 0 else "−"

## Summary tables (Markdown in Jupyter)

def show_table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(value) for value in row) + " |" for row in rows]
    try:
        from IPython.display import Markdown, display
        display(Markdown("\n".join(lines)))
    except ImportError:
        print("\n".join(lines))
    return rows

def summarize_topology(topology_dir="9.0_topology"):
    rows = []
    for entry in extract_topology_manifest(topology_dir):
        directory = entry["path"]
        parity, wcc = extract_parity(directory), extract_wcc(directory)
        if parity:
            invariant = f"ν = {parity['conditional_fu_kane_nu']} (parity)"
        elif wcc:
            invariant = ", ".join(f"N={n}: Z₂ = {m['z2']} (WCC)" for n, m in wcc["manifolds"].items())
        else:
            invariant = "not validated"
        fermi = extract_fermi_topology(directory)
        gaps, _ = extract_gap_extrema(directory)
        insulating = all(gap["valence_max_ev"] < fermi < gap["conduction_min_ev"] for gap in gaps)
        tr = extract_json(os.path.join(directory, "tr_evidence.json"))
        spin = extract_json(os.path.join(directory, "spin_screen", "spin_screen_summary.json"))
        rows.append([structure_labels[entry["id"]], entry["directory"], entry["expected_nelect"], invariant,
                     ", ".join(format_gap(gap["direct_gap_ev"]) for gap in gaps),
                     "insulator" if insulating else "metal",
                     f"{tr['magnetisation_density']['max_abs_m_muB_per_A3']:.1e}" if tr else "pending",
                     ("collapsed" if spin["all_seeds_collapsed"] else "MOMENT KEPT") if spin else "pending"])
    header = ["Structure", "Directory", "NELECT", "Conditional invariant", "Min. direct gap E(N+1)−E(N)",
              "Neutral E_F", "max \\|m(r)\\| (μB/Å³)", "Spin seeds"]
    return show_table(header, rows)

def summarize_parity(topology_dir="9.0_topology", validation=False):
    # Kramers-pair parities at the four TRIM, or (validation=True) the float64 wavefunction checks
    rows = []
    for entry in extract_topology_manifest(topology_dir):
        if not entry["inversion_expected"]:
            continue
        report = extract_parity(entry["path"])
        label = structure_labels[entry["id"]]
        if report is None or "wavefunction_validation" not in report:
            rows.append([label, "not validated"])
            continue
        checks = report["wavefunction_validation"]
        if validation:
            rows.append([label, f"{max(c['lowdin_inversion_singular_value_max_error'] for c in checks.values()):.1e}",
                         f"{max(c['plain_metric_opposite_parity_max_overlap'] for c in checks.values()):.1e}",
                         f"{max(c['inversion_residual_max'] for c in checks.values()):.1e}",
                         f"{max(c['plain_metric_same_parity_max_overlap'] for c in checks.values()):.1e}",
                         ", ".join(f"{v:.3g}" for v in report["irrep_plain_metric_orthogonality_messages"]) or "none"])
            continue
        cells = []
        for name in trim_points:
            pairs = " ".join(format_sign(p) for block in report["trims"][name]["blocks"]
                             for p in block["kramers_pair_parities_unordered_within_block"])
            above = checks[name]["next_block_above_N_parity"]
            cells.append(f"{pairs} \\| {format_sign(above) if above else '?'} ({format_gap(checks[name]['next_block_above_N_gap_ev'])})")
        rows.append([label, report["occupied_spinor_bands"], *cells,
                     " ".join(format_sign(report["trims"][n]["delta"]) for n in trim_points), report["conditional_fu_kane_nu"]])
    if validation:
        header = ["Structure", "Löwdin unitarity error", "Opposite-parity overlap", "Inversion residual",
                  "Same-parity overlap (PAW metric)", "IrRep plain-metric message"]
    else:
        header = ["Structure", "N", "Γ", "X", "Y", "M", "δ (Γ X Y M)", "ν"]
    return show_table(header, rows)

def summarize_gaps(topology_dir="9.0_topology"):
    rows = []
    for entry in extract_topology_manifest(topology_dir):
        fermi = extract_fermi_topology(entry["path"])
        gaps, _ = extract_gap_extrema(entry["path"])
        for gap in gaps:
            rows.append([structure_labels[entry["id"]], gap["N"], format_gap(gap["direct_gap_ev"]),
                         "(%.5f, %.5f)" % tuple(gap["direct_gap_k"][:2]), gap["direct_gap_stage"],
                         f"{gap['valence_max_ev'] - fermi:+.3f}", f"{gap['conduction_min_ev'] - fermi:+.3f}",
                         f"{gap['indirect_gap_ev']:+.3f}"])
    header = ["Structure", "N", "Min. direct gap", "k (fractional)", "Stage", "max E_N − E_F (eV)",
              "min E_N+1 − E_F (eV)", "Indirect gap (eV)"]
    return show_table(header, rows)

def summarize_gap_zoom(topology_dir="9.0_topology"):
    rows = []
    for entry in extract_topology_manifest(topology_dir):
        zoom = extract_json(os.path.join(entry["path"], "gap_refine_3_summary.json"))
        if zoom is None:
            continue
        for index, basin in enumerate(zoom["basins"]):
            final = basin["levels"][-1]
            rows.append([structure_labels[entry["id"]], "%d: (%.5f, %.5f)" % (index + 1, *basin["start"]), len(basin["levels"]),
                         " → ".join(f"{1000*level['min_direct_gap_ev']:.3f}" for level in basin["levels"]),
                         "(%.7f, %.7f)" % tuple(final["min_k"]), f"{final['spacing']:.1e}",
                         "yes" if final["min_interior"] else "no", f"{1000*final['slope_bound_ev']:+.3f}"])
    header = ["Structure", "Basin (start k)", "Levels", "Min. gap per level (meV)", "Final k", "Final spacing",
              "Interior", "Slope bound (meV)"]
    return show_table(header, rows)

def summarize_wcc(topology_dir="9.0_topology"):
    rows = []
    for entry in extract_topology_manifest(topology_dir):
        status = extract_wcc(entry["path"])
        if status is None:
            continue
        for bands, manifold in status["manifolds"].items():
            gap = min(stage["minimum_direct_gap_ev"] for stage in manifold["sampled_gaps"].values())
            rows.append([structure_labels[entry["id"]], bands, manifold["status"], manifold["z2"], format_gap(gap),
                         len(manifold["line_positions"]), manifold["topology_certified"]])
    header = ["Structure", "Bands", "Status", "Conditional Z₂", "Min. sampled direct gap", "Wilson loops", "Certified"]
    return show_table(header, rows)

def summarize_time_reversal(topology_dir="9.0_topology", spin_screen=False):
    # Converged-state evidence (tr_evidence.json), or (spin_screen=True) the magnetic-seed SCFs
    rows = []
    for entry in extract_topology_manifest(topology_dir):
        label = structure_labels[entry["id"]]
        if spin_screen:
            spin = extract_json(os.path.join(entry["path"], "spin_screen", "spin_screen_summary.json"))
            if spin is None:
                rows.append([label, "pending", "", "", "", ""])
                continue
            for seed, result in spin["seeds"].items():
                sites = result["site_moments_muB"] or [float("nan")]
                rows.append([label, seed, " ".join(f"{v:g}" for v in result["initial_magmom_muB"]),
                             f"{result['mag_muB']:.4f}", f"{max(abs(v) for v in sites):.4f}",
                             "yes" if result["moment_collapsed"] else "NO"])
            continue
        tr = extract_json(os.path.join(entry["path"], "tr_evidence.json"))
        if tr is None:
            continue
        density, minus_k = tr["magnetisation_density"], tr.get("scf_minus_k")
        rows.append([label, f"{density['max_abs_m_muB_per_A3']:.1e}", f"{density['integral_abs_m_muB']:.1e}",
                     f"{density['paw_occupancy_max_abs_re_m']:.1e} / {density['paw_occupancy_max_abs_im_m']:.1e}",
                     f"{tr['trim_kramers_max_splitting_ev']:.1e}",
                     f"{minus_k['max_abs_dE_ev']:.1e} ({minus_k['pairs']} pairs)" if minus_k else "—"])
    if spin_screen:
        header = ["Structure", "Seed", "Initial moments (μB)", "Final total (μB)", "Final max \\|site\\| (μB)", "Collapsed"]
    else:
        header = ["Structure", "max \\|m(r)\\| (μB/Å³)", "∫\\|m\\| (μB)", "One-centre m: max \\|Re\\| / max \\|Im\\|",
                  "TRIM Kramers splitting (eV)", "max \\|E(k)−E(−k)\\| (eV)"]
    return show_table(header, rows)

## Plotting

def create_matters_topology(matters_list):
    # Ensure input is a list of lists: [label, structure directory, colour family, second colour family]
    if isinstance(matters_list, list) and matters_list and not any(isinstance(i, list) for i in matters_list):
        source_data = matters_list[:]
        matters_list.clear()
        matters_list.append(source_data)
    matters = []
    for current_matter in matters_list:
        label, directory, *optional = current_matter
        colors = [optional[0] if len(optional) > 0 and optional[0] else "Blue",
                  optional[1] if len(optional) > 1 and optional[1] else "Orange"]
        entry = extract_topology_entry(directory)
        matters.append([label, directory, entry, extract_subspaces(entry), [color_sampling(c)[1] for c in colors]])
    return matters

def plot_topology_bands(suptitle, matters_list=None, eigen_range=None, legend_loc=True):
    # Help information
    help_info = """
    Usage: plot_topology_bands
        arg[0]: suptitle;
        arg[1]: matters list, [[label, structure directory, colour family, second colour family], ...];
        arg[2]: energy window relative to E_F, e.g. (-11, 5);
        arg[3]: figure legend (True/False);
    Upper panels: SOC bands with the lowest-N subspace coloured and Kramers-pair parities at the TRIM;
    lower panels: direct gap E_(N+1) - E_N along the path (log scale) and the minimum over the whole BZ.
    """
    if suptitle in ["help", "Help"]:
        print(help_info)
        return

    matters = create_matters_topology(matters_list)
    columns = min(3, len(matters))
    rows = math.ceil(len(matters) / columns)

    # Figure settings
    fig_setting = canvas_setting(6 * columns, 6.4 * rows)
    params = fig_setting[2]; plt.rcParams.update(params)
    fig = plt.figure(figsize=fig_setting[0], dpi=fig_setting[1])
    outer = gridspec.GridSpec(rows, columns, figure=fig, left=0.07, right=0.98, top=0.93, bottom=0.08, hspace=0.3, wspace=0.26)

    # Colors calling
    fermi_color = color_sampling("Violet")
    annotate_color = color_sampling("Grey")

    # Title
    fig.suptitle(f"{suptitle}", fontsize=fig_setting[3][0], y=0.985)

    # Data calling and plotting
    energy_window = eigen_range if eigen_range is not None else (-11, 5)
    for index, matter in enumerate(matters):
        label, directory, entry, subspaces, colors = matter
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[index // columns, index % columns],
                                                 height_ratios=[3.2, 1.3], hspace=0.06)
        ax_bands = fig.add_subplot(inner[0])
        ax_gap = fig.add_subplot(inner[1], sharex=ax_bands)
        kpath, eigenvalues, kpoints, positions, labels = extract_band_path(directory)
        eigenvalues = eigenvalues - extract_fermi_topology(directory)

        # Subspace bands, higher bands and Fermi energy
        lower = 0
        for color, subspace in zip(colors, subspaces):
            for band in range(lower, subspace):
                ax_bands.plot(kpath, eigenvalues[:, band], c=color, lw=1.8, zorder=4)
            lower = subspace
        for band in range(lower, eigenvalues.shape[1]):
            ax_bands.plot(kpath, eigenvalues[:, band], c=annotate_color[2], lw=1.2, zorder=3)
        ax_bands.axhline(y=0, color=fermi_color[0], alpha=0.8, linestyle="--", zorder=2)

        # Kramers-pair parities at the TRIM on the path
        parity = extract_parity(directory)
        if parity:
            annotate_parities(ax_bands, kpath, eigenvalues, kpoints, parity, annotate_color[0])

        # Direct gap along the path and the BZ minimum
        gaps, _ = extract_gap_extrema(directory)
        for order, (color, subspace, gap) in enumerate(zip(colors, subspaces, gaps)):
            ax_gap.semilogy(kpath, np.maximum(eigenvalues[:, subspace] - eigenvalues[:, subspace - 1], 1e-6), c=color, lw=1.6, zorder=4)
            ax_gap.axhline(y=gap["direct_gap_ev"], color=color, lw=1.0, linestyle=":", zorder=3)
            ax_gap.text(kpath[-1] * (0.98 if order == 0 else 0.02), gap["direct_gap_ev"] * 1.5, f"BZ min {format_gap(gap['direct_gap_ev'])}",
                        ha="right" if order == 0 else "left", va="bottom", fontsize=11, color=annotate_color[0])

        # High symmetry path
        for k_loc in positions[1:-1]:
            ax_bands.axvline(x=k_loc, color=annotate_color[1], linestyle="--", alpha=0.8, zorder=1)
            ax_gap.axvline(x=k_loc, color=annotate_color[1], linestyle="--", alpha=0.8, zorder=1)

        # Subtitle, axes and ranges
        if parity:
            subtitle = f"{label}: ν = {parity['conditional_fu_kane_nu']}"
        else:
            status = extract_wcc(directory)
            subtitle = f"{label}: Z₂ = " + "/".join(str(m["z2"]) for m in status["manifolds"].values()) + " (WCC)" if status else label
        ax_bands.set_title(subtitle, fontsize=fig_setting[3][1])
        ax_bands.set_ylim(energy_window[0], energy_window[1])
        ax_bands.set_xlim(kpath[0], kpath[-1])
        ax_gap.set_ylim(3e-4, 30)
        ax_gap.set_xticks(positions)
        ax_gap.set_xticklabels(labels)
        ax_bands.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True, labelbottom=False)
        ax_gap.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)
        if index % columns == 0:
            ax_bands.set_ylabel(r"$E-E_\mathrm{F}$ (eV)")
            ax_gap.set_ylabel(r"$E_{N+1}-E_N$ (eV)")

    # Legend
    if legend_loc:
        handles = [Line2D([], [], c=color_sampling("Blue")[1], lw=1.8, label="lowest-N subspace (bands 1…N)"),
                   Line2D([], [], c=color_sampling("Orange")[1], lw=1.8, label="1H-β bands 5–6 (lowest-6 subspace)"),
                   Line2D([], [], c=annotate_color[2], lw=1.2, label="bands above N"),
                   Line2D([], [], c=fermi_color[0], lw=1.0, linestyle="--", label="Fermi energy")]
        fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))

def annotate_parities(ax, kpath, eigenvalues, kpoints, parity, text_color):
    # Parities of the Kramers pairs of the lowest-N subspace; nearly degenerate pairs share one label
    for index in range(len(kpoints)):
        name = identify_trim(kpoints[index])
        if name is None or (index > 0 and abs(kpath[index] - kpath[index - 1]) < 1e-12):
            continue
        pairs = [p for block in parity["trims"][name]["blocks"] for p in block["kramers_pair_parities_unordered_within_block"]]
        groups = []
        for pair, value in enumerate(pairs):
            energy = eigenvalues[index, 2 * pair:2 * pair + 2].mean()
            if groups and energy - groups[-1][0] < 0.45:
                groups[-1][1].append(value)
            else:
                groups.append([energy, [value]])
        right_end = kpath[index] > 0.97 * kpath[-1]
        for energy, values in groups:
            ax.text(kpath[index] + (-0.025 if right_end else 0.025) * kpath[-1], energy, ",".join(format_sign(v) for v in values),
                    ha="right" if right_end else "left", va="center", fontsize=12, color=text_color, zorder=5,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))

def plot_direct_gap_maps(suptitle, matters_list=None, gap_range=None):
    # Help information
    help_info = """
    Usage: plot_direct_gap_maps
        arg[0]: suptitle;
        arg[1]: matters list, [[label, structure directory], ...];
        arg[2]: colour range of E_(N+1) - E_N in eV, default (1e-2, 10), log scale;
    Direct gap over the sampled reciprocal cell (105x105 SCF grid unfolded with the point group and k -> -k);
    crosses mark the refined minima and their symmetry images, circles the TRIM.
    """
    if suptitle in ["help", "Help"]:
        print(help_info)
        return

    matters = create_matters_topology(matters_list)
    panels = [(matter, subspace) for matter in matters for subspace in matter[3]]
    columns = min(4, len(panels))
    rows = math.ceil(len(panels) / columns)
    gap_limits = gap_range if gap_range is not None else (1e-2, 10)
    norm = LogNorm(vmin=gap_limits[0], vmax=gap_limits[1])

    # Figure settings
    fig_setting = canvas_setting(5.2 * columns + 1.2, 5.4 * rows)
    params = fig_setting[2]; plt.rcParams.update(params)
    fig, axes = plt.subplots(rows, columns, figsize=fig_setting[0], dpi=fig_setting[1], squeeze=False)

    # Colors calling
    annotate_color = color_sampling("Grey")
    minimum_color = color_sampling("Red")[1]

    # Title
    fig.suptitle(f"{suptitle}", fontsize=fig_setting[3][0], y=1.00)

    # Data calling and plotting
    image = None
    for ax, (matter, subspace) in zip(axes.flat, panels):
        label, directory = matter[0], matter[1]
        grid, mesh = extract_direct_gap_grid(directory, subspace)
        shift = [m // 2 for m in mesh]
        values = np.roll(grid, shift, axis=(0, 1))
        edges = [(np.arange(m + 1) - shift[i] - 0.5) / m for i, m in enumerate(mesh)]
        frac_1, frac_2 = np.meshgrid(edges[0], edges[1], indexing="ij")
        reciprocal = extract_reciprocal_2d(directory)
        kx = frac_1 * reciprocal[0, 0] + frac_2 * reciprocal[1, 0]
        ky = frac_1 * reciprocal[0, 1] + frac_2 * reciprocal[1, 1]
        image = ax.pcolormesh(kx, ky, np.clip(values, gap_limits[0] * 0.5, None), cmap=gap_colormap, norm=norm,
                              shading="flat", rasterized=True)
        for name, trim in trim_points.items():
            point = np.array(trim) @ reciprocal
            ax.plot(*point, marker="o", ms=6, mfc="white", mec=annotate_color[0], mew=0.9, zorder=4)
            if name == "Gamma":
                ax.annotate("Γ", point, xytext=(5, 5), textcoords="offset points", fontsize=13, color=annotate_color[0])
        zoom = extract_json(os.path.join(directory, "gap_refine_3_summary.json"))
        gaps, _ = extract_gap_extrema(directory)
        best = [gap for gap in gaps if gap["N"] == subspace][0]
        minima = [basin["levels"][-1]["min_k"] for basin in zoom["basins"]] if zoom else [best["direct_gap_k"][:2]]
        for kpoint in minima:
            points = extract_symmetry_images(directory, kpoint) @ reciprocal
            ax.plot(points[:, 0], points[:, 1], linestyle="none", marker="x", ms=10, mew=2.0, color=minimum_color, zorder=5)
        ax.set_title(f"{label}, N = {subspace}\nmin {format_gap(best['direct_gap_ev'])}", fontsize=fig_setting[3][1] - 4)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$k_x$ (Å$^{-1}$)")
        ax.set_ylabel(r"$k_y$ (Å$^{-1}$)")
        ax.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)
    for ax in list(axes.flat)[len(panels):]:
        ax.axis("off")

    # Colorbar
    fig.subplots_adjust(left=0.06, right=0.88, top=0.88, bottom=0.07, hspace=0.38, wspace=0.42)
    colorbar = fig.colorbar(image, ax=axes, shrink=0.8, pad=0.03, extend="min")
    colorbar.set_label(r"$E_{N+1}-E_N$ (eV)")

def plot_wcc(suptitle, matters_list=None):
    # Help information
    help_info = """
    Usage: plot_wcc
        arg[0]: suptitle;
        arg[1]: matters list, [[label, structure directory, colour family, second colour family], ...];
    Wannier charge centres of each Wilson-loop subspace on the surface [t, s/2, 0] (Z2Pack).
    """
    if suptitle in ["help", "Help"]:
        print(help_info)
        return

    matters = create_matters_topology(matters_list)
    panels = [(matter, color, bands, manifold) for matter in matters for status in [extract_wcc(matter[1])] if status
              for color, (bands, manifold) in zip(matter[4], status["manifolds"].items())]

    # Figure settings
    fig_setting = canvas_setting(6 * len(panels), 5.2)
    params = fig_setting[2]; plt.rcParams.update(params)
    fig, axes = plt.subplots(1, len(panels), figsize=fig_setting[0], dpi=fig_setting[1], squeeze=False)

    # Title
    fig.suptitle(f"{suptitle}", fontsize=fig_setting[3][0], y=1.00)

    # Data calling and plotting
    for ax, (matter, color, bands, manifold) in zip(axes.flat, panels):
        for line_position, centres in zip(manifold["line_positions"], manifold["wcc"]):
            ax.plot([line_position / 2] * len(centres), centres, linestyle="none", marker="o", ms=6, c=color, zorder=4)
        ax.set_title(f"Lowest {bands} bands: Z₂ = {manifold['z2']}", fontsize=fig_setting[3][1])
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$k_2$ (units of $\mathbf{b}_2$)")
        ax.set_ylabel(r"WCC along $\mathbf{a}_1$")
        ax.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)

    plt.tight_layout()

def plot_gap_zoom(title, matters_list=None):
    # Help information
    help_info = """
    Usage: plot_gap_zoom
        arg[0]: title;
        arg[1]: matters list, [[label, structure directory, colour family], ...];
    Minimum direct gap of every adaptive-zoom level against the sampling radius; the vertical whisker reaches
    the heuristic local bound (minimum - largest neighbour slope x sampling radius) where that bound is positive.
    """
    if title in ["help", "Help"]:
        print(help_info)
        return

    matters = create_matters_topology(matters_list)

    # Figure settings
    fig_setting = canvas_setting(10, 6)
    params = fig_setting[2]; plt.rcParams.update(params)
    plt.figure(figsize=fig_setting[0], dpi=fig_setting[1])
    plt.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)

    # Data calling and plotting
    markers = ["o", "s", "^", "D"]
    radii = []
    for matter in matters:
        label, directory, color = matter[0], matter[1], matter[4][0]
        zoom = extract_json(os.path.join(directory, "gap_refine_3_summary.json"))
        if zoom is None:
            continue
        for index, basin in enumerate(zoom["basins"]):
            radius = np.array([level["sampling_radius_inverse_angstrom"] for level in basin["levels"]])
            gap = np.array([1000 * level["min_direct_gap_ev"] for level in basin["levels"]])
            bound = np.array([1000 * level["slope_bound_ev"] for level in basin["levels"]])
            radii.extend(radius)
            name = label + (f" basin {index + 1}" if len(zoom["basins"]) > 1 else "")
            plt.plot(radius, gap, c=color, lw=1.8, marker=markers[index], ms=8, label=f"{name}: {gap[-1]:.3f} meV", zorder=4)
            positive = bound > 0
            plt.vlines(radius[positive], bound[positive], gap[positive], color=color, lw=1.2, zorder=3)
            plt.plot(radius[positive], bound[positive], linestyle="none", marker="_", ms=12, c=color, zorder=3)

    # Title and axes
    plt.title(f"{title}")
    plt.xscale("log")
    plt.xlim(max(radii) * 1.6, min(radii) / 1.6)
    plt.ylim(0, None)
    plt.xlabel(r"Sampling radius (Å$^{-1}$), finer $\rightarrow$")
    plt.ylabel(r"Min. $E_{N+1}-E_N$ in patch (meV)")

    # Legend
    plt.legend(loc="upper right", frameon=False)
    plt.tight_layout()
