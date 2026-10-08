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
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, LogNorm, Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
from matplotlib.patches import Circle, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D

from vmatplot.output_settings import color_sampling, canvas_setting, figure_version, LINE_WIDTH, FERMI_STYLE
from vmatplot.topology_figures import topology_canvas, STYLES

import matplotlib as mpl

mpl.rcParams["lines.solid_capstyle"] = "round"
mpl.rcParams["lines.dash_capstyle"]  = "round"
mpl.rcParams["lines.solid_joinstyle"] = "round"
mpl.rcParams["lines.dash_joinstyle"]  = "round"

# Four 2D TRIM in the reciprocal basis of the DFT cell
trim_points = {"Gamma": (0.0, 0.0), "X": (0.5, 0.0), "Y": (0.0, 0.5), "M": (0.5, 0.5)}
structure_labels = {"alpha": "α", "beta": "β (P-1, superseded)", "beta_p3m1": "β", "st": "trilayer cubic", "alpha_2h": "2H-α", "beta_1h": "1H-β", "beta_2h": "2H-β"}
# Manifest entries left out of the tables: the P-1 β cell is an artefact of an 11×11-k relaxation (P-3m1 at 27×27 k)
superseded_ids = {"beta"}

# Blue[1], Orange[1] and Cyan[1] pass the colour-blind all-pairs check (worst ΔE 15.4, tritan 6.1);
# Orange is below 3:1 on white, so it is always direct-labelled.
# White to navy; every RGB component is a multiple of five.
heatmap_rgb = np.array([(255, 255, 255), (210, 235, 245), (155, 210, 230),
                        (85, 165, 205), (30, 120, 175), (10, 60, 125), (5, 30, 85)])


def make_heatmap_colormap(norm, ticks, reverse=False):
    """Small colour jumps at major ticks, with a gradient inside each interval."""
    positions = np.asarray(norm(ticks))
    boundaries = np.r_[0., positions[(positions > 0) & (positions < 1)], 1.]
    x = np.linspace(0, 1, 4096)
    interval = np.clip(np.searchsorted(boundaries, x, side="right") - 1, 0, len(boundaries) - 2)
    low, high = boundaries[interval], boundaries[interval + 1]
    # Compress each interval slightly, leaving a visible colour step at its boundary.
    start = np.where(low == 0, 0, low + .08 * (high - low))
    end = np.where(high == 1, 1, high - .08 * (high - low))
    progress = start + (x - low) / (high - low) * (end - start)
    anchors = heatmap_rgb[::-1] if reverse else heatmap_rgb
    rgb = np.column_stack([np.interp(progress, np.linspace(0, 1, len(anchors)), channel)
                           for channel in anchors.T])
    return ListedColormap(rgb / 255, name="tick_aligned_steps")


def style_heatmap_colorbar(colorbar, ticks):
    colorbar.ax.set_label("<colorbar>")
    colorbar.set_ticks(ticks)
    colorbar.ax.minorticks_off()
    colorbar.ax.tick_params(which="both", direction="in")
    colorbar.ax.grid(False)

## Data extraction

def extract_json(file_path):
    if file_path is None or not os.path.isfile(file_path):
        return None
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)

def extract_latest(directory, pattern):
    matches = sorted(glob.glob(os.path.join(directory, pattern)))
    return matches[-1] if matches else None

def extract_topology_manifest(topology_dir, superseded=False):
    # Campaign manifest entries, each with its structure directory (superseded entries only on request)
    manifest = extract_json(os.path.join(topology_dir, "manifest.json"))
    entries = [entry for entry in manifest["structures"] if superseded or entry["id"] not in superseded_ids]
    for entry in entries:
        entry["path"] = os.path.join(topology_dir, entry["directory"])
    return entries

def extract_topology_entry(directory):
    # Manifest entry of one structure directory, e.g. "9.0_topology/a-Beryllene"
    directory = os.path.normpath(directory)
    for entry in extract_topology_manifest(os.path.dirname(directory), superseded=True):
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
    # Distorted hexagonal cells (neither hexagonal nor square): the TRIM on the path are named as in the parity tables,
    # X = (1/2, 0), Y = (0, 1/2), M = (1/2, 1/2); the KPOINTS labels follow the hexagonal parent (M = (1/2, 0))
    reciprocal = extract_reciprocal_2d(directory)
    lengths = np.linalg.norm(reciprocal, axis=1)
    cosine = abs(reciprocal[0] @ reciprocal[1]) / lengths.prod()
    distorted = abs(lengths[0] - lengths[1]) > 1e-3 * lengths[0] or min(abs(cosine - 0.5), cosine) > 1e-3
    positions, labels = [], []
    for segment in range(len(kpoints) // per_segment):
        for index, label in ((segment * per_segment, ends[2 * segment]), ((segment + 1) * per_segment - 1, ends[2 * segment + 1])):
            trim = identify_trim(kpoints[index]) if distorted else None
            label = {"Gamma": "Γ"}.get(trim, trim) if trim else label
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

def extract_eigenvalue_grid(directory, bands):
    # Eigenvalues of the selected (0-based) bands on the full Gamma-centred SCF grid, shape (M1, M2, len(bands))
    with open(os.path.join(directory, "scf", "KPOINTS"), "r", encoding="utf-8") as f:
        mesh = [int(v) for v in f.readlines()[3].split()[:2]]
    with h5py.File(os.path.join(directory, "scf", "vaspout.h5"), "r") as f:
        eigenvalues = f["/results/electron_eigenvalues/eigenvalues"][()][0][:, bands]
        kpoints = f["/results/electron_eigenvalues/kpoint_coords"][()][:, :2]
    grid = np.full(mesh + [len(bands)], np.nan)
    disagreement = 0.0
    for rotation in extract_point_group_2d(directory):
        for sign in (1, -1):
            images = sign * kpoints @ rotation * mesh
            if np.abs(images - np.rint(images)).max() > 1e-4:
                raise ValueError(f"Symmetry image off the SCF grid in {directory}")
            index = np.rint(images).astype(int) % mesh
            known = ~np.isnan(grid[index[:, 0], index[:, 1], 0])
            if known.any():
                disagreement = max(disagreement, float(np.abs(grid[index[known, 0], index[known, 1]] - eigenvalues[known]).max()))
            grid[index[:, 0], index[:, 1]] = eigenvalues
    if np.isnan(grid).any() or disagreement > 1e-4:
        raise ValueError(f"Unfolding incomplete or inconsistent ({disagreement:.2e} eV) in {directory}")
    return grid, mesh

def extract_lindhard_susceptibility(directory, sigma=0.05, window=1.0):
    # Constant-matrix-element static susceptibility chi0(q) and nesting function xi(q) on the SCF q grid.
    # chi0(q) = 1/N sum_{k,n,m} [f(e_nk) - f(e_m,k+q)] / (e_m,k+q - e_nk), Gaussian occupations of width sigma;
    # xi(q)   = 1/N sum_{k,n,m} d(e_nk - E_F) d(e_m,k+q - E_F). Bands with states within +/- window of E_F are used.
    from scipy.special import erfc
    fermi = extract_fermi_topology(directory)
    with h5py.File(os.path.join(directory, "scf", "vaspout.h5"), "r") as f:
        all_bands = f["/results/electron_eigenvalues/eigenvalues"][()][0]
    bands = [n for n in range(all_bands.shape[1])
             if all_bands[:, n].min() < fermi + window and all_bands[:, n].max() > fermi - window]
    grid, mesh = extract_eigenvalue_grid(directory, bands)
    x = (grid - fermi) / sigma
    occupation = 0.5 * erfc(x)
    delta = np.exp(-x ** 2) / (sigma * np.sqrt(np.pi))
    chi0, xi = np.zeros(mesh), np.zeros(mesh)
    count = mesh[0] * mesh[1]
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            shifted = np.roll(grid, (-i, -j), axis=(0, 1))
            f_shift = np.roll(occupation, (-i, -j), axis=(0, 1))
            d_shift = np.roll(delta, (-i, -j), axis=(0, 1))
            de = shifted[:, :, None, :] - grid[:, :, :, None]
            df = occupation[:, :, :, None] - f_shift[:, :, None, :]
            small = np.abs(de) < 1e-6
            ratio = np.where(small, delta[:, :, :, None], df / np.where(small, 1.0, de))
            chi0[i, j] = ratio.sum() / count
            xi[i, j] = (delta[:, :, :, None] * d_shift[:, :, None, :]).sum() / count
    return {"chi0": chi0, "xi": xi, "mesh": mesh, "bands_one_based": [bands[0] + 1, bands[-1] + 1],
            "fermi": fermi, "sigma": sigma}

def extract_symmetry_images(directory, kpoint):
    images = {tuple(np.round(((sign * np.asarray(kpoint[:2]) @ rotation) + 0.5) % 1 - 0.5, 9))
              for rotation in extract_point_group_2d(directory) for sign in (1, -1)}
    return np.array(sorted(images))

def extract_eigenval(file_path):
    # Fractional k, weights and band energies (eV) of a non-spin-polarised EIGENVAL
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    nk, nb = (int(v) for v in lines[5].split()[1:3])
    kpoints, weights, energies = [], [], []
    for index in range(7, 7 + nk * (nb + 2), nb + 2):
        values = lines[index].split()
        kpoints.append([float(v) for v in values[:3]])
        weights.append(float(values[3]))
        energies.append([float(lines[index + 1 + band].split()[1]) for band in range(nb)])
    return np.array(kpoints), np.array(weights), np.array(energies)

def extract_hse_bands(directory):
    # HSE06 (EIGENVAL) and same-cell PBE (EIGENVAL.pbe) of the zero-weight chunks band_<nn>_<from>-<to>_<i> in
    # 4.4_bandstructure_hse/<structure>; edge chunks (band_<nn>_edge-*) and the weighted scf IBZ only enter the extrema
    with open(os.path.join(directory, "scf", "POSCAR"), "r", encoding="utf-8") as f:
        lines = f.readlines()
    lattice = np.array([[float(v) for v in line.split()[:3]] for line in lines[2:5]]) * float(lines[1].split()[0])
    reciprocal = (2 * np.pi * np.linalg.inv(lattice).T)[:2, :2]
    with open(os.path.join(directory, "scf", "OUTCAR"), "r", encoding="utf-8", errors="replace") as f:
        occupied = int(round(float(next(line for line in f if "NELECT" in line).split()[2]))) // 2
    chunks = sorted(glob.glob(os.path.join(directory, "band_*")))
    results = {"occupied_bands": occupied}
    for functional, file_name in (("HSE06", "EIGENVAL"), ("PBE", "EIGENVAL.pbe")):
        path_k, path_e, all_k, all_e, ticks, segment = [], [], [], [], [], None
        for chunk in chunks:
            kpoints, weights, energies = extract_eigenval(os.path.join(chunk, file_name))
            zero = weights == 0
            all_k.append(kpoints[zero]); all_e.append(energies[zero])
            name = os.path.basename(chunk).split("_")[2]
            if name.startswith("edge"):
                continue
            if name != segment:  # first chunk of a segment: tick at its start point
                ticks.append((sum(len(k) for k in path_k), name.split("-")[0]))
                segment = name
            path_k.append(kpoints[zero]); path_e.append(energies[zero])
        ticks.append((sum(len(k) for k in path_k) - 1, segment.split("-")[1]))
        kpoints, weights, energies = extract_eigenval(os.path.join(directory, "scf", file_name))
        all_k.append(kpoints); all_e.append(energies)
        path_k, path_e = np.concatenate(path_k)[:, :2], np.concatenate(path_e)
        all_k, all_e = np.concatenate(all_k)[:, :2], np.concatenate(all_e)
        steps = np.linalg.norm(np.diff(path_k @ reciprocal, axis=0), axis=1)
        kpath = np.concatenate([[0.0], np.cumsum(steps)])
        valence, conduction = all_e[:, occupied - 1], all_e[:, occupied]
        direct = conduction - valence
        results[functional] = {"kpath": kpath, "energies": path_e, "kpoints": path_k,
                               "ticks": [(kpath[i], label) for i, label in ticks],
                               "vbm": (valence.max(), all_k[valence.argmax()]),
                               "cbm": (conduction.min(), all_k[conduction.argmin()]),
                               "direct": (direct.min(), all_k[direct.argmin()]), "points": len(all_k)}
    return results

def extract_hse_parity(directory):
    # hse_parity.json written by 9.0_topology/tools/hse_parity.py into 4.4_bandstructure_hse/<structure>/scf
    return extract_json(os.path.join(directory, "scf", "hse_parity.json"))

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

def extract_incar_tag(incar_path, tag):
    # Active (uncommented) value of one INCAR tag, or None
    if not os.path.isfile(incar_path):
        return None
    with open(incar_path, "r", encoding="utf-8") as f:
        for line in f:
            text = line.split("#")[0].strip()
            if text.upper().startswith(tag.upper()) and "=" in text:
                return text.split("=", 1)[1].strip()
    return None

def summarize_phonon_stability(matters_list):
    # Phonopy finite-displacement runs: [label, phonopy parent directory, optional symprec], one row per run
    # A symprec of 1e-3 restores the hexagonal group of cells that deviate from it by ~1e-4 Å; with 1e-5 phonopy
    # finds Cmmm / C2/m and breaks the C6 equivalence of supercell images in the interpolation.
    import phonopy
    rows = []
    for current_matter in matters_list:
        label, directory, *optional = current_matter
        symprec = optional[0] if optional else 1e-5
        displacements = sorted(glob.glob(os.path.join(directory, "disp-*")))
        incar = os.path.join(displacements[0], "INCAR") if displacements else os.path.join(directory, "INCAR")
        ivdw, sigma = extract_incar_tag(incar, "IVDW"), extract_incar_tag(incar, "SIGMA")
        kpoints = os.path.join(os.path.dirname(incar), "KPOINTS")
        mesh_k = "×".join(open(kpoints, encoding="utf-8").readlines()[3].split()[:2]) if os.path.isfile(kpoints) else "?"
        setting = f"PBE{'+D3' if ivdw else ''}, SIGMA {sigma}, k {mesh_k}"
        if not os.path.isfile(os.path.join(directory, "FORCE_SETS")):
            done = sum(os.path.isfile(os.path.join(d, "OUTCAR")) and "General timing" in open(os.path.join(d, "OUTCAR"), errors="replace").read()
                       for d in displacements)
            rows.append([label, directory, setting, "", f"pending ({done}/{len(displacements)} displacements finished)", "", "", ""])
            continue
        yaml = os.path.join(directory, "phonopy_disp.yaml")
        yaml = yaml if os.path.isfile(yaml) else os.path.join(directory, "phonopy.yaml")
        ph = phonopy.load(yaml, force_sets_filename=os.path.join(directory, "FORCE_SETS"), produce_fc=True, log_level=0,
                          symprec=symprec)
        first = ph.dataset["first_atoms"]
        pairs = [(np.array(a["forces"]) + np.array(b["forces"])) / 2 for i, a in enumerate(first) for b in first[i+1:]
                 if a["number"] == b["number"] and np.allclose(a["displacement"], -np.array(b["displacement"]), atol=1e-8)]
        residual = f"{np.abs(np.mean(pairs, axis=0)).max():.1e}" if pairs else "—"
        mesh = [int(round(v)) for v in np.diag(ph.supercell_matrix)[:2]]
        commensurate = [(i / mesh[0], j / mesh[1], 0) for i in range(mesh[0]) for j in range(mesh[1])]
        ph.run_qpoints(commensurate); f_comm = ph.get_qpoints_dict()["frequencies"]
        fine = [(i / 60, j / 60, 0) for i in range(60) for j in range(60)]
        ph.run_qpoints(fine); f_fine = ph.get_qpoints_dict()["frequencies"]
        index = int(np.argmin(f_fine.min(axis=1)))
        q_min = np.round(((np.array(fine[index][:2]) + 0.5) % 1 - 0.5), 3)
        rows.append([label, directory, setting, ph.symmetry.get_international_table(), residual, f"{f_comm.min():+.3f}",
                     f"{f_fine.min():+.3f} at q = ({q_min[0]:g}, {q_min[1]:g})", f"{f_fine.max():.1f}"])
    header = ["Structure", "Directory", "Forces", "Space group", "Residual force (eV/Å)", "Min. ν, commensurate q (THz)",
              "Min. ν, 60×60 grid (THz)", "Max. ν (THz)"]
    return show_table(header, rows)

def summarize_frozen_phonon(directory):
    # Frozen-phonon energies E0(A) - E0(A0) from the OSZICAR of each amplitude folder
    rows, reference = [], None
    for folder in sorted(glob.glob(os.path.join(directory, "A*"))):
        oszicar = os.path.join(folder, "OSZICAR")
        final = [line for line in open(oszicar, encoding="utf-8") if " E0= " in line] if os.path.isfile(oszicar) else []
        if not final:
            rows.append([os.path.basename(folder), "pending", ""])
            continue
        energy = float(final[-1].split("E0=")[1].split()[0])
        reference = energy if reference is None else reference
        rows.append([os.path.basename(folder), f"{energy:.6f}", f"{1000 * (energy - reference):+.2f}"])
    return show_table(["Amplitude", "E0 (eV)", "E0 − E0(A0) (meV)"], rows)

def extract_poscar_cartesian(file_path):
    # Lattice (Å), element symbol of every atom and Cartesian positions (Å) of a VASP5 POSCAR in direct coordinates
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    lattice = np.array([[float(v) for v in line.split()[:3]] for line in lines[2:5]]) * float(lines[1].split()[0])
    symbols = [s for s, n in zip(lines[5].split(), lines[6].split()) for _ in range(int(n))]
    start = 9 if lines[7].strip()[0] in "Ss" else 8
    fractional = np.array([[float(v) for v in line.split()[:3]] for line in lines[start:start + len(symbols)]])
    return lattice, symbols, fractional

def extract_frozen_phonon_fit(directory):
    # Fit E(Q) - E(0) = alpha Q^2 + beta Q^4 (meV per supercell) in mass-weighted amplitude Q (amu^1/2 Å) over the A* folders
    from phonopy.structure.atoms import atom_data, symbol_map
    lattice, symbols, reference = extract_poscar_cartesian(os.path.join(directory, "A0", "POSCAR"))
    masses = np.array([atom_data[symbol_map[s]][3] for s in symbols])
    energies, q2, umax = [], [], []
    for folder in sorted(glob.glob(os.path.join(directory, "A*"))):
        final = [line for line in open(os.path.join(folder, "OSZICAR"), encoding="utf-8") if " E0= " in line]
        displacement = extract_poscar_cartesian(os.path.join(folder, "POSCAR"))[2] - reference
        displacement = (displacement - np.round(displacement)) @ lattice
        energies.append(float(final[-1].split("E0=")[1].split()[0]))
        q2.append(float((masses * (displacement ** 2).sum(axis=1)).sum()))
        umax.append(float(np.sqrt((displacement ** 2).sum(axis=1)).max()))
    energies = 1000 * (np.array(energies) - energies[int(np.argmin(q2))])
    q2, umax = np.array(q2), np.array(umax)
    order = np.argsort(q2)[1:]
    matrix = np.vstack([q2[order], q2[order] ** 2]).T
    (alpha, beta), *_ = np.linalg.lstsq(matrix, energies[order], rcond=None)
    residual = np.abs(energies[order] - matrix @ np.array([alpha, beta])).max()
    return {"alpha": alpha, "beta": beta, "residual": residual, "energies": energies[order], "q2": q2[order], "umax": umax[order]}

def summarize_frozen_phonon_fit(matters_list):
    # Fit E(Q) - E(0) = alpha Q^2 + beta Q^4 in mass-weighted amplitude Q (amu^1/2 Å) over the A* folders:
    # matters [label, directory]; harmonic nu = sqrt(2 alpha) / 2 pi (imaginary shown negative); well depth alpha^2 / 4 beta
    rows = []
    for label, directory in matters_list:
        fit = extract_frozen_phonon_fit(directory)
        alpha, beta, residual, energies, q2, umax = (fit[key] for key in ("alpha", "beta", "residual", "energies", "q2", "umax"))
        omega2 = 2 * alpha / 0.1036427  # meV per amu Å^2 (rad/ps)^2
        nu = np.sign(omega2) * np.sqrt(abs(omega2)) / (2 * np.pi)
        well = (f"{alpha ** 2 / (4 * beta):.3f} at max \\|u\\| = {np.sqrt(-alpha / (2 * beta) / q2[0]) * umax[0]:.3f} Å"
                if alpha < 0 else "none")
        kpoints = open(os.path.join(directory, "A0", "KPOINTS"), encoding="utf-8").readlines()[3].split()[:2]
        rows.append([label, extract_incar_tag(os.path.join(directory, "A0", "INCAR"), "SIGMA"), "×".join(kpoints),
                     ", ".join(f"{e:+.3f} ({u:.3f} Å)" for e, u in zip(energies, umax)),
                     f"{nu:+.2f}", well, f"{residual:.3f}"])
    header = ["Series", "SIGMA (eV)", "Supercell k", "E − E(A0) in meV (max \\|u\\|)", "Harmonic ν (THz, − = imaginary)",
              "Double-well depth (meV / cell)", "Max. fit residual (meV)"]
    return show_table(header, rows)

def summarize_frozen_phonon_quantum(matters_list, temperatures=(0, 300)):
    # One-mode quantum treatment of the quartic fit V(Q) = alpha Q^2 + beta Q^4 (kinetic energy Q'^2 / 2, Q mass-weighted):
    # exact eigenstates on a finite-difference grid, and the self-consistent harmonic frequency of the single mode,
    # c Omega^2 = 2 alpha + 12 beta <Q^2>, <Q^2> = hbar / (2 c Omega) coth(hbar Omega / 2 kT). Coupling to the other modes is
    # left out, so this indicates the size of the quantum fluctuations; it is not a substitute for a full SSCHA.
    from scipy.linalg import eigh_tridiagonal
    from scipy.optimize import brentq
    hbar, c, kB, thz = 0.6582119569, 0.1036427, 0.08617333262, 4.135667696  # meV ps, meV/(amu Å^2 ps^-2), meV/K, meV/THz
    kinetic = hbar ** 2 / c / 2
    grid = np.linspace(-3.0, 3.0, 12001); step = grid[1] - grid[0]
    rows = []
    for label, directory in matters_list:
        fit = extract_frozen_phonon_fit(directory)
        alpha, beta = fit["alpha"], fit["beta"]
        levels, states = eigh_tridiagonal(alpha * grid ** 2 + beta * grid ** 4 + 2 * kinetic / step ** 2,
                                          -kinetic / step ** 2 * np.ones(len(grid) - 1), select="i", select_range=(0, 1))
        ground = states[:, 0] ** 2 / step
        if ground[0] > 1e-12 * ground.max():
            raise ValueError(f"Ground state of {label} reaches the edge of the Q grid")
        harmonic = np.sign(alpha) * np.sqrt(2 * abs(alpha) / c) / (2 * np.pi)
        depth = alpha ** 2 / (4 * beta) if alpha < 0 else 0.0
        q2 = lambda omega, t: hbar / (2 * omega * c) * (1 / np.tanh(hbar * omega / (2 * kB * t)) if t > 0 else 1.0)
        scha = [brentq(lambda omega: c * omega ** 2 - 2 * alpha - 12 * beta * q2(omega, t), 1e-3, 500) / (2 * np.pi)
                for t in temperatures]
        rows.append([label, f"{harmonic:+.2f}", f"{depth:.3f}", f"{levels[0] + depth:.2f}", f"{(levels[1] - levels[0]) / thz:.2f}",
                     *[f"{nu:.2f}" for nu in scha],
                     f"{np.sqrt((ground * grid ** 2).sum() * step):.3f} / {np.sqrt(max(-alpha / (2 * beta), 0)):.3f}"])
    header = ["Series", "Harmonic ν (THz)", "Well depth (meV)", "E₀ above the well bottom (meV)", "Exact E₁ − E₀ (THz)",
              *[f"One-mode SCHA ν at {t:g} K (THz)" for t in temperatures], "√⟨Q²⟩₀ / Q at the minimum (amu^½ Å)"]
    return show_table(header, rows)

def extract_final_energy(directory):
    # Final E0 (eV) of a VASP run from its OSZICAR
    with open(os.path.join(directory, "OSZICAR"), "r", encoding="utf-8") as f:
        return float([line for line in f if " E0= " in line][-1].split("E0=")[1].split()[0])

def extract_phonon_free_energy(directory, temperatures, mesh=40):
    # Harmonic vibrational free energy per cell (eV) at each temperature (0 K = ZPE), symprec as recorded in phonopy.yaml;
    # imaginary and zero modes are left out
    import phonopy, yaml
    with open(os.path.join(directory, "phonopy.yaml"), "r", encoding="utf-8") as f:
        symprec = float(yaml.safe_load(f)["phonopy"]["symmetry_tolerance"])
    ph = phonopy.load(os.path.join(directory, "phonopy_disp.yaml"), force_sets_filename=os.path.join(directory, "FORCE_SETS"),
                      produce_fc=True, log_level=0, symprec=symprec)
    ph.run_mesh([mesh, mesh, 1])
    data = ph.get_mesh_dict()
    energies = 4.135667696e-3 * np.where(data["frequencies"] > 0, data["frequencies"], 0.0)  # THz -> eV
    weights = data["weights"] / data["weights"].sum()
    result = []
    for temperature in temperatures:
        free = 0.5 * energies
        if temperature > 0:
            kt = 8.617333262e-5 * temperature
            free = free + kt * np.log1p(-np.exp(-np.where(energies > 0, energies, np.inf) / kt))
        result.append(float((np.where(energies > 0, free, 0.0) * weights[:, None]).sum()))
    return result

def extract_h2_zero_point(directory):
    # Zero-point energy (eV) of H2 from the stretching mode of an IBRION = 5 run (rotational residues excluded)
    with open(os.path.join(directory, "OUTCAR"), "r", encoding="utf-8", errors="replace") as f:
        modes = [float(line.split()[-2]) for line in f if " f  =" in line and "meV" in line]
    return max(modes) / 2000

def summarize_hydrogenation_energetics(matters_list, h2_dir, h2_freq_dir, temperature=298.15, delta_mu_h2=-0.31605):
    # Adsorption energy per H against 1/2 H2: electronic, + ZPE, and the harmonic free energy at temperature and 1 bar.
    # matters [label, relaxed slab, slab phonons, relaxed base, base phonons, number of H];
    # delta_mu_h2 = H(T) - H(0) - T S of H2 at 1 bar (JANAF, 298.15 K: 8.467 kJ/mol, 130.680 J/mol/K)
    e_h2 = extract_final_energy(h2_dir)
    z_h2 = extract_h2_zero_point(h2_freq_dir)
    rows = []
    for label, slab, slab_ph, base, base_ph, count in matters_list:
        energy = extract_final_energy(slab) - extract_final_energy(base) - count / 2 * e_h2
        f_slab, f_base = extract_phonon_free_energy(slab_ph, [0, temperature]), extract_phonon_free_energy(base_ph, [0, temperature])
        zpe = energy + f_slab[0] - f_base[0] - count / 2 * z_h2
        free = energy + f_slab[1] - f_base[1] - count / 2 * (z_h2 + delta_mu_h2)
        pressure = np.exp(2 * free / count / (8.617333262e-5 * temperature))
        rows.append([label, count, f"{energy / count:+.3f}", f"{zpe / count:+.3f}", f"{free / count:+.3f}", f"{pressure:.1e}"])
    header = ["Structure", "H per cell", "ΔE per H (eV)", "ΔE + ΔZPE per H (eV)", f"ΔG({temperature:g} K, 1 bar) per H (eV)",
              f"Equilibrium p(H₂) at {temperature:g} K (bar)"]
    return show_table(header, rows)

def summarize_polymorph_energies(matters_list):
    # Relaxed total energy per formula unit relative to the first entry: [label, relaxed directory, formula units]
    import spglib
    rows, reference = [], None
    for label, directory, units in matters_list:
        energy = extract_final_energy(directory) / units
        reference = energy if reference is None else reference
        lattice, symbols, fractional = extract_poscar_cartesian(os.path.join(directory, "CONTCAR"))
        numbers = [sorted(set(symbols)).index(s) for s in symbols]
        rows.append([label, directory, spglib.get_spacegroup((lattice, fractional, numbers), symprec=1e-3),
                     f"{energy:.6f}", f"{1000 * (energy - reference):+.1f}"])
    return show_table(["Structure", "Directory", "Space group (1e-3 Å)", "E0 per formula unit (eV)", "Relative (meV)"], rows)

def summarize_lindhard(matters_list, sigmas=(0.02, 0.05, 0.10)):
    # Rank of chosen q points in chi0(q) and xi(q): matters [label, structure directory, (q1, q2)]
    rows = []
    for label, directory, qpoint in matters_list:
        for sigma in sigmas:
            result = extract_lindhard_susceptibility(directory, sigma=sigma)
            mesh = result["mesh"]
            index = tuple(int(round(v * m)) % m for v, m in zip(qpoint, mesh))
            mask = np.ones(mesh, bool); mask[0, 0] = False
            cells = []
            for name in ("chi0", "xi"):
                grid = result[name]
                peak = np.unravel_index(np.argmax(np.where(mask, grid, -np.inf)), mesh)
                rank = int((grid[mask] > grid[index]).sum()) + 1
                cells += [f"{grid[index]:.3f}", f"{rank}/{mask.sum()}",
                          f"{grid[peak]:.3f} at ({(peak[0] / mesh[0] + 0.5) % 1 - 0.5:.3f}, {(peak[1] / mesh[1] + 0.5) % 1 - 0.5:.3f})"]
            rows.append([label, f"({qpoint[0]:g}, {qpoint[1]:g})", f"{sigma:.2f}", *cells])
    header = ["Structure", "q", "σ (eV)", "χ₀(q)", "χ₀ rank", "χ₀ max (q ≠ 0)", "ξ(q)", "ξ rank", "ξ max (q ≠ 0)"]
    return show_table(header, rows)

def summarize_hse_parity(matters_list):
    # Spinless HSE06 parities at the TRIM against the PBE+SOC screen: [label, 4.4_bandstructure_hse/<structure>]
    rows = []
    for label, directory in matters_list:
        parity = extract_hse_parity(directory)
        if parity is None:
            rows.append([label, "pending", "", "", "", "", ""])
            continue
        trims = parity["trims"]
        hse = " ".join(str(trims[name]["odd_bands"]) for name in trim_points)
        pbe = " ".join(str(parity["reference"]["odd_kramers_pairs"][name]) for name in trim_points)
        # Smallest gap above band n at a TRIM; a boundary inside an orbital doublet is split by SOC only
        split = [name for name in trim_points if trims[name]["boundary_inside_degenerate_block"]]
        name = min(trim_points, key=lambda n: trims[n]["gap_above_degenerate_block_ev"] if n in split else trims[n]["gap_above_n_ev"])
        shown = {"Gamma": "Γ"}.get(name, name)
        if split:
            critical = f"{'/'.join(split)}: SOC-split doublet (parity {trims[split[0]]['next_block_parity']:+d})"
        else:
            critical = f"{shown}: {trims[name]['gap_above_n_ev']:.3f} / {trims[name]['pbe_gap_above_n_ev']:.3f}"
        rows.append([label, 2 * parity["spinless_bands"], hse, pbe, parity["conditional_fu_kane_nu"],
                     parity["reference"]["conditional_fu_kane_nu"], critical])
    header = ["Structure", "N", "Odd pairs Γ X Y M (HSE06)", "Odd pairs Γ X Y M (PBE+SOC)", "ν (HSE06)", "ν (PBE+SOC)",
              "Smallest E_N+1 − E_N at a TRIM, HSE06 / PBE (eV)"]
    return show_table(header, rows)

def summarize_hse_gap(matters_list):
    # Band edges from all HSE06 chunks: [label, 4.4_bandstructure_hse/<structure>, optional 9.0_topology/<structure>]
    rows = []
    for current_matter in matters_list:
        label, directory, *optional = current_matter
        bands = extract_hse_bands(directory)
        for functional in ("PBE", "HSE06"):
            result = bands[functional]
            gap = result["cbm"][0] - result["vbm"][0]
            rows.append([label, f"{functional}, no SOC, 20 Å, 12×12 + {result['points'] - 43} path/edge k",
                         "(%.4f, %.4f)" % tuple(result["vbm"][1]), "(%.4f, %.4f)" % tuple(result["cbm"][1]),
                         f"{gap:.3f}", f"{result['direct'][0]:.3f} at ({result['direct'][1][0]:.4f}, {result['direct'][1][1]:.4f})"])
        if optional:
            reference = extract_gap_extrema(optional[0])[0][0]
            rows.append([label, "PBE+SOC, 40 Å, 105×105 + refinements (9.0_topology)",
                         "(%.4f, %.4f)" % tuple(reference["valence_max_k"][:2]), "(%.4f, %.4f)" % tuple(reference["conduction_min_k"][:2]),
                         f"{reference['indirect_gap_ev']:.3f}",
                         f"{reference['direct_gap_ev']:.3f} at ({reference['direct_gap_k'][0]:.4f}, {reference['direct_gap_k'][1]:.4f})"])
    header = ["Structure", "Calculation", "VBM k", "CBM k", "Indirect gap (eV)", "Min. direct gap (eV)"]
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
    Article version: up to three structures per row (three structures give a 1x3 row with the legend below);
    thesis version (figure_version("thesis")): two columns, so three structures give a 2x2 grid with the legend
    in the fourth panel.
    """
    if suptitle in ["help", "Help"]:
        print(help_info)
        return

    matters = create_matters_topology(matters_list)
    thesis = figure_version() == "thesis"
    columns = min(2 if thesis else 3, len(matters))
    rows = math.ceil(len(matters) / columns)
    legend_cell = (rows - 1, columns - 1) if legend_loc and rows * columns > len(matters) else None

    # Figure settings (margins in inches: suptitle above, legend or tick labels below)
    top_space = .65 if thesis else 1.15
    bottom_space = 1.15 if legend_loc and legend_cell is None else .5
    fig_setting = topology_canvas("bands", figure_version())
    width, height = fig_setting[0]
    params = fig_setting[2]; plt.rcParams.update(params)
    fig = plt.figure(figsize=fig_setting[0], dpi=fig_setting[1])
    outer = gridspec.GridSpec(rows, columns, figure=fig, left=(.8 if thesis else 1.0) / width, right=1 - .2 / width,
                              top=1 - top_space / height, bottom=bottom_space / height, hspace=.16 if thesis else .15, wspace=.20)

    # Colors calling
    fermi_color = color_sampling("Grey")
    annotate_color = color_sampling("Grey")

    # Title
    fig.suptitle(f"{suptitle}", fontsize=fig_setting[3][0], y=1 - 0.25 / height, va="top")

    # Data calling and plotting
    energy_window = eigen_range if eigen_range is not None else (-11, 5)
    for index, matter in enumerate(matters):
        label, directory, entry, subspaces, colors = matter
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[index // columns, index % columns],
                                                 height_ratios=[3.2, 1.3], hspace=0.06)
        ax_bands = fig.add_subplot(inner[0])
        ax_bands._topology_band_panel = True
        ax_gap = fig.add_subplot(inner[1], sharex=ax_bands)
        kpath, eigenvalues, kpoints, positions, labels = extract_band_path(directory)
        eigenvalues = eigenvalues - extract_fermi_topology(directory)

        # Subspace bands, higher bands and Fermi energy
        lower = 0
        for color, subspace in zip(colors, subspaces):
            for band in range(lower, subspace):
                ax_bands.plot(kpath, eigenvalues[:, band], c=color, lw=LINE_WIDTH, zorder=4)
            lower = subspace
        for band in range(lower, eigenvalues.shape[1]):
            ax_bands.plot(kpath, eigenvalues[:, band], c=annotate_color[2], lw=LINE_WIDTH, zorder=3)
        ax_bands.axhline(y=0, **FERMI_STYLE, zorder=2)
        ax_bands.set_ylim(*energy_window)

        # Kramers-pair parities at the TRIM on the path
        parity = extract_parity(directory)
        if parity:
            annotate_parities(ax_bands, kpath, eigenvalues, kpoints, parity, annotate_color[0])

        # Direct gap along the path and the BZ minimum
        gaps, _ = extract_gap_extrema(directory)
        gap_handles = []
        for color, subspace in zip(colors, subspaces):
            gap = next(g for g in gaps if g["N"] == subspace)
            ax_gap.semilogy(kpath, np.maximum(eigenvalues[:, subspace] - eigenvalues[:, subspace - 1], 1e-6), c=color, lw=LINE_WIDTH, zorder=4)
            ax_gap.axhline(y=gap["direct_gap_ev"], color=color, lw=LINE_WIDTH, linestyle=":", zorder=3)
            gap_handles.append(Line2D([], [], color=color, lw=LINE_WIDTH,
                                      label=f"N = {subspace}: {format_gap(gap['direct_gap_ev'])}"))
        ax_gap.legend(handles=gap_handles, loc="lower right", frameon=True, fancybox=True,
                      facecolor="white", framealpha=.9, borderpad=.35, labelspacing=.3,
                      handlelength=1.2, fontsize=STYLES[figure_version()]["note"])


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
        subtitle = f"({chr(97+index)}) {subtitle}"
        if thesis:
            subtitle = subtitle.replace(": Z₂", ":\nZ₂")
        ax_bands.set_title(subtitle, fontsize=fig_setting[3][1])
        ax_bands.set_ylim(energy_window[0], energy_window[1])
        for annotation in ax_bands.texts:
            annotation.set_visible(energy_window[0] <= annotation.get_position()[1] <= energy_window[1])
        ax_bands.set_xlim(kpath[0], kpath[-1])
        ax_gap.set_ylim(3e-4, 30)
        ax_gap.set_xticks(positions)
        ax_gap.set_xticklabels(labels)
        ax_bands.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True, labelbottom=False)
        ax_gap.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)
        if index % columns == 0:
            ax_bands.set_ylabel(r"$E-E_\mathrm{F}$ (eV)")
            ax_gap.set_ylabel(r"$E_{N+1}-E_N$ (eV)")

    # Legend: below the row (article) or in the free panel of the grid (thesis)
    if legend_loc:
        handles = [Line2D([], [], c=color_sampling("Blue")[1], lw=LINE_WIDTH, label="lowest-N subspace (bands 1…N)")]
        for label, directory, entry, subspaces, colors in matters:
            if len(subspaces) > 1:
                separator = "\n" if legend_cell is not None else " "
                handles.append(Line2D([], [], c=colors[1], lw=LINE_WIDTH,
                                      label=f"{label} bands {subspaces[0] + 1}–{subspaces[1]}{separator}(with blue: lowest-{subspaces[1]} subspace)"))
                break
        handles += [Line2D([], [], c=annotate_color[2], lw=LINE_WIDTH, label="bands above selected subspace(s)"),
                    Line2D([], [], c=annotate_color[0], lw=LINE_WIDTH, linestyle=":", label=("sampled BZ minimum of\n" + r"$E_{N+1}-E_N$") if legend_cell is not None else r"sampled BZ minimum of $E_{N+1}-E_N$"),
                    Line2D([], [], **FERMI_STYLE, label="Fermi energy")]
        if legend_cell is not None:
            ax_legend = fig.add_subplot(outer[legend_cell[0], legend_cell[1]])
            ax_legend.axis("off")
            ax_legend.legend(handles=handles, loc="center", frameon=True, fancybox=True, facecolor="white", framealpha=.9, fontsize=params["legend.fontsize"], borderpad=.45, labelspacing=.6)
        else:
            fig.legend(handles=handles, loc="lower center", ncol=3, frameon=True, fancybox=True, facecolor="white", framealpha=.9, bbox_to_anchor=(.5, .015), borderpad=.4, labelspacing=.4)

    # Use one label position for both the energy and direct-gap axes.
    for ax in fig.axes:
        if ax.get_ylabel():
            ax.yaxis.set_label_coords(-.16 if thesis else -.15, .5)

def annotate_parities(ax, kpath, eigenvalues, kpoints, parity, text_color):
    # Centre the actual glyph outlines, including mixed signs, in each circle.
    for index in range(len(kpoints)):
        name = identify_trim(kpoints[index])
        if name is None or (index > 0 and abs(kpath[index] - kpath[index - 1]) < 1e-12):
            continue
        pairs = [p for block in parity["trims"][name]["blocks"] for p in block["kramers_pair_parities_unordered_within_block"]]
        groups = []
        for pair, value in enumerate(pairs):
            energy = eigenvalues[index, 2 * pair:2 * pair + 2].mean()
            if groups and energy - groups[-1][0] < .45:
                groups[-1][1].append(value)
            else:
                groups.append([energy, [value]])
        for energy, values in groups:
            if not ax.get_ylim()[0] <= energy <= ax.get_ylim()[1]:
                continue
            signs = ",".join(format_sign(v) for v in values)
            glyph = TextPath((0, 0), signs, size=STYLES[figure_version()]["note"],
                             prop=FontProperties(family="serif"))
            bounds = glyph.get_extents()
            side = max(bounds.width, bounds.height) + 5
            drawing = DrawingArea(side, side, 0, 0)
            drawing.add_artist(Circle((side/2, side/2), side/2, fc="white", ec=text_color, lw=.8, alpha=.9))
            transform = Affine2D().translate(side/2-(bounds.x0+bounds.x1)/2, side/2-(bounds.y0+bounds.y1)/2)
            drawing.add_artist(PathPatch(glyph.transformed(transform), fc=text_color, lw=0))
            offset = (-1 if kpath[index] > .97*kpath[-1] else 1) * (side/2+2)
            annotation = AnnotationBbox(drawing, (kpath[index], energy), xybox=(offset, 0),
                                       boxcoords="offset points", frameon=False, pad=0, zorder=8)
            annotation._parity_signs = signs
            ax.add_artist(annotation)

def plot_direct_gap_maps(suptitle, matters_list=None, gap_range=None):
    """Direct gaps in reciprocal fractions: k = k1*b1 + k2*b2; minima are sampled values."""
    matters = create_matters_topology(matters_list)
    panels = [(matter, subspace) for matter in matters for subspace in matter[3]]
    thesis = figure_version() == "thesis"
    columns = min(2 if thesis else 4, len(panels))
    rows = math.ceil(len(panels) / columns)
    free_cells = rows * columns - len(panels)
    gap_limits = gap_range if gap_range is not None else (1e-2, 10)
    norm = LogNorm(vmin=gap_limits[0], vmax=gap_limits[1])
    ticks = np.logspace(np.ceil(np.log10(gap_limits[0])), np.floor(np.log10(gap_limits[1])),
                        int(np.floor(np.log10(gap_limits[1])) - np.ceil(np.log10(gap_limits[0]))) + 1)
    cmap = make_heatmap_colormap(norm, ticks, reverse=True)
    settings = topology_canvas("gap_maps", figure_version())
    width, height = settings[0]
    if not thesis:
        width = round(width * columns / 4)
    plt.rcParams.update(settings[2])
    fig, axes = plt.subplots(rows, columns, figsize=(width, height), dpi=settings[1], squeeze=False)
    fig.suptitle(suptitle, fontsize=settings[3][0], y=1-.16/height, va="top")
    grey = color_sampling("Grey")[0]
    minimum_color = "#FF4600"
    vector = r"\vec" if thesis else r"\mathbf"
    for index, (ax, (matter, subspace)) in enumerate(zip(axes.flat, panels)):
        label, directory = matter[:2]
        grid, mesh = extract_direct_gap_grid(directory, subspace)
        shift = [m // 2 for m in mesh]
        values = np.roll(grid, shift, axis=(0, 1))
        edges = [(np.arange(m+1)-shift[i]-.5)/m for i, m in enumerate(mesh)]
        k1, k2 = np.meshgrid(*edges, indexing="ij")
        image = ax.pcolormesh(k1, k2, np.clip(values, gap_limits[0]*.5, None), cmap=cmap,
                             norm=norm, shading="flat", rasterized=True)
        for name, trim in trim_points.items():
            ax.plot(*trim, marker="o", ms=6, mfc="white", mec=grey, mew=.9, zorder=4, clip_on=False)
            if name == "Gamma":
                ax.annotate("Γ", trim, xytext=(5, 5), textcoords="offset points", color=grey)
        gaps, _ = extract_gap_extrema(directory)
        best = next(g for g in gaps if g["N"] == subspace)
        zoom = extract_json(os.path.join(directory, "gap_refine_3_summary.json"))
        minima = [b["levels"][-1]["min_k"] for b in zoom["basins"]] if zoom else [best["direct_gap_k"][:2]]
        for point in minima:
            points = extract_symmetry_images(directory, point)
            ax.plot(points[:, 0], points[:, 1], ls="none", marker=r"$\odot$", ms=13,
                    color=minimum_color, zorder=6, clip_on=False)
        system = f"({chr(97+index)}) {label}"
        note = f"N = {subspace}: sampled min {format_gap(best['direct_gap_ev'])}"
        if thesis:
            box = dict(boxstyle="round", fc="white", ec=plt.rcParams["legend.edgecolor"], alpha=.75)
            ax.text(.025, .975, system, transform=ax.transAxes, va="top", zorder=8,
                    fontsize=settings[3][1], bbox=box)
            ax.text(.025, .025, note, transform=ax.transAxes, va="bottom", zorder=8,
                    fontsize=STYLES["thesis"]["note"], bbox=box)
        else:
            ax.set_title(system + "\n" + note, fontsize=settings[3][1])
        extent = .70 if thesis else .54
        ax.set(xlim=(-extent, extent), ylim=(-extent, extent), aspect="equal")
        ax.set_xticks([-.5, 0, .5]); ax.set_yticks([-.5, 0, .5])
        ax.set_xlabel(r"$k_1$ (units of $" + vector + r"{b}_1$)")
        if not thesis or index % columns == 0:
            ax.set_ylabel(r"$k_2$ (units of $" + vector + r"{b}_2$)")
        if thesis and index + columns < len(panels):
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)
    for ax in list(axes.flat)[len(panels):]:
        ax.axis("off")
    handles = [Line2D([], [], ls="none", marker="o", ms=6, mfc="white", mec=grey, mew=.9, label="TRIM"),
               Line2D([], [], ls="none", marker=r"$\odot$", ms=13, color=minimum_color,
                      label="refined minima and their symmetry images")]
    fig.subplots_adjust(left=(.75 if thesis else .95)/width,
                        right=1-(.18 if free_cells else (1.0 if thesis else 1.3))/width,
                        bottom=(.6 if free_cells else 1.05)/height,
                        top=1-(.48 if thesis else 1.1)/height,
                        wspace=.10 if thesis else .26, hspace=.12 if thesis else .45)
    if free_cells:
        ax_free = list(axes.flat)[len(panels)]
        cax = ax_free.inset_axes([.10, .15, .055, .72])
        colorbar = fig.colorbar(image, cax=cax, extend="both", extendfrac=.04)
        handles[1].set_label("refined minima\nand their\nsymmetry images")
        ax_free.legend(handles=handles, loc="center left", bbox_to_anchor=(.45, .5),
                       frameon=True, fancybox=True, facecolor="white", framealpha=.9,
                       handlelength=1., handletextpad=.5, borderpad=.4)
    else:
        fig.canvas.draw()
        boxes = [ax.get_position() for ax in axes.flat]
        bottom, top = min(b.y0 for b in boxes), max(b.y1 for b in boxes)
        right = max(b.x1 for b in boxes)
        colorbar = fig.colorbar(image, cax=fig.add_axes([right+.20/width, bottom, .13/width, top-bottom]), extend="both", extendfrac=.04)
        fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(.5, .015),
                   frameon=True, fancybox=True, facecolor="white", framealpha=.9)
    style_heatmap_colorbar(colorbar, ticks)
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
    fig_setting = topology_canvas("wcc", figure_version())
    params = fig_setting[2]; plt.rcParams.update(params)
    fig, axes = plt.subplots(1, len(panels), figsize=fig_setting[0], dpi=fig_setting[1], squeeze=False)

    # Title
    fig.suptitle(f"{suptitle}", fontsize=fig_setting[3][0], y=.97)

    # Data calling and plotting
    for index, (ax, (matter, color, bands, manifold)) in enumerate(zip(axes.flat, panels)):
        for line_position, centres in zip(manifold["line_positions"], manifold["wcc"]):
            ax.plot([line_position / 2] * len(centres), centres, linestyle="none", marker="o", ms=6, c=color, zorder=4)
        separator = "\n" if figure_version() == "thesis" else " "
        ax.set_title(f"({chr(97+index)}) Lowest {bands} bands:{separator}Z₂ = {manifold['z2']}", fontsize=fig_setting[3][1])
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$k_2$ (units of $\vec{b}_2$)" if figure_version() == "thesis" else r"$k_2$ (units of $\mathbf{b}_2$)")
        ax.set_ylabel(r"WCC along $\vec{a}_1$" if figure_version() == "thesis" else r"WCC along $\mathbf{a}_1$")
        ax.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)

    fig.subplots_adjust(left=.1, right=.98, bottom=.17, top=.84, wspace=.32)

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
    fig_setting = topology_canvas("gap_zoom", figure_version())
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
            # Legend without values: the minimum over all samples is reported with the sampled gap maps and band figures
            plt.plot(radius, gap, c=color, lw=LINE_WIDTH, marker=markers[index], ms=8, label=name, zorder=4)
            positive = bound > 0
            plt.vlines(radius[positive], bound[positive], gap[positive], color=color, lw=LINE_WIDTH, zorder=3)
            plt.plot(radius[positive], bound[positive], linestyle="none", marker="_", ms=12, c=color, zorder=3)

    # Title and axes
    plt.title(f"{title}")
    plt.xscale("log")
    plt.xlim(max(radii) * 1.6, min(radii) / 1.6)
    plt.ylim(0, None)
    plt.xlabel(r"Sampling radius (Å$^{-1}$), finer $\rightarrow$")
    plt.ylabel(r"Min. $E_{N+1}-E_N$ in patch (meV)")

    # Legend (lower left: the curves stay above 0.5 meV and the whiskers reach down at fine radii)
    plt.legend(loc="lower left", frameon=True, fancybox=True, facecolor="white", framealpha=.9)
    plt.tight_layout()


def plot_lindhard_susceptibility(suptitle, matters_list=None, sigma=0.05):
    # Help information
    help_info = """
    Usage: plot_lindhard_susceptibility
        arg[0]: suptitle;
        arg[1]: matters list, [[label, structure directory, (q1, q2)], ...];
        arg[2]: Gaussian width sigma in eV (default 0.05);
    Constant-matrix-element static susceptibility chi0(q) and nesting function xi(q) on the SCF q grid;
    orange-red circle-dot markers mark the chosen q and its symmetry images.
    """
    if suptitle in ["help", "Help"]:
        print(help_info)
        return

    # Keep the Cartesian reciprocal-cell geometry and place each independent colour bar beside its panel.
    ratio = max(np.ptp(corners[:, 0]) / np.ptp(corners[:, 1]) for matter in matters_list
                for corners in [np.array([[0.5, 0.5], [0.5, -0.5], [-0.5, 0.5], [-0.5, -0.5]]) @ extract_reciprocal_2d(matter[1])])
    panel_height, panel_gap = 4.8, 2.0
    panel_width = panel_height * ratio
    width, height = 2 * panel_width + 4.1, 6.1 * len(matters_list)
    fig_setting = canvas_setting(width, height)
    plt.rcParams.update(fig_setting[2])
    fig = plt.figure(figsize=fig_setting[0], dpi=fig_setting[1])
    fig.suptitle(suptitle, fontsize=fig_setting[3][0], y=1 - .10 / height)
    annotate_color = color_sampling("Grey")

    for row, (label, directory, qpoint) in enumerate(matters_list):
        result = extract_lindhard_susceptibility(directory, sigma=sigma)
        mesh = result["mesh"]
        shift = [m // 2 for m in mesh]
        edges = [(np.arange(m + 1) - shift[i] - 0.5) / m for i, m in enumerate(mesh)]
        frac_1, frac_2 = np.meshgrid(edges[0], edges[1], indexing="ij")
        reciprocal = extract_reciprocal_2d(directory)
        qx = frac_1 * reciprocal[0, 0] + frac_2 * reciprocal[1, 0]
        qy = frac_1 * reciprocal[0, 1] + frac_2 * reciprocal[1, 1]
        points = extract_symmetry_images(directory, qpoint) @ reciprocal
        bottom = .65 + (len(matters_list) - row - 1) * 6.1
        for column, (name, title, unit) in enumerate((("chi0", r"$\chi_0(\mathbf{q})$", r"eV$^{-1}$"),
                                                   ("xi", r"$\xi(\mathbf{q})$", r"eV$^{-2}$"))):
            if figure_version() == "thesis":
                title = title.replace(r"\mathbf{q}", r"\vec{q}\,")
            left = .80 + column * (panel_width + panel_gap)
            ax = fig.add_axes([left / width, bottom / height, panel_width / width, panel_height / height])
            values = np.roll(result[name], shift, axis=(0, 1))
            upper = np.delete(result[name].ravel(), 0).max()  # q = 0 (intraband self-nesting) excluded from the colour scale
            clipped = np.minimum(values, upper)
            norm = Normalize(vmin=clipped.min(), vmax=upper)
            ticks = MaxNLocator(nbins=5).tick_values(norm.vmin, norm.vmax)
            ticks = ticks[(ticks >= norm.vmin) & (ticks <= norm.vmax)]
            image = ax.pcolormesh(qx, qy, clipped, cmap=make_heatmap_colormap(norm, ticks), norm=norm,
                                 shading="flat", rasterized=True)
            ax.plot(points[:, 0], points[:, 1], linestyle="none", marker=r"$\odot$", ms=13,
                    color="#FF4600", zorder=5)
            ax.plot(0, 0, marker="o", ms=6, mfc="white", mec=annotate_color[0], mew=.9, zorder=4)
            ax.text(.04, .97, f"{title}\nσ = {sigma:g} eV", transform=ax.transAxes, ha="left", va="top",
                    fontsize=fig_setting[3][1] - 4,
                    bbox=dict(boxstyle="round", facecolor="white", edgecolor="grey", alpha=.75), zorder=6)
            ax.set_aspect("equal")
            ax.set_xlabel(r"$q_x$ (Å$^{-1}$)")
            ax.set_ylabel(r"$q_y$ (Å$^{-1}$)")
            ax.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)
            bar_left = left + panel_width + .20
            cax = fig.add_axes([bar_left / width, bottom / height, .16 / width, panel_height / height])
            colorbar = fig.colorbar(image, cax=cax, extend="both", extendfrac=.04, ticks=ticks)
            colorbar.set_label(f"{title} ({unit})")
            style_heatmap_colorbar(colorbar, ticks)


def plot_hse_bands(title, matters_list=None, eigen_range=None, legend_loc=True):
    # Help information
    help_info = """
    Usage: plot_hse_bands
        arg[0]: title;
        arg[1]: matters list, [label, 4.4_bandstructure_hse/<structure>, PBE line-mode directory or None,
                HSE colour family, PBE colour family];
        arg[2]: energy window relative to each VBM, e.g. (-8, 8);
        arg[3]: legend (True/False);
    HSE06 bands (no SOC) from the zero-weight chunks and PBE bands, each aligned at its own VBM: PBE from the given
    line-mode run (e.g. 4.0_bandstructure/<structure>) or, with None, from the PBE step of the HSE06 jobs.
    The shaded band marks the HSE06 gap; HSE06 band edges are taken over all chunks and the scf mesh.
    """
    if title in ["help", "Help"]:
        print(help_info)
        return

    label, directory, *optional = matters_list
    pbe_dir = optional[0] if len(optional) > 0 else None
    colors = {"HSE06": color_sampling(optional[1] if len(optional) > 1 and optional[1] else "Orange")[1],
              "PBE": color_sampling(optional[2] if len(optional) > 2 and optional[2] else "Blue")[1]}
    bands = extract_hse_bands(directory)
    occupied = bands["occupied_bands"]
    curves = {"HSE06": bands["HSE06"]}
    if pbe_dir is None:
        curves["PBE"] = bands["PBE"]
        positions = [position for position, _ in bands["HSE06"]["ticks"]]
        labels = [{"Gamma": r"$\Gamma$"}.get(name, name) for _, name in bands["HSE06"]["ticks"]]
    else:
        from vmatplot.bandstructure import extract_kpath, extract_eigenvalues_bands_nonpolarized, kpoints_path_lists
        energies = np.array(extract_eigenvalues_bands_nonpolarized(pbe_dir)).T
        valence, conduction = energies[:, occupied - 1], energies[:, occupied]
        curves["PBE"] = {"kpath": np.array(extract_kpath(pbe_dir)), "energies": energies,
                         "vbm": (valence.max(), None), "cbm": (conduction.min(), None)}
        positions, labels = kpoints_path_lists(pbe_dir)

    # Figure settings
    fig_setting = canvas_setting()
    params = fig_setting[2]; plt.rcParams.update(params)
    plt.figure(figsize=fig_setting[0], dpi=fig_setting[1])
    plt.tick_params(direction="in", which="both", top=True, right=True, bottom=True, left=True)

    # Colors calling
    vbm_color = color_sampling("Violet")
    annotate_color = color_sampling("Grey")

    # Data calling and plotting
    width = plt.rcParams["lines.linewidth"]
    for functional in ("PBE", "HSE06"):
        result = curves[functional]
        energies = result["energies"] - result["vbm"][0]
        gap = result["cbm"][0] - result["vbm"][0]
        for band in range(energies.shape[1]):
            plt.plot(result["kpath"], energies[:, band], c=colors[functional], lw=width, linestyle="solid",
                     zorder=4 if functional == "HSE06" else 3, label=f"{functional}: indirect gap {gap:.2f} eV" if band == 0 else None)
    hse = curves["HSE06"]
    plt.axhspan(0, hse["cbm"][0] - hse["vbm"][0], color=colors["HSE06"], alpha=0.12, lw=0, zorder=1,
                label="HSE06 indirect gap (shaded)")
    # Each functional is aligned at its own VBM; zero is not the Fermi level.
    plt.axhline(y=0, color=vbm_color[0], lw=width, alpha=0.8, linestyle="--", zorder=2)

    # High symmetry path
    for k_loc in positions[1:-1]:
        plt.axvline(x=k_loc, color=annotate_color[1], lw=width, linestyle="--", alpha=0.8, zorder=1)
    plt.xticks(positions, labels)

    # Title, axes and ranges
    plt.title(f"{title}")
    plt.xlim(positions[0], positions[-1])
    plt.ylim(*(eigen_range if eigen_range is not None else (-8, 8)))
    plt.ylabel(r"$E-E_\mathrm{VBM}$ (eV)")

    # Legend
    if legend_loc:
        plt.legend(loc="lower right", frameon=True, fancybox=True, facecolor="white", framealpha=.75)
    plt.tight_layout()
