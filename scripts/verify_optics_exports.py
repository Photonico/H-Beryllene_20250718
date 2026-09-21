#!/usr/bin/env python3
"""Independently compare archived optical PDF curves with the historical model.

The PDF's own tick strokes and numeric labels calibrate each axis. No dielectric
array is used to infer the coordinate transform. Historical Matplotlib output
simplified curves, so this checks retained, visible sample vertices and reports
clipped samples separately. Raw inputs and archived PDFs are never modified.

Requires numpy, h5py and PyMuPDF. Run from any working directory::

    python scripts/verify_optics_exports.py
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess

import h5py
import numpy as np
import pymupdf


ROOT = Path(__file__).resolve().parents[1]
LEGACY_REVISION = "f863ab02ebf93d5804f5552967ea750571778335"
FIGURES = {
    "dielectric": ("3.1_dielectric.pdf", "fig3.13_dielectric.pdf"),
    "absorption": ("4.1_absorption.pdf", "fig3.14a_absorption.pdf"),
    "refractive": ("4.2_refractive.pdf", "fig3.14b_refractive.pdf"),
    "extinction": ("4.3_extinction.pdf", "fig3.15_extinction.pdf"),
    "reflectivity": ("4.4_reflectivity.pdf", "S3.16_reflectivity.pdf"),
    "energy-loss": ("4.5_energy-loss.pdf", "S3.17_energy-loss.pdf"),
}


def historical_arrays():
    """Reimplement the archived model, without importing the repaired library."""
    arrays, sources = {}, []
    for name, thickness in [("a", 3.96), ("b", 0.0473771641975619 + 3.96)]:
        source = ROOT / f"6.0_dielectric_selection/{name}-Beryllene/vaspout.h5"
        with h5py.File(source, "r") as handle:
            energy = handle["results/linear_response/energies_dielectric_function"][()]
            tensor = handle["results/linear_response/density_density_dielectric_function"][()]
        mask = (energy >= 0) & (energy <= 12)
        energy = energy[mask]
        sources.append({"path": str(source.relative_to(ROOT)),
                        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                        "legacy_effective_thickness_angstrom": thickness})
        for component, index in [("xx", 0), ("zz", 2)]:
            eps1 = 1 + (tensor[index, index, mask, 0] - 1) * 40 / thickness
            eps2 = tensor[index, index, mask, 1] * 40 / thickness
            magnitude = np.sqrt(eps1**2 + eps2**2)
            n = np.sqrt((magnitude + eps1) / 2)
            k = np.sqrt((magnitude - eps1) / 2)
            arrays[name, component] = {
                "energy": energy, "real": eps1, "imag": eps2,
                "absorption": 2 * (energy / 4.135667662e-15) * k / 2.99792458e17,
                "refractive": n, "extinction": k,
                "reflectivity": ((n - 1)**2 + k**2) / ((n + 1)**2 + k**2),
                "energy-loss": eps2 / (eps1**2 + eps2**2),
            }
    return arrays, sources


def numeric_spans(page):
    result = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                label = span["text"].strip().replace("−", "-")
                if re.fullmatch(r"[-+]?\d+(?:\.\d+)?", label):
                    result.append((float(label), pymupdf.Rect(span["bbox"]), label))
    return result


def tick_transform(ticks):
    """Fit only self-consistent PDF tick labels, including vector-minus cases.

    Some legacy PDFs store the minus sign as a drawing, outside text extraction.
    A consensus among at least three tick labels avoids treating -5 as +5.
    All accepted/rejected ticks are recorded, making the calibration reviewable.
    """
    if len(ticks) < 3:
        raise ValueError(f"Too few independently readable ticks: {ticks}")
    values, coordinates = np.array([(t[0], t[1]) for t in ticks]).T
    best = None
    for i, j in itertools.combinations(range(len(ticks)), 2):
        if values[i] == values[j]:
            continue
        slope = (coordinates[j] - coordinates[i]) / (values[j] - values[i])
        intercept = coordinates[i] - slope * values[i]
        good = np.abs(coordinates - (slope * values + intercept)) < 0.02
        score = int(good.sum())
        if best is None or score > best[0]:
            best = score, good
    if best is None or best[0] < 3:
        raise ValueError(f"Inconsistent PDF ticks: {ticks}")
    good = best[1]
    slope, intercept = np.polyfit(values[good], coordinates[good], 1)
    return float(slope), float(intercept), {
        "accepted_ticks": [[t[0], t[1]] for t, use in zip(ticks, good) if use],
        "rejected_text_only_ticks": [[t[0], t[1]] for t, use in zip(ticks, good) if not use],
        "max_tick_residual_pt": float(np.max(np.abs(coordinates[good] - (slope * values[good] + intercept)))),
    }


def calibrate_axis(rect, drawings, spans):
    xticks, yticks = [], []
    for drawing in drawings:
        items = drawing["items"]
        if len(items) != 1 or items[0][0] != "l" or drawing.get("color") != (0, 0, 0):
            continue
        _, p, q = items[0]
        if abs(p.x - q.x) < 0.001 and 1 < abs(p.y - q.y) < 8 and abs(max(p.y, q.y) - rect.y1) < 0.01:
            candidates = [(v, box, label) for v, box, label in spans
                          if abs((box.x0 + box.x1) / 2 - p.x) < 3 and abs(box.y0 - rect.y1) < 12]
            if len(candidates) == 1:
                xticks.append((candidates[0][0], float(p.x)))
        if abs(p.y - q.y) < 0.001 and 1 < abs(p.x - q.x) < 8 and abs(min(p.x, q.x) - rect.x0) < 0.01:
            candidates = [(v, box, label) for v, box, label in spans
                          if 0 < rect.x0 - box.x1 < 20 and abs((box.y0 + box.y1) / 2 - p.y) < 10]
            if len(candidates) == 1:
                yticks.append((candidates[0][0], float(p.y)))
    sx, tx, xinfo = tick_transform(xticks)
    sy, ty, yinfo = tick_transform(yticks)
    return sx, tx, sy, ty, {"x": xinfo, "y": yinfo,
                            "page_transform": [sx, tx, sy, ty],
                            "axis_rect_pt": list(rect)}


def curve_points(drawing):
    points = []
    for item in drawing["items"]:
        if item[0] != "l":
            raise ValueError("Expected a polyline for an optical curve")
        if not points or tuple(points[-1]) != tuple(item[1]):
            points.append(item[1])
        points.append(item[2])
    return np.array([(point.x, point.y) for point in points])


def verify_legacy_pdf(pdf_bytes, property_name, arrays, tolerance_pt):
    document = pymupdf.open(stream=pdf_bytes, filetype="pdf")
    if len(document) != 1:
        raise ValueError("Expected one-page legacy figure")
    page = document[0]
    drawings = page.get_drawings()
    spans = numeric_spans(page)
    axes = [drawing for drawing in drawings if drawing["type"] == "f"
            and len(drawing["items"]) == 1 and drawing["items"][0][0] == "re"
            and drawing["rect"].width > 100 and drawing["rect"].height > 100
            and drawing["rect"] != page.rect]
    if len(axes) != (4 if property_name == "dielectric" else 2):
        raise ValueError(f"Unexpected axis count {len(axes)}")
    result = []
    for index, axis in enumerate(axes):
        rect = axis["rect"]
        stop = axes[index + 1]["seqno"] if index + 1 < len(axes) else float("inf")
        axis_drawings = [d for d in drawings if axis["seqno"] < d["seqno"] < stop]
        sx, tx, sy, ty, calibration = calibrate_axis(rect, axis_drawings, spans)
        component = "xx" if index % 2 == 0 else "zz"
        quantity = ("real" if index < 2 else "imag") if property_name == "dielectric" else property_name
        curves = [d for d in axis_drawings if d["type"] == "s" and len(d["items"]) > 20 and d["rect"].width > 100]
        if len(curves) != 2:
            raise ValueError(f"Expected two curves, got {len(curves)}")
        curve_results = []
        for drawing in curves:
            # The historical palettes uniquely identify alpha (blue) and beta (orange).
            material = "a" if drawing["color"][2] > drawing["color"][0] else "b"
            source = arrays[material, component]
            expected = np.column_stack([sx * source["energy"] + tx, sy * source[quantity] + ty])
            points = curve_points(drawing)
            inside = ((points[:, 0] >= rect.x0 - 1e-4) & (points[:, 0] <= rect.x1 + 1e-4)
                      & (points[:, 1] >= rect.y0 - 1e-4) & (points[:, 1] <= rect.y1 + 1e-4))
            visible = points[inside]
            # Match the retained vertex by its x sample, never by fitting coordinates.
            sample_indices = np.abs(expected[:, 0, None] - visible[None, :, 0]).argmin(axis=0)
            distances = np.linalg.norm(expected[sample_indices] - visible, axis=1)
            maximum = float(distances.max())
            clipped = ((expected[:, 1] < rect.y0) | (expected[:, 1] > rect.y1))
            curve_results.append({"material": material, "component": component, "quantity": quantity,
                                  "retained_visible_vertices": len(visible),
                                  "retained_outside_vertices": int((~inside).sum()),
                                  "raw_samples_clipped_by_legacy_axis": int(clipped.sum()),
                                  "max_page_residual_pt": maximum,
                                  "passed": maximum <= tolerance_pt})
        result.append({"component": component, "quantity": quantity,
                       "calibration": calibration, "curves": curve_results})
    return result


def run_legacy(reference, output, tolerance_pt):
    arrays, sources = historical_arrays()
    revision = subprocess.check_output(["git", "rev-parse", reference], cwd=ROOT, text=True).strip()
    files = []
    for quantity, filenames in FIGURES.items():
        for folder, filename in zip(("figures_for_thesis", "figures_for_publication"), filenames):
            path = f"{folder}/{filename}"
            content = subprocess.check_output(["git", "show", f"{revision}:{path}"], cwd=ROOT)
            axes = verify_legacy_pdf(content, quantity, arrays, tolerance_pt)
            files.append({"path": path, "sha256": hashlib.sha256(content).hexdigest(), "axes": axes})
    curves = [curve for file in files for axis in file["axes"] for curve in axis["curves"]]
    audit = {"mode": "legacy_pdf_vs_historical_density_density_model", "git_revision": revision,
             "calibration_method": "PDF axis rectangles, tick strokes and numeric tick labels only; no fitting to data arrays",
             "model": "40 Angstrom supercell; alpha thickness 3.96; beta thickness 4.0073771641975619; absorption uses E/h without 2*pi",
             "scope": "Retained visible PDF vertices only, since legacy curves are simplified and some loss peaks clipped. This is not a physical DFT error estimate.",
             "tolerance_pt": tolerance_pt, "input_sources": sources, "files": files,
             "summary": {"file_count": len(files), "curve_count": len(curves),
                         "retained_visible_vertices": sum(c["retained_visible_vertices"] for c in curves),
                         "max_page_residual_pt": max(c["max_page_residual_pt"] for c in curves),
                         "passed": all(c["passed"] for c in curves)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit["summary"], indent=2))
    return audit["summary"]["passed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-ref", default=LEGACY_REVISION,
                        help="Git revision containing the original PDFs")
    parser.add_argument("--legacy-output", type=Path, default=ROOT / ".agent/optics_audit/legacy_pdf_verification.json")
    parser.add_argument("--tolerance-pt", type=float, default=0.01)
    args = parser.parse_args()
    passed = run_legacy(args.legacy_ref, args.legacy_output, args.tolerance_pt)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
