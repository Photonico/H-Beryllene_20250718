#!/usr/bin/env python3
"""Validate and archive one VASP overlap calculation before Z2Pack reads it."""

import argparse
import json
import math
import re
import shutil
import sys
import traceback
import uuid
from datetime import datetime, timezone
from pathlib import Path


def number(value):
    return float(value.replace("D", "E").replace("d", "e"))


def read_eigenval(path):
    """Read a noncollinear/non-spin-polarized EIGENVAL without dependencies."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 6:
        raise ValueError(f"Truncated EIGENVAL: {path}")
    nelect, nkpts, nbands = map(int, lines[5].split())
    points, energies = [], []
    cursor = 6
    for _ in range(nkpts):
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        fields = lines[cursor].split()
        if len(fields) != 4:
            raise ValueError("Invalid EIGENVAL k-point record")
        points.append([number(x) for x in fields[:3]])
        cursor += 1
        row = []
        for band in range(1, nbands + 1):
            fields = lines[cursor].split()
            if len(fields) != 3 or int(fields[0]) != band:
                raise ValueError("Expected ordered spinor band records in EIGENVAL")
            row.append(number(fields[1]))
            cursor += 1
        if any(b < a - 1e-7 for a, b in zip(row, row[1:])):
            raise ValueError("EIGENVAL band energies are not ordered")
        energies.append(row)
    if not all(math.isfinite(x) for row in points + energies for x in row):
        raise ValueError("Non-finite EIGENVAL data")
    return {"nelect": nelect, "nkpts": nkpts, "nbands": nbands,
            "points": points, "energies": energies}


def gap_summary(eigenval, bands):
    if not 0 < bands < eigenval["nbands"]:
        raise ValueError("Selected manifold lacks a next band")
    values = [row[bands] - row[bands - 1] for row in eigenval["energies"]]
    index = min(range(len(values)), key=values.__getitem__)
    return {
        "bands": bands,
        "minimum_direct_gap_ev": values[index],
        "minimum_direct_gap_kpoint": eigenval["points"][index],
        "indirect_gap_ev": min(row[bands] for row in eigenval["energies"])
        - max(row[bands - 1] for row in eigenval["energies"]),
        "sampled_kpoints": eigenval["nkpts"],
    }


def check_outcar(path, *, fixed_density=False, expected_bands=None):
    text = Path(path).read_text(errors="replace")
    if "General timing and accounting" not in text:
        raise ValueError("VASP did not reach its normal termination summary")
    if "aborting loop because EDIFF is reached" not in text:
        raise ValueError("VASP electronic convergence was not recorded")
    for tag in ("LSORBIT", "LNONCOLLINEAR"):
        if not re.search(rf"\b{tag}\s*=\s*T\b", text):
            raise ValueError(f"OUTCAR does not confirm {tag}=T")
    if fixed_density and not re.search(r"\bICHARG\s*=\s*11\b", text):
        raise ValueError("WCC line did not use fixed SOC charge density")
    if expected_bands is not None:
        found = re.findall(r"\bNBANDS\s*=\s*(\d+)", text)
        if not found or int(found[-1]) != expected_bands:
            raise ValueError(f"VASP actual NBANDS differs from {expected_bands}")


def read_nnkp(path):
    lines = Path(path).read_text().splitlines()

    def block(name):
        start = next(i for i, line in enumerate(lines)
                     if line.strip().lower() == f"begin {name}")
        end = next(i for i in range(start + 1, len(lines))
                   if lines[i].strip().lower() == f"end {name}")
        return [line.strip() for line in lines[start + 1:end] if line.strip()]

    data = block("kpoints")
    nkpts = int(data[0])
    points = [[number(x) for x in line.split()] for line in data[1:]]
    if len(points) != nkpts or any(len(row) != 3 for row in points):
        raise ValueError("Invalid nnkp kpoints block")
    data = block("nnkpts")
    neighbors = int(data[0])
    links = [tuple(map(int, line.split())) for line in data[1:]]
    if len(links) != nkpts * neighbors or any(len(row) != 5 for row in links):
        raise ValueError("Invalid nnkp neighbor block")
    return points, neighbors, links


def equivalent(a, b, tolerance=2e-7):
    return all(abs(x - y - round(x - y)) < tolerance for x, y in zip(a, b))


def read_mmn(path, expected_bands, expected_points):
    links = []
    minimum_norm, maximum_norm = float("inf"), 0.0
    with Path(path).open() as stream:
        stream.readline()
        dimensions = tuple(map(int, stream.readline().split()))
        if dimensions != (expected_bands, expected_points, 1):
            raise ValueError(f"Unexpected mmn dimensions {dimensions}")
        for _ in range(expected_points):
            link = tuple(map(int, stream.readline().split()))
            if len(link) != 5:
                raise ValueError("Malformed mmn neighbor record")
            links.append(link)
            norm = 0.0
            for _ in range(expected_bands ** 2):
                fields = stream.readline().split()
                if len(fields) != 2:
                    raise ValueError("Truncated overlap matrix")
                real, imag = map(number, fields)
                if not math.isfinite(real) or not math.isfinite(imag):
                    raise ValueError("Non-finite overlap matrix")
                norm += real * real + imag * imag
            minimum_norm = min(minimum_norm, norm)
            maximum_norm = max(maximum_norm, norm)
        if any(line.strip() for line in stream):
            raise ValueError("Unexpected extra overlap records")
    if minimum_norm <= 1e-12:
        raise ValueError("Overlap matrix is numerically zero")
    return links, {"minimum_frobenius_norm_squared": minimum_norm,
                   "maximum_frobenius_norm_squared": maximum_norm}


def next_archive(parent):
    if not parent.is_dir():
        raise ValueError("Line archive directory must already exist")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    for _ in range(10):
        target = parent / f"line_{stamp}_{uuid.uuid4().hex[:8]}"
        try:
            target.mkdir()
            return target
        except FileExistsError:
            continue
    raise RuntimeError("Cannot allocate a unique line archive")


def quarantine_mmn():
    if Path("wannier90.mmn").exists():
        Path("wannier90.mmn").rename(f"wannier90.rejected.{uuid.uuid4().hex}.mmn")


def validate(args):
    if args.vasp_exit:
        raise ValueError(f"VASP returned status {args.vasp_exit}")
    requested = json.loads(Path("requested_loop.json").read_text())
    points = requested["points"]
    nkpts = len(points) - 1
    if nkpts < 4 or not equivalent(points[0], points[-1]):
        raise ValueError("Invalid requested closed loop")
    eig = read_eigenval("EIGENVAL")
    if (eig["nelect"], eig["nbands"], eig["nkpts"]) != (5, args.nbands, nkpts):
        raise ValueError("EIGENVAL electron, band or k-point count mismatch")
    check_outcar("OUTCAR", fixed_density=True, expected_bands=args.nbands)
    nnpoints, neighbors, nnlinks = read_nnkp("wannier90.nnkp")
    if len(nnpoints) != nkpts or neighbors != 1:
        raise ValueError("Expected one nearest neighbor per requested point")
    for index in range(nkpts):
        if not equivalent(eig["points"][index], points[index]):
            raise ValueError(f"EIGENVAL point {index + 1} is out of order")
        if not equivalent(nnpoints[index], eig["points"][index]):
            raise ValueError(f"nnkp point {index + 1} differs from VASP")
    for index, link in enumerate(nnlinks):
        source, target, *shift = link
        expected_target = (index + 1) % nkpts + 1
        if (source, target) != (index + 1, expected_target):
            raise ValueError("nnkp links do not follow the requested loop order")
        actual_step = [b + g - a for a, b, g in
                       zip(nnpoints[source - 1], nnpoints[target - 1], shift)]
        intended_step = [b - a for a, b in zip(points[index], points[index + 1])]
        if any(abs(a - b) > 2e-7 for a, b in zip(actual_step, intended_step)):
            raise ValueError("nnkp reciprocal translation follows a wrong periodic image")
    mmnlinks, norms = read_mmn("wannier90.mmn", args.bands, nkpts)
    if mmnlinks != nnlinks:
        raise ValueError("mmn matrix ordering differs from nnkp ordering")
    gaps = gap_summary(eig, args.bands)
    if gaps["minimum_direct_gap_ev"] <= args.threshold:
        raise ValueError(f"Manifold is not isolated on this line: {gaps}")
    return {"status": "validated", "nelect": 5, "vasp_nbands": args.nbands,
            "wcc_count": args.bands, "gap": gaps, "overlap_norms": norms}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--bands", type=int, choices=(4, 6), required=True)
    parser.add_argument("--nbands", type=int, default=28)
    parser.add_argument("--threshold", type=float, default=1e-4)
    parser.add_argument("--vasp-exit", type=int, required=True)
    args = parser.parse_args()
    destination = next_archive(args.archive.resolve())
    failure = None
    try:
        report = validate(args)
    except Exception as exc:
        failure = exc
        report = {"status": "failed", "error": str(exc),
                  "traceback": traceback.format_exc()}
    report["recorded_utc"] = datetime.now(timezone.utc).isoformat()
    report["archive"] = str(destination)
    try:
        # Never copy CHGCAR repeatedly; the immutable input copy is retained.
        for name in (
            "requested_loop.json", "INCAR", "KPOINTS", "OUTCAR", "OSZICAR",
            "EIGENVAL", "vasp.log", "wannier90.win", "wannier90.wout",
            "wannier90.werr", "wannier90.nnkp", "wannier90.mmn",
            "wannier90.eig",
        ):
            source = Path(name)
            if source.is_file():
                shutil.copy2(source, destination / name)
        with (destination / "validation.json").open("x") as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")
    except Exception:
        failure = RuntimeError("Unable to preserve WCC line evidence")
        traceback.print_exc()
    if failure:
        # Prevent an unchecked .mmn from being consumed even if a caller
        # neglects the subprocess exit status. This is private scratch only.
        quarantine_mmn()
        print(f"WCC line rejected: {failure}; evidence: {destination}", file=sys.stderr)
        return 1
    print(json.dumps(report, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        quarantine_mmn()
        traceback.print_exc()
        sys.exit(1)
