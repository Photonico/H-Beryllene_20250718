#!/usr/bin/env python3
"""Validate a completed static SOC VASP run and summarize sampled band gaps.

Only the supplied run directory receives new, exclusively-created reports.
No VASP calculation, input modification, or package installation is performed.
Python 3.9+; standard library only, with optional installed spglib support.
Example: python -B check_topology_outputs.py NEW_RUN_DIR --expected-nelect 4
"""

import argparse
import hashlib
import importlib
import itertools
import json
import math
from pathlib import Path
import re
import sys
from datetime import datetime, timezone

sys.dont_write_bytecode = True


class InvalidRun(ValueError):
    pass


def number(value):
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise InvalidRun("Nonfinite numeric value: " + value)
    return result


def integer(value, label):
    result = number(str(value))
    if result != int(result):
        raise InvalidRun(label + " must be an integer")
    return int(result)


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def read_outcar(path):
    result = {"timing_footer": False, "electronic_convergence": False}
    tags = ("LSORBIT", "LNONCOLLINEAR", "NELECT", "NKPTS", "NBANDS",
            "NELM", "NSW", "NIONS", "ENCUT", "ISYM", "ICHARG", "ISPIN")
    patterns = {tag: re.compile(r"\b" + tag + r"\s*=\s*([^\s;]+)") for tag in tags}
    last_iteration_line = last_convergence_line = footer_line = -1
    with path.open(errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            if "vasp." in line and "version" not in result:
                result["version"] = line.strip()
            for tag, pattern in patterns.items():
                match = pattern.search(line)
                if match:
                    result[tag] = match.group(1)
            match = re.search(r"Iteration\s+(\d+)\s*\(\s*(\d+)\s*\)", line)
            if match:
                result["final_ionic_iteration"] = int(match.group(1))
                result["final_electronic_iteration"] = int(match.group(2))
                last_iteration_line = line_number
            if "aborting loop because EDIFF is reached" in line:
                last_convergence_line = line_number
            if "General timing and accounting" in line:
                footer_line = line_number
            if "number of electron" in line and "magnetization" in line:
                result["last_magnetization_line"] = line.strip()
                magnetic_fields = line.split("magnetization", 1)[1].split()
                if len(magnetic_fields) in (1, 3):
                    result["last_total_magnetization_muB"] = [number(v) for v in magnetic_fields]
    result["timing_footer"] = footer_line > last_iteration_line >= 0
    result["electronic_convergence"] = (
        footer_line > last_convergence_line > last_iteration_line >= 0
    )
    result["evidence_lines"] = {"last_iteration": last_iteration_line,
                                "ediff_reached": last_convergence_line,
                                "timing_footer": footer_line}
    errors = []
    for tag in ("LSORBIT", "LNONCOLLINEAR"):
        if result.get(tag, "").strip(".").upper() not in ("T", "TRUE"):
            errors.append(tag + " is not confirmed true in OUTCAR")
    if not result["timing_footer"]:
        errors.append("OUTCAR lacks a timing footer after the final iteration")
    if not result["electronic_convergence"]:
        errors.append("Final electronic iteration lacks a subsequent EDIFF convergence marker")
    for tag in ("NELECT", "NKPTS", "NBANDS", "NELM", "NSW"):
        if tag not in result:
            errors.append("Missing OUTCAR parameter " + tag)
    if "NELM" in result and "final_electronic_iteration" in result:
        if result["final_electronic_iteration"] >= integer(result["NELM"], "NELM"):
            errors.append("Final electronic iteration reached/exceeded NELM")
    if "NSW" in result and integer(result["NSW"], "NSW") != 0:
        errors.append("Expected a fixed-geometry run with NSW=0; no ionic convergence is assumed")
    return result, errors


def determinant(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def inverse(a):
    det = determinant(a)
    if abs(det) < 1e-12:
        raise InvalidRun("Singular POSCAR lattice")
    return [[(
        a[(j + 1) % 3][(i + 1) % 3] * a[(j + 2) % 3][(i + 2) % 3]
        - a[(j + 1) % 3][(i + 2) % 3] * a[(j + 2) % 3][(i + 1) % 3]
    ) / det for j in range(3)] for i in range(3)]


def row_product(vector, matrix):
    return [sum(vector[i] * matrix[i][j] for i in range(3)) for j in range(3)]


def read_poscar(path):
    lines = path.read_text().splitlines()
    if len(lines) < 8:
        raise InvalidRun("POSCAR is incomplete")
    scale = [number(v) for v in lines[1].split()]
    lattice = [[number(v) for v in line.split()[:3]] for line in lines[2:5]]
    if any(len(row) != 3 for row in lattice):
        raise InvalidRun("POSCAR lattice must have three components per vector")
    if len(scale) == 1 and scale[0] != 0:
        factor = scale[0]
        if factor < 0:
            if abs(determinant(lattice)) < 1e-12:
                raise InvalidRun("Singular unscaled POSCAR lattice")
            factor = (-factor / abs(determinant(lattice))) ** (1 / 3)
        factors = [factor] * 3
    elif len(scale) == 3 and all(v > 0 for v in scale):
        factors = scale
    else:
        raise InvalidRun("Unsupported/invalid POSCAR scale")
    lattice = [[v * factors[j] for j, v in enumerate(row)] for row in lattice]
    inv = inverse(lattice)
    species = lines[5].split()
    if all(re.fullmatch(r"\d+", s) for s in species):
        raise InvalidRun("VASP4 POSCAR lacks species labels; provide a VASP5-format copy")
    counts = [integer(v, "species count") for v in lines[6].split()]
    if len(species) != len(counts) or any(c <= 0 for c in counts):
        raise InvalidRun("Invalid POSCAR species/counts")
    index = 7
    if lines[index].lstrip().lower().startswith("s"):
        index += 1
    mode = lines[index].lstrip().lower()
    if not mode or mode[0] not in ("d", "c", "k"):
        raise InvalidRun("Unknown POSCAR coordinate convention")
    index += 1
    positions = []
    for line in lines[index:index + sum(counts)]:
        coords = [number(v) for v in line.split()[:3]]
        if len(coords) != 3:
            raise InvalidRun("Incomplete POSCAR coordinate row")
        if mode[0] != "d":
            coords = row_product([coords[j] * factors[j] for j in range(3)], inv)
        positions.append([v % 1 for v in coords])
    if len(positions) != sum(counts):
        raise InvalidRun("POSCAR has fewer positions than declared atoms")
    numbers = [i + 1 for i, count in enumerate(counts) for _ in range(count)]
    return {"lattice": lattice, "fractional_positions": positions,
            "species": species, "counts": counts, "numbers": numbers}


def periodic_distance(diff, lattice, inv):
    initial = [v - round(v) for v in diff]
    radius = math.sqrt(sum(v * v for v in row_product(initial, lattice)))
    bounds = [radius * math.sqrt(sum(inv[i][j] ** 2 for i in range(3)))
              for j in range(3)]
    ranges = [range(math.ceil(-diff[j] - bounds[j] - 1e-10),
                    math.floor(-diff[j] + bounds[j] + 1e-10) + 1) for j in range(3)]
    if math.prod(len(r) for r in ranges) > 100000:
        raise InvalidRun("Extremely skew lattice exceeds safe periodic-image search size")
    best = radius * radius
    for shift in itertools.product(*ranges):
        cart = row_product([diff[j] + shift[j] for j in range(3)], lattice)
        best = min(best, sum(v * v for v in cart))
    return math.sqrt(max(0, best))


def perfect_matching(distances, tolerance):
    matched = {}

    def visit(i, seen):
        for j, distance in enumerate(distances[i]):
            if distance > tolerance or j in seen:
                continue
            seen.add(j)
            if j not in matched or visit(matched[j], seen):
                matched[j] = i
                return True
        return False

    return all(visit(i, set()) for i in range(len(distances)))


def inversion_check(structure, tolerance, use_spglib):
    positions, labels = structure["fractional_positions"], structure["numbers"]
    lattice, inv = structure["lattice"], inverse(structure["lattice"])
    groups = [[i for i, label in enumerate(labels) if label == kind]
              for kind in sorted(set(labels))]
    anchor_group = min(groups, key=len)
    anchor = anchor_group[0]
    best = None
    for partner in anchor_group:
        tau = [(positions[anchor][j] + positions[partner][j]) % 1 for j in range(3)]
        matrices = []
        for group in groups:
            matrices.append([[periodic_distance(
                [tau[k] - positions[i][k] - positions[j][k] for k in range(3)], lattice, inv)
                for j in group] for i in group])
        thresholds = sorted(set(v for matrix in matrices for row in matrix for v in row))
        lo, hi = 0, len(thresholds) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if all(perfect_matching(m, thresholds[mid] + 1e-12) for m in matrices):
                hi = mid
            else:
                lo = mid + 1
        residual = thresholds[lo]
        if best is None or residual < best[0]:
            best = (residual, tau)
    result = {"method": "species-preserving periodic bijection; atom-pinned inversion centers",
              "symprec_angstrom": tolerance, "candidate_max_residual_angstrom": best[0],
              "candidate_translation_fractional": best[1],
              "has_inversion": best[0] <= tolerance,
              "borderline_center_fit": tolerance < best[0] <= 2 * tolerance,
              "note": "Geometry does not establish electronic time-reversal symmetry."}
    if use_spglib:
        try:
            spglib = importlib.import_module("spglib")
            symmetry = spglib.get_symmetry((lattice, positions, labels), symprec=tolerance)
            if symmetry is None:
                raise ValueError("spglib returned no symmetry dataset")
            translations = [list(map(float, tau)) for rot, tau in
                            zip(symmetry["rotations"], symmetry["translations"])
                            if all(int(rot[i][j]) == (-1 if i == j else 0)
                                   for i in range(3) for j in range(3))]
            result["spglib"] = {"version": spglib.__version__,
                                 "inversion_translations_fractional": translations}
            result["has_inversion"] = bool(translations)
        except Exception as error:
            # A broken optional installation must not disable the standard-library check.
            result["spglib_unavailable"] = str(error)
    return result


def eigenval_scan(path, expected_nelect, manifolds, gap_tolerance):
    with path.open() as stream:
        header = [stream.readline() for _ in range(6)]
        if not all(header):
            raise InvalidRun("EIGENVAL header is incomplete")
        if integer(header[0].split()[-1], "EIGENVAL ISPIN") != 1:
            raise InvalidRun("Expected one spinor eigenvalue column, not collinear ISPIN=2")
        nelect, nkpts, nbands = [integer(v, "EIGENVAL header") for v in header[5].split()]
        if nelect != expected_nelect or nkpts <= 0 or nbands <= 0:
            raise InvalidRun("EIGENVAL NELECT differs from reference, or dimensions are invalid")
        selections = sorted(set([expected_nelect] + manifolds))
        if expected_nelect % 2:
            selections = sorted(set(selections + [expected_nelect - 1, expected_nelect + 1]))
        if any(n < 1 or n >= nbands for n in selections):
            raise InvalidRun("Each requested manifold needs 1 <= N < actual NBANDS")
        statistics = {n: {"N": n, "direct_gap_ev": math.inf,
                          "valence_max_ev": -math.inf, "conduction_min_ev": math.inf}
                      for n in selections}
        trim_found = set()
        weight_sum = 0.0
        line_number = 6

        def nonblank():
            nonlocal line_number
            while True:
                line = stream.readline()
                line_number += 1
                if not line:
                    raise InvalidRun("Truncated EIGENVAL near line " + str(line_number))
                if line.strip():
                    return line.split()

        for index in range(nkpts):
            fields = nonblank()
            if len(fields) != 4:
                raise InvalidRun("Malformed EIGENVAL k-point row at line " + str(line_number))
            kpoint = [number(v) for v in fields[:3]]
            weight = number(fields[3])
            if weight < -1e-12:
                raise InvalidRun("Negative k-point weight")
            weight_sum += weight
            reduced = [v % 1 for v in kpoint]
            if min(reduced[2], 1 - reduced[2]) < 1e-6:
                if all(abs(2 * v - round(2 * v)) < 2e-6 for v in reduced[:2]):
                    trim_found.add(tuple(int(round(2 * v)) % 2 for v in reduced[:2]))
            energies = []
            for band in range(1, nbands + 1):
                row = nonblank()
                if len(row) != 3 or integer(row[0], "band index") != band:
                    raise InvalidRun("Expected spinor band/energy/occupation row at line " + str(line_number))
                energy, occupation = number(row[1]), number(row[2])
                if occupation < -1e-6 or occupation > 1 + 1e-6:
                    raise InvalidRun("Spinor occupation outside [0,1]")
                if energies and energy < energies[-1] - 1e-6:
                    raise InvalidRun("EIGENVAL bands are not ordered by energy")
                energies.append(energy)
            for n, stats in statistics.items():
                lower, upper = energies[n - 1], energies[n]
                if upper - lower < stats["direct_gap_ev"]:
                    stats.update(direct_gap_ev=upper - lower, direct_gap_k=kpoint,
                                 direct_gap_k_index_1based=index + 1)
                if lower > stats["valence_max_ev"]:
                    stats.update(valence_max_ev=lower, valence_max_k=kpoint)
                if upper < stats["conduction_min_ev"]:
                    stats.update(conduction_min_ev=upper, conduction_min_k=kpoint)
        if any(line.strip() for line in stream):
            raise InvalidRun("Unexpected trailing nonblank EIGENVAL content")
    for n, stats in statistics.items():
        stats["indirect_gap_ev"] = stats["conduction_min_ev"] - stats["valence_max_ev"]
        stats["time_reversal_compatible_dimension"] = n % 2 == 0
        stats["neutral_filling_split"] = n == expected_nelect
        stats["sampled_direct_gap_above_tolerance"] = stats["direct_gap_ev"] > gap_tolerance
        stats["sampled_indirect_gap_above_tolerance"] = stats["indirect_gap_ev"] > gap_tolerance
    return {"NELECT": nelect, "NKPTS": nkpts, "NBANDS": nbands,
            "weight_sum": weight_sum, "manifolds": list(statistics.values()),
            "trim_present_fractional": [[i / 2, j / 2, 0] for i, j in sorted(trim_found)],
            "all_four_2d_trim_present": len(trim_found) == 4,
            "gap_tolerance_ev": gap_tolerance,
            "scope": "Only supplied k points; no continuous-BZ gap/convergence or Z2 certification.",
            "odd_neutral_filling": bool(nelect % 2)}


def manifest_entry(path, structure_id):
    data = json.loads(path.read_text())
    collection = data.get("structures", data.get("systems", data))
    if structure_id:
        if isinstance(collection, list):
            matches = [item for item in collection if item.get("id", item.get("structure_id", item.get("name"))) == structure_id]
            if len(matches) != 1:
                raise InvalidRun("Manifest structure ID must match exactly one entry")
            entry = matches[0]
        elif (isinstance(collection, dict)
              and collection.get("id", collection.get("structure_id")) == structure_id):
            entry = collection
        elif isinstance(collection, dict) and structure_id in collection:
            entry = collection[structure_id]
        else:
            raise InvalidRun("Cannot locate structure ID in manifest")
    else:
        entry = data
    for key in ("expected_nelect", "NELECT", "nelect"):
        if isinstance(entry, dict) and key in entry:
            return entry, integer(entry[key], key)
    raise InvalidRun("Manifest entry lacks expected_nelect/NELECT/nelect")


def text_summary(report):
    lines = ["PASS" if report["validation_passed"] else "FAIL", "Run: " + report["run_directory"]]
    lines += ["ERROR: " + error for error in report["errors"]]
    if "magnetism_warning" in report:
        lines.append("WARNING: " + report["magnetism_warning"])
    if "eigenvalues" in report:
        eigen = report["eigenvalues"]
        for stats in eigen["manifolds"]:
            lines.append("N={N}: sampled direct gap={direct_gap_ev:.9g} eV; "
                         "sampled indirect gap={indirect_gap_ev:.9g} eV".format(**stats))
        lines.append("All four 2D TRIM explicitly present: " + str(eigen["all_four_2d_trim_present"]))
        if eigen["odd_neutral_filling"]:
            lines.append("Odd neutral filling: no ordinary TR-symmetric primitive-cell band insulator; even-N manifolds are separate candidates.")
        lines.append(eigen["scope"])
    if "inversion" in report:
        lines.append("Structural inversion: " + str(report["inversion"]["has_inversion"]))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--expected-nelect", type=int)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--structure-id")
    parser.add_argument("--manifold", type=int, action="append", default=[])
    parser.add_argument("--symprec", type=float, default=1e-5, help="Structural tolerance in Angstrom")
    parser.add_argument("--gap-tolerance", type=float, default=1e-6, help="Reporting threshold in eV, not a convergence claim")
    parser.add_argument("--no-spglib", action="store_true")
    parser.add_argument("--require-wavecar", action="store_true")
    args = parser.parse_args()
    run = args.run_directory.resolve()
    if not run.is_dir():
        parser.error("run_directory must already exist and be the NEW calculation directory")
    if (not math.isfinite(args.symprec) or not math.isfinite(args.gap_tolerance)
            or args.symprec <= 0 or args.gap_tolerance < 0):
        parser.error("symprec must be positive and gap-tolerance nonnegative")
    targets = [run / "topology_validation.json", run / "topology_validation.txt"]
    if any(path.exists() or path.is_symlink() for path in targets):
        parser.error("Refusing to overwrite an existing validation report")
    report = {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
              "run_directory": str(run), "validation_passed": False, "errors": [],
              "topology_certified": False, "provenance": {}}
    try:
        expected = args.expected_nelect
        if args.manifest:
            entry, manifest_expected = manifest_entry(args.manifest, args.structure_id)
            if expected is not None and expected != manifest_expected:
                raise InvalidRun("CLI and manifest electron counts disagree")
            expected = manifest_expected
            report["manifest"] = {"path": str(args.manifest.resolve()), **fingerprint(args.manifest)}
            report["manifest_entry"] = entry
        if expected is None or expected <= 0:
            raise InvalidRun("Provide a positive expected electron count from the chosen reference")
        report["expected_nelect"] = expected
        for name in ("OUTCAR", "EIGENVAL", "POSCAR", "INCAR", "KPOINTS", "POTCAR"):
            path = run / name
            if not path.is_file() or path.stat().st_size == 0:
                raise InvalidRun("Missing/empty required file: " + name)
            report["provenance"][name] = fingerprint(path)
        report["provenance"]["validator"] = fingerprint(Path(__file__))
        metadata, errors = read_outcar(run / "OUTCAR")
        report["outcar"] = metadata
        report["errors"].extend(errors)
        report["time_reversal_electronic_state_certified"] = False
        if any(abs(v) > 1e-4 for v in metadata.get("last_total_magnetization_muB", [])):
            report["magnetism_warning"] = "Nonzero total magnetization: do not apply the ordinary TR Z2 criterion without checking magnetic symmetry."
        if "NELECT" in metadata and abs(number(metadata["NELECT"]) - expected) > 1e-6:
            report["errors"].append("OUTCAR NELECT differs from the chosen reference")
        if args.require_wavecar:
            wavecar = run / "WAVECAR"
            if not wavecar.is_file() or wavecar.stat().st_size == 0:
                report["errors"].append("Required WAVECAR is absent/empty")
            else:
                report["wavecar_bytes"] = wavecar.stat().st_size
        structure = read_poscar(run / "POSCAR")
        report["structure"] = structure
        entry = report.get("manifest_entry", {})
        if "atom_count" in entry and integer(entry["atom_count"], "manifest atom_count") != sum(structure["counts"]):
            report["errors"].append("POSCAR atom count differs from the manifest")
        if "NIONS" in metadata and integer(metadata["NIONS"], "NIONS") != sum(structure["counts"]):
            report["errors"].append("POSCAR atom count differs from OUTCAR")
        report["inversion"] = inversion_check(structure, args.symprec, not args.no_spglib)
        if isinstance(entry.get("inversion_expected"), bool):
            if entry["inversion_expected"] != report["inversion"]["has_inversion"]:
                report["errors"].append("POSCAR inversion differs from the manifest expectation")
        if report["errors"]:
            report["analysis_skipped"] = "Do not analyze eigenvalues from a run that failed completion/SOC checks."
        else:
            eigen = eigenval_scan(run / "EIGENVAL", expected, args.manifold, args.gap_tolerance)
            report["eigenvalues"] = eigen
            for tag in ("NKPTS", "NBANDS"):
                if integer(metadata[tag], tag) != eigen[tag]:
                    report["errors"].append("OUTCAR/EIGENVAL " + tag + " mismatch")
        report["validation_passed"] = not report["errors"]
    except (InvalidRun, OSError, ValueError, IndexError, KeyError, TypeError) as error:
        report["errors"].append(type(error).__name__ + ": " + str(error))
    summary = text_summary(report)
    # Exclusive creation also protects against a race after the initial existence check.
    try:
        with targets[0].open("x") as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")
        with targets[1].open("x") as stream:
            stream.write(summary)
    except OSError as error:
        print("Could not exclusively create validation reports: " + str(error), file=sys.stderr)
        return 2
    print(summary, end="")
    return 0 if report["validation_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
