#!/usr/bin/env python3
"""Fail-closed four-TRIM parity screening with exactly IrRep 2.1.3.

Examples: python -B parity_analysis.py CAMPAIGN_ROOT
          python -B parity_analysis.py NEW_SOC_RUN --expected-nelect 4
Creates a unique analysis directory inside NEW_SOC_RUN. Original VASP files
are only read. A parity product alone does not establish a globally isolated
manifold, electronic time reversal, or a quantum spin Hall insulator.
"""

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import uuid

sys.dont_write_bytecode = True
from check_topology_outputs import (InvalidRun, fingerprint, integer,
                                    inversion_check, inverse, periodic_distance,
                                    read_outcar, read_poscar)

TRIMS = {"Gamma": [0.0, 0.0, 0.0], "X": [0.5, 0.0, 0.0],
         "Y": [0.0, 0.5, 0.0], "M": [0.5, 0.5, 0.0]}


def require(condition, message):
    if not condition:
        raise InvalidRun(message)


def complex_tree(real, imaginary):
    if isinstance(real, list):
        require(isinstance(imaginary, list) and len(real) == len(imaginary),
                "Malformed complex array")
        return [complex_tree(a, b) for a, b in zip(real, imaginary)]
    return complex(real, imaginary)


def monty_object(value):
    """Decode the numeric Monty JSON forms written by pinned IrRep."""
    if value.get("@module") == "numpy" and value.get("@class") == "array":
        data = value["data"]
        if "complex" in value.get("dtype", ""):
            require(len(data) == 2, "Unsupported complex ndarray encoding")
            return complex_tree(data[0], data[1])
        return data
    if value.get("@module") == "builtins" and value.get("@class") == "complex":
        return complex(value["real"], value["imag"])
    return value


def load_irrep(path):
    with path.open() as stream:
        return json.load(stream, object_hook=monty_object)


def trim_name(k):
    require(len(k) == 3, "Malformed k point")
    names = [name for name, target in TRIMS.items()
             if max(abs(float(a) - b - round(float(a) - b))
                    for a, b in zip(k, target)) < 1e-7]
    require(len(names) == 1, "Unexpected WAVECAR k point: " + str(k))
    return names[0]


def points_by_name(subspace, require_all=True):
    result = {}
    for point in subspace["k points"]:
        name = trim_name(point["k"])
        require(name not in result, "Duplicate TRIM " + name)
        result[name] = point
    if require_all:
        require(set(result) == set(TRIMS), "Need exactly all four 2D TRIM")
    return result


def select_inversion(data, structure, geometry, tolerance):
    require(data["spacegroup"]["spinor"] is True, "IrRep did not use spinors")
    require(data["spacegroup"]["magnetic"] is False,
            "Unexpected magnetic-space-group mode; this runner assumes a nonmagnetic model")
    candidates = []
    lattice = structure["lattice"]
    for index, symmetry in data["spacegroup"]["symmetries"].items():
        rotation = symmetry["rotation matrix"]
        if not all(abs(rotation[i][j] - (-1 if i == j else 0)) < 1e-8
                   for i in range(3) for j in range(3)):
            continue
        tau = symmetry["translation"]
        difference = [float(tau[i]) - geometry["candidate_translation_fractional"][i]
                      for i in range(3)]
        if periodic_distance(difference, lattice, inverse(lattice)) <= tolerance:
            candidates.append((int(index), [float(t) for t in tau]))
    require(candidates, "No IrRep inversion matches the independently checked POSCAR center")
    return min(candidates)


def extract_trace_parities(data, inversion_index, occupied, trace_tolerance):
    require(data.get("separated by symmetry") is False,
            "Expected an unseparated IrRep character report")
    entries = data["characters and irreps"]
    require(len(entries) == 1, "Unexpected original IrRep subspaces")
    subspace = entries[0]["subspace"]
    points = points_by_name(subspace)
    results = {}
    for name, point in points.items():
        energies = point["energies_raw"]
        dimensions = [integer(v, "degenerate-block dimension") for v in point["dimensions"]]
        require(len(energies) == occupied and sum(dimensions) == occupied,
                name + ": incomplete occupied manifold")
        require(all(math.isfinite(float(e)) for e in energies), "Nonfinite energies")
        syms = [int(i) for i in point["symmetries"]]
        require(inversion_index in syms, name + ": inversion missing from little group")
        column = syms.index(inversion_index)
        characters = point["characters"]  # DFT cell, including actual inversion translation.
        require(len(characters) == len(dimensions), "Character/block dimension mismatch")
        require(len(point["energies_mean"]) == len(dimensions), "Energy/block mismatch")
        pairs, blocks, band_start = [], [], 1
        for block_index, (dimension, row) in enumerate(zip(dimensions, characters)):
            require(dimension > 0 and dimension % 2 == 0,
                    name + ": energy block does not contain complete Kramers pairs")
            require(len(row) == len(syms), "Character/symmetry dimension mismatch")
            trace = complex(row[column])
            require(math.isfinite(trace.real) and math.isfinite(trace.imag)
                    and abs(trace.imag) <= trace_tolerance, "Inversion trace is not real")
            odd_states_float = (dimension - trace.real) / 2
            odd_states = round(odd_states_float)
            require(abs(odd_states_float - odd_states) <= trace_tolerance
                    and 0 <= odd_states <= dimension and odd_states % 2 == 0,
                    name + ": inversion trace cannot represent complete even/odd Kramers pairs")
            signs = [-1] * (odd_states // 2) + [1] * ((dimension - odd_states) // 2)
            pairs.extend(signs)
            blocks.append({"bands_one_based": [band_start, band_start + dimension - 1],
                           "energy_mean_ev": float(point["energies_mean"][block_index]),
                           "dimension": dimension, "inversion_trace_real": trace.real,
                           "inversion_trace_imag": trace.imag,
                           "kramers_pair_parities_unordered_within_block": signs})
            band_start += dimension
        require(len(pairs) == occupied // 2, "Incomplete Kramers-pair count")
        results[name] = {"k_fractional": TRIMS[name], "blocks": blocks,
                         "pair_parities_by_energy_block": pairs,
                         "odd_kramers_pairs": pairs.count(-1), "delta": math.prod(pairs)}
    return results, subspace


def check_separated(data, inversion_index, expected, occupied, tolerance):
    require(data.get("separated by symmetry") is True, "Missing symmetry separation")
    require(data["separating symmetries"] == [inversion_index], "Wrong separating symmetry")
    counts = {name: {-1: 0, 1: 0} for name in TRIMS}
    for entry in data["characters and irreps"]:
        eigenvalues = entry["symmetry eigenvalues"]
        require(len(eigenvalues) == 1, "Unexpected number of separating eigenvalues")
        value = complex(eigenvalues[0])
        sign = 1 if value.real >= 0 else -1
        require(abs(value - sign) < tolerance, "Separated inversion eigenvalue is not +/-1")
        for name, point in points_by_name(entry["subspace"], require_all=False).items():
            count = sum(integer(d, "separated dimension") for d in point["dimensions"])
            require(count == len(point["energies_raw"]) and count % 2 == 0,
                    name + ": separated parity subspace has incomplete Kramers pairs")
            counts[name][sign] += count
    for name, count in counts.items():
        require(sum(count.values()) == occupied, name + ": separation lost occupied states")
        require(count[-1] == 2 * expected[name]["odd_kramers_pairs"],
                name + ": parity trace and explicit symmetry separation disagree")
    return {name: {"odd_states": count[-1], "even_states": count[1]}
            for name, count in counts.items()}


def run_irrep(directory, run, occupied, degeneracy, separation=None):
    directory.mkdir()
    command = [sys.executable, "-B", "-c", "from irrep.cli import cli; cli()",
               "-code=vasp", "-spinor", "-fPOS=" + str(run / "POSCAR"),
               "-fWAV=" + str(run / "WAVECAR"), "-IBstart=1", "-IBend=" + str(occupied),
               "-degenThresh=" + str(degeneracy), "-json_file=irrep.json", "-v"]
    if separation is not None:
        command.extend(["-isymsep=" + str(separation), "-groupKramers"])
    with (directory / "command.json").open("x") as stream:
        json.dump(command, stream, indent=2)
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    with (directory / "irrep.log").open("x") as log:
        completed = subprocess.run(command, cwd=str(directory), env=environment,
                                   stdout=log, stderr=subprocess.STDOUT, check=False)
    require(completed.returncode == 0, "IrRep failed; see " + str(directory / "irrep.log"))
    log_text = (directory / "irrep.log").read_text(errors="replace")
    require(not re.search(r"WARNING|orthogonality.*>|non.?unitary|non-integer", log_text,
                          re.IGNORECASE), "IrRep numerical warning requires review: " + str(directory))
    return load_irrep(directory / "irrep.json")


def campaign_main(campaign, args):
    manifest_path = campaign / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    entries = manifest["structures"]
    require(isinstance(entries, list) and entries, "Campaign manifest needs a structures list")
    selected = [entry for entry in entries if entry.get("inversion_expected") is True]
    require(selected, "Campaign has no inversion-expected structures")
    ids = [entry["id"] for entry in selected]
    require(len(set(ids)) == len(ids), "Duplicate structure IDs")
    results = []
    for entry in selected:
        sid = entry["id"]
        require(isinstance(sid, str) and re.fullmatch(r"[A-Za-z0-9_-]+", sid),
                "Unsafe/invalid structure ID")
        run = campaign / sid / "trim"
        require(run.resolve().is_relative_to(campaign), "Run directory escapes campaign")
        validation = json.loads((run / "topology_validation.json").read_text())
        require(validation.get("manifest_entry", {}).get("id") == sid,
                "Validator manifest structure ID mismatch")
        require(validation.get("expected_nelect") == entry["expected_nelect"],
                "Validator electron count differs from campaign manifest")
        code = main([str(run), "--expected-nelect", str(entry["expected_nelect"]),
                     "--symprec", str(args.symprec), "--degen-thresh", str(args.degen_thresh),
                     "--trace-tolerance", str(args.trace_tolerance)])
        results.append({"id": sid, "run": str(run), "exit_code": code})
    path = campaign / ("parity_batch_" + uuid.uuid4().hex + ".json")
    with path.open("x") as stream:
        json.dump({"manifest": fingerprint(manifest_path), "results": results,
                   "topology_certified": False}, stream, indent=2)
    print(str(path))
    return 0 if all(item["exit_code"] == 0 for item in results) else 2


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--expected-nelect", type=int)
    parser.add_argument("--symprec", type=float, default=1e-5)
    parser.add_argument("--degen-thresh", type=float, default=1e-5)
    parser.add_argument("--trace-tolerance", type=float, default=1e-3)
    args = parser.parse_args(argv)
    run = args.run_dir.resolve(strict=True)
    output = None
    try:
        if args.expected_nelect is None:
            return campaign_main(run, args)
        require(args.expected_nelect in (2, 4, 6), "Expected project manifolds have 2, 4 or 6 spinor bands")
        require(all(0 < v <= 1e-2 for v in (args.symprec, args.degen_thresh, args.trace_tolerance)),
                "Numerical tolerances must be positive and <=1e-2")
        for name in ("WAVECAR", "POSCAR", "OUTCAR"):
            require((run / name).is_file() and (run / name).stat().st_size > 0,
                    "Missing/nonempty input required: " + name)
        validation = json.loads((run / "topology_validation.json").read_text())
        require(validation.get("validation_passed") is True,
                "Prerequisite topology_validation.json did not pass")
        require(validation.get("expected_nelect") == args.expected_nelect,
                "Prerequisite validator electron count disagrees")
        for name in ("POSCAR", "OUTCAR"):
            require(validation["provenance"][name] == fingerprint(run / name),
                    "Validated file changed since validation: " + name)
        require(importlib.metadata.version("irrep") == "2.1.3", "This runner requires exactly IrRep 2.1.3")
        metadata, errors = read_outcar(run / "OUTCAR")
        require(not errors, "; ".join(errors))
        require(float(metadata["NELECT"]) == args.expected_nelect, "NELECT differs from manifest expectation")
        require(integer(metadata["NKPTS"], "NKPTS") == 4, "Expected exactly four explicit TRIM")
        require(integer(metadata["NBANDS"], "NBANDS") > args.expected_nelect,
                "Need at least one extra band to inspect the manifold boundary")
        structure = read_poscar(run / "POSCAR")
        geometry = inversion_check(structure, args.symprec, False)
        require(geometry["has_inversion"], "POSCAR does not pass structural inversion check")
        output = run / ("parity-analysis-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                        + "-" + uuid.uuid4().hex[:8])
        output.mkdir()
        source_fingerprints = {name: fingerprint(run / name) for name in ("POSCAR", "OUTCAR", "WAVECAR")}
        raw = run_irrep(output / "characters", run, args.expected_nelect, args.degen_thresh)
        operation, tau = select_inversion(raw, structure, geometry, args.symprec)
        parities, subspace = extract_trace_parities(raw, operation, args.expected_nelect,
                                                   args.trace_tolerance)
        gap = float(subspace["Minimal direct gap (eV)"])
        require(math.isfinite(gap) and gap > args.degen_thresh,
                "Selected manifold touches upper bands at a sampled TRIM")
        separated = run_irrep(output / "separated", run, args.expected_nelect,
                              args.degen_thresh, operation)
        counts = check_separated(separated, operation, parities, args.expected_nelect,
                                 args.trace_tolerance)
        require(source_fingerprints == {name: fingerprint(run / name) for name in source_fingerprints},
                "Source files changed during parity analysis")
        product = math.prod(point["delta"] for point in parities.values())
        report = {"status": "parity_screen_passed_global_topology_not_yet_certified",
                  "irrep_version": "2.1.3", "source_run": str(run),
                  "inputs": source_fingerprints, "outcar": metadata,
                  "geometry": geometry, "occupied_spinor_bands": args.expected_nelect,
                  "inversion_operation_index": operation,
                  "inversion_translation_fractional": tau,
                  "inversion_center_fractional": [t / 2 for t in tau],
                  "character_convention": "DFT-cell characters for actual {R=-I|tau}; no reference-cell phase substitution",
                  "trims": parities, "separated_state_counts": counts,
                  "minimum_direct_gap_at_four_TRIM_ev": gap,
                  "delta_product": product, "conditional_fu_kane_nu": 0 if product == 1 else 1,
                  "not_verified": ["Electronic time-reversal symmetry of the converged state",
                                   "Isolation of this fixed-band manifold over the full 2D Brillouin zone",
                                   "Global indirect gap and insulating Fermi level"],
                  "note": "Kramers pairs inside a larger exactly degenerate block have no unique band ordering."}
        with (output / "parity_summary.json").open("x") as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
        print(str(output / "parity_summary.json"))
        return 0
    except Exception as error:
        if output is not None:
            with (output / "PARITY_NOT_VALIDATED.json").open("x") as stream:
                json.dump({"status": "failed_closed", "error": str(error)}, stream, indent=2)
        print("Parity not validated: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
