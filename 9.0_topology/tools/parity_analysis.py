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


def read_wavecar(path):
    """Plane-wave coefficients of a noncollinear VASP WAVECAR (all k, all bands)."""
    import numpy as np
    with path.open("rb") as stream:
        reclen, nspin, rtag = (int(v) for v in np.fromfile(stream, np.float64, 3))
        require(nspin == 1 and rtag in (45200, 45210), "Unsupported WAVECAR header")
        dtype = np.complex64 if rtag == 45200 else np.complex128
        stream.seek(reclen)
        header = np.fromfile(stream, np.float64, 12)
        nk, nb, encut = int(header[0]), int(header[1]), float(header[2])
        lattice = header[3:12].reshape(3, 3)
        points, record = [], 2
        for _ in range(nk):
            stream.seek(record * reclen)
            info = np.fromfile(stream, np.float64, 4 + 3 * nb)
            npw, k = int(info[0]), info[1:4]
            coefficients = np.empty((nb, npw), dtype)
            for band in range(nb):
                stream.seek((record + 1 + band) * reclen)
                coefficients[band] = np.fromfile(stream, dtype, npw)
            points.append((k, info[4:].reshape(nb, 3)[:, 0], coefficients))
            record += nb + 1
    return rtag, encut, lattice, points


def wavecar_gvectors(lattice, k, encut, coefficients_per_band):
    """G vectors in WAVECAR order (x fastest); the count must match the file."""
    import numpy as np
    reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
    kmax2 = encut / (13.605826 * 0.529177249 ** 2)  # VASP HSQDTM = RYTOEV*AUTOA^2
    bound = np.ceil(np.sqrt(kmax2) * np.linalg.norm(lattice, axis=1) / (2 * np.pi)).astype(int) + 1
    axes = [np.concatenate([np.arange(0, b + 1), np.arange(-b, 0)]) for b in bound]
    gz, gy, gx = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
    g = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    kg = (g + k) @ reciprocal
    g = g[np.einsum("ij,ij->i", kg, kg) < kmax2]
    require(2 * len(g) == coefficients_per_band,
            "WAVECAR plane-wave count does not match the reconstructed G sphere")
    return g


def independent_wavefunction_check(run, occupied, tau, degeneracy, tolerance,
                                   residual_tolerance, expected):
    """Double-precision inversion representation from the WAVECAR itself.

    PAW pseudo-wavefunctions are orthonormal in the S metric, not the plain
    plane-wave metric of IrRep's orthogonality message, and IrRep normalises
    the single-precision WAVECAR coefficients in single precision. S commutes
    with inversion, so genuine parity mixing appears as an inversion residual,
    an opposite-parity plain overlap, or a non-unitary inversion matrix on the
    Lowdin-orthonormalised N-band subspace; all are required below tolerance.
    Same-parity plain overlaps (the PAW metric) are only reported.
    """
    import numpy as np
    rtag, encut, lattice, points = read_wavecar(run / "WAVECAR")
    tau = np.array(tau, float)
    results = {}
    for k, energies, raw in points:
        name = trim_name(k)
        require(name not in results, "Duplicate TRIM " + name)
        g = wavecar_gvectors(lattice, k, encut, raw.shape[1])
        ng = len(g)
        # Bands 1..N plus the next complete degenerate block(s) up to N+2.
        top = occupied + 2
        while top < len(energies) and energies[top] - energies[top - 1] <= degeneracy:
            top += 1
        coefficients = raw[:top].astype(np.complex128)
        raw_norm = np.einsum("ij,ij->i", coefficients.conj(), coefficients).real
        coefficients /= np.sqrt(raw_norm)[:, None]
        # {-1|tau}: coefficient at G' = -G - 2k is c(G) exp(2 pi i (k+G).tau).
        shift = np.rint(2 * k).astype(int)
        index = {tuple(v): i for i, v in enumerate(g)}
        try:
            target = np.array([index[tuple(-v - shift)] for v in g])
        except KeyError:
            raise InvalidRun(name + ": G sphere is not closed under inversion")
        phase = np.exp(2j * np.pi * ((g + k) @ tau))
        transformed = np.zeros_like(coefficients)
        for spin in (0, 1):
            transformed[:, spin * ng + target] = coefficients[:, spin * ng:(spin + 1) * ng] * phase
        overlap = coefficients.conj() @ coefficients.T
        values, vectors = np.linalg.eigh(overlap)
        lowdin = vectors @ np.diag(values ** -0.5) @ vectors.conj().T
        inversion_all = lowdin.conj().T @ (coefficients.conj() @ transformed.T) @ lowdin
        values, vectors = np.linalg.eigh(overlap[:occupied, :occupied])
        lowdin = vectors @ np.diag(values ** -0.5) @ vectors.conj().T
        inversion = lowdin.conj().T @ (coefficients[:occupied].conj() @ transformed[:occupied].T) @ lowdin
        blocks, start = [], 0
        for stop in range(1, top + 1):
            if stop == top or energies[stop] - energies[stop - 1] > degeneracy:
                blocks.append((start, stop))
                start = stop
        require(any(b[1] == occupied for b in blocks), name + ": degenerate block crosses band N")
        parity = [0] * top
        odd_states, eigenvalue_error = 0, 0.0
        off_block = inversion.copy()
        for b1, b2 in blocks:
            matrix = inversion[b1:b2, b1:b2] if b2 <= occupied else inversion_all[b1:b2, b1:b2]
            eigenvalues = np.linalg.eigvals(matrix)
            signs = np.where(eigenvalues.real >= 0, 1, -1)
            if abs(signs.sum()) == b2 - b1:
                parity[b1:b2] = [int(signs[0])] * (b2 - b1)
            if b2 <= occupied:
                off_block[b1:b2, b1:b2] = 0
                eigenvalue_error = max(eigenvalue_error, float(np.abs(eigenvalues - signs).max()))
                odd_states += int((signs < 0).sum())
        singular = np.linalg.svd(inversion, compute_uv=False)
        same = opposite = 0.0
        for m in range(top):
            for n in range(top):
                if m != n and parity[m] and parity[n]:
                    if parity[m] == parity[n]:
                        same = max(same, abs(overlap[m, n]))
                    else:
                        opposite = max(opposite, abs(overlap[m, n]))
        residual = [float(np.linalg.norm(transformed[n] - parity[n] * coefficients[n]))
                    for n in range(top) if parity[n]]
        above = [b for b in blocks if b[0] == occupied][0]
        record = {
            "wavecar_single_precision": rtag == 45200,
            "bands_checked_one_based": [1, top],
            "raw_pseudo_norm_range": [float(raw_norm[:occupied].min()), float(raw_norm[:occupied].max())],
            "plain_metric_same_parity_max_overlap": float(same),
            "plain_metric_opposite_parity_max_overlap": float(opposite),
            "inversion_residual_max": max(residual),
            "bands_without_definite_parity_one_based": [n + 1 for n in range(top) if not parity[n]],
            "lowdin_inversion_singular_value_max_error": float(np.abs(singular - 1).max()),
            "lowdin_inversion_off_block_max": float(np.abs(off_block).max()),
            "lowdin_inversion_eigenvalue_max_error": eigenvalue_error,
            "odd_states": odd_states,
            "next_block_above_N_bands_one_based": [above[0] + 1, above[1]],
            "next_block_above_N_gap_ev": float(energies[above[0]] - energies[occupied - 1]),
            "next_block_above_N_parity": parity[above[0]] or None,
        }
        for key in ("plain_metric_opposite_parity_max_overlap",
                    "lowdin_inversion_singular_value_max_error",
                    "lowdin_inversion_off_block_max",
                    "lowdin_inversion_eigenvalue_max_error"):
            require(record[key] <= tolerance, name + ": " + key + " = %.3g exceeds %.1g" % (record[key], tolerance))
        require(record["inversion_residual_max"] <= residual_tolerance,
                name + ": inversion residual %.3g exceeds %.1g" % (record["inversion_residual_max"], residual_tolerance))
        require(all(parity[:occupied]), name + ": a band-N degenerate block has no definite parity")
        require(odd_states == 2 * expected[name]["odd_kramers_pairs"],
                name + ": WAVECAR inversion matrix and IrRep trace parity disagree")
        results[name] = record
    require(set(results) == set(TRIMS), "WAVECAR must contain exactly the four 2D TRIM")
    return results


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
    # IrRep 2.1.3 prints these two informational warnings for the deliberately
    # chosen DFT-cell, characters-only mode. They are not numerical failures.
    numerical_log = "\n".join(line for line in log_text.splitlines()
                              if not line.startswith((
                                  "Warning: transformation to the convenctional unit cell",
                                  "Warning: -kpnames not specified. Only traces of")))
    # Kpoint.Separate prints max|<psi|psi>-1| of single-precision PAW pseudo-
    # wavefunctions in the plain metric. It is recorded here and replaced by
    # independent_wavefunction_check, which the caller must run.
    orthogonality = [float(value) for value in re.findall(
        r"^orthogonality \(largest of diag\. <psi_nk\|psi_mk>\):\s*(\S+) > 1e-5\s*$",
        numerical_log, re.MULTILINE)]
    numerical_log = re.sub(r"^orthogonality \(largest of diag\. <psi_nk\|psi_mk>\):\s*\S+ > 1e-5\s*$",
                           "", numerical_log, flags=re.MULTILINE)
    require(not re.search(r"WARNING|orthogonal|non.?unitary|non-integer", numerical_log,
                          re.IGNORECASE), "IrRep numerical warning requires review: " + str(directory))
    return load_irrep(directory / "irrep.json"), orthogonality


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
        run = campaign / entry.get("directory", sid) / "trim"
        require(run.resolve().is_relative_to(campaign), "Run directory escapes campaign")
        validation = json.loads((run / "topology_validation.json").read_text())
        require(validation.get("manifest_entry", {}).get("id") == sid,
                "Validator manifest structure ID mismatch")
        require(validation.get("expected_nelect") == entry["expected_nelect"],
                "Validator electron count differs from campaign manifest")
        code = main([str(run), "--expected-nelect", str(entry["expected_nelect"]),
                     "--symprec", str(args.symprec), "--degen-thresh", str(args.degen_thresh),
                     "--trace-tolerance", str(args.trace_tolerance),
                     "--wavefunction-tolerance", str(args.wavefunction_tolerance),
                     "--residual-tolerance", str(args.residual_tolerance)])
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
    parser.add_argument("--wavefunction-tolerance", type=float, default=1e-6)
    parser.add_argument("--residual-tolerance", type=float, default=1e-5)
    args = parser.parse_args(argv)
    run = args.run_dir.resolve(strict=True)
    output = None
    try:
        if args.expected_nelect is None:
            return campaign_main(run, args)
        require(args.expected_nelect in (2, 4, 6), "Expected project manifolds have 2, 4 or 6 spinor bands")
        require(all(0 < v <= 1e-2 for v in (args.symprec, args.degen_thresh, args.trace_tolerance)),
                "Numerical tolerances must be positive and <=1e-2")
        require(0 < args.wavefunction_tolerance <= 1e-5 and 0 < args.residual_tolerance <= 1e-5,
                "Wavefunction tolerances must be positive and no looser than IrRep's 1e-5")
        for name in ("WAVECAR", "POSCAR", "OUTCAR"):
            require((run / name).is_file() and (run / name).stat().st_size > 0,
                    "Missing/nonempty input required: " + name)
        validation = json.loads((run / "topology_validation.json").read_text())
        require(validation.get("validation_passed") is True,
                "Prerequisite topology_validation.json did not pass")
        require(validation.get("expected_nelect") == args.expected_nelect,
                "Prerequisite validator electron count disagrees")
        for name in ("POSCAR", "OUTCAR", "EIGENVAL"):
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
        raw, _ = run_irrep(output / "characters", run, args.expected_nelect, args.degen_thresh)
        operation, tau = select_inversion(raw, structure, geometry, args.symprec)
        parities, subspace = extract_trace_parities(raw, operation, args.expected_nelect,
                                                   args.trace_tolerance)
        # IrRep measures against the mean of the final degenerate block.
        # Use the actual boundary E_(N+1)-E_N from the validated eigenvalues.
        eigenvalues = validation["eigenvalues"]
        require(eigenvalues.get("all_four_2d_trim_present") is True,
                "Validated eigenvalues do not contain all four 2D TRIM")
        boundaries = [item for item in eigenvalues["manifolds"]
                      if item["N"] == args.expected_nelect]
        require(len(boundaries) == 1, "Validated target band manifold is missing")
        gap = float(boundaries[0]["direct_gap_ev"])
        irrep_gap = float(subspace["Minimal direct gap (eV)"])
        require(math.isfinite(gap) and gap > args.degen_thresh,
                "Selected manifold touches upper bands at a sampled TRIM")
        separated, orthogonality = run_irrep(output / "separated", run, args.expected_nelect,
                                             args.degen_thresh, operation)
        counts = check_separated(separated, operation, parities, args.expected_nelect,
                                 args.trace_tolerance)
        wavefunctions = independent_wavefunction_check(run, args.expected_nelect, tau,
                                                       args.degen_thresh, args.wavefunction_tolerance,
                                                       args.residual_tolerance, parities)
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
                  "irrep_plain_metric_orthogonality_messages": orthogonality,
                  "wavefunction_tolerance": args.wavefunction_tolerance,
                  "inversion_residual_tolerance": args.residual_tolerance,
                  "wavefunction_validation": wavefunctions,
                  "minimum_direct_gap_at_four_TRIM_ev": gap,
                  "irrep_gap_to_degenerate_block_mean_ev": irrep_gap,
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
