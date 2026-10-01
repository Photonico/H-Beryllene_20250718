#!/usr/bin/env python3
"""Spinless inversion parities at the four 2D TRIM from a scalar (no-SOC) WAVECAR.

Example: python -B hse_parity.py ../../4.4_bandstructure_hse/b-Beryllene_bb/scf --bands 3 \
                 --reference ../b-Beryllene_bb/trim
Used for the HSE06 check in 4.4_bandstructure_hse: the WAVECAR holds the HSE06 IBZ of a
Gamma-centred mesh, EIGENVAL.pbe the PBE step in the same cell. Without SOC every spinless
band is one Kramers pair, so the lowest n = N/2 bands give the Fu-Kane delta of the lowest-N
spinor subspace. TRIM missing from the IBZ are mapped from an equivalent TRIM with the space-
group operation and its lattice-translation phase. The same float64 checks as
parity_analysis.independent_wavefunction_check are required; the result is written to
<run>/hse_parity.json. Like the SOC screen, this is not a topology certificate.
"""

import argparse
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from check_topology_outputs import InvalidRun, fingerprint, read_poscar
from parity_analysis import TRIMS, read_wavecar, require, wavecar_gvectors


def trim_of(k):
    names = [name for name, target in TRIMS.items()
             if max(abs(a - b - round(a - b)) for a, b in zip(k, target)) < 1e-7]
    return names[0] if names else None


def read_eigenval(path):
    """k points and band energies of a non-spin-polarised EIGENVAL."""
    lines = path.read_text().splitlines()
    nk, nb = (int(v) for v in lines[5].split()[1:3])
    points, index = [], 7
    for _ in range(nk):
        k = [float(v) for v in lines[index].split()[:3]]
        energies = [float(lines[index + 1 + b].split()[1]) for b in range(nb)]
        points.append((k, energies))
        index += nb + 2
    return points


def inversion_operation(structure, reference_tau, symprec):
    import numpy as np
    import spglib
    cell = (structure["lattice"], structure["fractional_positions"], structure["numbers"])
    dataset = spglib.get_symmetry(cell, symprec=symprec)
    operations = list(zip(dataset["rotations"], dataset["translations"]))
    candidates = [np.array(t, float) for r, t in operations if (r == -np.eye(3, dtype=int)).all()]
    require(candidates, "No inversion operation in the POSCAR")
    if reference_tau is not None:
        # Same in-plane inversion centre as the SOC parity run; z is irrelevant for k_z = 0.
        def offset(t):
            d = np.array(t[:2]) - np.array(reference_tau[:2])
            return np.abs(d - np.rint(d)).max()
        candidates.sort(key=offset)
        require(offset(candidates[0]) < 1e-5, "No inversion with the in-plane centre of the reference run")
    return candidates[0], operations


def scalar_parities(k, energies, raw, lattice, encut, bands, tau, degeneracy, tolerance,
                    residual_tolerance):
    """float64 inversion representation of bands 1..n+2 at one TRIM (one spinor component)."""
    import numpy as np
    g = wavecar_gvectors(lattice, k, encut, raw.shape[1], components=1)
    top = bands + 2
    while top < len(energies) and energies[top] - energies[top - 1] <= degeneracy:
        top += 1
    coefficients = raw[:top].astype(np.complex128)
    raw_norm = np.einsum("ij,ij->i", coefficients.conj(), coefficients).real
    coefficients /= np.sqrt(raw_norm)[:, None]
    shift = np.rint(2 * np.array(k)).astype(int)
    index = {tuple(v): i for i, v in enumerate(g)}
    try:
        target = np.array([index[tuple(-v - shift)] for v in g])
    except KeyError:
        raise InvalidRun(str(k) + ": G sphere is not closed under inversion")
    transformed = np.zeros_like(coefficients)
    transformed[:, target] = coefficients * np.exp(2j * np.pi * ((g + k) @ tau))
    overlap = coefficients.conj() @ coefficients.T
    values, vectors = np.linalg.eigh(overlap)
    lowdin = vectors @ np.diag(values ** -0.5) @ vectors.conj().T
    inversion = lowdin.conj().T @ (coefficients.conj() @ transformed.T) @ lowdin
    blocks, start = [], 0
    for stop in range(1, top + 1):
        if stop == top or energies[stop] - energies[stop - 1] > degeneracy:
            blocks.append((start, stop))
            start = stop
    # Without SOC band n may sit inside an orbitally degenerate block that only SOC splits. SOC
    # commutes with inversion, so a block of uniform parity still fixes the parities up to band n.
    limit = [b2 for b1, b2 in blocks if b1 < bands <= b2][0]
    parity, eigenvalue_error = [0] * top, 0.0
    off_block = inversion.copy()
    for b1, b2 in blocks:
        eigenvalues = np.linalg.eigvals(inversion[b1:b2, b1:b2])
        signs = np.where(eigenvalues.real >= 0, 1, -1)
        if abs(signs.sum()) == b2 - b1:
            parity[b1:b2] = [int(signs[0])] * (b2 - b1)
        off_block[b1:b2, b1:b2] = 0
        if b2 <= limit:
            eigenvalue_error = max(eigenvalue_error, float(np.abs(eigenvalues - signs).max()))
    opposite = max([abs(overlap[m, n]) for m in range(top) for n in range(top)
                    if m != n and parity[m] * parity[n] == -1] or [0.0])
    residual = [float(np.linalg.norm(transformed[n] - parity[n] * coefficients[n]))
                for n in range(top) if parity[n]]
    above = [b for b in blocks if b[0] == limit][0]
    record = {"bands_checked_one_based": [1, top],
              "plain_metric_opposite_parity_max_overlap": float(opposite),
              "inversion_residual_max": max(residual),
              "lowdin_inversion_singular_value_max_error":
                  float(np.abs(np.linalg.svd(inversion[:limit, :limit], compute_uv=False) - 1).max()),
              "lowdin_inversion_off_block_max": float(np.abs(off_block[:limit, :limit]).max()),
              "lowdin_inversion_eigenvalue_max_error": eigenvalue_error}
    for key in ("plain_metric_opposite_parity_max_overlap", "lowdin_inversion_singular_value_max_error",
                "lowdin_inversion_off_block_max", "lowdin_inversion_eigenvalue_max_error"):
        require(record[key] <= tolerance, str(k) + ": " + key + " = %.3g exceeds %.1g" % (record[key], tolerance))
    require(record["inversion_residual_max"] <= residual_tolerance,
            str(k) + ": inversion residual %.3g exceeds %.1g" % (record["inversion_residual_max"], residual_tolerance))
    require(all(parity[:limit]), str(k) + ": a degenerate block up to band n has no definite parity")
    record.update({"parities": parity[:bands], "boundary_inside_degenerate_block": limit > bands,
                   "next_block_bands_one_based": [above[0] + 1, above[1]],
                   "next_block_parity": parity[above[0]] or None,
                   "gap_above_n_ev": float(energies[bands] - energies[bands - 1]),
                   "gap_above_degenerate_block_ev": float(energies[limit] - energies[limit - 1])})
    return record


def map_trim(name, available, operations, tau):
    """Parity factor of an equivalent TRIM: I' = g^-1 I g = {1|L} I with L = R^-1(tau - 2t) - tau."""
    import numpy as np
    for source in available:
        k = np.array(TRIMS[source])
        for rotation, translation in operations:
            inverse = np.linalg.inv(rotation)
            if trim_of(k @ inverse) != name:
                continue
            lattice_vector = inverse @ (tau - 2 * translation) - tau
            require(np.abs(lattice_vector - np.rint(lattice_vector)).max() < 1e-5,
                    "Conjugated inversion is not a lattice translate")
            factor = int(np.rint(np.cos(2 * np.pi * k @ np.rint(lattice_vector))))
            return source, factor
    raise InvalidRun(name + " is neither in the IBZ nor equivalent to an IBZ TRIM")


def main(argv=None):
    import numpy as np
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--bands", type=int, required=True, help="spinless bands n = N/2")
    parser.add_argument("--reference", type=Path, help="SOC trim run with a parity-analysis-*/parity_summary.json")
    parser.add_argument("--symprec", type=float, default=1e-5)
    parser.add_argument("--degen-thresh", type=float, default=1e-5)
    parser.add_argument("--wavefunction-tolerance", type=float, default=1e-6)
    parser.add_argument("--residual-tolerance", type=float, default=1e-5)
    args = parser.parse_args(argv)
    run = args.run_dir.resolve(strict=True)
    reference = None
    if args.reference is not None:
        summaries = sorted(args.reference.resolve(strict=True).glob("parity-analysis-*/parity_summary.json"))
        require(summaries, "Reference run has no parity_summary.json")
        reference = json.loads(summaries[-1].read_text())
        require(reference["occupied_spinor_bands"] == 2 * args.bands, "Reference N differs from 2n")
    structure = read_poscar(run / "POSCAR")
    tau, operations = inversion_operation(structure, reference and reference["inversion_translation_fractional"],
                                          args.symprec)
    rtag, encut, lattice, points = read_wavecar(run / "WAVECAR")
    pbe = {trim_of(k): energies for k, energies in read_eigenval(run / "EIGENVAL.pbe") if trim_of(k)}
    trims = {}
    for k, energies, raw in points:
        name = trim_of(k)
        if name is None:
            continue
        require(name not in trims, "Duplicate TRIM " + name)
        record = scalar_parities(np.array(k), energies, raw, lattice, encut, args.bands, tau,
                                 args.degen_thresh, args.wavefunction_tolerance, args.residual_tolerance)
        record.update({"source": "IBZ", "hse_energies_ev": [float(e) for e in energies[:args.bands + 2]],
                       "pbe_energies_ev": pbe[name][:args.bands + 2],
                       "pbe_gap_above_n_ev": pbe[name][args.bands] - pbe[name][args.bands - 1]})
        trims[name] = record
    require("Gamma" in trims, "Gamma missing from the WAVECAR")
    for name in TRIMS:
        if name not in trims:
            source, factor = map_trim(name, [n for n in TRIMS if n in trims], operations, tau)
            record = dict(trims[source], source="mapped from " + source, lattice_phase=factor)
            record["parities"] = [factor * p for p in record["parities"]]
            if record["next_block_parity"]:
                record["next_block_parity"] *= factor
            trims[name] = record
    for record in trims.values():
        record["odd_bands"] = sum(p < 0 for p in record["parities"])
        record["delta"] = math.prod(record["parities"])
    product = math.prod(record["delta"] for record in trims.values())
    report = {"run": str(run), "inputs": {name: fingerprint(run / name) for name in ("POSCAR", "WAVECAR", "EIGENVAL.pbe")},
              "wavecar_single_precision": rtag == 45200, "spinless_bands": args.bands,
              "inversion_translation_fractional": [float(t) for t in tau],
              "trims": {name: trims[name] for name in TRIMS},
              "delta_product": product, "conditional_fu_kane_nu": 0 if product == 1 else 1,
              "topology_certified": False}
    if reference is not None:
        odd = {name: point["odd_kramers_pairs"] for name, point in reference["trims"].items()}
        report["reference_pair_parities"] = {name: point["pair_parities_by_energy_block"]
                                             for name, point in reference["trims"].items()}
        report["reference"] = {"run": reference["source_run"], "conditional_fu_kane_nu": reference["conditional_fu_kane_nu"],
                               "odd_kramers_pairs": odd,
                               "same_odd_counts": all(odd[name] == trims[name]["odd_bands"] for name in TRIMS)}
    with (run / "hse_parity.json").open("w") as stream:
        json.dump(report, stream, indent=1)
    print(json.dumps({"nu": report["conditional_fu_kane_nu"],
                      "odd": {name: trims[name]["odd_bands"] for name in TRIMS},
                      "reference": report.get("reference", {}).get("odd_kramers_pairs"),
                      "gap_hse": {name: round(trims[name]["gap_above_n_ev"], 4) for name in TRIMS},
                      "gap_pbe": {name: round(trims[name]["pbe_gap_above_n_ev"], 4) for name in TRIMS},
                      "next_parity": {name: trims[name]["next_block_parity"] for name in TRIMS}}))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except InvalidRun as error:
        print("FAILED: " + str(error), file=sys.stderr)
        sys.exit(2)
