#!/usr/bin/env python3
"""Conditionally evaluate isolated lowest-4/6 SOC subspaces of neutral 1H-beta."""

import argparse
import fcntl
import hashlib
import json
import os
import shlex
import shutil
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from wcc_line_check import check_outcar, gap_summary, read_eigenval


OWNER = "beryllene-direct-wcc-v1"
THRESHOLD = 1e-4
NBANDS = 28


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_new(path, data):
    with path.open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def to_native(value):
    if isinstance(value, dict):
        return {str(key): to_native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_native(item) for item in value]
    if hasattr(value, "tolist"):
        return to_native(value.tolist())
    if hasattr(value, "item"):
        return value.item()
    return value


def read_incar(path):
    values = {}
    for line in path.read_text().splitlines():
        line = line.split("!", 1)[0].split("#", 1)[0]
        for field in line.split(";"):
            if "=" in field:
                key, value = field.split("=", 1)
                values[key.strip().upper()] = value.strip()
    return values


def incar_text(source, bands):
    values = read_incar(source)
    # Retain all physical settings from the validated SOC density calculation.
    for tag in ("ISPIN", "NPAR", "KSPACING", "KGAMMA", "LOCPROJ", "WANNIER90_WIN"):
        values.pop(tag, None)
    values.update({
        "SYSTEM": f"1H-beta isolated lowest-{bands} spinor-band WCC",
        "ISTART": "0", "ICHARG": "11", "NSW": "0", "IBRION": "-1",
        "LSORBIT": ".TRUE.", "LNONCOLLINEAR": ".TRUE.", "MAGMOM": "9*0",
        "ISYM": "-1", "SAXIS": "0 0 1", "ALGO": "Normal", "EDIFF": "1E-8",
        "NELM": "200", "NELMDL": "-5", "NBANDS": str(NBANDS),
        "NCORE": "6", "KPAR": "1", "LWAVE": ".FALSE.", "LCHARG": ".FALSE.",
        "LWANNIER90": ".TRUE.", "LWANNIER90_RUN": ".FALSE.",
        "LWRITE_MMN_AMN": ".TRUE.", "LWRITE_UNK": ".FALSE.",
        "NUM_WANN": str(bands), "LOPTICS": ".FALSE.", "LEPSILON": ".FALSE.",
    })
    return "\n".join(f"{key} = {value}" for key, value in values.items()) + "\n"


def assert_converged(result):
    report = to_native(result.convergence_report)
    for scope, names in {"line": ("PosCheck",),
                         "surface": ("GapCheck", "MoveCheck")}.items():
        for name in names:
            entry = report.get(scope, {}).get(name)
            if (not isinstance(entry, dict) or not entry.get("PASSED")
                    or entry.get("FAILED") or entry.get("MISSING")):
                raise RuntimeError(f"WCC convergence failed: {scope}/{name}: {entry}")
    return report


def prepare_manifold(root, structure, checker, bands, owner):
    destination = root / f"lowest_{bands}"
    try:
        destination.mkdir()
    except FileExistsError:
        if not (destination / "input_ready.json").is_file():
            raise RuntimeError(f"Existing incomplete/unowned manifold directory: {destination}")
        saved = json.loads((destination / "input_ready.json").read_text())
        if saved.get("owner") != owner:
            raise RuntimeError("Existing WCC inputs differ; use a fresh run directory")
        for name, digest in saved.get("input_sha256", {}).items():
            source = destination / "input" / name
            if not source.is_file() or sha256(source) != digest:
                raise RuntimeError(f"WCC input copy changed: {name}")
        if len(saved.get("input_sha256", {})) != 6:
            raise RuntimeError("Incomplete WCC input provenance")
        return destination
    inputs = destination / "input"
    inputs.mkdir()
    (destination / "lines").mkdir()
    for name in ("POSCAR", "POTCAR", "CHGCAR"):
        shutil.copy2(structure / "scf" / name, inputs / name)
    shutil.copy2(checker, inputs / "wcc_line_check.py")
    (inputs / "INCAR").write_text(incar_text(structure / "scf" / "INCAR", bands))
    (inputs / "wannier90.win").write_text(
        f"num_wann = {bands}\nnum_bands = {bands}\nspinors = true\n"
        f"exclude_bands = {bands + 1}-{NBANDS}\npostproc_setup = true\n"
        "skip_B1_tests = true\nnum_iter = 0\n"
    )
    write_json_new(destination / "input_ready.json", {
        "owner": owner,
        "input_sha256": {path.name: sha256(path) for path in sorted(inputs.iterdir())},
    })
    return destination


def run_manifold(destination, tools, bands, attempt):
    import z2pack

    if (destination / "success.json").is_file():
        return json.loads((destination / "success.json").read_text())
    inputs = destination / "input"
    checkpoint = destination / "checkpoint.json"
    executable = tools / "vasp" / "bin" / "vasp_ncl"
    python = tools / "venv" / "bin" / "python"
    check = " ".join(shlex.quote(str(x)) for x in (
        python, "wcc_line_check.py", "--archive", destination / "lines",
        "--bands", bands, "--nbands", NBANDS, "--threshold", THRESHOLD,
    ))
    command = (f"mpirun -np 42 {shlex.quote(str(executable))} > vasp.log 2>&1; "
               f"vasp_status=$?; {check} --vasp-exit \"$vasp_status\"")

    def loop_request(points):
        return json.dumps({"points": [[float(x) for x in point] for point in points]},
                          allow_nan=False) + "\n"

    system = z2pack.fp.System(
        input_files=[str(inputs / name) for name in (
            "CHGCAR", "INCAR", "POSCAR", "POTCAR", "wannier90.win", "wcc_line_check.py")],
        kpt_fct=[z2pack.fp.kpoint.vasp, z2pack.fp.kpoint.wannier90_full, loop_request],
        kpt_path=["KPOINTS", "wannier90.win", "requested_loop.json"],
        command=command, executable="/bin/bash",
        build_folder=str(destination / "private_build"),
        mmn_path="wannier90.mmn", num_wcc=bands,
    )
    result = z2pack.surface.run(
        system=system, surface=lambda s, t: [t, s / 2, 0],
        pos_tol=0.005, gap_tol=0.3, move_tol=0.3, num_lines=21,
        min_neighbour_dist=0.001, iterator=range(12, 101, 4),
        save_file=str(checkpoint), serializer="json",
        load=checkpoint.exists(), load_quiet=False,
    )
    report = to_native(result.convergence_report)
    write_json_new(destination / f"convergence_{attempt}.json", report)
    assert_converged(result)
    if any(len(line) != bands for line in result.wcc):
        raise RuntimeError("Unexpected WCC count")
    invariant = int(z2pack.invariant.z2(result, check_kramers_pairs=True))
    summary = {
        "status": "converged_sampled_subspace", "bands": bands, "z2": invariant,
        "interpretation": "Isolated lowest-even-band subspace; not a Fermi-level insulating invariant.",
        "limits": "Finite k-point sampling does not prove a nonzero direct gap everywhere.",
        "gap_threshold_ev": THRESHOLD, "convergence": report,
        "surface": "[t, s/2, 0]", "line_positions": to_native(result.t),
        "wcc": to_native(result.wcc), "z2pack_version": z2pack.__version__,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_new(destination / "success.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure_dir", type=Path)
    parser.add_argument("tools_dir", type=Path)
    parser.add_argument("--nelect", type=int, choices=(5,), default=5)
    args = parser.parse_args()
    structure, tools = args.structure_dir.resolve(), args.tools_dir.resolve()
    python = tools / "venv" / "bin" / "python"
    # Execute the driver in the dedicated, validated tool environment.
    if not python.is_file():
        raise RuntimeError("Missing WCC Python environment")
    if Path(sys.prefix).resolve() != (tools / "venv").resolve():
        os.execv(str(python), [str(python), str(Path(__file__).resolve()),
                             str(structure), str(tools), "--nelect", str(args.nelect)])
    checker = Path(__file__).resolve().with_name("wcc_line_check.py")
    sources = [structure / "scf" / name for name in
               ("POSCAR", "POTCAR", "INCAR", "CHGCAR", "EIGENVAL", "OUTCAR")]
    sources += [structure / "trim" / name for name in ("EIGENVAL", "OUTCAR")]
    sources += [checker, Path(__file__).resolve()]
    for source in sources:
        if not source.is_file() or source.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty input: {source}")
    poscar = (structure / "scf" / "POSCAR").read_text().splitlines()
    if sum(map(int, poscar[6].split())) != 3:
        raise RuntimeError("Expected three atoms in the 1H-beta POSCAR")
    owner = {"owner": OWNER, "structure": str(structure),
             "source_sha256": {str(path): sha256(path) for path in sources}}
    root = structure / "wcc_direct"
    try:
        root.mkdir()
        write_json_new(root / "owner.json", owner)
    except FileExistsError:
        if not (root / "owner.json").is_file():
            raise RuntimeError("Refusing to use an existing unowned WCC directory")
        if json.loads((root / "owner.json").read_text()) != owner:
            raise RuntimeError("Source inputs changed; refusing to reuse WCC results")
    lock = (root / "driver.lock").open("a")
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("Another WCC driver is already running")
    attempt = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    status = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "nelect": 5, "threshold_ev": THRESHOLD, "manifolds": {}}
    try:
        eigenvals = {}
        for stage in ("scf", "trim"):
            check_outcar(structure / stage / "OUTCAR")
            eigenvals[stage] = read_eigenval(structure / stage / "EIGENVAL")
            if eigenvals[stage]["nelect"] != 5:
                raise RuntimeError("This driver only handles neutral five-electron 1H-beta")
        eligible = []
        for bands in (4, 6):
            gaps = {stage: gap_summary(data, bands) for stage, data in eigenvals.items()}
            allowed = all(value["minimum_direct_gap_ev"] > THRESHOLD for value in gaps.values())
            status["manifolds"][str(bands)] = {
                "status": "eligible" if allowed else "scientifically_skipped",
                "sampled_gaps": gaps,
            }
            if allowed:
                eligible.append(bands)
        write_json_new(root / f"screening_{attempt}.json", status)
        if not eligible:
            status["status"] = "scientifically_skipped_no_isolated_sampled_manifold"
            write_json_new(root / f"status_{attempt}.json", status)
            print(json.dumps(status, indent=2, allow_nan=False))
            return 0
        ready = tools / "READY"
        if not ready.is_file() or not ready.read_text().strip():
            raise RuntimeError("Missing nonempty TOOLS_DIR/READY validation flag")
        executable = tools / "vasp" / "bin" / "vasp_ncl"
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise RuntimeError("Missing executable VASP/Wannier90 SOC build")
        tool_owner = dict(owner, tools_dir=str(tools),
                          tools_ready_sha256=sha256(ready),
                          vasp_sha256=sha256(executable))
        failures = []
        for bands in eligible:
            try:
                destination = prepare_manifold(root, structure, checker, bands, tool_owner)
                result = run_manifold(destination, tools, bands, attempt)
                status["manifolds"][str(bands)].update(result)
            except Exception as exc:
                failures.append(bands)
                status["manifolds"][str(bands)].update(
                    status="failed", error=str(exc), traceback=traceback.format_exc())
                print(f"Lowest-{bands} WCC failed: {exc}", file=sys.stderr)
        status["status"] = "failed" if failures else "completed_conditional_screening"
        write_json_new(root / f"status_{attempt}.json", status)
        print(json.dumps(to_native(status), indent=2, allow_nan=False))
        return 1 if failures else 0
    except Exception as exc:
        status.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        write_json_new(root / f"failure_{attempt}.json", status)
        raise
    finally:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()


if __name__ == "__main__":
    sys.exit(main())
