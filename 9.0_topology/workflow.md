# SOC topology calculation workflow

Calculation setup, submission records, monitoring commands, and interpretation notes are maintained here. The [results notebook](../9.0_topology.ipynb) reads existing output files only. Paths and Python snippets below are relative to the repository root unless explicitly absolute. Run the Python snippets in order from the repository root.

For tool implementation details, see [tools/README.md](tools/README.md).


# Beryllene SOC topology: submitted calculation campaign

Submitted UTC: **20260922T055824Z**. Repository: `/cmt2/lniu6305/H-Beryllene_20250718`.

All paths use the current `9.0_topology` directory. `20260922_draft` is preserved as collaborator-provided completed legacy calculations; its manuscript topology claims are not treated as validated for the current reference geometries.

## Structures and analysis methods

| Structure | Directory under 9.0_topology | Neutral spinor filling | Geometric inversion | Method |
|---|---|---:|---|---|
| alpha | `a-Beryllene` | 2 | Yes | Fu–Kane parity |
| beta | `b-Beryllene` | 4 | Yes | Fu–Kane parity |
| st | `c-Beryllene_trilayer` | 6 | Yes | Fu–Kane parity |
| alpha_2h | `a-Beryllene_hh` | 4 | Yes | Fu–Kane parity |
| beta_1h | `b-Beryllene_b` | 5 | No | Direct SOC-DFT Wilson loop / WCC, lowest 4 and 6 bands |
| beta_2h | `b-Beryllene_bb` | 6 | Yes | Fu–Kane parity |

Geometry is taken from the existing frozen `6.0_dielectric_selection` CONTCARs. The existing `4.0_bandstructure` k paths and PAW potentials are reused. Inversion was independently checked on each complete Be/H input structure at 1e-5 angstrom tolerance, allowing shifted inversion centers. Electronic time reversal still requires checking the converged state.

## Submitted PBS jobs

| Task | PBS job ID | afterok dependencies | Script |
|---|---|---|---|
| tools | `41847.headnode` | None | `9.0_topology/tools/build_tools.pbs` |
| alpha | `41848.headnode` | None | `9.0_topology/a-Beryllene/run.pbs` |
| beta | `41849.headnode` | None | `9.0_topology/b-Beryllene/run.pbs` |
| st | `41850.headnode` | None | `9.0_topology/c-Beryllene_trilayer/run.pbs` |
| alpha_2h | `41851.headnode` | None | `9.0_topology/a-Beryllene_hh/run.pbs` |
| beta_1h | `41852.headnode` | None | `9.0_topology/b-Beryllene_b/run.pbs` |
| beta_2h | `41853.headnode` | None | `9.0_topology/b-Beryllene_bb/run.pbs` |
| parity | `41854.headnode` | 41848.headnode, 41849.headnode, 41850.headnode, 41851.headnode, 41853.headnode | `9.0_topology/tools/parity.pbs` |
| wcc | `41855.headnode` | 41847.headnode, 41852.headnode | `9.0_topology/tools/wcc.pbs` |

Submission means accepted by PBS, not calculation completion. The authoritative new records are `9.0_topology/tools/ACTIVE_SUBMISSION.json` and `9.0_topology/tools/submitted_jobs_20260922T055824Z.jsonl`. Legacy submission records are not used.

## Calculation sequence and parameters

Each material pipeline runs SOC SCF on a Gamma-centered 105 x 105 x 1 grid; explicit Gamma/X/Y/M TRIM; the original band path; then two local direct-gap refinement stages around sampled minima. Parameters: PBE, ENCUT 600 eV, Accurate precision, EDIFF 1e-8 eV, SIGMA 0.01 eV, 28 spinor bands, zero initial moments, existing slab/dipole setup, fixed structures. Site VASP: `/cmt2/ocon2505/VASP/vasp.6.5.0/bin/vasp_ncl`. A scalar-relativistic CHGCAR supplies the initial density; its WAVECAR is not reused for SOC.

For the lowest N=2M spinor bands, every stage reports min_k(E_(N+1)-E_N), min_k E_(N+1)-max_k E_N, and the minimizing k point. A discrete grid and local refinement provide numerical evidence, not a rigorous all-k proof. Very small gaps or unstable minima remain unresolved pending convergence study. Negative indirect gap must not be described as an insulating Fermi-level gap.

Fu–Kane analysis uses IrRep 2.1.3 with SOC spinors and four distinct TRIM, retaining the inversion operation, degenerate blocks/Kramers-pair parities, delta_i and conditional Z2. The two square-lattice X-type TRIM are retained separately. The parity job does not depend on Wannier90 compilation.

1H-beta has five neutral electrons per primitive cell: lowest-4/lowest-6 manifolds are conditional isolated subspaces, not the ordinary insulating occupied manifold. WCC runs only when all sampled stages keep direct separation above 0.1 meV. It uses direct SOC-DFT overlaps through the VASP–Wannier90 interface and Z2Pack; it does not claim an already-fitted Wannier Hamiltonian. Loop-position, gap and movement convergence are retained.

## Software and next review

Existing Python environment: IrRep 2.1.3, Z2Pack 2.2.1, spglib 2.5.0, NumPy 1.26.4 and SciPy 1.13.1. The tools job builds downloaded Wannier90 3.1.0 and a private VASP 6.5.0 SOC executable linked to it in `/cmt2/lniu6305/H-Beryllene_20250718/9.0_topology/packages/soc_topology`. WCC waits for both that build and the 1H-beta pipeline. No edge-state jobs are included at this stage.

After the jobs finish, use the read-only Python snippets below and the root results notebook to collect completion, gap, parity and WCC reports. Inspect failed prerequisites, magnetism, local-gap convergence and matching band indices before making final topology claims.


## Queue and pipeline checks

```python
from pathlib import Path
import json, subprocess
repo = Path.cwd()
if not (repo / "9.0_topology").is_dir():
    repo = Path("/cmt2/lniu6305/H-Beryllene_20250718")
topology = repo / "9.0_topology"
submission = json.loads((topology / "tools/ACTIVE_SUBMISSION.json").read_text())
print(subprocess.run(["qstat", "-x", *[j["id"] for j in submission["jobs"]]], capture_output=True, text=True).stdout)
for entry in json.loads((topology / "manifest.json").read_text())["structures"]:
    folder = topology / entry["directory"]
    print(entry["id"], "pipeline complete:", (folder / "PIPELINE_COMPLETE").exists())
```

## Reading the generated reports

```python
for entry in json.loads((topology / "manifest.json").read_text())["structures"]:
    folder = topology / entry["directory"]
    print("\n", entry["id"])
    for stage in ("scf", "trim", "bands", "gap_refine_1", "gap_refine_2"):
        report = folder / stage / "topology_validation.json"
        if report.exists():
            data = json.loads(report.read_text())
            print(stage, "valid:", data.get("validation_passed"))
            for gap in data.get("eigenvalues", {}).get("manifolds", []):
                print("  N=", gap["N"], "direct eV=", gap["direct_gap_ev"], "indirect eV=", gap["indirect_gap_ev"], "k=", gap["direct_gap_k"])
    for pattern in ("trim/parity-analysis-*/parity_summary.json", "trim/parity-analysis-*/PARITY_NOT_VALIDATED.json", "wcc_direct/status_*.json"):
        for report in sorted(folder.glob(pattern)):
            print(report.relative_to(repo))
            print(json.dumps(json.loads(report.read_text()), indent=2))
```

## Previous notebook notes (preserved)


# Topology

## references

Coexistence of superconductivity and topological aspects in beryllenes

## 2026-09-23: continue SOC topology calculations

The three beta SCF jobs hit the 200 GB memory limit before electronic iterations. Their SCF KPAR is now 1 instead of 7, retaining 42 MPI ranks, NCORE=6, 200 GB, and all physical settings. The preparation template also uses KPAR=1. Successful alpha, st, and alpha_2h calculations are retained.

K-point review: all six SCF meshes remain Gamma-centered 105x105x1 with zero shift. Separate reciprocal-coordinate TRIM calculations include (0,0,0), (1/2,0,0), (0,1/2,0), and (1/2,1/2,0). Band paths retain the original reference coordinates; for distorted beta cells the nominal K label is a reference-path label, not an exact high-symmetry assignment. Local refinement and WCC reciprocal coordinates are unchanged. Finite sampling does not certify a gap everywhere.

WCC screening now includes bands and verifies completed stage validation and input fingerprints. Reported subspace Z2 values remain conditional on electronic time-reversal symmetry; they are not automatically Fermi-level insulating invariants.

Submitted jobs: beta=42143.headnode, beta_1h=42144.headnode, beta_2h=42145.headnode, parity=42146.headnode, wcc=42147.headnode. Parity waits for beta and beta_2h; WCC waits for beta_1h. No backup or recovery script is required; changes are managed in Git.
