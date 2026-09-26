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
    for pattern in ("trim/parity-analysis-*/parity_summary.json", "trim/parity-analysis-*/PARITY_NOT_VALIDATED.json",
                    "wcc_direct_v3/status_*.json", "tr_evidence.json", "gap_refine_3_summary.json",
                    "spin_screen/spin_screen_summary.json"):
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

## 2026-09-26: parity diagnosis, time-reversal evidence, isolation and spin screen

**WCC (1H-β).** The final WCC job 42175 finished with exit 0 after three fixes: VASP writes PEAD-type overlaps only with NCORE=1; the k-point connectivity check must compare Wannier90 coordinates modulo reciprocal lattice vectors (Wannier90 folds them into the first zone); Z2Pack checkpoints need `serializer="auto"`. Results are in `b-Beryllene_b/wcc_direct_v3/`: the lowest-4 and lowest-6 subspaces converge with conditional Z2 = 0 and 0; status `completed_conditional_screening`, `topology_certified: false`.

**Parity job 42146 (exit 2) — diagnosis.** Only α passed; β, ST, 2H-α and 2H-β failed on IrRep's line `orthogonality (largest of diag. <psi_nk|psi_mk>): X > 1e-5`, printed by `Kpoint.Separate` (a check the IrRep source marks "Rm once tests are fixed"). The logged values are 1.1e-5–1.2e-5 (β), 1.3e-5–6.5e-3 (ST), 3.6e-5–4.4e-5 (2H-α) and 2.9e-3–1.55e-2 (2H-β). Two sources were identified; neither is poor wavefunction quality:

1. *Single-precision normalisation.* WAVECAR stores complex64 coefficients (RTAG 45200) and IrRep normalises and forms ⟨ψ|ψ⟩ in complex64 over ≈1.1×10⁴ coefficients sorted by |k+G|. Its diagonal error (4.4e-6 to 4.5e-5; α passed only because its maximum was 9.5e-6) vanishes (≤1e-14) when the raw coefficients are normalised in float64.
2. *PAW metric.* Pseudo-wavefunctions are orthonormal in the S metric S = 1 + Σ|p⟩q⟨p|, not the plain plane-wave metric. The remaining off-diagonal overlaps (up to 2.6e-2 over bands 1..N+2) occur only between same-parity bands at different energies; opposite-parity overlaps are ≤4.0e-8. An independent reconstruction of S from the POTCAR projectors and augmentation charges makes all 28 bands S-orthonormal to ≤1.3e-7 (diagonal) and ≤6.0e-8 (off-diagonal). IrRep's own `symm_matrix` uses the Gram right-inverse and is therefore correct in the non-orthogonal basis.

TRIM NSCF runs converged in 7–8 Davidson steps (final dE ≤ 6.7e-9 eV); Kramers pairs of bands 1..N+2 are split by ≤3.8e-6 eV; the TRIM coordinates are the four 2D TRIM of spglib-primitive cells; the lowest N bands equal NELECT (Be 2s² and H 1s valence only); IrRep uses the geometric inversion {−1|(0,0,0.2)}.

**Change to `tools/parity_analysis.py` (no threshold relaxed).** The regex no longer fails on the plain-metric orthogonality line; the values are recorded in `irrep_plain_metric_orthogonality_messages`, and every other IrRep warning still fails closed. `independent_wavefunction_check` rebuilds the inversion matrix from the WAVECAR in float64 (G sphere with VASP's HSQDTM = RYTOEV·AUTOA²) and requires, at every TRIM: Löwdin-orthonormalised inversion matrix on bands 1..N unitary, block-diagonal and with eigenvalues ±1 to 1e-6; opposite-parity plain overlap over bands 1..N+2 ≤1e-6; per-band inversion residual ‖Iψ−pψ‖ ≤1e-5 (IrRep's own scale); and odd-state counts identical to IrRep's trace parities. The rerun (2026-09-26, no VASP) passed for all five centrosymmetric structures: `parity_batch_57b2275a8b154208a7148a579825f89d.json`. Worst values: unitarity 2.5e-11, opposite-parity overlap 4.0e-8, residual 6.4e-6. Conditional ν: α 1, β 0, ST 1, 2H-α 0, 2H-β 0.

**Electronic time reversal (existing outputs, `tools/tr_evidence.py` → `<structure>/tr_evidence.json`).** The converged SOC states are TR-invariant within numerical precision: max|m(r)| ≤ 2.6e-7 μB/Å³ (≤3.4e-7 of max ρ), ∫|m| ≤ 8.9e-7 μB, TR-odd one-centre Re(m) ≤ 3e-7 against the TR-even spin–orbit Im(m) of 7–9e-5, Kramers splitting of bands 1..N+2 at the TRIM ≤ 3.8e-6 eV. For 1H-β (C1, SCF on the full unsymmetrised 105×105 grid) E_n(k) = E_n(−k) holds to 4.0e-7 eV for bands 1–8 over 5512 ±k pairs. Every SCF started from m = 0, and for the five centrosymmetric structures ISYM=2 removed inversion-odd m by construction, so these data do not test stability against magnetic order.

**Spin-polarised stability screen (`tools/spin_screen.py`, `tools/spin_screen.pbs`).** Collinear ISPIN=2 SCFs (vasp_std, no SOC) on the same geometry, 105×105 mesh, cutoff, smearing and dipole setup, seeded with ferromagnetic moments (Be 1 μB, H 0.5 μB), a layer-alternating Be pattern when there are two or more Be, and for ST an additional inversion-odd (+, 0, −) pattern. A seed "collapses" when the total and every LORBIT=11 site moment fall below 1e-3 μB. Supercell magnetic orders are not sampled.

**Isolation of the lowest-N subspace.** The 105×105 SCF grid was unfolded with the symmetry VASP used (all IBZ weights reproduced) and every direct-gap basin below 0.3 eV located. The earlier refinements only followed the global minimum of each stage; three basins were unresolved: β N=4 B1 near (0.3234, 0.3619) and B2 near (0.328, 0.314), and ST N=6 on the Γ–M diagonal near (0.3194, 0.3194), never refined, with a Dirac-cone fit consistent with a near-zero gap. `tools/refine_gap_targeted.py` zooms into explicit basins with fixed-density SOC NSCF (ICHARG=11, ISYM=−1), 9×9 patches re-centred on each basin's minimum and shrunk 4× only when that minimum is interior; stages `gap_refine_3a…`, summary `gap_refine_3_summary.json`. α at K (1.15 meV, fit-converged) receives the same zoom for a local slope bound. 2H-α, 2H-β and 1H-β have eV-scale margins and need nothing further.

**Fermi level.** Only 2H-α has E_F in a global gap (indirect +4.84 eV). α, β, ST and 2H-β have indirect overlaps of −4.03, −3.52, −3.75 and −6.15 eV; 1H-β has odd filling (bands 5 and 6 each half occupied). Their parity/WCC indices describe isolated band subspaces of metals.

**Submitted 2026-09-26 (cmt queue, 24 cores each so that jobs backfill the free slot on cmt02; a 168-core request would wait for a whole node).**

| Task | PBS job | Script |
|---|---|---|
| β zoom (B1, B2; 5 levels) | 43036 | `qsub -N refine_beta -v STRUCTURE=b-Beryllene,CENTERS='0.32338:0.36190:0.000595;0.3283:0.3135:0.00476',LEVELS=5 tools/refine_gap_targeted.pbs` |
| ST zoom (diagonal basin; 5 levels) | 43037 | `CENTERS='0.31936:0.31936:0.00238',LEVELS=5` |
| α zoom at K (4 levels) | 43038 | `CENTERS='0.333333333333:0.333333333333:0.000595',LEVELS=4` |
| Spin screens α, β, ST, 2H-α, 2H-β, 1H-β | 43030–43035 | `qsub -N spin_<id> -v STRUCTURE=<dir> tools/spin_screen.pbs` |

The six first-attempt spin jobs (43021–43026, 42 cores) were deleted before starting and resubmitted at 24 cores.

**Results notebook (2026-09-26).** `../9.0_topology.ipynb` only calls `vmatplot/band_topology.py`: `extract_*` readers for the campaign reports, `summarize_*` Markdown tables and `plot_topology_bands`, `plot_direct_gap_maps`, `plot_gap_zoom`, `plot_wcc`, whose figures are saved to `figures/9_topology/9.1–9.4_*.pdf`. The module only reads files; `plot_*("help")` prints the argument order.

