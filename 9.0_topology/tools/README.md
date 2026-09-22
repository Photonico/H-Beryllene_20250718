# Beryllene SOC topology campaign

This campaign freezes the six final structures in `6.0_dielectric_selection` and reuses their PBE PAW potentials and scalar charge densities. It does not modify any earlier calculation. Legacy beta `_h` / `_hh` names mean bridge `_b` / `_bb` respectively.

| ID | Reference directory | Neutral spinor filling | Geometric inversion | Planned analysis |
|---|---|---:|---|---|
| alpha | a-Beryllene | 2 | Yes | Four-TRIM parity and sampled gaps |
| beta | b-Beryllene | 4 | Yes | Four-TRIM parity and sampled gaps |
| st | c3-Beryllene | 6 | Yes | Four-TRIM parity and sampled gaps |
| alpha_2h | a-Beryllene_hh | 4 | Yes | Four-TRIM parity and sampled gaps |
| beta_1h | b-Beryllene_b | 5 | No | Sampled gaps; conditional lowest-4/6 Wilson loops |
| beta_2h | b-Beryllene_bb | 6 | Yes | Four-TRIM parity and sampled gaps |

## Calculations

Each structure has one PBS pipeline: static SOC SCF on a Gamma 105 x 105 x 1 mesh; fixed-density SOC at all four explicit 2D TRIM; fixed-density SOC bands on the existing `4.0_bandstructure` path. Inputs use ENCUT 600 eV, Accurate precision, EDIFF 1e-8, Gaussian width 0.01 eV, 28 spinor bands, zero initial spin moment, and the reference dipole correction. Geometry and potentials are copied with source SHA256 fingerprints. No relaxation, optics, hybrid functional, or empirical dispersion is added. Scalar WAVECARs are not reused. Large dense-mesh WAVECARs are not written; the small explicit-TRIM WAVECAR is retained.

The five centrosymmetric systems use IrRep 2.1.3 to extract parity from the newly generated spinor WAVECARs. This avoids the legacy irvsp/VASP-output patch. Geometric inversion does not itself certify electronic time reversal. Pair consistency and the four actual TRIM must pass before a conditional Fu-Kane invariant is reported.

The odd electron count of neutral beta_1h excludes an ordinary fully occupied, time-reversal-symmetric band-insulator interpretation. Lowest-4 and lowest-6 subspaces are studied only when their sampled direct separation exceeds 0.1 meV. Z2Pack checks every added loop, rejects a closing direct gap, checks overlap ordering/dimensions, Kramers pairing and WCC convergence, and retains the evidence. No edge-state calculation is needed before a reliable nontrivial isolated subspace is identified.

Sampled gaps are not a rigorous all-k proof. Tiny gaps require targeted convergence checks before publication. Metallic indirect overlap must not be described as a QSH insulating gap.

## Private tools and dependencies

`manifest.json` records a separate private tools prefix outside Git. Wannier90 3.1.0 is selected because the VASP 6.5 interface supports its established serial library; a fresh copy of the licensed site VASP source is built there with the existing Intel toolchain. The original VASP installation remains untouched. `READY` confirms compilation, dynamic linkage and Python imports; the first WCC loop performs the actual interface/overlap validation. A failed prerequisite stops dependent PBS jobs.

`submitted_jobs.jsonl` records each accepted job immediately. `SUBMISSION_COMPLETE` means all nine intended PBS jobs were accepted, not that their calculations finished. Scheduler logs and per-stage validation reports distinguish completion, numerical failure and scientific ineligibility. Submission refuses a preexisting record to prevent duplicate jobs after an interruption.

## Git and data policy

Only new campaign scripts, documentation, manifests and small input files are staged. Existing runtime-output and licensed-POTCAR ignore rules are retained. For the relocated soc_topology directory, only individual files exceeding the size threshold are newly ignored. Repository paths larger than 99,000,000 bytes are added to ignore rules; already tracked historical content is preserved. Preexisting unrelated modifications are excluded from this commit.

## Reference documentation

- https://vasp.at/wiki/Makefile.include
- https://vasp.at/wiki/LWANNIER90
- https://vasp.at/wiki/NUM_WANN
- https://github.com/wannier-developers/wannier90/tree/v3.1.0
- https://github.com/irreducible-representations/irrep
- https://z2pack.greschd.ch/en/latest/reference/fp.html
- https://z2pack.greschd.ch/en/latest/reference/surface.html

## Toolchain location (2026-09-22)

The complete toolchain now lives in `9.0_topology/packages/soc_topology/`. Current scripts and manifests use this location. The old `/cmt2/lniu6305/Packages/soc_topology_20260922_1435` is a compatibility symlink for already-submitted PBS copies; retain it until those jobs finish. Historical submission records and build logs retain their original paths.

This directory is not ignored as a whole. The scan found one file larger than 99,000,000 bytes: wannier90-3.1.0.tar.gz, which is individually listed in .gitignore. Repeat the scan when adding files; Git cannot automatically ignore files by size. Copy ignored files separately or download them again. Linux binaries require site Intel/MPI/HDF5 libraries. On another machine, recreate the Python environment and rebuild with the tools scripts and locally available licensed VASP sources. The material jobs still use the site VASP executable; WCC uses the private build.
