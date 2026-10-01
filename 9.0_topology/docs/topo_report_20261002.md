# SOC band topology and stability of beryllene and hydrogenated beryllenes

Report, 2026-10-02. Update of [`topo_report_20261001.md`](topo_report_20261001.md). Results are displayed in [`9.0_topology.ipynb`](../../9.0_topology.ipynb); the analysis code is `vmatplot/band_topology.py`.

> **Status: conditional screening, not certified.** Every Z₂ value belongs to a fixed lowest-N spinor-band subspace and is conditional on electronic time-reversal symmetry and on isolation of that subspace over the whole 2D Brillouin zone (BZ). All reports carry `topology_certified: false`.

## What changed since 2026-10-01

1. **β-beryllene recomputed in the P-3m1 geometry.** The topology, HSE06, PBE phonons, bands, PDoS and dielectric run are now redone. The conclusions do not change: β is a metal with a trivial conditional ν = 0, and it is dynamically stable. Its smallest subspace gap is now an SOC gap pinned at K (Section 1).
2. **Hydrogenation thermodynamics with zero-point energy and free energy.** Only 2H-α stays stable against H₂ release at room temperature and 1 bar. 1H-β becomes endothermic once ZPE is included. 2H-β needs about 90 bar of H₂ (Section 2).
3. **2H-α is not the lowest BeH₂ monolayer.** A puckered square layer with tetrahedral Be is 38 meV per BeH₂ lower (Section 3).

All jobs 44069–44079 and 44093–44094 ended with exit 0 and converged SCFs/relaxations; the PDoS rerun 44095 was still running when this was written.

## 1. β-beryllene, P-3m1

Geometry (`2.0_geometry_optimization/b-Beryllene_p3m1`, 27 × 27 k): a = 2.16672 Å, buckling 1.9095 Å, Be–Be 2.2828 Å, E0 = −6.6577664 eV. Every later stage reuses the INCAR and KPOINTS of the old β directory next to it, with only POSCAR changed. The sibling directories are named `b-Beryllene_p3m1` and the old P-1 runs are kept.

| Quantity | P-3m1 (new) | P-1 (superseded) |
|---|---|---|
| Conditional ν (N = 4) | 0 | 0 |
| Kramers-pair parities Γ / X / Y / M | +− / +− / +− / +− | same |
| Min. direct gap E₅ − E₄ | **0.898 meV at K** (zoom to 6 × 10⁻⁵ Å⁻¹, slope bound +0.65 meV) | 0.854 meV near K |
| Position of that gap | +1.60 eV above E_F | +1.58 eV |
| Indirect gap (band overlap) | −3.51 eV (metal) | −3.52 eV |
| TR evidence: max\|m\| / Kramers splitting | 3.0 × 10⁻⁸ μB/Å³ / 5.2 × 10⁻⁸ eV | 8.9 × 10⁻⁸ μB/Å³ / 1.7 × 10⁻⁷ eV |
| Magnetic seeds (FM, AFM) | both collapse | both collapse |
| HSE06 parities (no SOC) | identical, ν = 0; Γ gap 3.471 eV (PBE 2.795 eV) | identical |
| PBE phonons | **stable everywhere**; residual force 5.7 × 10⁻³ eV/Å; max 21.6 THz | stable at commensurate q; residual 4.5 × 10⁻² eV/Å; ZA artefact |

- **The K gap is a Dirac point opened by SOC.** In P-3m1 the minimum direct gap sits exactly at K. The adaptive zoom (job 44093) keeps the minimum at K through four levels and ends with a positive local slope bound. In the P-1 cell it had moved off K, because the distortion broke C₃.
- **Fail-closed checks.** The parity screen (`9.0_topology/b-Beryllene_p3m1/trim/parity-analysis-*`) passes every float64 tolerance. The IrRep plain-metric message is 1.9 × 10⁻⁵ (single precision plus the PAW metric, as before).
- **Band and DOS reruns.** The β band and PDoS runs were repeated with NBANDS = 42 (jobs 44094, 44095), because VASP's default band count depends on the number of MPI ranks. The 24-rank runs had only 6 bands, whereas the old 42-rank runs had 42.
- **Notebook and figures.** `9.0_topology.ipynb` now uses `b-Beryllene_p3m1` for β, and keeps the P-1 rows as "β (P-1, superseded)" for comparison. Figures 9.1–9.3 and 9.11 show β in P-3m1; the 9.11 phonon path is Γ–M–K–Γ.
- **Not yet repointed:** the optics and thesis notebooks (`6.0_dielectric_function_convergence`, `6.1_dielectric_function`, `7.1`, `7.2`, `0.0`, `0.1`) still read `*/b-Beryllene`. The new dielectric run is in `6.0_dielectric_selection/b-Beryllene_p3m1`.

## 2. Hydrogenation thermodynamics

Adsorption energy per H against ½H₂, computed from:
- **Electronic energies:** PBE relaxations at 27 × 27 k, with β in P-3m1.
- **Zero-point energies:** harmonic PBE phonons on a 40 × 40 q mesh, imaginary modes left out. For H₂ the stretching mode is 4341 cm⁻¹, giving ZPE 0.269 eV (`1.0_convergence/H2_molecule_freq`, job 44077).
- **Free energy at 298.15 K and 1 bar:** the vibrational free energy of the layers, plus H(T) − H(0) − TS of H₂ from JANAF (−0.316 eV per H₂).

| Structure | ΔE | ΔE + ΔZPE | ΔG(298 K, 1 bar) | Equilibrium p(H₂) at 298 K |
|---|---|---|---|---|
| 2H-α | −0.331 eV | −0.243 eV | **−0.084 eV** | 1 × 10⁻³ bar |
| 1H-β | −0.072 eV | **+0.012 eV** | +0.166 eV | 4 × 10⁵ bar |
| 2H-β | −0.173 eV | −0.090 eV | +0.058 eV | 9 × 10¹ bar |

**Implications:**
- Zero-point energy raises every adsorption energy by 0.08–0.09 eV per H, because the Be–H modes outweigh the H₂ stretch.
- **2H-α** is thermodynamically stable against H₂ release at ambient conditions.
- **2H-β** is stable only under H₂ pressure (~90 bar at room temperature) or at low temperature.
- **1H-β** is endothermic even at 0 K once ZPE is included. In vacuum it can exist only as a kinetically trapped phase.
- These are harmonic estimates. The imaginary 2H-β pocket (one branch, 1% of the BZ) is left out of its free energy; its contribution is below 0.1 meV per H.
- The manuscript should state the ZPE- and temperature-corrected values, not only ΔE.

## 3. BeH₂ monolayer polymorphs

Two square seeds were relaxed with the same settings as 2H-α (`2.1_geometry_optimization_with_hydrogen_a`, jobs 44078, 44079):

| Structure | Space group | Be coordination | E per BeH₂ relative to 2H-α |
|---|---|---|---|
| 2H-α (H in opposite hollows) | P-3m1 | 6 H at 1.599 Å | 0 |
| 2H-α, H in aligned hollows | Cmmm | — | +741 meV |
| BeH₂ square, puckered | **P-4m2** | 4 H at 1.468 Å (tetrahedral, each H bridging two Be) | **−38 meV** |
| BeH₂ square, planar | P4/mmm | 4 H at 1.479 Å (planar) | +1233 meV |

- 2H-α is therefore metastable with respect to the puckered square layer, whose tetrahedral Be–H–Be bridging resembles bulk BeH₂. Converting 2H-α into it requires rearranging the Be lattice from triangular to square, so 2H-α can be kinetically stable as the product of hydrogenating α-beryllene. The manuscript should say so.
- Only two seeds were tried; this is not a global structure search.
- The puckered layer's own dynamical stability was not computed.

## 4. Notebooks and figures updated

The new data are now used in notebooks 0.0–7.2 and `9.0_topology.ipynb`. Every replaced PDF was kept next to the new one with the suffix `_old`: 59 in `figures/`, `figures_for_publication/` and `figures_for_thesis_ch3/`, and 8 in `figures/9_topology/`.

| Notebook | Change | Replaced figures |
|---|---|---|
| 0.0, 0.1 thesis demos | β dielectric from `6.0_dielectric_selection/b-Beryllene_p3m1`; thesis figures now saved to `figures_for_thesis_ch3/` (the old target `figures_for_thesis/` no longer exists) | fig3.13–3.15, S3.16, S3.17; thesis 3.1, 4.1–4.5 |
| 2.0, 2.1 geometry | β from `b-Beryllene_p3m1` (a = 2.16672 Å) | — (printed output only) |
| 3.0 phonons | α, β, ST from the PBE runs in `3.8_phonon_dispersion_pbe` (old: VASP finite differences with PBE+D3); paths as published (Γ–K–M–Γ, ST Γ–M–X–Γ) | S3.13b, S3.13d, fig3.4; figures 3.11, 3.12, 3.21 |
| 3.5 selected phonons | 2H-α, 1H-β, 2H-β from `3.7_phonon_dispersion_with_hydrogen_pbe` (old: PBE+D3 `3.5_*` and `*_ollie` runs); 1H-β and 2H-β on Γ–X–M–Γ–Y–M′–Γ, which shows the 2H-β soft mode | fig3.5, fig3.6a/b (+ `_alt`); figures 3.31–3.33 |
| 4.0, 4.1 band gaps | β → P-3m1; `extract_bandgap_outcar` fixed (it matched the band index anywhere in a line and returned meaningless values, also in the committed outputs); 4.1 shows the HSE06 gap table | — |
| 4.2 bands | β → P-3m1 | S3.15d; figure 4.12 |
| 4.3 bands with H | 2H-α: PBE (blue, solid) and HSE06 (orange, solid), each aligned at its VBM | fig3.9; figure 4.31 |
| 5.0, 5.2, 5.4 DOS/PDoS | β → P-3m1 (NBANDS 42) | figures 5.10, 5.11, 5.12 (PDoS) |
| 6.1 dielectric | β → P-3m1 | fig3.16, fig3.18; figures 6.10, 6.11, 6.20, 6.22 |
| 7.1, 7.2 optics | β → P-3m1 | fig3.17a/b, fig3.19a/b, S3.18a–c, S3.19a–c; figures 7.11–7.15, 7.21–7.25 |
| 9.0 topology | HSE06 orange, PBE blue (9.8); phonon paths as above (9.5, 9.10–9.12) | 9.1–9.3, 9.5 (`_old` = last committed version), 9.8, 9.10, 9.12, 9.11 (`_old` = the P-1 β) |

**Other changes:**
- Paths whose case differs from the directory names (e.g. `4.0_bandstructure/b-beryllene`) were corrected on executable lines, so the notebooks also run on a case-sensitive file system. Savefig names now match the existing PDF names exactly.
- Figures whose data did not change were restored byte for byte after execution.

**Not changed:**
- The bulk hcp and bcc phonons (S3.13a, S3.13c) still use PBE+D3; they were not recomputed.
- `3.5_phonon_dispersion_with_hydrogen_phonopy.ipynb` (working figures 3.31–3.33 `_phonopy`) still shows the old PBE+D3 phonopy runs.
- `fig3.5_phonon_a_hh_alt.pdf`, whose savefig is commented out, is the old version.
- The 1H-β and 2H-β band structures (fig3.10) come from the `*_ollie` runs. Their geometry matches the references, and D3 does not change eigenvalues.
- The β dielectric convergence tests in 6.0 (S3.7, S3.8) still use the P-1 geometry.
- 4.1 cells 2–3 point to `4.0_bandstructure/b-Beryllene_h`, `b-Beryllene_hh`, which contain no OUTCAR (as before).

## 5. Pending work

| Item | Status |
|---|---|
| Repoint the optics and thesis notebooks to `b-Beryllene_p3m1` | not done (user notebooks) |
| Anharmonic phonons of 2H-β (SSCHA) | not submitted: python-sscha and QE are not installed, and the 5 × 5 supercell with the dense k mesh the soft mode needs costs ~150 core-h per configuration, ~10⁴–10⁵ core-h in total |
| Phonons of the puckered BeH₂ layer | optional; only if the polymorph is discussed beyond an energy comparison |

The limitations of the topology screening are unchanged from the 2026-10-01 report.

## Files

| Item | Location |
|---|---|
| β P-3m1 chain | `*/b-Beryllene_p3m1` in 2.0, 2.4, 3.8, 4.0, 4.4, 5.0, 6.0_dielectric_selection, 9.0_topology |
| H₂ frequency | `1.0_convergence/H2_molecule_freq` |
| BeH₂ polymorphs | `2.1_geometry_optimization_with_hydrogen_a/BeH2_square_*`, `BeH2_square_note` |
| Energetics tables | `summarize_hydrogenation_energetics`, `summarize_polymorph_energies` in `vmatplot/band_topology.py` |
| Updated values | `readme.md` (β section, ZPE / ΔG sub-bullets, phonon lines) |
