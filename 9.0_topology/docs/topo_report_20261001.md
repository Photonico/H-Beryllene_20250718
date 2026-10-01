# SOC band topology and stability of beryllene and hydrogenated beryllenes

> **Partly superseded by [`topo_report_20261002.md`](topo_report_20261002.md):** β recomputed in P-3m1 (Section 4 here), hydrogenation energies with ZPE and free energy, BeH₂ polymorphs.

Report, 2026-10-01. Update of [`topo_report_20260928.md`](topo_report_20260928.md); methods and references are in [`topo_report_20260926.md`](topo_report_20260926.md). Results are displayed in [`9.0_topology.ipynb`](../../9.0_topology.ipynb); the analysis code is `vmatplot/band_topology.py`.

> **Status: conditional screening, not certified.** Every Z₂ value belongs to a fixed lowest-N spinor-band subspace and is conditional on electronic time-reversal symmetry and on isolation of that subspace over the whole 2D Brillouin zone (BZ). Five of the six structures are metallic at neutral filling, so their indices are not Fermi-level quantum-spin-Hall invariants. All reports carry `topology_certified: false`.

## What changed since 2026-09-28

1. **HSE06 leaves the band order unchanged.** All parities at the four TRIM are identical to PBE+SOC for the five centrosymmetric structures, so every conditional ν stands. The HSE06 gap of 2H-α is 6.00 eV (PBE 4.84 eV) (Section 2).
2. **PBE phonons are complete.** α, ST and 2H-α are stable over the whole BZ. 2H-β has a harmonically imaginary Be mode at q = (0.4, 0.6) (−1.65 THz), but its double well is only 0.02–0.03 meV per 20-atom cell, so no static distortion is expected (Section 3).
3. **The β-beryllene geometry is superseded.** Its P-1 cell was an artefact of the 11 × 11 k mesh used in the relaxation. With 27 × 27 k, β relaxes to the hexagonal P-3m1 buckled honeycomb, 4.5 meV per cell lower. Every β result, here and in the earlier DFT chain, refers to the old geometry (Section 4). Recomputing them needs a decision.
4. **Analysis fix for hexagonal cells.** Phonopy's default tolerance (1e-5) misses the hexagonal symmetry of α and 2H-α, whose cells deviate from it by about 1e-4 Å. This created an artificial imaginary pocket in α. These two cells are now analysed with 1e-3 (Section 3).

## 1. Summary

| Structure | Conditional Z₂ (PBE+SOC) | HSE06 check (no SOC) | Neutral Fermi level | Dynamical stability (PBE, no D3) |
|---|---|---|---|---|
| α | ν = 1 (N = 2) | same parities, ν = 1 | metal | stable (P6/mmm) |
| β (P-1, superseded) | ν = 0 (N = 4) | same parities, ν = 0 | metal | stable at commensurate q; P-3m1 not yet computed |
| ST | ν = 1 (N = 6) | same parities, ν = 1; boundary at M is a SOC-split orbital doublet | metal | stable |
| 2H-α | ν = 0 (N = 4) | same parities, ν = 0; gap 6.00 eV | insulator | stable |
| 1H-β | Z₂ = 0 / 0 (WCC, N = 4 / 6) | — (no inversion) | metal (odd filling) | stable at commensurate q; ZA artefact near Γ |
| 2H-β | ν = 0 (N = 6) | same parities, ν = 0; M gap 0.171 eV | metal | imaginary Be mode at q = (0.4, 0.6), −1.65 THz; double well 0.02–0.03 meV per 20 atoms |

## 2. HSE06 check

**Setup** (`4.4_bandstructure_hse`, jobs 43630–43644, all exit 0):
- HSE06 (HFSCREEN 0.2, AEXX 0.25), full 12 × 12 exchange q mesh (no NKRED), PRECFOCK = Normal, no SOC.
- 20 Å cell with dipole correction, Γ-centred 12 × 12 k mesh containing all four TRIM.
- Each job first runs PBE in the same cell, then restarts HSE06 from it (ALGO = Damped). All SCFs reached EDIFF = 1e-6 within 12–21 steps.

**Parities.** `9.0_topology/tools/hse_parity.py` builds the float64 inversion matrix from the HSE06 WAVECAR. It applies the same fail-closed tolerances as the SOC screen, and every one is met by a wide margin:

| Check | Largest value | Tolerance |
|---|---|---|
| Opposite-parity overlap | 1.8 × 10⁻⁸ | 10⁻⁶ |
| Off-block element of the Löwdin inversion matrix | 2.5 × 10⁻⁸ | 10⁻⁶ |
| Inversion residual | 1.2 × 10⁻⁷ | 10⁻⁵ |

Without SOC every spinless band is one Kramers pair, so the lowest N/2 bands give the Fu–Kane δ of the lowest-N subspace. TRIM missing from the IBZ (Y of α, ST and 2H-α) are mapped from an equivalent TRIM with the lattice-translation phase of the conjugated inversion.

| Structure | Parities Γ / X / Y / M, bands 1…N/2 (HSE06) | Same as PBE+SOC | ν HSE06 / PBE+SOC | E(N+1) − E(N) at the critical TRIM, HSE06 / PBE (20 Å) / PBE+SOC (40 Å) |
|---|---|---|---|---|
| α | + / − / − / − | yes | 1 / 1 | Γ: 6.209 / 5.169 / 5.103 eV |
| β (P-1) | +− / +− / +− / +− | yes | 0 / 0 | Γ: 3.471 / 2.808 / 2.810 eV |
| ST | +−+ / +−+ / +−+ / −−+ | yes | 1 / 1 | M: 0 / 0 / 0.0013 eV (SOC-split doublet; next gap 1.98 eV) |
| 2H-α | +− / −+ / −+ / −+ | yes | 0 / 0 | X: 6.847 / 5.644 / 5.645 eV |
| 2H-β | +−+ / −++ / +−− / −+− | yes | 0 / 0 | M: 0.171 / 0.178 / 0.177 eV |

- **No inversion is created or removed.** HSE06 widens the TRIM gaps of α, β, ST and 2H-α by 0.4–1.4 eV. In 2H-β the gaps at X and Y (+0.05 eV) and the critical gap at M (0.178 → 0.171 eV) barely move, and the M gap stays open.
- **ST at M.** The N = 6 boundary lies inside an orbital doublet that only SOC splits (1.3 meV in PBE+SOC). SOC commutes with inversion and the doublet has a single parity (+), so δ_M does not depend on which SOC-split pair lies lower. ST's ν = 1 therefore does not hinge on the SOC ordering at M. The SOC gap still sets the isolation of the subspace.
- **Cell and mesh are adequate.** The PBE step in the 20 Å cell reproduces the 40 Å / 105 × 105 PBE+SOC gaps at the TRIM:

  | Structure | Agreement |
  |---|---|
  | 2H-α | within 1 meV |
  | 2H-β (M gap) | within 0.7 meV |
  | β, ST, 2H-β (other TRIM) | within 8 meV |
  | α | within 65 meV; a metal on a coarse SCF mesh, but its gaps exceed 5 eV |

  No denser mesh or larger vacuum is needed for these conclusions.
- **Limits.** HSE06 is applied without SOC and only at the TRIM. The meV spin–orbit gaps that isolate the α and ST subspaces elsewhere in the BZ are not tested by it.

**2H-α band gap.** The band edges are taken over 132 k points: the 12 × 12 mesh, a 71-point Γ–K–M–Γ path, and 9 points around each PBE extremum. These were split into 10 independent HSE06 jobs.

| Calculation | VBM k | CBM k | Indirect gap | Min. direct gap |
|---|---|---|---|---|
| PBE, 20 Å (same jobs) | (0.2546, 0.2546) | (0.4268, 0) | 4.842 eV | 4.925 eV |
| **HSE06**, 20 Å | (0.2546, 0.2546) | (0.4268, 0) | **5.998 eV** | **6.089 eV** |
| PBE+SOC, 40 Å, full-BZ search (9.0_topology) | (0.2546, 0.2546) | (0.4268, 0) | 4.842 eV | 4.925 eV |

- The band edges do not move under HSE06. Both extrema stay at the centre of their sampling windows.
- The gap opens by 1.16 eV. This is the scissor shift to use if the PBE optical spectra of 2H-α are to be corrected.
- Figure: `figures/9_topology/9.8_hse_bands_a_hh.pdf`.

## 3. Phonons with consistent settings

**Settings:**
- PBE without D3, SIGMA 0.05, ±0.01 Å displacements.
- Hydrogenated structures in `3.7_phonon_dispersion_with_hydrogen_pbe`; pristine ones in `3.8_phonon_dispersion_pbe` (5 × 5 × 1, ST 4 × 4 × 1).
- k: 13 × 13 for the 5 × 5 metals, 16 × 16 for ST (65 × 65 / 64 × 64 unit-cell equivalent), 7 × 7 for the insulator 2H-α.
- The pristine runs replace the earlier VASP finite-difference runs in `3.0_phonon_dispersion`, which used PBE+D3.

| Structure | Space group (phonopy) | Residual force (eV/Å) | Min. ν, commensurate q | Min. ν, 60 × 60 grid | Max. ν |
|---|---|---|---|---|---|
| α | P6/mmm (tolerance 1e-3) | — (single displacement; zero by symmetry) | +0.000 | +0.000 | 24.4 THz |
| β (P-1 geometry) | P-1 | 4.5 × 10⁻² | −0.000 | −0.30 THz near Γ (ZA artefact) | 22.1 THz |
| ST | P4/mmm | 2.7 × 10⁻² | −0.000 | −0.000 | 21.5 THz |
| 2H-α | P-3m1 (tolerance 1e-3) | 4.6 × 10⁻⁴ | −0.000 | −0.000 | 44.4 THz |
| 1H-β | P1 | 9.9 × 10⁻³ | −0.000 | −1.45 THz near Γ (ZA artefact) | 48.6 THz |
| 2H-β | P-1 | 6.9 × 10⁻³ | **−1.649 THz at (0.4, 0.6)** | −1.76 THz | 50.5 THz |

**Symmetry tolerance.** The α and 2H-α cells deviate from hexagonal by about 1e-4 Å; for α the lattice angle is 120.003°. At phonopy's default 1e-5 tolerance the cells are found as Cmmm / C2/m.
- The interpolation then weights equidistant supercell images unequally and breaks the C6 equivalence of the Γ–M directions. For α this gave +1.47 THz along b₁ and −1.28 THz along b₁ − b₂ (a 7% "pocket").
- The old PBE+D3 force constants of `3.0_phonon_dispersion/a-Beryllene` show the same artefact (−0.92 THz) at 1e-5 and none at 1e-3.
- With 1e-3 (P6/mmm, P-3m1) both structures are stable everywhere.
- No recomputation was needed: for α the single displacement is paired with its inversion image by symmetry.

**β residual force.** The residual is 0.045 eV/Å, compared with < 0.01 eV/Å in the other metals. It signals the geometry problem of Section 4.

**2H-β soft mode.**
- **Location.** It sits at the commensurate q = (0.4, 0.6) (POSCAR reciprocal basis). This is equivalent to −(−0.4, 0.4), on the Γ–M′ segment of figure 9.9. Its value is exact for the 5 × 5 supercell and does not depend on interpolation.
- **Same mode as before.** The eigenvector is identical to that of the old PBE+D3 run (overlap 1.000). That run gave −0.646 THz, so D3 stiffens the mode. The next three modes at this q agree within 0.07 THz.
- **Character.** 99% of the eigenvector weight is on Be.
- **Extent.** The imaginary region covers 1.1% of the BZ around ±(0.4, −0.4).

**Frozen phonon** (`3.6_phonon_soft_mode_check/b-Beryllene_bb`, jobs 43606–43614, all exit 0). The structure is modulated along this eigenvector, and E = αQ² + βQ⁴ is fitted in the mass-weighted amplitude Q:

| SIGMA | Supercell k (unit-cell equivalent) | E − E(0) at max\|u\| = 0.025 / 0.049 / 0.098 Å | Harmonic ν | Double-well depth (per 20 atoms) | Fit residual |
|---|---|---|---|---|---|
| 0.02 eV | 16 × 97 (97 × 97) | +0.015 / −0.007 / +10.03 meV | (1.75i THz) | (0.19 meV) | 0.16 meV: not converged in k |
| 0.05 eV | 11 × 65 (65 × 65) | −0.013 / +0.578 / +12.45 meV | 1.02i THz | 0.021 meV at 0.019 Å | 0.001 meV |
| 0.10 eV | 11 × 65 (65 × 65) | −0.012 / +0.495 / +11.93 meV | 1.13i THz | 0.032 meV at 0.022 Å | 0.018 meV |

- **The harmonic term is negative** (SIGMA 0.05 and 0.10 eV), in agreement with the finite-displacement value. That value is larger, 1.65i THz, because it probes smaller amplitudes (0.01 Å). Undistorted 2H-β is therefore a saddle point of the PBE energy surface, as the 09-28 report did not yet establish.
- **The instability is negligible in energy.** The double well is ≤ 0.03 meV per 20-atom cell, about 1–1.5 μeV per atom. The quartic term dominates beyond 0.02 Å. The zero-point energy of a ~1 THz mode is ħω/2 ≈ 2 meV, two orders of magnitude larger. No static charge-density-wave distortion is expected at any temperature, and the symmetric structure is the appropriate one, but the mode is strongly anharmonic.
- **No sign of an electronically driven softening.** The harmonic term does not soften from SIGMA 0.10 to 0.05 eV, and χ₀(q) has no peak at this q (09-28 report, Section 5). This is consistent with a q-dependent electron–phonon matrix element or a lattice effect, not with Fermi-surface nesting. The SIGMA 0.02 eV series is too noisy in k to resolve the harmonic term.
- **Consequence for electron–phonon work.** A harmonic DFPT/EPW calculation will reproduce this imaginary pocket in 2H-β, with or without D3. λ from such a calculation is ill-defined there. An anharmonic treatment (e.g. SSCHA or temperature-dependent effective potentials), or at least an explicit estimate of the pocket's contribution, is required before a T_c of 2H-β can be quoted.

**Figures:**

| Figure | Structure | Shows |
|---|---|---|
| `9.9_phonon_b_bb_pbe.pdf` | 2H-β | soft mode just after M′ |
| `9.10_phonon_a_pbe.pdf` | α | dispersion |
| `9.11_phonon_b_pbe.pdf` | β | dispersion (P-1 geometry) |
| `9.12_phonon_st_pbe.pdf` | ST | dispersion |

## 4. β-beryllene geometry

- **Original relaxation.** `2.0_geometry_optimization/b-Beryllene` used 11 × 11 k (α 27 × 27, ST 25 × 25, others ≥ 23 × 23) and ended in a P-1 cell (the starting cell was already distorted): |a₁| = 2.147 Å, |a₂| = 2.186 Å, angle 119.80°.
- **Static check at 27 × 27 k** (`b-Beryllene_k27_check`, job 43605). The same geometry has forces up to 0.066 eV/Å and an in-plane shear stress of 1.9 kB. The phonon supercell (65 × 65 equivalent) independently gives 0.045 eV/Å. Both are far above the relaxation criterion of 10⁻³ eV/Å.
- **Re-relaxation at 27 × 27 k** (`b-Beryllene_k27_relax`, job 44066), same INCAR. The cell became hexagonal, and E0 dropped by 4.5 meV per cell (2.3 meV per atom). The run ended with a ZBRENT error after the energy had been constant to 10⁻⁷ eV for nine ionic steps.
- **Final symmetrised relaxation** (`b-Beryllene_p3m1`, job 44067, exit 0, converged):

  | Quantity | Value |
  |---|---|
  | Space group | **P-3m1** (low-buckled honeycomb, D3d) |
  | a | 2.16672 Å |
  | Buckling | 1.9095 Å |
  | Be–Be | 2.2828 Å |
  | E0 | −6.6577664 eV, vs −6.6532613 eV for the P-1 geometry at the same k |

**Consequences:**
- Every β result refers to a geometry that is not an equilibrium structure at converged k. This covers the topology screening (ν = 0, minimum direct gap 0.854 meV), the HSE06 check, the PBE phonons, and the band, DOS and optical results of the earlier chain (`6.0_dielectric_selection/b-Beryllene` and everything built on it). The higher symmetry (D3d with C3) can change the band degeneracies, so the β topology result cannot simply be carried over.
- 1H-β and 2H-β were relaxed separately, with hydrogen, at 27 × 27 k, so their geometries are not affected. Their energies relative to pristine β do change, because `readme.md` uses the 11 × 11 energy of β (−6.6470385 eV). At 27 × 27 k the β energy is −6.6532613 eV in the P-1 geometry (k-mesh effect −6.2 meV) and −6.6577664 eV in P-3m1 (relaxation −4.5 meV).
  - Adsorption per H vs ½H₂: 1H-β −0.083 → **−0.073 eV**; 2H-β −0.178 → **−0.173 eV**.
  - β cohesive energy: −6.570 → **−6.581 eV** per cell (−3.290 eV per atom).
  - `readme.md` has not been changed yet.
- **Recommended recomputation** (not yet submitted; needs a decision because it redefines β in the manuscript):
  1. β topology chain (`9.0_topology/b-Beryllene`: SOC SCF, bands, TRIM parities, gap refinement, TR evidence, spin screen);
  2. PBE phonons of β (P-3m1);
  3. HSE06 check of β;
  4. the earlier β band, DOS and optics chain.

  Items 1–3 are small (2-atom cell).

## 5. Consequences

- **Manuscript and thesis.**
  - The conditional Z₂ values stand under HSE06 for α, ST, 2H-α and 2H-β. Keep the wording "Z₂-nontrivial isolated band subspace of a metal" for α and ST.
  - β must be redescribed as the P-3m1 buckled honeycomb, and its results replaced once recomputed.
  - Report the 2H-α gap as 4.84 eV (PBE) and 6.00 eV (HSE06).
  - Report 2H-β as harmonically unstable at q = (0.4, 0.6) with a negligible double well, i.e. stabilised by anharmonicity, not as "stable".
- **Superconductivity collaborator** (`10_superconductivity/note_superconductivity`, updated):
  - Use the P-3m1 geometry for β.
  - Expect an imaginary harmonic pocket at q ≈ (0.4, 0.6) in 2H-β. The undistorted structure is still the right one; its λ needs an anharmonic treatment or an explicit analysis of the pocket.

## 6. Pending work and limitations

| Item | Status |
|---|---|
| β recomputation in the P-3m1 geometry (Section 4) | decision pending |
| Anharmonic phonons of 2H-β (e.g. SSCHA) | not started; needed for a quantitative λ / T_c of 2H-β |

Remaining limitations of the topology screening:
1. No continuous-BZ gap proof.
2. No cutoff or density convergence study of the meV gaps.
3. HSE06 only without SOC and at the TRIM.
4. Supercell magnetic orders are not sampled.
5. No edge-state calculation.

## Files and figures

| Item | Location |
|---|---|
| Results notebook | `9.0_topology.ipynb` |
| HSE06 jobs, parities (`scf/hse_parity.json`) | `4.4_bandstructure_hse/`, note in `4.4_bandstructure_hse/note` |
| HSE06 parity tool | `9.0_topology/tools/hse_parity.py` |
| PBE phonons | `3.7_phonon_dispersion_with_hydrogen_pbe/`, `3.8_phonon_dispersion_pbe/` |
| Frozen phonon | `3.6_phonon_soft_mode_check/b-Beryllene_bb/` |
| β geometry | `2.0_geometry_optimization/b-Beryllene_k27_check`, `b-Beryllene_k27_relax`, `b-Beryllene_p3m1` |
| Figures | `figures/9_topology/9.8_hse_bands_a_hh.pdf`, `9.9–9.12_phonon_*_pbe.pdf` |
| Phonon audit | `3.5_phonon_dispersion_with_hydrogen_selected/audit_20260927.md` |
