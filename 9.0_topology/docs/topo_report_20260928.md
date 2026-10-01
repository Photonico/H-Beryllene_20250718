# SOC band topology and stability of beryllene and hydrogenated beryllenes

> **Partly superseded by [`topo_report_20261001.md`](topo_report_20261001.md):** the β geometry (P-1 → P-3m1), the 2H-β soft mode (harmonically imaginary, −1.65 THz, with a negligible double well), the HSE06 check and the pristine phonons.

Report, 2026-09-28. Update of [`topo_report_20260926.md`](topo_report_20260926.md), which keeps the full methods, parity tables and references. Results are displayed in [`9.0_topology.ipynb`](../../9.0_topology.ipynb); the analysis code is `vmatplot/band_topology.py`.

> **Status: conditional screening, not certified.** Every Z₂ value belongs to a fixed lowest-N spinor-band subspace and is conditional on electronic time-reversal (TR) symmetry and on isolation of that subspace over the whole 2D Brillouin zone (BZ). Five of the six structures are metallic at neutral filling, so their indices are not Fermi-level quantum-spin-Hall invariants. All reports carry `topology_certified: false`.

## What changed since 2026-09-26

1. **Magnetic stability checked.** All 11 spin-polarised seeds collapse to the nonmagnetic state (Section 2).
2. **Hydrogen reference energies corrected.** The isolated H atom had been computed without spin polarisation, and adsorption was referenced to atoms instead of H₂. Corrected adsorption energies are 10–40 times smaller (Section 3.1).
3. **Phonons of the hydrogenated structures audited and recomputed.** The selected force calculations used PBE+D3 on PBE geometries, and the band path missed a 2H-β soft mode. 2H-α and 1H-β are recomputed with consistent PBE settings; 2H-β is running (Section 3.2).

## 1. Summary

| Structure | Conditional Z₂ | Min. direct gap E(N+1)−E(N) | Neutral Fermi level | Magnetic seeds | Dynamical stability (PBE) |
|---|:--:|---|---|---|---|
| α | ν = 1 (N = 2) | 1.127 meV near K, +2.72 eV above E_F | metal (−4.03 eV) | collapse | not recomputed in this campaign |
| β | ν = 0 (N = 4) | 0.854 meV, +1.58 eV above E_F | metal (−3.52 eV) | collapse | not recomputed |
| ST | ν = 1 (N = 6) | 0.735 meV on Γ–M, +1.19 eV above E_F | metal (−3.75 eV) | collapse | not recomputed |
| 2H-α | ν = 0 (N = 4) | 4.925 eV | insulator (+4.84 eV) | collapse | **stable** |
| 1H-β | Z₂ = 0 / 0 (WCC, N = 4 / 6) | 0.929 / 0.415 eV | metal (odd filling) | collapse | **stable** at all commensurate q; ZA artefact near Γ |
| 2H-β | ν = 0 (N = 6) | 0.177 eV | metal (−6.15 eV) | collapse | nearly zero-frequency mode at q = (0.4, 0.6); full PBE run in progress |

Values in parentheses after "metal"/"insulator" are the sampled indirect gaps (negative = band overlap).

**Topology conclusions are unchanged.** Only 2H-α has a Fermi-level classification, and it is trivial. The ν = 1 of α and ST describes isolated band subspaces of metals, protected by spin–orbit gaps of about 1 meV that lie 1.2–2.7 eV above E_F. These band inversions are therefore decoupled from the Fermi surface and cannot be linked to superconductivity.

## 2. Electronic time reversal

**Converged SOC states** (unchanged, `<structure>/tr_evidence.json`):
- max|m(r)| ≤ 2.6 × 10⁻⁷ μB/Å³.
- Kramers splitting at the TRIM ≤ 3.8 × 10⁻⁶ eV.
- For 1H-β, E(k) = E(−k) to 4.0 × 10⁻⁷ eV over 5512 ±k pairs.

**Stability against magnetic order** (jobs 43030–43035, all exit 0, `<structure>/spin_screen/`):

| Structure | Seeds | Result |
|---|---|---|
| α | FM | collapses to m = 0 |
| β | FM; AFM (inversion-odd) | both collapse |
| ST | FM; (+, −, +); inversion-odd (+, 0, −) | all collapse |
| 2H-α | FM | collapses |
| 1H-β | FM; AFM | both collapse |
| 2H-β | FM; AFM (inversion-odd) | both collapse |

- These are collinear ISPIN = 2 SCFs on the same geometry, 105 × 105 mesh, smearing and dipole setup.
- Every final total and site moment is below 10⁻³ μB.
- All seeds of one structure end at the same energy (spread ≤ 7 × 10⁻⁷ eV).
- This supports a nonmagnetic, TR-symmetric ground state against ferromagnetic and in-cell antiferromagnetic order at the PBE level. Supercell (q ≠ 0) magnetic orders were not sampled.

## 3. Corrections to the earlier DFT results

### 3.1 Hydrogen reference energies

`1.0_convergence/single_Hydrogen` used ISPIN = 1, giving E_H = −0.0142 eV for a non-spin-polarised H atom. The readme energetics used it for both cohesive and adsorption energies.

**New reference calculations**, same ENCUT, PREC and functional as the slabs, in a 20 Å box:
- **H atom:** spin polarised (ISPIN = 2), E_H = −1.11742 eV, moment 1.000 μB (`1.0_convergence/single_Hydrogen_spin`).
- **H₂ molecule:** E_H₂ = −6.77276 eV, relaxed H–H distance 0.750 Å (`1.0_convergence/H2_molecule`).
- **Check:** these give a PBE H₂ binding energy of 4.538 eV.

| Structure | Adsorption per H vs ½H₂ (new) | Old "absorption" (vs non-spin-polarised H atoms) | Cohesive energy new / old (eV per formula unit) |
|---|---|---|---|
| 2H-α | **−0.331 eV** | −7.405 eV for 2 H | −8.133 / −10.339 |
| 1H-β | **−0.083 eV** | −3.455 eV for 1 H | −8.922 / −10.026 |
| 2H-β | **−0.178 eV** | −7.100 eV for 2 H | −11.464 / −13.671 |

**Implications:**
- Hydrogenation is still exothermic, but only weakly. 1H-β is close to thermoneutral with respect to H₂.
- `readme.md` has been corrected.
- The article and thesis text on the Mac must be updated as well.

### 3.2 Phonons of the hydrogenated structures

The full audit is in `3.5_phonon_dispersion_with_hydrogen_selected/audit_20260927.md`.

**Problems in the earlier phonon results:**
1. **Functional mismatch.** The force calculations used PBE+D3 (IVDW = 12) and SIGMA 0.1 on geometries relaxed with PBE without D3 and SIGMA 0.05. The structures were therefore not at equilibrium for the force field used: residual forces were 0.014–0.063 eV/Å.
2. **Incomplete band path.** The hexagonal Γ–K–M–Γ path was used for C2/m, P1 and P-1 cells, so inequivalent directions were never plotted. It missed a 2H-β mode at −0.646 THz at the commensurate q = (0.4, 0.6) (POSCAR reciprocal basis). The notebook therefore reported only −0.013 THz at Γ.
   The mode appears only with the dense supercell k mesh: the otherwise identical run with 5 × 5 k (25 × 25 unit-cell equivalent, `b-Beryllene_bb_alt2`) has no imaginary frequency at any commensurate q (−0.000 THz), whereas 13 × 13 k (65 × 65) gives −0.646 THz at q = (0.4, 0.6). It also survives removal of the residual forces. This k sensitivity points to an electronic origin rather than a numerical artefact. It is not Fermi-surface nesting, however (Section 5): the bare susceptibility has no peak at this q.
3. **Interpolation artefact.** The 1H-β imaginary pocket is not a real instability: all commensurate q are stable, and the pocket lies at |q| ≤ 0.12 on the flexural (ZA) branch.

**Recomputation** (`3.7_phonon_dispersion_with_hydrogen_pbe`):
- PBE without D3, SIGMA 0.05, the 6.0 reference geometries.
- 5 × 5 × 1 supercell, ±0.01 Å displacements.
- k meshes: 7 × 7 for 2H-α (insulator); 13 × 13 for the metals, equivalent to a 65 × 65 unit-cell mesh.
- New band path Γ–X–M–Γ–Y–M′–Γ, which crosses all inequivalent TRIM.

| Run | Forces | Residual force (eV/Å) | Min. ν, commensurate 5 × 5 q (THz) | Min. ν, 60 × 60 grid (THz) | Max. ν (THz) |
|---|---|---|---|---|---|
| 2H-α, old | PBE+D3, σ 0.1 | 1.4 × 10⁻² | −0.000 | −0.000 (Γ) | 44.5 |
| **2H-α, new** | PBE, σ 0.05 | 4.6 × 10⁻⁴ | −0.000 | −0.000 (Γ) | 44.4 |
| 1H-β, old | PBE+D3, σ 0.1 | 5.2 × 10⁻² | −0.000 | −1.43 near Γ | 48.5 |
| **1H-β, new** | PBE, σ 0.05 | 9.9 × 10⁻³ | −0.000 | −1.45 near Γ (ZA artefact) | 48.6 |
| 2H-β, old | PBE+D3, σ 0.1 | 4.7 × 10⁻² | **−0.646** at (0.4, 0.6) | −0.75 | 50.5 |
| 2H-β, new | PBE, σ 0.05 | running (array 43597[]) | — | — | — |

**2H-α.** Stable everywhere. The frequencies agree with the old run to within 0.05 THz at Γ, so D3 hardly affects this insulator.

**1H-β.**
- Stable at every commensurate q.
- The imaginary pocket near Γ remains after force-constant symmetrisation. It is the known finite-supercell failure of the ZA branch; a clean dispersion needs rotational-invariance (Huang) constraints or a larger supercell.
- The inconsistency did matter here: some modes shift by several THz. At q = (0.4, 0.2), for example, one branch moves from 19.41 THz (old) to 15.46 THz (new).
- The residual force (9.9 × 10⁻³ eV/Å) is larger than for 2H-α. Its source is the k-mesh difference between the relaxation (27 × 27) and the phonon runs (65 × 65 equivalent) in this metal; the ± displacement pairs cancel it to second order.

**2H-β soft mode, frozen phonon** (`3.6_phonon_soft_mode_check/b-Beryllene_bb`): 20-atom supercell A1 = 5a1, A2 = a1 + a2, displaced along the soft eigenvector; PBE, no D3, SIGMA 0.05, unit-cell-equivalent 65 × 65 k.

| Amplitude | Max. displacement | E0 − E0(A0) |
|---|---|---|
| A0 | 0 | 0 |
| A1 | 0.048 Å | +0.58 meV |
| A2 | 0.097 Å | +12.45 meV |

- The energy rises, so the undistorted 2H-β does not lower its energy by a charge-density-wave distortion along this mode.
- The rise is strongly anharmonic: the ratio is 21 instead of 4. A fit E = aA² + bA⁴ gives a quadratic term close to zero (slightly negative); any double well would be only ~0.02 meV deep, far below k_BT at 1 K.
- The mode is therefore nearly zero-frequency in the harmonic sense. A finite-displacement phonon calculation may show it as a small imaginary frequency, and an electron–phonon calculation must converge it with care.
- The full PBE phonon set for 2H-β (12 displacements, 250 GB each) is running. A first attempt with 100 GB was killed for lack of memory.

## 4. Consequences

- **Manuscript and thesis.** Keep the conditional wording of the 2026-09-26 report (§5 there): "Z₂-nontrivial isolated band subspace of a metal", never "topological insulator". Replace the hydrogen adsorption and cohesive energies. Report phonon stability from the PBE recomputation, on paths through all inequivalent TRIM, stating that the ZA pocket of 1H-β is an interpolation artefact.
- **Superconductivity collaborator** (brief in `note_superconductivity`):
  - Use the 6.0 geometries with PBE without D3, and treat the structures as nonmagnetic without SOC.
  - Converge the 2H-β mode at q = (0.4, 0.6), which is nearly zero-frequency and may dominate λ.
  - Do not claim topological superconductivity: the band inversions lie 1.2–2.7 eV above E_F. This corrects the "coexistence of superconductivity and topological aspects" reading of Li et al., Mater. Today Phys. 38, 101257 (2023), and the collaborator draft, which called α-beryllene a topological insulator.
  - 2H-α is an insulator and is not superconducting unless doped.

## 5. Novelty assessment and next steps

Ranked by weight. Most results so far are corrections rather than discoveries.

1. **Most promising, not yet established: the nearly zero-frequency 2H-β mode at q = (0.4, 0.6).**
   - It is k-mesh sensitive (Section 3.2) and strongly anharmonic (frozen phonon: +0.58 / +12.45 meV, ratio 21 instead of 4). Hydrogenated β-beryllene therefore sits at the edge of a lattice (charge-density-wave-type) instability.
   - If this mode carries a large electron–phonon coupling λ, the story "hydrogenation drives β-beryllene to the edge of a lattice instability and enhances superconductivity" has real physical content, comparable to strongly anharmonic modes in other phonon-mediated superconductors.
   - Checks:
     - (a) **Done: no nesting.** The constant-matrix-element bare susceptibility χ₀(q) and nesting function ξ(q) were computed from the unfolded 105 × 105 SOC SCF eigenvalues (bands 5–8). At q = (0.4, 0.6), χ₀ ranks 8700–8900 of 11 024 q points and ξ ranks 3800–8400, for σ = 0.02, 0.05 and 0.10 eV. χ₀ peaks instead at small q ≈ (0.11, 0.10). The softening is therefore not nesting-driven. Together with its k sensitivity, it most likely comes from a q-dependent electron–phonon matrix element (as for the charge-density wave of NbSe₂, where nesting is also weak) or from the lattice itself. Only the collaborator's electron–phonon calculation can tell these apart. Figure: `figures/9_topology/9.7_lindhard_b_bb.pdf`.
     - (b) **Running:** frozen-phonon energies at σ = 0.02 and 0.10 eV, plus a smaller amplitude (0.024 Å) at σ = 0.05 eV to fix the sign of the harmonic term (jobs 43606–43614, `3.6_phonon_soft_mode_check/b-Beryllene_bb`). A CDW-type mode softens further as σ decreases.
     - (c) **Running:** full PBE phonons (array 43597[]). The collaborator's electron–phonon calculation (λ_q at this q) remains the decisive test.
2. **Solid, but a correction: topology and superconductivity are decoupled.** α and ST are metals. Their Z₂ band inversions lie 1.2–2.7 eV above E_F and are protected only by ~1 meV spin–orbit gaps, so they are unrelated to the Fermi surface that sets T_c (Section 4). This suits one section of a combined paper.
3. **Worth stating, limited novelty.**
   - Full hydrogenation turns α-beryllene into a trivial wide-gap insulator (2H-α, PBE gap 4.84 eV), analogous to graphene → graphane.
   - The corrected hydrogenation energies are small: −0.08 to −0.33 eV per H vs ½H₂, i.e. −0.66 eV per H₂ for 2H-α. 2H-α has the BeH₂ composition (≈ 18 wt% H), so a remark on reversible hydrogen uptake is possible. However, BeH₂ is well studied and Be is toxic, so this should not be a selling point.
   - The Z₂ classification of the hydrogenated structures is trivial in every case and serves completeness.
4. **Not physics novelty.** The IrRep orthogonality false alarm (single-precision normalisation plus the PAW metric) and the corrected phonon settings belong in the methods or supplement.

**Suggested paper focus:** hydrogenation-tuned superconductivity and the near-unstable 2H-β mode, with the topology section as a clarification of earlier claims. Because nesting is ruled out, the mode is interesting only if its electron–phonon coupling turns out large. This should be settled before it is made a headline result.

## 6. Pending work and limitations

| Item | Status |
|---|---|
| 2H-β PBE phonons (`3.7_phonon_dispersion_with_hydrogen_pbe/b-Beryllene_bb`, array 43597[]) | running; dispersion and stability table will be added to the notebook |
| Pristine α, β, ST phonons | not re-audited in this campaign (no D3 check done) |
| β-beryllene relaxation used an 11 × 11 k mesh (others 27 × 27) | recorded in `readme.md`; a denser-k re-relaxation is advisable because β is metallic |

Remaining limitations of the topology screening (unchanged):
1. No continuous-BZ gap proof.
2. No cutoff or density convergence study of the meV gaps.
3. PBE only; hybrid functionals could reorder bands separated by about 1 meV (α, ST) or 0.18 eV (2H-β at M).
4. Supercell magnetic orders are not sampled.
5. No edge-state calculation.

## Files and figures

| Item | Location |
|---|---|
| Results notebook | `9.0_topology.ipynb` |
| Topology and stability figures | `figures/9_topology/9.1–9.6_*.pdf` (9.5 and 9.6: PBE phonons of 2H-α and 1H-β) |
| Topology methods, job records | `9.0_topology/workflow.md`, `9.0_topology/tools/README.md` |
| Phonon audit | `3.5_phonon_dispersion_with_hydrogen_selected/audit_20260927.md` |
| Frozen phonon | `3.6_phonon_soft_mode_check/b-Beryllene_bb/` |
| PBE phonons | `3.7_phonon_dispersion_with_hydrogen_pbe/` |
| Reference energies | `1.0_convergence/single_Hydrogen_spin`, `1.0_convergence/H2_molecule`, `readme.md` |
| Collaborator brief | `note_superconductivity` |
