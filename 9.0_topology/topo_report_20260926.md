# SOC band topology of beryllene and hydrogenated beryllenes

Report, 2026-09-26. Supersedes the topology claims of the collaborator draft `20260922_draft/topo_draft.md` for the current reference geometries. Results and figures are read by `../9.0_topology.ipynb`; methods, job records and troubleshooting are in [`workflow.md`](workflow.md).

> **Status: conditional screening, not certified.** Every Z₂ value below belongs to a fixed lowest-N spinor-band subspace and is conditional on electronic time-reversal (TR) symmetry and on isolation of that subspace over the whole 2D Brillouin zone (BZ). Five of the six structures are metallic at neutral filling, so their indices are not Fermi-level quantum-spin-Hall (QSH) invariants. All reports carry `topology_certified: false`.

## 1. Summary

| Structure | Space group (1e-5 Å) | N | Method | Conditional Z₂ | Min. sampled direct gap E(N+1)−E(N) | Neutral Fermi level |
|---|---|---:|---|:--:|---|---|
| α | Cmmm (P6/mmm at 1e-2 Å) | 2 | Fu–Kane parity | **ν = 1** | 1.127 meV at the K-point Dirac point (0.33333, 0.33333) | metal (indirect −4.03 eV) |
| β | P-1 | 4 | Fu–Kane parity | ν = 0 | 0.854 meV near (0.3284, 0.3137); 0.874 meV near (0.3234, 0.3619) | metal (−3.52 eV) |
| ST (square trilayer) | P4/mmm | 6 | Fu–Kane parity | **ν = 1** | 0.735 meV at (0.31935, 0.31935) and its C₄ images; 1.282 meV at M | metal (−3.75 eV) |
| 2H-α | C2/m (P-3m1 at 1e-2 Å) | 4 | Fu–Kane parity | ν = 0 | 4.925 eV | **insulator (+4.84 eV)** |
| 2H-β | P-1 | 6 | Fu–Kane parity | ν = 0 | 0.177 eV | metal (−6.15 eV) |
| 1H-β | P1 (no inversion) | 4 / 6 | Wilson loop (WCC) | Z₂ = 0 / 0 | 0.929 eV / 0.415 eV | metal (odd filling, 5 e) |

N is the number of spinor bands in the subspace (N = NELECT for the even-electron structures). "Indirect" is min E(N+1) − max E(N) over all sampled k. A negative value means bands N and N+1 both cross E_F.

- **The only insulating case is 2H-α, and it is trivial.** Its gap is +4.84 eV in PBE and its ν = 0 is a genuine Fermi-level classification, still conditional on TR stability (Section 4.3).
- **α and ST are ν = 1, but only as isolated subspaces of metals.** The lowest two (α) and six (ST) bands form Z₂-nontrivial subspaces, separated from the next band everywhere sampled by spin–orbit direct gaps of 1.13 meV (α) and 0.73 meV (ST). They are not QSH insulators.
- **β, 2H-β and 1H-β are trivial** in the analysed subspaces.

## 2. Computational methods

First-principles calculations used VASP 6.5.0 [1] with PAW potentials [2,3] (PAW_PBE Be 06Sep2000, 2s² valence; H 15Jun2001) and the PBE functional [4].

- **Geometries.** The six reference structures were taken unchanged from the project's frozen `6.0_dielectric_selection` CONTCARs. They were not re-relaxed and no dispersion correction was applied in the topology runs.
- **Cell and electrostatics.** c = 40 Å, with a dipole correction along z (IDIPOL = 3) [11].
- **Electronic settings.**
  - Plane-wave cutoff 600 eV, PREC = Accurate, LASPH, EDIFF = 10⁻⁸ eV.
  - Gaussian smearing of 0.01 eV and 28 spinor bands.
  - Noncollinear SOC with SAXIS = (0, 0, 1) and all initial magnetic moments set to zero.
- **Self-consistent step.** A Γ-centred 105 × 105 × 1 mesh (ISYM = 2).
- **Fixed-density stages.** All later runs reused that charge density (ICHARG = 11, ISYM = −1, no symmetry imposed):
  - the four 2D TRIM Γ, X = (½, 0), Y = (0, ½) and M = (½, ½) of each primitive cell;
  - the existing band path (384 k);
  - two local refinement levels around the sampled extrema;
  - an adaptive zoom into every unresolved direct-gap basin (Section 4.4).

**Parity eigenvalues and ν.** For the five centrosymmetric structures, parities were obtained with IrRep 2.1.3 [5] directly from the SOC spinor WAVECAR. IrRep evaluates the characters of the actual inversion {−1 | (0, 0, 0.2)} (centre at z = 4 Å) in the DFT cell, then separates the states explicitly by inversion eigenvalue. The Z₂ index follows the Fu–Kane criterion [6]:

$$
\delta_i=\prod_{m=1}^{N/2}\xi_{2m}(\Lambda_i),\qquad (-1)^{\nu}=\prod_{i=1}^{4}\delta_i ,
$$

Here ξ₂ₘ(Λᵢ) = ±1 is the common parity of the m-th Kramers pair at TRIM Λᵢ. All four TRIM are computed explicitly, because in the current Cmmm and P-1 cells X, Y and M are not all equivalent by symmetry. For α they correspond to the three M points of the ideal hexagonal zone.

**Independent parity check.** Every parity is cross-checked by an independent double-precision reconstruction of the inversion matrix from the WAVECAR (Section 4.1).

**1H-β (no inversion).** Z₂ was computed from the evolution of Wannier charge centres (WCC) [7,8]. Overlaps came directly from SOC-DFT through the VASP–Wannier90 3.1.0 interface [9] and were fed to Z2Pack 2.2.1 on the half-BZ surface [t, s/2, 0] (23 Wilson loops). The line-position, gap and movement convergence checks all passed.

## 3. Parity analysis and conditional topological classification

**Table 1.** Inversion parity of each Kramers pair at the four TRIM, ordered from band 1 upward. The entry after "|" is the next pair above band N, with the direct gap E(N+1) − E(N) at that TRIM in parentheses. δᵢ is the product over the N/2 pairs.

| Structure | N | Γ | X | Y | M | δ (Γ, X, Y, M) | (−1)^ν | ν |
|---|---:|---|---|---|---|---|:--:|:--:|
| α | 2 | + \| − (5.10 eV) | − \| + (6.37 eV) | − \| + (6.37 eV) | − \| + (6.37 eV) | +, −, −, − | −1 | 1 |
| β | 4 | +, − \| + (2.81 eV) | +, − \| − (5.24 eV) | +, − \| − (5.95 eV) | +, − \| − (5.52 eV) | −, −, −, − | +1 | 0 |
| ST | 6 | +, −, + \| − (2.22 eV) | +, −, + \| − (5.10 eV) | +, −, + \| − (5.10 eV) | −, −, + \| + (1.28 meV) | −, −, −, + | −1 | 1 |
| 2H-α | 4 | +, − \| + (9.79 eV) | −, + \| − (5.64 eV) | −, + \| − (5.64 eV) | −, + \| − (5.64 eV) | −, −, −, − | +1 | 0 |
| 2H-β | 6 | +, −, + \| − (8.37 eV) | −, +, + \| − (1.30 eV) | +, −, − \| + (1.30 eV) | −, +, − \| + (0.177 eV) | −, −, +, + | +1 | 0 |

- **α.** The single occupied Kramers pair is even at Γ and odd at the three zone-boundary TRIM. This gives (−1)^ν = δ_Γ δ_X δ_Y δ_M = −1: the ordering of a Kane–Mele-type [10] band inversion. The two bands are separated from band 3 everywhere sampled. The smallest separation is the SOC gap at the Dirac point at K (1.127 meV).
- **β.** δ = −1 at all four TRIM, so ν = 0.
- **ST.** δ_M = +1 against δ = −1 at Γ, X and Y, so ν = 1. At M the pair above band 6 (bands 7–8, 1.28 meV higher) has the same parity (+) as bands 5–6. Swapping that near-degenerate ordering therefore cannot change δ_M or ν. The global minimum of the direct gap (0.735 meV) lies at a SOC-gapped crossing on the Γ–M diagonal, not at a TRIM.
- **2H-α.** Trivial, with every gap at the TRIM in the eV range.
- **2H-β.** Trivial. At every TRIM the pair above band 6 has the opposite parity. ν therefore depends on the level ordering, most weakly at M, where that pair lies 0.177 eV above.

**1H-β (Wilson loops).** The lowest-4 and lowest-6 subspaces both give converged Z₂ = 0. The minimum sampled direct gaps are 0.929475 eV and 0.41504 eV (gap_refine_2). With five electrons the Fermi level lies inside the half-filled band-5/6 Kramers manifold. Both subspaces are isolated sub-manifolds of a metal, not an occupied manifold.

## 4. Validation of the conditions behind the indices

### 4.1 Wavefunction quality: the IrRep orthogonality flags

Parity job 42146 accepted α and marked β, ST, 2H-α and 2H-β `PARITY_NOT_VALIDATED`. The only trigger was IrRep's diagnostic line `orthogonality (largest of diag. <psi_nk|psi_mk>): X > 1e-5`. It has two separate sources, and neither indicates poor wavefunctions.

1. **Single-precision normalisation.**
   - WAVECAR stores complex64 coefficients. IrRep normalises them and forms ⟨ψ|ψ⟩ in complex64 over about 1.1 × 10⁴ plane waves.
   - The resulting diagonal error is 4.4 × 10⁻⁶ to 4.5 × 10⁻⁵. It disappears (≤ 10⁻¹⁴) when the same coefficients are normalised in float64.
   - α passed only because its error happened to be 9.5 × 10⁻⁶.
2. **PAW metric.**
   - PAW pseudo-wavefunctions are orthonormal in S = 1 + Σ|p⟩q⟨p|, not in the plain plane-wave metric.
   - The remaining off-diagonal overlaps (up to 1.55 × 10⁻² in 2H-β) occur only between same-parity bands at different energies.
   - Overlaps between opposite-parity bands are ≤ 4.0 × 10⁻⁸, as they must be because S commutes with inversion.
   - An independent reconstruction of S from the POTCAR projectors makes all 28 bands S-orthonormal to ≤ 1.3 × 10⁻⁷.

**Convergence and consistency.**
- The TRIM runs converged in 7–8 Davidson steps (final ΔE ≤ 6.7 × 10⁻⁹ eV).
- Eigenvalues of bands 1..N+2 agree across stages to 1 µeV.
- IrRep's own inversion matrices use the Gram inverse of the non-orthogonal basis: their eigenvalues are ±1 to 10⁻¹⁵, with off-diagonal elements ≤ 3.8 × 10⁻⁷.

**Revised validation and rerun.** The parity validation therefore no longer fails on this plain-metric message; the message is recorded. It is replaced by stricter float64 tests at every TRIM:

| Test | Requirement | Worst value found |
|---|---|---|
| Löwdin-orthonormalised inversion matrix on bands 1..N | unitary, block-diagonal and ±1 to 10⁻⁶ | unitarity 2.5 × 10⁻¹¹ |
| Opposite-parity overlap over bands 1..N+2 | ≤ 10⁻⁶ | 4.0 × 10⁻⁸ |
| Per-band inversion residual ‖Iψ − pψ‖ | ≤ 10⁻⁵ | 6.4 × 10⁻⁶ |
| Odd-state counts | identical to IrRep's | identical |

All other IrRep warnings still fail closed, and no tolerance was relaxed. The rerun needed no new VASP calculation and passed for all five structures (`parity_batch_57b2275a8b154208a7148a579825f89d.json`).

### 4.2 Band and k-point selection

- **Cells and k-points.** spglib confirms that all cells are primitive (1, 2, 3, 3 and 4 atoms). The WAVECAR k-points are exactly the four 2D TRIM of those cells.
- **Bands.** Only Be 2s² and H 1s electrons are in the valence, so the lowest N bands equal NELECT; no semicore states are involved.
- **Inversion operation.** The operation used by IrRep coincides with the geometric inversion found independently (residual ≤ 3 × 10⁻¹³ Å).

### 4.3 Electronic time reversal

**The converged SOC states are TR-invariant within numerical precision** (`<structure>/tr_evidence.json`).

| Quantity | Largest value found |
|---|---|
| Magnetisation density max\|m(r)\| | ≤ 2.6 × 10⁻⁷ μB/Å³ (≤ 3.4 × 10⁻⁷ of the maximum charge density) |
| ∫\|m\| dV | ≤ 8.9 × 10⁻⁷ μB |
| TR-odd real part of the one-centre PAW magnetisation | ≤ 3 × 10⁻⁷ |
| Kramers splitting of bands 1..N+2 at the TRIM | ≤ 3.8 × 10⁻⁶ eV |

- The ≤ 3 × 10⁻⁷ TR-odd part compares with a TR-even spin–orbit part of 7–9 × 10⁻⁵.
- **1H-β is the most direct test.** It has no inversion and its SCF ran on the full, unsymmetrised 105 × 105 grid. There E_n(k) = E_n(−k) holds to 4.0 × 10⁻⁷ eV for bands 1–8 over 5512 ±k pairs, while spin–orbit splittings at the same k reach 7 × 10⁻⁴ eV.

**These data cannot test whether a magnetic state has lower energy.** Every SCF started from m = 0, and for the centrosymmetric structures ISYM = 2 removed inversion-odd magnetisation by construction. The instability is most plausible for 1H-β, whose neutral Fermi level sits in a half-filled Kramers band, although its DOS at E_F is modest (0.36 states/eV/cell).

**Stability screen submitted (jobs 43030–43035, results pending).**
- Collinear spin-polarised SCFs on the same geometry, mesh and smearing.
- Seeds: ferromagnetic (Be 1 μB, H 0.5 μB); layer-alternating Be (inversion-odd for β, 2H-β and 1H-β); and for ST also the inversion-odd (+, 0, −) pattern.
- Collapse of every seed would support, but not prove, a nonmagnetic ground state; supercell orders are not sampled.

### 4.4 Isolation of the fixed-band subspace over the 2D BZ

The 105 × 105 SCF grid was unfolded with the symmetry VASP actually used, which reproduces every IBZ weight. Every local minimum of E(N+1) − E(N) below 0.3 eV was then examined.

- **Earlier refinements.** They followed only the global minimum of each stage, which left three basins unresolved.
- **Adaptive zoom.** Each basin was refined with fixed-density SOC NSCF 9 × 9 patches, re-centred on the minimum and shrunk fourfold whenever that minimum was interior.
- **Slope bound.** Each level reports the minimum sampled gap minus (largest neighbour slope) × (sampling radius). This is a heuristic local bound: a positive value supports a nonzero gap but is not a proof.

| Structure/N | Basin | Levels | Minimum direct gap | Final spacing (fractional) | Slope bound | Status |
|---|---|---|---|---|---|---|
| β/4 | B1, start (0.32338, 0.36190) | 1.013 → 1.013 → 0.877 → 0.874 → 0.874 meV | 0.8737 meV at (0.3233725, 0.3619238) | 5.8 × 10⁻⁷ | +0.869 meV | converged, interior |
| β/4 | B2, start (0.3283, 0.3135) | 6.09 → 4.33 → 1.17 → 0.882 → 0.854 meV | 0.8541 meV at (0.3284209, 0.3137278) | 4.7 × 10⁻⁶ | +0.781 meV | converged, interior |
| ST/6 | M (TRIM) | 13.8 (grid) → 1.282 → 1.282 meV | 1.282 meV | — | pinned at M | converged |
| ST/6 | diagonal D, start (0.31936, 0.31936) | 0.865 → 0.865 → 0.865 → 0.735 → 0.735 meV | 0.7348 meV at (0.3193507, 0.3193507), ×4 by C₄ | 2.3 × 10⁻⁶ | +0.715 meV | converged, interior (never refined before) |
| α/2 | K = (⅓, ⅓) | 1.156 → 1.156 → 1.156 → 1.127 meV | 1.1272 meV at (0.3333287, 0.3333310) | 2.3 × 10⁻⁶ | +1.068 meV | converged, interior |
| 2H-α/4, 2H-β/6, 1H-β/4, 1H-β/6 | all basins | — | 4.925 eV, 0.177 eV, 0.929 eV, 0.415 eV | — | whole-grid Lipschitz bound > 0 | isolated |

**Corrections.** The earlier β minimum of 8.57 meV was not converged: β's lowest four bands are separated from band 5 by two SOC-gapped Dirac-like points with gaps of 0.854 and 0.874 meV. ST's minimum is 0.735 meV on the Γ–M diagonal, not 1.28 meV at M. Both values are consistent with the ~1 meV spin–orbit scale of Be. Both subspaces remain isolated in the sampled sense, but by less than 1 meV. The α gap of 1.127 meV sits within 5 × 10⁻⁶ (fractional) of K, as expected for its nearly hexagonal Cmmm cell.

**Scope of the gap claims.** These are sampled, fixed-density PBE gaps. They do not prove a gap on the continuous BZ, and they come with no ENCUT or charge-density convergence claim. Across stages, SCF and fixed-density eigenvalues differ by a nearly rigid 0.1–1 meV offset, while direct gaps change by ≤ 30 µeV. Gaps of about 1 meV are also far below the accuracy of PBE band ordering, so the ν = 1 subspaces of α and ST should be described as SOC-gapped band inversions of semimetal/metal character. Their nontrivial character could be sensitive to the functional.

### 4.5 Fermi level

Only 2H-α has E_F inside a global gap: max E₄ lies 0.065 eV below E_F and min E₅ lies 4.78 eV above it.

| Structure | Bands crossing E_F | Indirect overlap |
|---|---|---|
| α | 1–4 | −4.03 eV |
| β | 1–6 | −3.52 eV |
| ST | 5–8 | −3.75 eV |
| 2H-β | 5–8 | −6.15 eV |
| 1H-β | 3–6 | — |

For 1H-β the neutral E_F lies 3.00 eV above min E₅ and 3.33 eV below max E₆, with a hole pocket of 83 meV in bands 3–4. Odd filling combined with TR excludes an insulator.

## 5. Comparison with the 2026-09-22 draft

The draft (legacy calculations in `20260922_draft/`) used PBE-D3, 520 eV, a 24 × 24 mesh, 15 Å vacuum, irvsp, and relaxed hexagonal cells: α in P6/mmm with a = 2.120 Å, β in P-3m1 with a = 2.143 Å. The current campaign uses the project's frozen reference geometries: α with a = 2.129 Å (hexagonal to 3 × 10⁻³°) and β distorted to P-1 (a = 2.147 Å, b = 2.186 Å, γ = 119.8°).

| Item in the draft | Current result | Required change |
|---|---|---|
| α: ν = 1, "non-trivial Z₂ topological insulator" (citing [12]) | ν = 1 for the lowest Kramers-pair subspace, but α is metallic (indirect −4.03 eV); SOC direct gap 1.13 meV at K | Replace "topological insulator" with a conditional Z₂-nontrivial isolated band subspace in a metal |
| β: ν = 0, "topologically trivial" | ν = 0 (δ = −1 at all four TRIM); metallic; subspace gap about 0.85 meV | Keep ν = 0 and state that it is conditional and metallic |
| ST: ν = 1, "nontrivial Z₂ topological phase" | ν = 1 (δ_M = +1); metallic; minimum direct gap 0.735 meV on the Γ–M diagonal | Same qualification as α |
| Table 1 lists state-level parities (e.g. α at M "+, −") | Kramers partners share a parity; listing per pair gives α M = (−) | Use pair parities at all four TRIM |
| Γ and M only, with δ_M³ or δ_X² equivalences | X, Y and M are not all symmetry-equivalent in the current cells; all four TRIM computed | Use the four-TRIM product |
| Reference [7] "Phys. Rev. B 76, 45302" | Phys. Rev. B 76, 045302 | Correct the article number |
| — | 2H-α: trivial insulator, ν = 0, gap 4.84 eV; 2H-β ν = 0; 1H-β WCC Z₂ = 0 / 0 | New results |

The legacy ν values for α, β and ST coincide with the current ones. The draft's insulating and certified wording does not follow from either set of calculations.

**Suggested wording** (keeps the qualifiers): *"Within PBE+SOC, the lowest Kramers-pair subspace of α-beryllene and the lowest six spinor bands of square-trilayer beryllene carry a conditional Fu–Kane index ν = 1. Both systems are metallic at neutral filling, and the nontrivial subspaces are separated from the next band only by spin–orbit gaps of about 1 meV (1.13 and 0.73 meV, sampled on dense and locally refined k meshes). The index therefore characterises a band inversion rather than a quantum-spin-Hall insulating state, and is conditional on electronic time-reversal symmetry and on the isolation of the subspace; it has not been certified as a global topological invariant."*

## 6. Pending calculations and remaining limitations

| Job | Purpose | Status (2026-09-26) |
|---|---|---|
| 43036 | β zoom into basins B1 and B2 | finished (Section 4.4) |
| 43037 | ST zoom into the unrefined diagonal basin | finished (Section 4.4) |
| 43038 | α zoom at K (local slope bound) | finished (Section 4.4) |
| 43030–43035 | collinear spin-polarised stability screens (α, β, ST, 2H-α, 2H-β, 1H-β) | queued |

Limitations that remain after these jobs:

1. No continuous-BZ gap proof.
2. No convergence study of the meV-scale gaps with respect to cutoff and density.
3. PBE only; hybrid functionals or GW could reorder bands separated by ≲ 1 meV (α, ST) or 0.18 eV (2H-β at M).
4. Supercell magnetic orders are not sampled.
5. No edge-state calculation, which is only meaningful after a subspace is established as robustly nontrivial.

## References

[1] G. Kresse and J. Furthmüller, Efficient iterative schemes for *ab initio* total-energy calculations using a plane-wave basis set, Phys. Rev. B **54**, 11169 (1996).

[2] P. E. Blöchl, Projector augmented-wave method, Phys. Rev. B **50**, 17953 (1994).

[3] G. Kresse and D. Joubert, From ultrasoft pseudopotentials to the projector augmented-wave method, Phys. Rev. B **59**, 1758 (1999).

[4] J. P. Perdew, K. Burke, and M. Ernzerhof, Generalized gradient approximation made simple, Phys. Rev. Lett. **77**, 3865 (1996).

[5] M. Iraola, J. L. Mañes, B. Bradlyn, M. K. Horton, T. Neupert, M. G. Vergniory, and S. S. Tsirkin, IrRep: Symmetry eigenvalues and irreducible representations of *ab initio* band structures, Comput. Phys. Commun. **272**, 108226 (2022).

[6] L. Fu and C. L. Kane, Topological insulators with inversion symmetry, Phys. Rev. B **76**, 045302 (2007).

[7] A. A. Soluyanov and D. Vanderbilt, Computing topological invariants without inversion symmetry, Phys. Rev. B **83**, 235401 (2011).

[8] D. Gresch, G. Autès, O. V. Yazyev, M. Troyer, D. Vanderbilt, B. A. Bernevig, and A. A. Soluyanov, Z2Pack: Numerical implementation of hybrid Wannier centers for identifying topological materials, Phys. Rev. B **95**, 075146 (2017).

[9] G. Pizzi *et al.*, Wannier90 as a community code: new features and applications, J. Phys.: Condens. Matter **32**, 165902 (2020).

[10] C. L. Kane and E. J. Mele, Z₂ topological order and the quantum spin Hall effect, Phys. Rev. Lett. **95**, 146802 (2005).

[11] J. Neugebauer and M. Scheffler, Phys. Rev. B **46**, 16067 (1992).

[12] J. Li, M. Guo, J. Si, L. Shi, X. Shi, J.-J. Ma, Q. Zhang, D. J. Singh, P.-F. Liu, and B.-T. Wang, Coexistence of superconductivity and topological aspects in beryllenes, Mater. Today Phys. **38**, 101257 (2023).
