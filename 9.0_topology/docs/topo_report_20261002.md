# SOC band topology and stability of beryllene and hydrogenated beryllenes

Report, 2026-10-02. Update of [`topo_report_20261001.md`](topo_report_20261001.md). Results are displayed in [`9.0_topology.ipynb`](../../9.0_topology.ipynb); the analysis code is `vmatplot/band_topology.py`.

> **Status: conditional screening, not certified.** Every Z₂ value belongs to a fixed lowest-N spinor-band subspace and is conditional on electronic time-reversal symmetry and on isolation of that subspace over the whole 2D Brillouin zone (BZ). All reports carry `topology_certified: false`.

## What changed since 2026-10-01

1. **β-beryllene recomputed in the P-3m1 geometry.** The topology, HSE06, PBE phonons, bands, PDoS and dielectric run are now redone. The conclusions do not change: β is a metal with a trivial conditional ν = 0, and it is dynamically stable. Its smallest subspace gap is now an SOC gap pinned at K (Section 1).
2. **Hydrogenation thermodynamics with zero-point energy and free energy.** Only 2H-α stays stable against H₂ release at room temperature and 1 bar. 1H-β becomes endothermic once ZPE is included. 2H-β needs about 90 bar of H₂ (Section 2).
3. **2H-α is not the lowest BeH₂ monolayer.** A puckered square layer with tetrahedral Be is 38 meV per BeH₂ lower (Section 3).
4. **One figure package (2026-10-04).** All figures are in `figures/`, as article (`<name>.pdf`) and thesis (`<name>_thesis.pdf`) versions; the energetics and phonon analyses moved from `9.0_topology.ipynb` to the notebooks of their topics (Section 6).

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

The new data are now used in notebooks 0.0–8.0 and `9.0_topology.ipynb`. Previous versions of every replaced PDF are in git (commit `c3acdcda` and earlier).

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
| 9.0 topology | HSE06 orange, PBE blue (9.8); phonon paths as above (9.5, 9.10–9.12) | 9.1–9.3, 9.5, 9.8, 9.10–9.12 |

**Other changes:**
- Paths whose case differs from the directory names (e.g. `4.0_bandstructure/b-beryllene`) were corrected on executable lines, so the notebooks also run on a case-sensitive file system. Savefig names now match the existing PDF names exactly.
- Figures whose data did not change were restored byte for byte after execution.

**Review fixes (after an independent check of the notebooks):**

| Notebook | Problem | Fix | Figures |
|---|---|---|---|
| 4.1 band gaps | 1H-β and 2H-β read `4.0_bandstructure/b-Beryllene_h`, `_hh` (no OUTCAR) | read the `*_ollie` band runs used for fig3.10 | — |
| 4.3 bands with H | the 2H-β cell saved into the 1H-β working file `4.32_BS_b_b.pdf` | 2H-β now saves to `4.33_BS_b_bb.pdf` | 4.32, 4.33 |
| 5.0 DOS | in the hexagonal panel the α-bulk and α-monolayer directories were swapped | swapped back | 5.11 |
| 5.1 DOS with H | pristine β still read the P-1 run | `5.0_PDoS/b-Beryllene_p3m1` | fig3.12a, 5.30 |
| 5.3 PDoS analysis | `b-beryllene_h`, `_hh` do not exist | `b-Beryllene_b`, `b-Beryllene_bb` | — |
| 5.4 PDoS | the β panel (y limit 1.2) clipped the P-3m1 peak of 1.25 states/eV | y limit 1.4 | 5.12 |
| 8.0 AIMD | `time_step = 1.0` although POTIM = 0.5 fs: the time axis was doubled and only the first half of the 5.5 ps runs was shown | `time_step = 0.5` | fig3.7a1, fig3.7b1, fig3.7c1; 8_aimd/* |
| 1.x, 2.2, 3.1–3.5, 5.1, 5.5, 8.0 | lowercase paths and savefig names that only work on a case-insensitive disk | case corrected | — |

All notebooks 0.0–9.0 now execute without errors on a case-sensitive file system.

**Recomputed after the check:**
- **Bulk phonons (S3.13a, S3.13c).** They used PBE+D3, and the "hcp" run `3.0_phonon_dispersion/a-Beryllium_3` is the bcc supercell: its POSCAR and DYNMAT are identical to `c0-Beryllium_3`, so S3.13a showed bcc force constants on the hexagonal path. They are recomputed with PBE in `3.8_phonon_dispersion_pbe/a-Beryllium` (hcp, 4 × 4 × 3, job array 44153[]) and `c0-Beryllium` (bcc, 4 × 4 × 4, job array 44154[]); see `3.8_phonon_dispersion_pbe/bulk_note`.
  - **bcc is harmonically unstable:** −3.90 THz at N, exact because N is commensurate with the supercell; the minimum on a 20³ mesh is −3.72 THz. bcc is the high-temperature phase of Be and is stabilised only by anharmonicity, so this is expected. The old S3.13c was plotted from 0 THz, which would hide the instability; the new figure starts at −5 THz.
  - **hcp is stable over the whole BZ:**
    - maximum 20.96 THz, the upper Γ optical mode;
    - Γ optical modes at 13.86 and 20.96 THz (the Raman E2g mode of Be is ≈ 13.7 THz, 456 cm⁻¹);
    - residual force 3.7 × 10⁻⁴ eV/Å.
  - **Figures:** S3.13a (hcp, Γ–K–M–Γ–A–H–L–A–Γ) and S3.13c (bcc, Γ–H–N–Γ–P–H|P–N, y from −5 THz) are regenerated with notebook 3.0. `3.0_phonon_dispersion_alt.ipynb` (S3.13c `_alt`) now uses the same PBE bcc run.
- **β dielectric convergence tests (S3.7, S3.8).** Redone in the P-3m1 geometry (`6.0_dielectric_function/b-Beryllene_p3m1_*`, jobs 44155–44163, all exit 0).
  - The legend shows the NBANDS that VASP actually used with 24 ranks: 72, 96, 144, 192 and 264.
  - The k series is converged at 75 × 75 and 105 × 105.
  - In P-3m1 the xx and yy components coincide, as required by the hexagonal symmetry.
- **`exported_figures/`.** 26 of its 34 PDFs are replaced by the current versions: the 24 optics figures that contain β, plus S3.7 and S3.8. The README and manifest record the refresh; its β numbers and validation arrays are marked as superseded.

**Not changed:**
- The 1H-β and 2H-β band structures (fig3.10) come from the `*_ollie` runs. Their geometry matches the references, and D3 does not change eigenvalues.

**Full regeneration (end of 2026-10-02):**
- All 35 notebooks and the two figure-merge notebooks (`figures/0_structure_selected/merge.ipynb`, `figures_for_publication/fig3.7/merge.ipynb`) were re-executed with the current data and library. None has an error output, and every figure in `figures/`, `figures_for_publication/` and `figures_for_thesis_ch3/` was regenerated.
- `3.0_phonon_dispersion_alt.ipynb` now plots the PBE bcc run with `vmatplot.phonon`; the module `vmatplot.phonon_alt` no longer exists.
- `3.5_phonon_dispersion_with_hydrogen_phonopy.ipynb` shows the PBE runs of `3.7_phonon_dispersion_with_hydrogen_pbe`.
- The AIMD panels in `figures_for_publication/fig3.7/` are synchronised with the corrected fig3.7a1/b1/c1.
- All 34 PDFs in `exported_figures/` match their sources.
- **Removed (history in git):**
  - `3.5_phonon_dispersion_with_hydrogen_selected_backup.ipynb`, which would have overwritten fig3.5/fig3.6 with old data if run;
  - the orphan PBE+D3 figures `S3.13c_phonon_bcc.pdf`, `3.20_phonon_c-Beryllium.pdf`, `3.32_phonon_b_b_alt.pdf`, `3.33_phonon_b_bb_alt1.pdf` and `3.33_phonon_b_bb_alt2.pdf`;
  - all `_old` copies.
- **Not regenerable here:** the structure renders (VESTA PNGs in `figures/0_structure_selected/`, `figures_for_publication/fig3.7/`) are made in VESTA by hand. The β render still shows the P-1 cell; it differs from P-3m1 by 2% in one lattice length and 0.2° in angle, which is invisible at figure scale, but it should be re-rendered from `2.0_geometry_optimization/b-Beryllene_p3m1/CONTCAR` for exactness.

## 5. Pending work

| Item | Status |
|---|---|
| Anharmonic phonons of 2H-β (SSCHA) | not submitted: python-sscha and QE are not installed, and the 5 × 5 supercell with the dense k mesh the soft mode needs costs ~150 core-h per configuration, ~10⁴–10⁵ core-h in total. One-mode estimate from the frozen-phonon fit (2026-10-05, `summarize_frozen_phonon_quantum`, notebook 3.5): about 4 THz at 0 K and 5.5–5.7 THz at 300 K for all three SIGMA; a static distortion would need α ≈ −45 to −50 meV amu⁻¹ Å⁻² (≈ 5i THz harmonic) against −2.1 to −6.2 from the fits. Mode–mode coupling is not included. |
| Phonons of the puckered BeH₂ layer | optional; only if the polymorph is discussed beyond an energy comparison |

The limitations of the topology screening are unchanged from the 2026-10-01 report.

## 6. Figure package (2026-10-04)

The figure folders `figures/<n>_<topic>/`, `figures_for_publication/`, `figures_for_thesis_ch3/` and `exported_figures/` are replaced by one flat folder, `figures/`. The file names in Section 4 refer to the old folders; the table at the end of this section maps them to the new names.

- **Names.** `<quantity>_<structure or set>.pdf`, for example `phonon_2H-beta-beryllene.pdf` or `optics_pristine_absorption.pdf`. The names carry no figure numbers; the article and the thesis number their figures themselves.
- **Two versions.** Every notebook that saves figures sets `figure_version("article")` in its prework cell. The article version is the default and is saved as `figures/<name>.pdf`. With `figure_version("thesis")` the same notebook saves `figures/<name>_thesis.pdf` in the thesis typography of `o-B14_20241024` (the rcParams of its thesis notebook): titles 24 pt, axis labels and ticks 18 pt, legend 12 pt, title pad 10 pt, TrueType fonts. This switch replaces the separate thesis blocks of `0.0_thesis_demo.ipynb`. Those blocks gave the previous thesis copies of the α/β optics (thesis 3.1, 4.1–4.5) larger legends (18 pt) and axis labels (20 pt); all thesis figures now share the `o-B14_20241024` values, which are set in one place, `thesis_params` in `vmatplot/output_settings.py`.
- **Layouts that depend on the version.**
  - `topology_bands_pristine`, `topology_bands_hydrogenated`: article 1 × 3 with the legend below; thesis 2 × 2 with the legend in the fourth panel.
  - `topology_direct-gap_pristine`: article 1 × 3 with the colour bar at the right; thesis 2 × 2 with the colour bar and the marker legend in the fourth panel.
  - `topology_direct-gap_hydrogenated` has four panels, because 1H-β has two subspaces (N = 4 and N = 6): article 1 × 4, thesis 2 × 2 with the colour bar at the right.
  - `optics_pristine_*`, `optics_hydrogenated_*` (xx, yy, zz): article 1 × 3, each panel with its legend; thesis 2 × 2 with the shared legend in the fourth panel, as for the optical rows of `o-B14_20241024`.
  - All other figures keep their layout in both versions; only the typography changes. As in `o-B14_20241024`, the two-panel rows and the 2 × 3 dielectric grids are not rearranged.
- **Text that would not fit.** The thesis sizes are applied to the article canvases. Where a title would then leave the canvas, or rotated tick labels would overlap (the k-mesh labels of the convergence plots), `save_figure` scales that text down until it fits (`fit_thesis_text` in `vmatplot/output_settings.py`); the long convergence titles end up at about 17–18 pt. The article version is saved unchanged.
- **Code.** `vmatplot/output_settings.py` has `figure_version()` (returns or sets the version) and `save_figure(name)` (writes the PDF without a creation date); `canvas_setting()` applies the thesis typography. `plot_topology_bands` and `plot_direct_gap_maps` (`vmatplot/band_topology.py`) and `plot_linear_optical_property` (`vmatplot/linear_optical_properties.py`) choose their grid from the version.
- **Inputs that are not computed** are in `figures/sources/`:
  - `structure_overview/` (bulk, pristine and hydrogenated grids: `structure_bulk`, `structure_pristine`, `structure_hydrogenated`) and `structure_single/` (top and side view of one hydrogenated layer: `structure_2H-alpha-beryllene`, `structure_1H-beta-beryllene`, `structure_2H-beta-beryllene`): VESTA files, renders and a merge notebook each. The merge notebooks have their own `figure_version` switch and write into `figures/`.
  - `structure_library/`: the full set of VESTA renders.
  - `brillouin_zone/`: GIMP source of the earlier raster sketch. The Brillouin-zone sketch is now drawn as vector graphics by `figures/brillouin_zone.ipynb` (`brillouin_zone.pdf`, `brillouin_zone_thesis.pdf`). Its right panel shows the path Γ–X–M–Γ–Y–M′–Γ of the 1H-β and 2H-β phonons and bands (Section 7).
- **Moved out of `9.0_topology.ipynb`**, which now shows only the topology screening:
  - hydrogenation thermodynamics and BeH₂ polymorphs → `2.2_geometry_information_with_hydrogen.ipynb`;
  - phonon stability tables → `3.0_phonon_dispersion.ipynb` (pristine) and `3.5_phonon_dispersion_with_hydrogen_selected.ipynb` (hydrogenated, with the 2H-β frozen phonon and the nesting analysis);
  - HSE06 gap of 2H-α → `4.1_bandgap_with_hydrogen.ipynb`.

  The PBE phonon and HSE06 band figures of 9.0 duplicated those of 3.0, 3.5 and 4.3 and were dropped. The tables show only the current runs: the superseded P-1 β cell and the PBE+D3 phonon runs are left out. The cubic trilayer is labelled "cubic trilayer" instead of "ST".
- **Check.** The 19 figure notebooks and the two merge notebooks ran in both versions without errors: 92 article PDFs, each with its `_thesis` twin. A full rerun of both versions reproduces all 184 PDFs byte for byte. Of the 144 earlier figure files that have a successor, 143 keep their page size; the exception is the nesting figure, whose panels were brought next to their colour bars.
- **Content check against the notebooks before the reorganisation.** The notebooks of the previous commit (3.0, 5.1, 5.4, 6.0, 6.1, 8.0 and the two merge notebooks) were run unchanged into a scratch folder: all 57 outputs match the new article figures pixel for pixel, so the reorganisation changed no figure content. Of the 137 earlier files with a successor, 92 are identical to the new article figures, `topology_gap-zoom` differs on purpose, and the remaining 44 differ because the committed PDFs (commit `c3acdcda`, the same files as `figures_old/`) predate the corrections of 2026-10-02 or were rendered on the Mac. Manuscripts that still use those files carry outdated figures, in particular the old S3.13a (hcp phonons computed in the bcc supercell), S3.13c, fig3.12a, S3.3–S3.12 and fig3.7a1–c1; use the files in `figures/` instead.
- **Printed size in the thesis.** The thesis text width is 15.2 cm (measured on the chapter pages in `.agent/manuscript_previews/`). With the insertion fractions of `o-B14_20241024` (0.6 `\textwidth` for 10 × 6 in plots, 0.8 for two-panel rows and 2 × 2 grids, 1.0 for the 24 × 12 in grids), tick labels print at 4.5–7.3 pt and legends at 3.0–5.1 pt, the same range as the `o-B14_20241024` thesis figures (4.3–7.2 pt and 2.9–4.8 pt). The parity labels of `topology_bands_*` stay at 12 pt in the thesis version: at 14 pt neighbouring labels at X of the cubic trilayer overlap.
- **Stale tools.** `.agent/sync_manuscripts.py` and `.agent/manuscript_sync.json` still copy from `exported_figures/`, which no longer exists.

| New name (`figures/`) | Publication / thesis name before | Working copy before |
|---|---|---|
| `dielectric_alpha-beta` | `fig3.13_dielectric`, thesis `3.1_dielectric` | — |
| `optics_alpha-beta_absorption` | `fig3.14a_absorption`, thesis `4.1_absorption` | — |
| `optics_alpha-beta_refractive` | `fig3.14b_refractive`, thesis `4.2_refractive` | — |
| `optics_alpha-beta_extinction` | `fig3.15_extinction`, thesis `4.3_extinction` | — |
| `optics_alpha-beta_reflectivity` | `S3.16_reflectivity`, thesis `4.4_reflectivity` | — |
| `optics_alpha-beta_energy-loss` | `S3.17_energy-loss`, thesis `4.5_energy-loss` | — |
| `convergence_energy_encut_g-beryllene` | — | `1_convergence/1.1_conv_g_encut` |
| `convergence_energy_kpoints_g-beryllene` | — | `1_convergence/1.2_conv_g` |
| `convergence_energy_kpoints_hcp-beryllium` | `S3.2a_conv_hcp` | `1_convergence/1.3_conv_a0` |
| `convergence_energy_kpoints_alpha-beryllene` | `S3.2b_conv_alpha` | `1_convergence/1.4_conv_a1` |
| `convergence_energy_kpoints_bcc-beryllium` | `S3.2c_conv_bcc` | `1_convergence/1.5_conv_c0` |
| `convergence_energy_kpoints_cubic-beryllene-monolayer` | — | `1_convergence/1.6_conv_c1` |
| `convergence_energy_kpoints_cubic-beryllene-bilayer` | — | `1_convergence/1.7_conv_c2` |
| `convergence_energy_kpoints_cubic-beryllene-trilayer` | `S3.2d_conv_cubic` | `1_convergence/1.8_conv_c3` |
| `convergence_cohesive_encut_g-beryllene` | — | `1_convergence/2.1_coh_g_encut` |
| `convergence_cohesive_kpoints_g-beryllene` | — | `1_convergence/2.2_coh_g_kpoints` |
| `convergence_cohesive_kpoints_hcp-beryllium` | `S3.1a_conv_hcp` | `1_convergence/2.3_coh_a0` |
| `convergence_cohesive_kpoints_alpha-beryllene` | `S3.1b_conv_alpha` | `1_convergence/2.4_coh_a1` |
| `convergence_cohesive_kpoints_bcc-beryllium` | `S3.1c_conv_bcc` | `1_convergence/2.5_coh_c0` |
| `convergence_cohesive_kpoints_cubic-beryllene-monolayer` | — | `1_convergence/2.6_coh_c1` |
| `convergence_cohesive_kpoints_cubic-beryllene-bilayer` | — | `1_convergence/2.7_coh_c2` |
| `convergence_cohesive_kpoints_cubic-beryllene-trilayer` | `S3.1d_conv_cubic` | `1_convergence/2.8_coh_c3` |
| `phonon_hcp-beryllium` | `S3.13a_phonon_hcp` | `3_phonon/3.10_phonon_a-Beryllium` |
| `phonon_alpha-beryllene` | `S3.13b_phonon_alpha` | `3_phonon/3.11_phonon_a-Beryllene`, `9_topology/9.10_phonon_a_pbe` |
| `phonon_beta-beryllene` | `S3.13d_phonon_beta` | `3_phonon/3.12_phonon_b-Beryllene`, `9_topology/9.11_phonon_b_pbe` |
| `phonon_bcc-beryllium` | `S3.13c_phonon_bcc_alt2` | `3_phonon/3.20_phonon_c-Beryllium_alt2` |
| `phonon_cubic-beryllene-trilayer` | `fig3.4_phonon_c-Beryllene` | `3_phonon/3.21_phonon_c-Beryllene_trilayer`, `9_topology/9.12_phonon_st_pbe` |
| `phonon_2H-alpha-beryllene` | `fig3.5_phonon_a_hh` | `3_phonon/3.31_phonon_a_hh`, `9_topology/9.5_phonon_a_hh_pbe` |
| `phonon_1H-beta-beryllene` | `fig3.6a_phonon_b_b_alt` | `9_topology/9.6_phonon_b_b_pbe` |
| `phonon_2H-beta-beryllene` | `fig3.6b_phonon_b_bb_alt` | `9_topology/9.9_phonon_b_bb_pbe` |
| `bands_hcp-beryllium` | `S3.15a_bs_hcp` | `4_bs/4.10_BS_a_bulk` |
| `bands_alpha-beryllene` | `S3.15b_bs_alpha` | `4_bs/4.11_BS_a` |
| `bands_beta-beryllene` | `S3.15d_bs_beta` | `4_bs/4.12_BS_b` |
| `bands_bcc-beryllium` | `S3.15c_bs_bcc` | `4_bs/4.20_BS_c_bulk` |
| `bands_cubic-beryllene-trilayer` | `fig3.8_bs_cubic-trilayer-Beryllene` | `4_bs/4.21_BS_c` |
| `bands_2H-alpha-beryllene` | `fig3.9_bs_a-Beryllene_hh` | `4_bs/4.31_BS_a_hh`, `9_topology/9.8_hse_bands_a_hh` |
| `bands_1H-beta-beryllene` | `fig3.10a_bs_b-Beryllene_b_alt` | — |
| `bands_2H-beta-beryllene` | `fig3.10b_bs_b-Beryllene_bb_alt` | — |
| `dos_overview` | — | `5_dos/5.10_DoS_overview` |
| `dos_hexagonal-family` | — | `5_dos/5.11_DoS_hex` |
| `dos_cubic-family` | `fig3.11_DoS_cubic` | `5_dos/5.12_DoS_cubic` |
| `dos_hydrogenated` | `fig3.12a_Total_DoS` | `5_dos/5.30_DoS_hydrogen` |
| `pdos_hcp-beryllium` | — | `5_pdos/5.10_a-Beryllium` |
| `pdos_alpha-beryllene` | — | `5_pdos/5.11_a-Beryllene` |
| `pdos_beta-beryllene` | — | `5_pdos/5.12_b-Beryllene` |
| `pdos_bcc-beryllium` | — | `5_pdos/5.20_c0-Beryllium` |
| `pdos_cubic-beryllene-trilayer` | — | `5_pdos/5.21_c3-Beryllene` |
| `pdos_2H-alpha-beryllene` | `fig3.12b_a_hh` | `5_pdos/5.31_a_hh` |
| `pdos_1H-beta-beryllene` | `fig3.12c_b_b` | `5_pdos/5.32_b_h` |
| `pdos_2H-beta-beryllene` | `fig3.12d_b_bb` | `5_pdos/5.33_b_hh` |
| `convergence_dielectric_nbands_hcp-beryllium` | `S3.3_dielec_hcp` | `6_dielectric_test/6.10_dielec_a-Beryllium_nbands` |
| `convergence_dielectric_kpoints_hcp-beryllium` | `S3.4_dielec_hcp` | `6_dielectric_test/6.10_dielec_a-Beryllium_kpoints` |
| `convergence_dielectric_nbands_alpha-beryllene` | `S3.5_dielec_alpha` | `6_dielectric_test/6.11_dielec_a-Beryllene_nbands` |
| `convergence_dielectric_kpoints_alpha-beryllene` | `S3.6_dielec_alpha` | `6_dielectric_test/6.11_dielec_a-Beryllene_kpoints` |
| `convergence_dielectric_nbands_beta-beryllene` | `S3.7_dielec_beta` | `6_dielectric_test/6.12_dielec_b-Beryllene_nbands` |
| `convergence_dielectric_kpoints_beta-beryllene` | `S3.8_dielec_beta` | `6_dielectric_test/6.12_dielec_b-Beryllene_kpoints` |
| `convergence_dielectric_nbands_bcc-beryllium` | `S3.9_dielec_bcc` | `6_dielectric_test/6.20_dielec_c-Beryllium_nbands` |
| `convergence_dielectric_kpoints_bcc-beryllium` | `S3.10_dielec_bcc` | `6_dielectric_test/6.20_dielec_c-Beryllium_kpoints` |
| `convergence_dielectric_nbands_cubic-beryllene-trilayer` | `S3.11_dielec_cubic` | `6_dielectric_test/6.21_dielec_c-Beryllene_trilayer_nbands` |
| `convergence_dielectric_kpoints_cubic-beryllene-trilayer` | `S3.12_dielec_cubic` | `6_dielectric_test/6.21_dielec_c-Beryllene_trilayer_kpoints` |
| `dielectric_pristine` | `fig3.16_dielec` | `6_dielectric/6.10_dielec` |
| `dielectric_hexagonal-family` | — | `6_dielectric/6.11_dielec_hex` |
| `dielectric_cubic-family` | — | `6_dielectric/6.12_dielec_cubic` |
| `dielectric_hydrogenated` | `fig3.18_dielec_H` | `6_dielectric/6.20_dielec_with_Hydrogen` |
| `dielectric_alpha-family` | — | `6_dielectric/6.21_dielec_a` |
| `dielectric_beta-family` | — | `6_dielectric/6.22_dielec_b` |
| `optics_pristine_absorption` | `fig3.17a_abs` | `7_optical/7.11_absorption` |
| `optics_pristine_refractive` | `S3.18a_refractive` | `7_optical/7.12_refractive` |
| `optics_pristine_extinction` | `S3.18c_extinction` | `7_optical/7.13_extinction` |
| `optics_pristine_reflectivity` | `S3.18b_reflectivity` | `7_optical/7.14_reflectivity` |
| `optics_pristine_energy-loss` | `fig3.17b_energy-loss` | `7_optical/7.15_energy-loss` |
| `optics_hydrogenated_absorption` | `fig3.19a_abs_H` | `7_optical/7.21_absorption_H` |
| `optics_hydrogenated_refractive` | `S3.19a_refractive_H` | `7_optical/7.22_refractive_H` |
| `optics_hydrogenated_extinction` | `S3.19c_extinction_H` | `7_optical/7.23_extinction_H` |
| `optics_hydrogenated_reflectivity` | `S3.19b_reflectivity_H` | `7_optical/7.24_reflectivity_H` |
| `optics_hydrogenated_energy-loss` | `fig3.19b_energy-loss_H` | `7_optical/7.25_energy-loss_H` |
| `aimd_2H-alpha-beryllene` | `fig3.7a1` | `8_aimd/a-Beryllene_hh` |
| `aimd_1H-beta-beryllene` | `fig3.7b1` | `8_aimd/b-Beryllene_b` |
| `aimd_2H-beta-beryllene` | `fig3.7c1_aimd` | `8_aimd/b-Beryllene_bb` |
| `structure_bulk` | `fig3.1_bulks_view` | — |
| `structure_pristine` | `fig3.2_layers_view` | — |
| `structure_hydrogenated` | `fig3.3_layers_H_view` | — |
| `structure_2H-alpha-beryllene` | `fig3.7a2` | — |
| `structure_1H-beta-beryllene` | `fig3.7b2` | — |
| `structure_2H-beta-beryllene` | `fig3.7c2` | — |
| `brillouin_zone` | `S3.14` | — |
| `topology_bands_pristine` + `topology_bands_hydrogenated` | — | `9_topology/9.1_topology_bands` |
| `topology_direct-gap_pristine` + `topology_direct-gap_hydrogenated` | — | `9_topology/9.2_direct_gap_maps` |
| `topology_gap-zoom` | — | `9_topology/9.3_gap_zoom` |
| `topology_wcc_1H-beta-beryllene` | — | `9_topology/9.4_wcc_beta_1h` |
| `nesting_2H-beta-beryllene` | — | `9_topology/9.7_lindhard_b_bb` |

## 7. Corrections (2026-10-06)

- **Brillouin-zone sketch.** The right panel of `brillouin_zone.pdf` now shows the path Γ–X–M–Γ–Y–M′–Γ that the phonons and the bands of 1H-β and 2H-β use, with the reciprocal vectors **b**₁ and **b**₂ (drawn at 0° and 60°): X = **b**₁/2, Y = **b**₂/2 and M′ = (**b**₂ − **b**₁)/2 lie on the zone boundary; M = (**b**₁ + **b**₂)/2 lies outside the first zone on the line Γ–K and is equivalent to M′. The earlier labels X, S, X₁, Y and the arrows **a**, **b** matched none of the computed paths.
- **Band paths of 1H-β and 2H-β.** The `*_ollie` paths were written for a cell with 60° between **a**₁ and **a**₂; in the actual cells (117.0° and 132.7°) their "K" is a general point inside the zone. New non-self-consistent runs on the phonon path: `4.0_bandstructure/b-Beryllene_b_trim`, `b-Beryllene_bb_trim` (jobs 45819, 45820, exit 0). Notebooks 4.1 and 4.3 read them; details in `4.0_bandstructure/note`.
- **Fermi level of the PBE band plots.** A non-self-consistent line-mode run obtains E-fermi from the path k-points only. The band plots used that value, which is off by 0.05–0.08 eV (hcp, α, β, cubic trilayer), 0.40 eV (bcc), 0.58 eV (1H-β) and 0.77 eV (2H-β). `plot_bandstructure` now takes the self-consistent run as the last entry of each `bs_list` (`extract_reference_fermi`, `vmatplot/bandstructure.py`), and notebook 4.2 and 4.3 pass `2.4_self_consistent/<structure>` and `2.5_self_consistent_with_hydrogen/<structure>`. The 2H-α plot is aligned at the VBM and is unchanged. The DOS and PDoS runs are non-self-consistent too, but on 105 × 105 (bulk 105³) meshes over the whole zone, so their E-fermi is reliable; the topology plots use the self-consistent value.
- **Topology band figure.** `extract_band_path` names the TRIM on the paths of the distorted cells as in the parity tables: the point (½, 0) of 1H-β and 2H-β is now X instead of M. K = (1/3, 1/3) is the zone corner of the hexagonal parent and lies close to, not at, the corner of these cells (for 2H-β 12% beyond it along Γ–M); the article caption says so.
- **β render.** `figures/sources/structure_overview/2b_b-beryllene.vesta` and its two PNGs now hold the P-3m1 cell; `structure_pristine` (both versions) is regenerated. The renders were made here with VESTA 3.5.8 (`VESTA-gtk3.tar.bz2`) under a user-space Xvfb unpacked from the Rocky 9 RPMs (`xorg-x11-server-Xvfb`, `libXfont2`, `libxkbfile`, `libXdmcp`, `libfontenc`, `xkbcomp`; the hard-coded `/usr/bin` of xkbcomp in the Xvfb binary patched to `./xkbin`). Same canvas as before (VESTA.ini window 1345 × 893 gives a 1031 × 618 canvas; `VESTA -open <file> -export_img scale=4 <png> -close`) and same views: side view as stored in the file with translation 1.785 and scale 4.300, top view along **c** with **a** horizontal and scale 5.160. Rendering the old P-1 file this way reproduces the old PNGs to a mean pixel difference of 1–2 (out of 765), so only the structure changed. `structure_library/` still holds the P-1 β renders; no figure uses them.
- **Article text** (`../H-Beryllene_article_2026`): α in-plane lattice constant 6.1% (not 6.0%) below hcp; the cubic trilayer, built from three bcc (001) layers, has an in-plane lattice constant 12.7% below the cubic lattice constant of bcc and 0.8% above its nearest-neighbour distance, and an interlayer distance 20.2% above the bcc (001) spacing (the old "reduced by 17%" is not reproduced by any reference). `optics_update.md` and its manifest now list the current file names.

## Files

| Item | Location |
|---|---|
| β P-3m1 chain | `*/b-Beryllene_p3m1` in 2.0, 2.4, 3.8, 4.0, 4.4, 5.0, 6.0_dielectric_selection, 9.0_topology |
| H₂ frequency | `1.0_convergence/H2_molecule_freq` |
| BeH₂ polymorphs | `2.1_geometry_optimization_with_hydrogen_a/BeH2_square_*`, `BeH2_square_note` |
| Energetics tables | `summarize_hydrogenation_energetics`, `summarize_polymorph_energies` in `vmatplot/band_topology.py` |
| Updated values | `readme.md` (β section, ZPE / ΔG sub-bullets, phonon lines) |
