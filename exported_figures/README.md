# Corrected optical figures and audit

This delivery contains the 34 figures requested in `.agent/ai_problem.md`: 33 regenerated PDFs and one independently checked, unchanged PDF. All twelve A figures and all twelve B figures were regenerated. Nine of the ten C convergence figures needed correction; `S3.3_dielec_hcp.pdf` is byte-identical to the archived original. The original filenames and the two source layouts are retained.

Use `figures_for_thesis/` for the six thesis-layout figures and `figures_for_publication/` for the remaining 28 figures. The former maps to the external Thesis `figures_ch3/` directory and the latter to `figures_proj3/`. `manifest.json` gives the exact per-file mapping, generating notebook and cell, status, and old/new SHA-256. The files here were subsequently synchronized to the local article and Thesis repositories on 2026-09-21, together with the corresponding TeX corrections; each repository contains an `optics_update_manifest.json` receipt.

The complete source changes are in the parent repository. This delivery was produced from base commit `f863ab02ebf93d5804f5552967ea750571778335` plus the current uncommitted corrections. No new commit or public archive release has been created. Source hashes, actual software versions and input hashes are recorded in `validation/`.

## Findings and corrections

The three central findings in the feedback are supported by the source, the original HDF5 arrays, the structures and the old PDF curves:

1. The beta thickness incorrectly added a fractional coordinate span to a length. The same error also affected cubic trilayer and all three hydrogenated structures.
2. The absorption call supplied ordinary frequency to an angular-frequency expression, losing exactly one factor of `2*pi` at fixed dielectric response.
3. Fixed y limits concealed real energy-loss peaks, notably the alpha zz maximum of about 8.878 at 10.637 eV. Additional clipping occurred in some convergence and comparison plots.

`vmatplot/structure.py` now reads the actual relaxed geometry for each structure, converts coordinates to lengths, projects the full cell height along the surface normal, and unwraps a slab across a periodic boundary. Bulk factors remain `(1, 1)`. The original surface-centre/radius convention is retained: Be 1.98 Angstrom, H 0.70 Angstrom; the radii of the outermost atoms define the two endpoints. No unique physical thickness is claimed. In particular, the H/H convention can make the effective 2H-alpha thickness smaller than the pristine-alpha value.

Corrected effective thicknesses are alpha 3.960000, beta 5.855087, cubic trilayer 6.966490, 2H-alpha 3.144512, 1H-beta 5.523544 and 2H-beta 5.192831 Angstrom. Full precision is retained in computation. All 58 selected/convergence HDF5 structures match their CONTCAR identities; the maximum coordinate discrepancy is about 1.78e-15 Angstrom.

`energy_to_frequency` still returns Hz. The absorption function explicitly accepts Hz and converts to rad/s once; an explicitly named rad/s alternative is also available. `frequency_to_energy` now multiplies by h. Its old inverse-conversion bug was confirmed, but no actual notebook output was found to call it. The main, merged and backup plotting APIs were checked. Their wavelength branches also had a separate error: they plotted a raw dielectric component instead of the requested derived quantity and mixed nm/eV window filtering. Those branches are now covered by independent numerical regression tests.

The six computation notebooks were run in cell order against the repaired library and their embedded plot outputs were refreshed. The formula notebook now has the correct reflectivity equation and explicit angular-frequency, imaginary-part and unit definitions. Export checks include the actual visible samples and interpolated x-window boundaries. New PDFs disable path simplification; the unchanged historical S3.3 PDF retains its original simplification and is checked against its own tick calibration.

## Numerical results

The four corrected absorption maxima over the available 0-12 eV samples are:

- Alpha xx: 0.021562 -> 0.135476 nm^-1.
- Alpha zz: 0.027854 -> 0.175010 nm^-1.
- Beta xx: 0.031402 -> 0.159298 nm^-1.
- Beta zz: 0.046811 -> 0.241453 nm^-1.

Alpha's other optical quantities are unchanged. Beta's `(epsilon1 - 1)` and `epsilon2` scale by about 0.6844266293, but n, k, R and loss were recomputed nonlinearly from the corrected epsilon. The beta xx absorption ratio is therefore not constrained to be `2*pi`.

Beta xx/zz refractive maxima are about 4.902274/4.662807; extinction maxima 3.149900/4.000426; reflectivity maxima 0.537446/0.766017. Alpha zz loss reaches 8.878109, while beta zz reaches 1.001770. Beta xx reaches 1.506761 at the final available sample near 11.9903 eV while still rising; this boundary value is not a resolved peak.

`validation/numeric_audit.json` records complete old/new maxima, sample positions, local extrema and model parameters for xx, yy and zz of all eight selected materials. `validation/corrected_spectra_0-12eV.npz` contains the actual unrounded arrays. Sampling is around 0.02 eV; the displayed numerical precision is not a physical accuracy estimate.

## File-by-file status

A: all six physical quantities in both layouts were regenerated from the original density-density response. Both layouts have identical numerical curves and x/y limits.

B: both dielectric comparison figures and all ten derived-property comparison figures were regenerated. Each actual xx, yy and zz component was used, with each structure's own factor; yy was never substituted by xx.

C: the ten convergence figures were individually audited:

- S3.3 hcp NBANDS: unchanged, including its original file hash.
- S3.4 hcp k points: expanded axes to display the original large low-energy response.
- S3.5 alpha NBANDS: expanded axes to include zz extrema.
- S3.6 alpha k points: expanded axes in its original 4-8 eV window.
- S3.7/S3.8 beta: corrected the thickness for every calculation and the selected reference.
- S3.9 bcc NBANDS: corrected the actual total-band label from 98 to 126.
- S3.10 bcc k points: expanded axes and removed a duplicate entry for the same 27-cubed calculation.
- S3.11/S3.12 cubic: corrected the thickness for every calculation and the selected reference.

The corresponding auxiliary PDFs under `figures/6_dielectric`, `figures/6_dielectric_test` and `figures/7_optical` were also refreshed. Across the source directories there are 60 checked exports, including the two aliases of the unchanged hcp NBANDS figure. The handoff folder intentionally contains only the requested 34 PDFs.

## Evidence and reproducibility

From the repository root, with Python and the packages in `requirements-optics.txt`:

```sh
python -m pytest -q tests
python scripts/verify_optics_exports.py
python scripts/regenerate_optics.py
python scripts/audit_optics.py
```

The active verification environment used for this delivery is `/tmp/h-beryllene-optics-venv/bin/python`. It is isolated from the system Python and is disposable. Package versions and Python version are saved in `validation/environment.json`.

The regeneration script executes the real notebook code cells, checks plotted ranges before saving, reopens exported PDF vector paths, and packages the existing filenames. It refuses to overwrite an unchanged C figure if its hash no longer matches Git. Baseline reproduction requires a separate checkout/extracted source directory containing the original `vmatplot` and `9.0_thesis_demo.ipynb`; `--baseline-source` writes to an isolated baseline folder and cannot overwrite this delivery. The old source and PDFs are recoverable from the base Git commit; the local old-model reproduction is under `.agent/optics_audit/baseline/`.

Evidence files include:

- `validation/legacy_pdf_verification.json`: all twelve original A PDFs, 56 curves and 13,102 visible retained vertices agree with the old model; maximum page residual about 0.00009831 pt. The PDF coordinate transform comes only from its ticks, not a fit to the arrays.
- `validation/figure_data.json.gz`: the arrays actually plotted, their structures, factors, component labels, limits and page coordinates for every source export. All regenerated vector paths were compared with these arrays. The preserved legacy S3.3 file is calibrated from its own ticks, with retained vertices checked within 0.01 pt and omitted samples within 0.25 pt of its simplified polyline.
- `validation/numeric_audit.json`: independently evaluates complex square roots, `alpha=4*pi*k/lambda_vacuum`, Fresnel R and `Im(-1/epsilon)` directly from HDF5, independently checks factors against structures, and binds arrays, source, inputs and all delivered PDFs by hashes.
- `validation/input_hashes_before.json` and `input_hashes_after.json`: unchanged raw calculation files.
- `validation/structure_h5_identity.json`: HDF5/CONTCAR identity checks.
- `validation/delivery_checks.json`: final test, visual-review, file-count and source-provenance checks.

The independent numerical audit covers 1,208 curve instances and 529,116 sample values across the 60 source exports. Its maximum absolute difference is approximately 2.04e-12. All six A layout pairs are numerically identical. These are post-processing/export checks, not estimates of physical DFT error.

## Manuscript context and remaining scope

The supplied 30-page `.agent/draft_2609.pdf` was reviewed alongside the arrays. See `draft_corrections.md` for page-specific English replacements, corrected SI thickness/band entries, and figure references. Two broad conclusions particularly need revision: cubic trilayer absorption is not enhanced relative to bcc throughout the visible region, and hydrogenation does not enhance every alpha polarization. The corrected cubic xx/yy absorption is below bulk bcc across the sampled 1.65-3.26 eV interval. The loss-peak ordering and several hydrogenated dielectric peaks also change.

The draft confirms the assigned-thickness normalization, but refers the numerical thickness convention to an SI that was not supplied. The original notebook comments, rather than the PDF alone, establish the retained Be/H endpoint-radius convention. The supplied draft PDF itself and all original VASP files have been preserved. After the user supplied the article and Thesis repository paths, the corresponding editable TeX and SI corrections were applied there; `draft_corrections.md` retains the original page-specific reasoning.

The actual output parameters for all eight selected materials have `WPLASMAI=0`, zero WPLASMA tensors and no added post-processing Drude term. The results remain the original independent-particle interband response. No new DFT calculation, electron-hole treatment, local-field correction, finite-sheet reflection model or lifetime calculation was introduced. Negative epsilon1, n<1, opposite-sign tensor components and scalar loss peaks alone do not establish propagating modes or lifetimes of a finite sheet. The [VASP LOPTICS documentation](https://vasp.at/wiki/index.php/LOPTICS) and [WPLASMAI documentation](https://vasp.at/wiki/WPLASMAI) explain the approximation and intraband option.

The hcp 45-cubed low-energy outlier is now visible; this correction does not prove low-frequency metallic convergence. Some selected reference curves change NBANDS and k mesh together and should not be described as additional strict one-variable convergence tests. Non-optical results, unrelated library copies, and remote/public archive publication are outside the completed local correction. No messages were sent to external collaborators.
