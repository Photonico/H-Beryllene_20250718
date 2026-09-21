# Optical corrections for draft_2609.pdf

Reviewed source: `.agent/draft_2609.pdf`, dated September 21, 2026, 30 pages.
SHA-256: `73dbd0338fabd9c3dfec3bd0a7f753edeae8d93f753efe1888f9fae371e7334d`.
Page numbers below are the one-based PDF pages and also match the printed page numbers.
The source PDF has not been changed. These are replacement suggestions for the editable manuscript.

The review covers the optical methods and interpretation on pages 6-8 and 17-26, with related statements on pages 1, 4 and 5. It does not validate the separate structural, adsorption-energy or topological claims. The supplied PDF does not include the Supplementary Information (SI).

Numerical evidence comes from `6.0_dielectric_selection/*/vaspout.h5`, specifically `results/linear_response/energies_dielectric_function` and `results/linear_response/density_density_dielectric_function`, and the associated relaxed `CONTCAR` files. Each xx, yy and zz component was evaluated separately. Values below use the corrected effective thickness, the original constants, and `alpha = 4*pi*k/lambda_vacuum`. Peak positions refer to available samples, generally spaced by about 0.02 eV; their displayed precision is not a physical uncertainty estimate.

## 1. Thickness convention: what the supplied draft establishes

Page 8 defines `d_SC` as the full supercell length normal to the layer and `d_2D` as the assigned structure thickness. Equations (13)-(14) use the same volume normalization for the displayed components:

```text
eps1_effective = 1 + (d_SC / d_2D) * (eps1_supercell - 1)
eps2_effective =     (d_SC / d_2D) * eps2_supercell
```

This agrees with the normalization preserved in the corrected project. Neither page 8 nor any other page in the supplied PDF specifies the atomic radii, how hydrogen termination determines the two surfaces, or the actual SI thickness values. Consequently, this PDF alone cannot independently establish the H/H endpoint convention.

The convention is documented by the original project notebooks: the separation of the outermost atomic centres plus the radii of the atoms at those two surfaces, with Be = 1.98 Angstrom and H = 0.70 Angstrom. The original factors explicitly add H+H for double-sided termination and H+Be for single-sided termination. The correction preserves that choice and restores the missing conversion from fractional coordinates to Angstrom. It does not substitute the extent of a union of atomic spheres.

The corrected thicknesses to carry into the SI are:

- Pristine alpha: 3.960000000000 Angstrom; Be/Be surfaces.
- Pristine beta: 5.855086567902 Angstrom; Be/Be surfaces.
- Cubic trilayer: 6.966489780213 Angstrom; Be/Be surfaces.
- 2H-alpha: 3.144512351252 Angstrom; H/H surfaces.
- 1H-beta: 5.523544103908 Angstrom; Be/H surfaces.
- 2H-beta: 5.192831416813 Angstrom; H/H surfaces.

The 2H-alpha effective thickness is smaller than the pristine alpha value because the selected endpoint radii change. This is a consequence of the retained thickness convention, not evidence that hydrogenation physically compresses the entire material into a smaller unique optical thickness. Absolute effective dielectric and absorption comparisons inherit this convention.

Suggested addition on page 8:

> We assign an effective optical thickness equal to the separation of the outermost atomic centres along the surface normal, plus the radii of the atoms at the two surfaces. The radii are 1.98 Angstrom for Be and 0.70 Angstrom for H. Thus, the surface-radius contribution is Be+Be for pristine structures, H+H for double-sided hydrogen termination, and Be+H for single-sided termination. The thickness is evaluated from the relaxed Cartesian geometry, with periodic wrapping removed. The resulting effective response is specific to this thickness convention.

## 2. Confirmed numerical or implementation corrections

### 2.1 Page 5: full supercell height is not vacuum thickness

Original phrase: "a vacuum region of 40 Angstrom is included."

The selected optical structures have a full normal cell height of 40 Angstrom. For a slab of nonzero extent, the empty region is smaller than this. The unrelated 15 Angstrom topological setup on page 6 was not audited here.

Replacement:

> For the two-dimensional optical calculations, the full supercell height normal to the layer is 40 Angstrom, including the slab and the surrounding vacuum.

### 2.2 Pages 6-7: photon energy and angular frequency must be distinguished

Original phrase on page 6: "omega denotes the photon energy."

Equation (8) then uses the angular-frequency formula with `c` in ordinary velocity units. The actual old Python call supplied ordinary frequency `E/h`, without the required `2*pi` conversion. This was an implementation error as well as ambiguous manuscript notation.

Recommended consistent convention: use `E` for photon energy throughout Equations (5)-(14), including the transition delta function in Equation (6), and explicitly define `omega = E/hbar` when needed. Replace Equation (8) by:

```text
alpha_lambda(E) = sqrt(2) * E / (hbar*c)
                 * sqrt(sqrt(eps1_lambda(E)^2 + eps2_lambda(E)^2)
                        - eps1_lambda(E))
               = 2*E*k_lambda(E)/(hbar*c).
```

Suggested text:

> Here E is the photon energy and omega = E/hbar is the angular frequency. With E in eV, hbar in eV s and c in nm/s, the absorption coefficient is reported in nm^-1.

Also label the absorption axes in Figures 10, 13 and 15 with `nm^-1`. The old figures omitted this unit.

For Equation (7), make the energy/frequency convention consistent and state the actual VASP parameter `CSHIFT = 0.1 eV`. Do not equate an unexplained additive `i*eta` in a squared-energy denominator with a quantity in eV. This notation needs an explicitly defined finite-broadening convention, or the usual principal-value zero-broadening relation followed by a separate description of the implemented numerical broadening.

### 2.3 Page 18: beta's opposite-sign dielectric interval

Original phrase: "approximately 1.2 eV to 2.2 eV".

For the corrected xx/zz pair, the low-energy samples with opposite signs extend from 1.323339 to 2.225616 eV. These are sampled interval endpoints, not interpolated zeros.

Replacement:

> In the corrected effective dielectric response, the xx and zz real components have opposite signs over sampled photon energies of approximately 1.32-2.23 eV.

The following claim that this alone establishes direction-dependent plasmonic modes should be replaced as described in Section 3 below.

### 2.4 Pages 18-19 and Figure 10: absorption magnitude

The corrected 0-12 eV sampled maxima are:

- Alpha xx: 0.135476 nm^-1 at approximately 6.46 eV.
- Alpha zz: 0.175010 nm^-1 at approximately 7.75 eV.
- Beta xx: 0.159298 nm^-1 at approximately 7.18 eV.
- Beta zz: 0.241453 nm^-1 at approximately 5.96 eV.

Suggested quantitative replacement for the broad magnitude comparison:

> The effective absorption is strongly anisotropic. Within the 0-12 eV window, the largest sampled coefficients for alpha-beryllene are approximately 0.135 and 0.175 nm^-1 in xx and zz, respectively; the corresponding beta-beryllene values are approximately 0.159 and 0.241 nm^-1. The beta out-of-plane maximum occurs near 5.96 eV.

Do not state that beta exceeds alpha at every energy; the curves cross. Alpha absorption changes by exactly `2*pi` at fixed dielectric response, whereas beta also changes because of the thickness correction.

### 2.5 Page 19: the n < 1 threshold is direction- and material-dependent

Original phrase: "For both structures, the refractive index drops below one above approximately 8 eV".

The first downward crossings are bracketed by the following adjacent corrected samples:

- Alpha xx: 7.4059-7.4260 eV.
- Alpha zz: 7.8073-7.8274 eV; it later crosses back above one around 11.70-11.72 eV.
- Beta xx: 7.3385-7.3586 eV.
- Beta zz: 6.3560-6.3761 eV.

Replacement:

> The first downward n = 1 crossing occurs near 7.42 eV for alpha xx, 7.82 eV for alpha zz, 7.35 eV for beta xx and 6.37 eV for beta zz. These crossings describe the dispersion of the effective optical constants and do not by themselves identify a plasmon resonance.

### 2.6 Page 19 and Figure 11: beta out-of-plane extinction

Original phrase: "then decreases rapidly toward the near-ultraviolet region."

The strongest corrected beta zz extinction peak is approximately 4.000426 at 5.955027 eV. The original direction-of-change description misses this peak.

Replacement:

> The beta out-of-plane extinction spectrum develops its strongest peak near 5.96 eV, where k is approximately 4.00, followed by weaker higher-energy features.

### 2.7 Page 20: visible and ultraviolet loss peaks were misranked

Original phrase: "with the visible-range peak being stronger."

This is contradicted by both the old and corrected beta xx arrays. In the corrected data, the visible local peak is approximately 0.449, while a high-energy local maximum is approximately 1.334 near 11.77 eV. The last available beta xx sample is approximately 1.507 at 11.99 eV and is still rising, so it is not a resolved peak.

Replacement:

> The beta in-plane loss spectrum contains a visible-range local feature of approximately 0.45 and larger high-energy features, including a local maximum of approximately 1.33 near 11.77 eV. It reaches approximately 1.51 at the upper sampled energy of 11.99 eV while still rising; the complete high-energy peak therefore lies outside the resolved window, if the rise continues.

Add to the out-of-plane comparison:

> The alpha out-of-plane spectrum has a dominant loss maximum of approximately 8.88 near 10.64 eV, substantially above the beta maximum of approximately 1.00 near 9.97 eV. The old common vertical limit concealed the alpha peak.

The source alpha zz maximum is unchanged by the thickness correction; the change in its visibility is a plotting correction.

### 2.8 Pages 21-22 and Figure 13: cubic trilayer absorption enhancement is not general

Original phrase: "larger absorption coefficients throughout the visible and near-ultraviolet regions".

This claim does not survive the cubic thickness correction. For an explicitly defined visible interval of 1.65-3.26 eV, the corrected cubic xx/yy absorption is below bulk bcc at every available cubic sample. For comparison, bulk bcc was linearly interpolated to the cubic energy grid. The respective visible maxima are approximately 0.060369 and 0.074598 nm^-1. Cubic zz exceeds bcc over only part of that interval. The statement therefore cannot be retained for all directions or the full interval.

Replacement:

> Reducing bulk bcc beryllium to a cubic trilayer changes the energy distribution and anisotropy of the effective interband response. With the corrected thickness normalization, the trilayer xx and yy absorption is lower than bulk bcc throughout the sampled 1.65-3.26 eV interval, whereas the zz response exceeds the bulk value over part of that interval. The strongest sampled trilayer absorption occurs near 5.40 eV in xx/yy and 5.25 eV in zz, with coefficients of approximately 0.148 and 0.166 nm^-1, respectively. These results show a redistribution of optical response rather than a universal absorption enhancement.

The earlier page 21 phrase "enhanced visible-region optical response" must likewise specify the actual quantity and direction. The corrected cubic epsilon2 maxima are approximately 19.75 at 3.74 eV in xx/yy and 27.32 at 2.71 eV in zz. These dielectric features should not be conflated with absorption maxima, which include the energy prefactor and the nonlinear relation to epsilon.

### 2.9 Page 23 and Figure 14: the 2H-alpha dielectric maximum

Original phrase: "values reaching approximately 60-70 near 5 eV."

The corrected 2H-alpha xx/yy real dielectric maxima are approximately 31.16 at 4.9066 eV. The original value was inflated by the erroneous thickness. The 2H-alpha zz real maximum is approximately 16.95 at 9.64 eV, so the statement also needs a component label.

Replacement:

> In the xx and yy components, 2H-alpha-beryllene has a real dielectric maximum of approximately 31.2 near 4.91 eV under the stated effective-thickness convention. Its strongest imaginary dielectric peak in these directions occurs near 5.59 eV, with a value of approximately 31.2. The out-of-plane response has a different spectral distribution, with its strongest imaginary peak near 9.77 eV.

### 2.10 Pages 23 and 25: hydrogenation of alpha does not enhance every optical direction

Original phrase: "substantially stronger optical absorption in the visible to near-ultraviolet region".

This requires a direction label. In the 1.65-3.26 eV interval, 2H-alpha xx/yy exceeds pristine alpha, but the absolute maximum is only approximately 0.001367 nm^-1. In zz, 2H-alpha is below pristine alpha at every sampled point in this interval; its maximum is approximately 0.000332 nm^-1 compared with approximately 0.000729 nm^-1 for pristine alpha interpolated onto the same grid. A similar directional distinction holds in the explicitly evaluated 3.26-4.13 eV interval.

Replacement:

> Hydrogenation of alpha-beryllene produces a strongly direction-dependent redistribution of the effective interband absorption. The xx/yy coefficients increase in the sampled visible range, but remain small compared with the much stronger ultraviolet features near 5.71 eV. The zz coefficient instead decreases over the sampled 1.65-4.13 eV interval, and its strongest peak shifts to approximately 9.83 eV.

The weak low-energy tails in a calculation with finite numerical broadening should not be described as demonstrating strong visible absorption or a new band-gap threshold.

### 2.11 Page 24: the 2H-beta features differ between xx, yy and zz

Original phrases: "negative values near 6 eV to 6.5 eV" and "a strong broad absorption peak appears around 6 eV to 6.5 eV."

The second phrase occurs in a discussion of epsilon2, not the derived absorption coefficient. Neither phrase describes all three components. Corrected 2H-beta epsilon2 maxima are approximately:

- xx: 14.96 near 6.59 eV.
- yy: 7.68 near 9.26 eV.
- zz: 14.97 near 8.24 eV.

The corresponding first epsilon1 sign-change brackets are 6.6307-6.6499 eV in xx, 9.3136-9.3328 eV in yy and 8.0871-8.1063 eV in zz. Derived absorption maxima occur near 6.75, 9.39 and 8.36 eV, respectively.

Replacement:

> Double-sided hydrogenation of beta-beryllene produces distinct dielectric spectra in all three Cartesian directions. The strongest imaginary dielectric features occur near 6.59 eV in xx, 9.26 eV in yy and 8.24 eV in zz. The first real-part sign changes occur near 6.64, 9.32 and 8.10 eV, respectively. Hydrogen coverage therefore changes the spectral distribution and the in-plane anisotropy; it does not produce a common peak or zero crossing for all directions.

Avoid a general monotonic claim that increasing coverage enhances screening: some corrected components and energy ranges decrease. Also avoid referring to these zero crossings as an "effective plasma frequency" without the additional mode analysis described below.

### 2.12 Pages 25-26: update the optical conclusion

Original claims requiring changes include "enhances visible-to-near-ultraviolet absorption" for the cubic trilayer and the unqualified enhancement following hydrogenation.

Suggested replacement optical paragraph:

> Dimensional reduction and hydrogen functionalization strongly redistribute the calculated effective interband optical response. Pristine alpha- and beta-beryllene exhibit pronounced directional differences. Cubic trilayer beryllene develops distinct low-energy dielectric features, but its corrected absorption does not exceed bulk bcc across all directions or throughout the visible range. Hydrogen functionalization changes the characteristic excitation energies and anisotropy, enhancing some components and reducing others. These comparisons use the stated effective-thickness convention and do not include a Drude intraband contribution.

The abstract on page 1 should also replace the general phrase "enhanced visible-to-ultraviolet absorption" with "direction-dependent redistribution of visible-to-ultraviolet optical response" unless a particular material, component and interval are named.

## 3. Methodological caveats: limits of the existing evidence

These are qualifications of the physical interpretation. They do not mean that a new DFT calculation has been performed or that every existing dielectric feature is erroneous.

### 3.1 Pages 6-8: disclose the actual optical calculation

The actual `vasprun.xml` parameters were checked for all eight selected structures, including compressed XML for bulk hcp. They contain `LOPTICS=T`, `ISPIN=1`, `LSORBIT=F`, `ISMEAR=0`, `SIGMA=0.01`, `CSHIFT=0.1`, `WPLASMAI=0` and nine zero WPLASMA components. The existing post-processing does not add a Drude term. This conclusion is based on output parameters, not an inferred default from a missing HDF5 key.

Suggested addition:

> The optical calculations use the long-wavelength independent-particle interband response with LOPTICS, Gaussian smearing SIGMA = 0.01 eV and complex broadening CSHIFT = 0.1 eV. Spin-orbit coupling is not included in these optical datasets. No Drude intraband contribution is included (WPLASMAI = 0), and no electron-hole or microscopic local-field correction is added. The resulting spectra do not represent the complete low-frequency response of the metallic phases.

Actual total NBANDS values are 66 for bulk hcp, 70 for bulk bcc, 133 for each pristine alpha/beta/cubic trilayer, and 140 for each of the three hydrogenated structures. These are total bands, not counts of unoccupied bands. Requested INCAR values can differ because VASP adjusts band counts for parallel execution. [VASP NBANDS](https://vasp.at/wiki/NBANDS), [LOPTICS](https://vasp.at/wiki/index.php/LOPTICS) and [WPLASMAI](https://vasp.at/wiki/WPLASMAI) document these distinctions.

### 3.2 Pages 7-8: effective response and reflectivity

Original phrase: "restore the intrinsic dielectric function".

The manuscript cites the volume-renormalization method of [Yang and Gao](https://pubs.rsc.org/en/content/articlelanding/2021/nr/d1nr04896a/unauth). The implemented correction retains its equations and the project's same normalization for xx, yy and zz. No alternate inverse-dielectric out-of-plane model has been substituted.

Suggested qualification:

> We report an effective dielectric response under the specified optical-thickness convention. Equation (11) gives the normal-incidence scalar reflectivity of a homogeneous semi-infinite medium with the selected effective optical constants. It is not a calculation of the full reflection of a finite anisotropic layer on a substrate.

Thickness ambiguity and electromagnetic geometry matter for atomically thin layers; this is supported by [PRB 106, 035431](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.035431) and [Merano, PRA 93, 013832](https://arxiv.org/abs/1509.04136). The repair should not silently change the physical model in an attempt to remove this limitation.

### 3.3 Pages 17-18: zero crossings, indefinite response and plasmon lifetime

Original phrases: "indicates two plasmon excitations", "signals direction-dependent plasmonic modes", and "produce finite plasmon lifetimes".

Negative epsilon1, opposite-sign tensor components, small epsilon2 and peaks in the plotted `Im(-1/epsilon)` can motivate an analysis of collective response. These scalar, zero-wavevector spectra do not establish propagating modes of the finite two-dimensional geometry or determine their lifetimes. The numerical broadening is an input, not a computed scattering lifetime.

Suggested replacement for the interpretive part of the page 17 paragraph:

> The effective dielectric spectra contain sign changes and loss features that differ between the two materials and directions. These features identify candidate energy ranges for further analysis of collective response. Establishing propagating modes in the finite two-dimensional geometry and their lifetimes requires a specified electromagnetic or wavevector-dependent response calculation and a physically justified treatment of damping; these are not determined here.

Suggested sentence following the corrected interval on page 18:

> This is an opposite-sign region of the assigned effective tensor; it is not by itself a demonstration of a propagating hyperbolic or plasmonic mode in an isolated finite layer.

### 3.4 Pages 19-25: avoid causal inferences from optical curves alone

Apply the same qualification to page 19's association of n < 1 with a plasmon resonance, page 20's identification of a decrease in epsilon2 with a "plasmonic optical response", page 21's "onset of ... plasmonic behaviour", page 22's "low-energy plasmonic peak", and pages 23-25's claims of strengthened collective oscillations or increased effective plasma frequency after adsorption.

Use "real-part sign change", "loss feature" or "interband spectral redistribution" where that is what the calculated quantity directly demonstrates. Assignments to particular surface bands or Be-H hybridized transitions require band- and matrix-element-resolved evidence. They can be presented as hypotheses consistent with separately established electronic-structure results, rather than conclusions established by peak magnitude alone.

The same wording restriction applies to the abstract (page 1), the final introduction sentence (page 4) and the conclusion (page 26). A suitable general phrase is "anisotropic dielectric and energy-loss response".

## 4. Figures and SI to replace alongside the text

Use the regenerated files in this delivery and preserve their existing names:

- Figure 9, page 18: `figures_for_publication/fig3.13_dielectric.pdf`.
- Figure 10, page 19: `fig3.14a_absorption.pdf` and `fig3.14b_refractive.pdf`.
- Figure 11, page 20: `fig3.15_extinction.pdf`.
- Figure 12, page 21: `fig3.16_dielec.pdf`.
- Figure 13, page 22: `fig3.17a_abs.pdf` and `fig3.17b_energy-loss.pdf`.
- Figure 14, page 24: `fig3.18_dielec_H.pdf`.
- Figure 15, page 25: `fig3.19a_abs_H.pdf` and `fig3.19b_energy-loss_H.pdf`.

The SI referenced on pages 5 and 8 is absent from the supplied PDF. Its thickness and total-NBANDS entries, the reflectivity/refractive/extinction comparisons, and the dielectric convergence figures therefore need the corresponding manuscript-side replacements. `manifest.json` records the delivered figure statuses. In particular, only the hcp NBANDS dielectric convergence figure S3.3 was verified as requiring no change; the other C-group figures needed thickness, axis-limit, duplicate-entry or band-label corrections.

The hcp k-point convergence dataset at 45^3 has a very large low-energy dielectric value that the old plot concealed. Displaying it honestly is a plotting correction; it does not establish that the low-frequency metallic spectrum is converged. The final selected curves also change NBANDS as well as k-point mesh in some convergence comparisons, so they are final-parameter references rather than additional strict one-parameter convergence points.

No editable manuscript or SI source was included with this review, and no external manuscript repository was changed. This file records the precise changes needed to bring that text into agreement with the repaired source results.

## Subsequent application to the editable manuscripts

After the article and Thesis repository paths were supplied, these optical corrections were applied to their editable TeX and SI files. The local repositories are `../H-Beryllene_article_2026` and `../PhD_thesis_20251216`; each now contains `optics_update.md` and `optics_update_manifest.json`. The original draft PDF remains unchanged. No local backups, commits or remote pushes were created as part of the file synchronization.
