from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path('/Users/lu/Repos/H-Beryllene_20250718')
ARTICLE=Path('/Users/lu/Repos/H-Beryllene_article_2026')
THESIS=Path('/Users/lu/Repos/PhD_thesis_20251216')
manifest=json.loads((ROOT/'exported_figures/manifest.json').read_text())
changes={}

def load(path):
    return changes.get(path,path.read_text())

def exact(path,old,new):
    s=load(path)
    assert s.count(old)==1,(str(path),old[:100],s.count(old))
    changes[path]=s.replace(old,new,1)

def span(path,start,end,new):
    s=load(path)
    assert s.count(start)==1,(path,start,s.count(start))
    a=s.index(start);b=s.index(end,a+len(start))
    changes[path]=s[:a]+new+s[b:]

# Keep frequency notation consistent with the absorption formula.
p=ARTICLE/'computational_methods.tex'
exact(p,r'a vacuum region of \qty{40}{\angstrom} is included.',r'the full supercell height normal to the layers is \qty{40}{\angstrom}, including the slab and the surrounding vacuum.')
exact(p,r'and $\omega$ denotes the photon energy.',r'and $\omega$ is the angular frequency, related to the photon energy by $E=\hbar\omega$.')
exact(p,r'\delta\!\left(E_{c\mathbf{k}}-E_{v\mathbf{k}}-\omega\right)',r'\delta\!\left(E_{c\mathbf{k}}-E_{v\mathbf{k}}-\hbar\omega\right)')
exact(p,r'The real part is obtained by applying the Kramers--Kronig transformation:',r'In the zero-broadening limit, the real part is related to the imaginary part by the Kramers--Kronig transformation:')
exact(p,"{\\omega'^{2}-\\omega^{2}+i\\eta}","{\\omega'^{2}-\\omega^{2}}")
exact(p,r'''$\mathcal{P}$ denotes the Cauchy principal value,
and $\eta$ is the complex shift,
with $\eta=0.1$ in the present calculations.''',r'''and $\mathcal{P}$ denotes the Cauchy principal value.
The numerical spectra are calculated in VASP 6.5.0 with \texttt{LOPTICS=T},
Gaussian occupation smearing of \qty{0.01}{\electronvolt} (\texttt{ISMEAR=0}, \texttt{SIGMA=0.01}),
and complex broadening \texttt{CSHIFT=0.1} in eV.
All eight optical datasets are spin-unpolarized and exclude spin-orbit coupling (\texttt{ISPIN=1}, \texttt{LSORBIT=F}).
They contain the independent-particle interband response without microscopic local-field or electron-hole corrections.
The actual output has \texttt{WPLASMAI=0}; no Drude term is added in the postprocessing.
These spectra therefore do not describe the complete low-frequency response of the metallic phases.''')
exact(p,r'$c$ is the speed of light in vacuum.',r'''$c$ is the speed of light in vacuum.
The absorption coefficient uses $\omega=E/\hbar=2\pi E/h$;
with $E$ in eV, $h$ in eV\,s and $c$ in nm/s, it is reported in $\mathrm{nm}^{-1}$.
Equation~\ref{eq:reflectivity} describes the scalar normal-incidence reflectivity of a homogeneous semi-infinite medium,
not the full reflection of a finite anisotropic layer on a substrate.''')
# Confirm the existing equation label before writing any file.
assert r'\label{eq:reflectivity}' in load(p)
exact(p,r'We therefore restore the intrinsic dielectric function of each two-dimensional structure using the dielectric-function renormalization method of Yang and Gao~\cite{yang2021method}:',r'We therefore obtain an effective dielectric response under the assigned-thickness convention using the dielectric-function renormalization method of Yang and Gao~\cite{yang2021method}:')
exact(p,r'are the intrinsic dielectric function of the two-dimensional structure and the dielectric function of the full supercell,',r'are the effective dielectric function assigned to the two-dimensional structure and the dielectric function of the full supercell,')
thickness=r'''We assign $d_{\mathrm{2D}}$ as the separation of the outermost atomic centres along the surface normal,
plus the radii of the atoms at those two surfaces: \qty{1.98}{\angstrom} for Be and \qty{0.70}{\angstrom} for H.
The relaxed coordinates are converted to Cartesian lengths with periodic wrapping removed.
The same optical volume normalization is applied to the $xx$, $yy$ and $zz$ components.
The absolute effective response depends on this thickness convention.
'''
exact(p,r'The values of $d_{\mathrm{2D}}$ and the values of \texttt{NBANDS} used in the dielectric function calculations are provided in the Supplementary Information.',thickness+r'The effective thicknesses and actual total numbers of bands are provided in the Supplementary Information.')

# Synchronize specific optical paragraphs, retaining each manuscript's figure blocks.
paths=[ARTICLE/'results.tex',THESIS/'src/project_3_main.tex']
alpha_intro=r'''For $\beta$-beryllene,
the in-plane real part contains zero crossings and a low-energy interval in which its sign differs from the out-of-plane component.
On the calculated grid, these components have opposite signs from approximately \qtyrange{1.32}{2.23}{\electronvolt};
these are sampled endpoints rather than interpolated zeros.
For $\alpha$-beryllene, the out-of-plane real part crosses zero near \qty{10.64}{\electronvolt}, where the imaginary part is comparatively small.
The sign changes and loss features identify candidate energy ranges for further analysis of collective response.
They do not by themselves establish propagating modes in a finite two-dimensional layer or determine their lifetimes.
The numerical broadening is an input to the calculation, not a calculated scattering lifetime.

'''
# The thesis already has a more detailed corrected alpha/beta interpretation; retain it.
s=load(ARTICLE/'results.tex');a=s.index('For $\\beta$-beryllene,',s.index(r'\label{3_optical_alpha_beta}'));b=s.index(r'\begin{figure}',a)
changes[ARTICLE/'results.tex']=s[:a]+alpha_intro+s[b:]
for p in paths:
    span(p,'The absorption coefficient and refractive index are shown in figure~\\ref{fig:proj3_alpha_beta_abs_refractive}.','The refractive index spectra follow',r'''The absorption coefficient and refractive index are shown in figure~\ref{fig:proj3_alpha_beta_abs_refractive}.
The effective absorption is strongly anisotropic.
Within the \qtyrange{0}{12}{\electronvolt} window, the largest sampled coefficients for $\alpha$-beryllene are approximately
\qty{0.135}{\per\nano\metre} in $xx$ and \qty{0.175}{\per\nano\metre} in $zz$;
the corresponding $\beta$-beryllene values are approximately \qty{0.159}{\per\nano\metre} and \qty{0.241}{\per\nano\metre}.
The strongest $\beta$ out-of-plane feature occurs near \qty{5.96}{\electronvolt}.
The spectra cross at some energies, so the larger maxima do not imply that $\beta$-beryllene absorbs more strongly at every energy.

''')
    old=r'''$\beta$-beryllene remains stronger over the low-energy range and then decreases rapidly toward the near-ultraviolet region.'''
    new=r'''$\beta$-beryllene develops its strongest extinction peak near \qty{5.96}{\electronvolt}, where $k\approx4.00$, followed by weaker higher-energy features.'''
    if old in load(p):exact(p,old,new)
    else:exact(p,r'$\beta$-beryllene develops its strongest extinction in the ultraviolet region.',new)
    span(p,'The energy-loss spectrum of $\\beta$-beryllene shows',r'\subsubsection{Bulk hcp, bulk bcc, and cubic trilayer beryllene}',r'''The in-plane energy-loss spectrum of $\beta$-beryllene has a visible-range local maximum of approximately 0.45
and a larger high-energy local maximum of approximately 1.33 near \qty{11.77}{\electronvolt}.
It reaches approximately 1.51 at the last available sample near \qty{11.99}{\electronvolt} while still rising;
this boundary value is not a resolved peak.
For out-of-plane response, $\alpha$-beryllene has a dominant loss maximum of approximately 8.88 near \qty{10.64}{\electronvolt},
substantially above the $\beta$ maximum of approximately 1.00 near \qty{9.97}{\electronvolt}.
These spectra describe the inverse dielectric response under the stated effective-thickness convention.

''')
    exact(p,'coinciding with the minimum in the real part and indicating a plasmonic optical response.','coinciding with a minimum in the real part and a change in the inverse dielectric response.')
    span(p,'For cubic trilayer beryllene,\nthe most prominent feature',r'\begin{figure}',r'''For cubic trilayer beryllene,
the dielectric response is redistributed and becomes direction dependent relative to bulk bcc beryllium.
The largest imaginary dielectric features occur near \qty{3.74}{\electronvolt} in $xx$ and $yy$, with a value of approximately 19.75,
and near \qty{2.71}{\electronvolt} in $zz$, with a value of approximately 27.32.
The real part also develops direction-dependent sign changes.
These dielectric features do not coincide in general with absorption maxima, which also depend on the photon energy and the nonlinear relation between the dielectric function and the extinction coefficient.
Assignments to specific surface bands or transitions would require band- and matrix-element-resolved evidence.

''')
    span(p,'The absorption coefficient and energy-loss spectra are shown in figure~\\ref{fig:proj3_bulk_cubic_abs_loss}.',r'\begin{figure}',r'''The absorption coefficient and energy-loss spectra are shown in figure~\ref{fig:proj3_bulk_cubic_abs_loss}.
With the corrected thickness normalization, cubic trilayer beryllene does not show a universal absorption enhancement relative to bulk bcc beryllium.
In the explicitly defined visible interval of \qtyrange{1.65}{3.26}{\electronvolt},
the trilayer $xx$ and $yy$ absorption is lower than bulk bcc at every available trilayer sample;
the bulk spectrum is linearly interpolated to the trilayer grid for this comparison.
The respective visible maxima are approximately \qty{0.0604}{\per\nano\metre} and \qty{0.0746}{\per\nano\metre}.
The trilayer $zz$ response exceeds the bulk value over only part of this interval.
Its strongest sampled absorption occurs near \qty{5.40}{\electronvolt} in $xx$ and $yy$,
and near \qty{5.25}{\electronvolt} in $zz$, with coefficients of approximately \qty{0.148}{\per\nano\metre} and \qty{0.166}{\per\nano\metre}, respectively.
These results show a redistribution of the effective interband response rather than enhanced absorption throughout the visible and near-ultraviolet regions.

The loss spectra also differ between directions and structures.
The low-energy $\beta$ in-plane feature and the high-energy $\alpha$ out-of-plane maximum are retained in the comparison.
The scalar loss features identify variations in the inverse dielectric response; they do not alone establish propagating plasmon modes.

''')
    start=load(p).index(r'\label{3_optical_hydrogenated_structures}')
    s=load(p);a=s.index('For $\\alpha$-beryllene,',start);b=s.index(r'\begin{figure}',a)
    changes[p]=s[:a]+r'''For $\alpha$-beryllene,
hydrogen adsorption changes the magnitude and spectral distribution of the effective dielectric response.
In $xx$ and $yy$, 2H-$\alpha$-beryllene has a real dielectric maximum of approximately 31.2 near \qty{4.91}{\electronvolt},
and its strongest imaginary dielectric feature is approximately 31.2 near \qty{5.59}{\electronvolt}.
The out-of-plane spectrum is different: its largest imaginary dielectric feature occurs near \qty{9.77}{\electronvolt}.
These comparisons retain the specified H/H surface-radius convention for the hydrogenated layer.
Changes in peak height therefore describe an effective response under that convention, not a thickness-independent measure of oscillator strength.

For $\beta$-beryllene,
single- and double-sided hydrogen adsorption redistribute the dielectric response and change its directional dependence.
Double-sided adsorption produces distinct spectra in all three Cartesian components.
The strongest imaginary dielectric features occur near \qty{6.59}{\electronvolt} in $xx$,
\qty{9.26}{\electronvolt} in $yy$ and \qty{8.24}{\electronvolt} in $zz$,
with values of approximately 14.96, 7.68 and 14.97, respectively.
The first real-part sign changes occur near \qty{6.64}{\electronvolt}, \qty{9.32}{\electronvolt} and \qty{8.10}{\electronvolt}, respectively.
Thus, increasing hydrogen coverage changes the spectral distribution and in-plane anisotropy;
it does not produce a common peak energy, a common zero crossing, or a monotonic enhancement in every component.
The imaginary dielectric peaks should also be distinguished from maxima of the derived absorption coefficient.

'''+s[b:]
    span(p,'The corresponding absorption coefficient and energy-loss spectra are shown in figure~\\ref{fig:proj3_hydrogenated_abs_loss}.',r'\begin{figure}',r'''The corresponding absorption coefficient and energy-loss spectra are shown in figure~\ref{fig:proj3_hydrogenated_abs_loss}.
For $\alpha$-beryllene, hydrogen adsorption produces a strongly direction-dependent redistribution of absorption.
The $xx$ and $yy$ coefficients increase over the sampled visible interval,
but the 2H-$\alpha$ visible maximum is only approximately \qty{0.00137}{\per\nano\metre}, much smaller than its ultraviolet maximum near \qty{5.71}{\electronvolt}.
The $zz$ coefficient instead decreases throughout the sampled \qtyrange{1.65}{4.13}{\electronvolt} interval,
and its strongest peak occurs near \qty{9.83}{\electronvolt}.
For 2H-$\beta$-beryllene, the strongest sampled absorption peaks occur near \qty{6.75}{\electronvolt} in $xx$,
\qty{9.39}{\electronvolt} in $yy$ and \qty{8.36}{\electronvolt} in $zz$.
The loss spectra likewise change with hydrogen coverage and direction.
These results establish changes in the calculated interband and inverse dielectric response;
weak low-energy tails under finite numerical broadening are not evidence of a sharply defined absorption onset or strong visible absorption.

''')
    span(p,'These results show that both dimensionality and hydrogen functionalization strongly modify the optical response of beryllene.',r'\section{' if 'project_3_main' in p.name else '\x00',r'''These results show that both dimensionality and hydrogen functionalization redistribute the effective interband optical response of beryllene.
The pristine phases exhibit pronounced anisotropy, while the cubic trilayer and hydrogenated phases enhance some components and reduce others.
The comparisons depend on the stated thickness convention and do not include a Drude intraband contribution.

''') if 'project_3_main' in p.name else None
    if p.name=='results.tex':
        s=load(p);a=s.index('These results show that both dimensionality and hydrogen functionalization strongly modify the optical response of beryllene.')
        changes[p]=s[:a]+r'''These results show that both dimensionality and hydrogen functionalization redistribute the effective interband optical response of beryllene.
The pristine phases exhibit pronounced anisotropy, while the cubic trilayer and hydrogenated phases enhance some components and reduce others.
The comparisons depend on the stated thickness convention and do not include a Drude intraband contribution.
'''
    # Units in all derived comparison captions, preserving existing line and figure structure.
    changes[p]=load(p).replace('Upper: absorption coefficient; Lower:',r'Upper: absorption coefficient ($\mathrm{nm}^{-1}$); Lower:')

# Correct article A index crossings without undoing the existing thesis version.
p=ARTICLE/'results.tex'
exact(p,r'''For both structures,
the refractive index drops below one above approximately \qty{8}{\electronvolt},
which occurs close to the negative region of the real part of the dielectric function.
This behaviour is consistent with the strong interaction between light and plasmonic oscillations near the plasmon resonance.''',r'''The first downward $n=1$ crossings occur near \qty{7.42}{\electronvolt} for $\alpha$ $xx$,
\qty{7.82}{\electronvolt} for $\alpha$ $zz$, \qty{7.35}{\electronvolt} for $\beta$ $xx$ and \qty{6.37}{\electronvolt} for $\beta$ $zz$.
The $\alpha$ out-of-plane index rises above one again near \qty{11.71}{\electronvolt}.
These crossings describe dispersion of the effective optical constants and do not by themselves identify a plasmon resonance.''')

# SI: thicknesses, total bands, and the limits of convergence evidence.
si_text=r'''The effective thicknesses $d_{\mathrm{2D}}$ used in the dielectric-function renormalization are
\qty{3.960}{\angstrom} for pristine $\alpha$-Beryllene,
\qty{5.855}{\angstrom} for pristine $\beta$-Beryllene,
\qty{6.966}{\angstrom} for cubic trilayer beryllene,
\qty{3.145}{\angstrom} for double-sided hydrogen-terminated $\alpha$-Beryllene,
\qty{5.524}{\angstrom} for single-sided hydrogen-terminated $\beta$-Beryllene,
and \qty{5.193}{\angstrom} for double-sided hydrogen-terminated $\beta$-Beryllene.
'''+thickness+r'''The endpoint-radius contribution is Be+Be for pristine structures, H+H for double-sided termination, and Be+H for single-sided termination.
Consequently, the assigned 2H-$\alpha$ thickness is smaller than the pristine value; this is a consequence of the chosen convention and is not a unique physical layer thickness.
Unrounded structure-derived values are used in all calculations.
The actual total numbers of bands are 140 for each hydrogenated structure,
133 for pristine $\alpha$- and $\beta$-beryllene and cubic trilayer beryllene,
70 for bulk bcc beryllium and 66 for bulk hcp beryllium.
These are total bands, not counts of unoccupied bands; the actual output may differ from the requested \texttt{NBANDS} because of parallel band allocation.

'''
convergence=r'''The plots retain the full range of every curve within the displayed energy window.
In particular, the bulk hcp \num{45} by \num{45} by \num{45} k-point calculation has a large low-energy dielectric response that is visible in figure~\ref{fig:proj3_hcp_dielectric_kpoint_convergence}.
These comparisons do not establish convergence of the complete low-frequency metallic response, which would also require an intraband treatment.
The selected reference may change both the k-point grid and the total band count relative to a convergence series;
it is a final-parameter reference rather than an additional strictly one-parameter convergence point.

'''
for p in [ARTICLE/'supplement.tex',THESIS/'src/project_3_SI.tex']:
    s=load(p)
    marker='The values of the thickness' if p.name=='supplement.tex' else 'The effective thicknesses'
    a=s.index(marker);b=s.index('and 66 for bulk hcp beryllium.',a)+len('and 66 for bulk hcp beryllium.')
    changes[p]=s[:a]+si_text.rstrip()+s[b:]
    changes[p]=load(p).replace('number of unoccupied bands','total number of bands').replace('numbers of unoccupied bands','total numbers of bands')
    if p.name=='supplement.tex':
        exact(p,'These results confirm that the computational settings adopted in the optical calculations are sufficient to converge the principal spectral features.',convergence.rstrip())
    else:
        exact(p,r'\subsection*{Dielectric-function convergence}'+'\n',r'\subsection*{Dielectric-function convergence}'+'\n\n'+convergence)

# Restrict claims consistently in the abstracts, introductions and conclusions.
optical_conclusion=r'''The calculated optical properties show a direction-dependent redistribution of the effective interband response after dimensional reduction and hydrogen functionalization.
Pristine $\alpha$- and $\beta$-beryllene exhibit pronounced anisotropy and distinct loss features.
Cubic trilayer absorption does not exceed bulk bcc throughout the visible range or in all directions,
and hydrogen adsorption enhances some components while reducing others.
These comparisons use the stated effective-thickness convention and exclude a Drude intraband contribution.

'''
span(ARTICLE/'conclusion.tex','The calculated optical properties demonstrate', 'Parity analysis of',optical_conclusion)
span(THESIS/'src/project_3_main.tex','The calculated optical properties demonstrate', 'Topological analysis based',optical_conclusion)
exact(ARTICLE/'abstract.tex','including enhanced visible-to-ultraviolet absorption, modified plasmonic behaviour,\nand direction-dependent dielectric responses associated with collective electronic excitations.', 'including a direction-dependent redistribution of visible-to-ultraviolet optical response\nand distinct dielectric and energy-loss features under the stated effective-thickness convention.')
exact(THESIS/'src/project_3_main.tex','including visible-to-ultraviolet absorption,\nmodified plasmonic behaviour,\nand direction-dependent dielectric response associated with collective electronic excitations.', 'including a direction-dependent redistribution of visible-to-ultraviolet optical response\nand distinct dielectric and energy-loss features under the stated effective-thickness convention.')
for p in [ARTICLE/'introduction.tex',THESIS/'src/project_3_main.tex']:
    exact(p,'together with plasmon-related features and visible-to-ultraviolet absorption.','together with anisotropic dielectric and energy-loss spectra and a redistribution of visible-to-ultraviolet optical response.')
p=THESIS/'src/project_3_main.tex'
exact(p,'For the pristine $\\alpha$- and $\\beta$-beryllene spectra discussed in section~\\ref{3_optical_alpha_beta},','For all eight selected optical datasets,')
exact(p,'The output contains 133 bands in total, while the input requests \\texttt{NBANDS=128}.','The pristine two-dimensional outputs contain 133 bands in total, while their inputs request \\texttt{NBANDS=128};\nthe actual total band counts for all systems are given in section~\\ref{3_supplementary}.')
exact(p,'The thickness values used for this correction are given in section~\\ref{3_supplementary}.',thickness+'The full normal supercell height for the two-dimensional optical calculations is \\qty{40}{\\angstrom}.\nThe thickness values used for this correction are given in section~\\ref{3_supplementary}.')
for p,old in [(THESIS/'src/abstract.tex','The dielectric response reveals strong anisotropy, visible-to-ultraviolet absorption, modified plasmonic behaviour, and direction-dependent screening.'),(THESIS/'src/conclusion.tex','The optical spectra show strong anisotropy, visible-to-ultraviolet absorption, modified plasmonic behaviour, and direction-dependent dielectric screening.')]:
    exact(p,old,'The beryllene optical spectra show strong anisotropy and direction-dependent redistribution of the effective interband and energy-loss response.\nNeither dimensional reduction nor hydrogenation enhances absorption in every direction or throughout the visible range.')

# Validate all destinations and the source hashes before the first mutation.
records=[]
for item in manifest:
    src=ROOT/'exported_figures'/item['export']
    assert hashlib.sha256(src.read_bytes()).hexdigest()==item['after_sha256']
    destinations=[THESIS/item['thesis_target']]
    if item['source'].startswith('figures_for_publication/'):
        destinations.append(ARTICLE/'figures_repo'/Path(item['source']).name)
    for dest in destinations:
        assert dest.is_file(),dest
        records.append({'source':str(src),'target':str(dest),'sha256':item['after_sha256'],
                        'previous_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
existing={}
for repo in [ARTICLE,THESIS]:
    existing[str(repo)]=subprocess.check_output(['git','status','--short'],cwd=repo,text=True)
untouched=[ARTICLE/'.DS_Store',THESIS/'thesis.bib',THESIS/'figures_ch4/figures_ch4.ipynb',THESIS/'src/fundamentals_c.tex',THESIS/'figures_ch4/4.1_band_subspace.pdf',THESIS/'figures_ch4/4.2_trim.pdf']
original_hash={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in untouched}
for p,s in changes.items():p.write_text(s)
for record in records:shutil.copyfile(record['source'],record['target'])
for p,h in original_hash.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
for r in records:assert hashlib.sha256(Path(r['target']).read_bytes()).hexdigest()==r['sha256']
audit={'initial_status':existing,'figure_copies':records,'text_files':[str(p) for p in changes],
       'preexisting_unrelated_files_unchanged':original_hash,'backups_created':False}
(ROOT/'.agent/manuscript_sync.json').write_text(json.dumps(audit,indent=2)+'\n')
print('Synced',len(records),'figure destinations and',len(changes),'TeX files; no backups.')
