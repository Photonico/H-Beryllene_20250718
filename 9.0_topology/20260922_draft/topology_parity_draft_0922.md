# Computational methods

First-principles calculations were performed within the framework of density functional theory (DFT) using the Vienna Ab initio Simulation Package (VASP) [1]. The electron–ion interactions were modeled using the projector-augmented wave (PAW) method [2,3], while exchange–correlation effects were approximated by the Perdew–Burke–Ernzerhof (PBE) generalized gradient approximation [4]. Dispersion interactions were incorporated using the Grimme D3 correction [5]. A plane-wave kinetic energy cutoff of 520 eV was used, along with a 24 × 24 × 1 $k$-point mesh. The self-consistent convergence criterion was set to $10^{-7}$ eV, and Gaussian smearing with a width of 0.02 eV was used. A vacuum distance of 15 Å was added perpendicular to the layer to reduce interactions between periodic images. The atomic coordinates and lattice parameters were relaxed until the force acting on each atom was less than 0.01 eV/Å. Spin–orbit coupling (SOC) was considered in the calculations of the band structure and parity.

Parity eigenvalues for SOC spinor wavefunctions were calculated using irvsp [6]. Because of the existence of spatial inversion and time-reversal symmetries in all three structures, the $Z_2$ index of each was calculated based on the parity eigenvalues of occupied Kramers pairs at the time-reversal-invariant momenta (TRIM) according to the Fu-Kane criterion [7]. The parity product for the TRIM $\Lambda_i$ is given by:

$$
\delta_i = \prod_{m=1}^{N} \xi_{2m}(\Lambda_i)
$$

where $\xi_{2m}(\Lambda_i) = \pm 1$ is the parity eigenvalue for one state of the $m$-th Kramers pair.

The topological index is then obtained from:

$$
(-1)^{\nu} = \prod_{i} \delta_i
$$

Here, $\nu = 1$ and $\nu = 0$ denote nontrivial and trivial band topology, respectively.

# Parity analysis and topological classification

Table 1 gives the results for the parity eigenvalues for the SOC cases along with their parity products. The three M points are equivalent by symmetry for α- and β-beryllene, which means the Fu-Kane product will be $(-1)^{\nu} = \delta_{\Gamma} \delta_{M}^{3} = \delta_{\Gamma} \delta_{M}$. However, for square-trilayer beryllene, the $Z_2$ invariant is found from the parity products at the TRIM points as $(-1)^{\nu} = \delta_{\Gamma} \delta_{M} \delta_{X}^{2} = \delta_{\Gamma} \delta_{M}$, where $\delta_{X}^{2}$ refers to the product for the two X TRIM points.

**Table 1.** Parity eigenvalues at the time-reversal invariant momenta (TRIM), parity products, $Z_2$ topological indices ($\nu$), and corresponding topological classifications for α-beryllene, β-beryllene, and square-trilayer (ST)-beryllene.

| Structure | Γ | M | Product | ν | Classification |
| --- | --- | --- | --- | --- | --- |
| α-beryllene | +, + | +, − | −1 | 1 | Nontrivial |
| β-beryllene | +, −, +, − | −, +, +, − | +1 | 0 | Trivial |
| ST-beryllene | −, +, +, +, −, + | −, +, −, +, −, + | −1 | 1 | Nontrivial |

Parity analysis reveals unique topological properties for each of the three beryllene allotropes. The parity products at Γ and M for α-beryllene have opposite signs, giving rise to a Fu-Kane parity product value of $(-1)^{\nu} = -1$ with a $Z_2$ index value of $\nu = 1$. Such values imply that α-beryllene is a non-trivial $Z_2$ topological insulator, as reported previously [8]. On the other hand, β-beryllene has positive parity products for both Γ and M; therefore, its Fu-Kane parity product value is given by $(-1)^{\nu} = +1$, with $\nu = 0$. It can be regarded as a topologically trivial material. Square-trilayer beryllene exhibits the same parity-product sign reversal between Γ and M as α-beryllene, leading to $(-1)^{\nu} = -1$ and $\nu = 1$. This confirms that square-trilayer beryllene is also a nontrivial $Z_2$ topological phase.

# References

[1] G. Kresse and J. Furthmüller, Efficient Iterative Schemes for *ab Initio* Total-Energy Calculations Using a Plane-Wave Basis Set, Phys. Rev. B **54**, 11169 (1996).

[2] P. E. Blöchl, Projector Augmented-Wave Method, Phys. Rev. B **50**, 17953 (1994).

[3] G. Kresse and D. Joubert, From Ultrasoft Pseudopotentials to the Projector Augmented-Wave Method, Phys. Rev. B **59**, 1758 (1999).

[4] J. P. Perdew, K. Burke, and M. Ernzerhof, Generalized Gradient Approximation Made Simple, Phys. Rev. Lett. **77**, 3865 (1996).

[5] Stefan Grimme, Jens Antony, Stephan Ehrlich, and Helge Krieg, A Consistent and Accurate *ab Initio* Parametrization of Density Functional Dispersion Correction (DFT-D) for the 94 Elements H-Pu, J. Chem. Phys. **132**, 154104 (2010).

[6] J. Gao, Q. Wu, C. Persson, and Z. Wang, Irvsp: To Obtain Irreducible Representations of Electronic States in the VASP, Comput. Phys. Commun. **261**, 107760 (2021).

[7] L. Fu and C. L. Kane, Topological Insulators with Inversion Symmetry, Phys. Rev. B **76**, 45302 (2007).

[8] J. Li, M. Guo, J. Si, L. Shi, X. Shi, J.-J. Ma, Q. Zhang, D. J. Singh, P.-F. Liu, and B.-T. Wang, Coexistence of Superconductivity and Topological Aspects in Beryllenes, Mater. Today Phys. **38**, 101257 (2023).
