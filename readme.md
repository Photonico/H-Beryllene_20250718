# H-beryllene

Investigate the structural, electronic, optical, and topological properties of H-beryllene.

## License

The source code in the folder `vmatplot` is released under the MIT License.

You are welcome to use, modify, distribute, and adapt the code for academic, educational,
or commercial purposes, provided that the original copyright and license notice are retained.

Unless otherwise stated, this license applies to the Python source code and related scripts in this repository.
Research data, manuscript text, figures, and third-party files may be subject to separate terms.

## systems information

* single atom energy
  * Hydrogen atom:      `E_H = -.11174239E+01`
    * spin polarized, 1 μB, 20 Å box
    * folder:           `1.0_convergence/single_Hydrogen_spin`
  * Beryllium atom:     `E_Be = -.38385171E-01` (closed-shell 2s², ISPIN = 1 is correct)
  * Hydrogen molecule:  `E_H2 = -.67727611E+01`
    * relaxed H–H 0.750 Å, 20 Å box
    * folder:           `1.0_convergence/H2_molecule`
  * Correction (2026-09-26):
    * the earlier       `E_H = -.14200286E-01` is a non-spin-polarized H atom
    * old folder:       `1.0_convergence/single_Hydrogen`
    * it must not be used; hydrogen adsorption is referenced to H2

* general information
  * energy cutoff:    `ENCUT = 600`
  * xenes k-points:   `105 105 1`
  * bulks k-points:   `105 105 105`
  * functional:       `GGA = PE` (PBE, no D3 in relaxations and energies)
  * PAW potentials:   `PAW_PBE Be 06Sep2000` (2s2), `PAW_PBE H 15Jun2001`
  * precision:        `PREC = Accurate`
    * relaxation:     `EDIFF = 1E-6`
    * SOC topology:   `EDIFF = 1E-8`
  * smearing:         `ISMEAR = 0`
    * relaxation:     `SIGMA = 0.05`
    * dielectric, SOC topology: `SIGMA = 0.01`
  * spin:             `ISPIN = 1`
    * nonmagnetic: FM and AFM seeds collapse (`9.0_topology/*/spin_screen`)
  * relaxation:       `ISIF = 3`, `EDIFFG = -1E-3`
    * xenes:          `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.`
  * relaxation k-points: `27 27 1`
    * β-beryllene xene: `11 11 1`
    * trilayer:       `25 25 1`
    * bulks:          `27 27 27`
  * slab setup:       `c = 40 Å`, `LDIPOL = .TRUE.`, `IDIPOL = 3`, `DIPOL` at the slab centre
  * reference geometries: `6.0_dielectric_selection/<structure>/CONTCAR`

* cohesive energy formula: `E_tot - m*E_Be - n*E_H` (spin-polarized atoms)
* adsorption energy formula: `(E_tot - E_base - n/2*E_H2)/n`
  * per H, relative to H2; negative = hydrogenation exothermic

### α-family

* α-Beryllium bulk
  * material energy:  `E_a_bulk = -.75339682E+01`
  * cohesive energy:  `E_a_bulk - 2*E_Be = -7.457197858`
  * lattice constant: `2.2664238526360245`
  * occupied states:  `2/66`
  * angle:            `120.000`
  * relaxation:       `2.0_geometry_optimization/a-Beryllium`
  * space group:      `P6_3/mmc` (hcp)
  * cohesive per atom: `-3.729` eV/atom

* α-beryllene xene
  * material energy:  `E_a_xene = -.29724647E+01`
  * cohesive energy:  `E_a_xene - 1*E_Be = -2.934079529`
  * lattice constant: `2.1289827591301242`
  * occupied states:  `2/133`
  * thickness:        `2*1.98`
  * angle:            `120.000`
  * relaxation:       `2.0_geometry_optimization/a-Beryllene`
  * space group:      `Cmmm` (`P6/mmm` within 1E-2 Å)
  * cohesive per atom: `-2.934` eV/atom
  * electronic (SOC): metal, indirect overlap `-4.03` eV
  * topology (SOC):   conditional `ν = 1` (lowest 2 bands)
    * min direct gap: `1.127` meV at K
    * position:       `+2.72` eV above E_F
  * topology folder:  `9.0_topology/a-Beryllene`

### α-family with Hydrogen adhesion

* α-beryllene double side adhesion
  * material energy:  `E_a_hh = -.10406248E+02`
  * cohesive energy:  `E_a_hh - 1*E_Be - 2*E_H = -8.133015029` (-2.711 eV/atom)
  * adsorption:       `(E_a_hh - E_a_xene - E_H2)/2 = -0.330511100` eV/H
  * lattice constant: `2.322032201315347`
  * H-Be bond length: `1.59940(0), 1.59941(0), (2:1)`
  * occupied states:  `2/140`
  * thickness:        `0.04361280878130169 + 2*0.70`
  * angle:            `120.000`
  * relaxation:       `2.1_geometry_optimization_with_hydrogen_a/a-Beryllene_hh_staggered`
  * space group:      `C2/m` (`P-3m1` within 1E-2 Å)
  * electronic (SOC): insulator, indirect gap `4.84` eV (PBE)
  * topology (SOC):   conditional `ν = 0` (lowest 4 bands, trivial insulator)
  * phonon:           stable at the 5x5 commensurate q (old PBE+D3 run)
    * PBE rerun:      `3.7_phonon_dispersion_with_hydrogen_pbe`
  * topology folder:  `9.0_topology/a-Beryllene_hh`

### β-family

* β-beryllene xene
  * material energy:  `E_b_xene = -.66470385E+01`
  * cohesive energy:  `E_b_xene - 2*E_Be = -6.570268158`
  * lattice constant: `2.1469504557842742`
  * occupied states:  `3/133`
  * thickness:        `0.0473771641975619 + 2*1.98`
  * angle:            `119.800`
  * lattice constant b: `2.1863886759357474`
  * relaxation:       `2.0_geometry_optimization/b-Beryllene`
  * space group:      `P-1` (distorted buckled honeycomb)
  * cohesive per atom: `-3.285` eV/atom
  * electronic (SOC): metal, indirect overlap `-3.52` eV
  * topology (SOC):   conditional `ν = 0` (lowest 4 bands)
    * min direct gap: `0.854` meV
    * position:       `+1.58` eV above E_F
  * topology folder:  `9.0_topology/b-Beryllene`

### β-family with Hydrogen adhesion

* β-beryllene single side adhesion
  * material energy:  `E_b_b = -.10116602E+02`
  * cohesive energy:  `E_b_b - 2*E_Be - 1*E_H = -8.922407758` (-2.974 eV/atom)
  * adsorption:       `E_b_b - E_b_xene - E_H2/2 = -0.083182950` eV/H
  * lattice constant: `2.322152643164112`
  * H-Be bond length: `1.45061(0), 1.45058(0), (1:1)`
  * occupied states:  `3/140`
  * thickness:        `0.07108860259770651 + 0.70 + 1.98`
  * angle:            `116.989`
  * lattice constant b: `2.107952248259575`
  * relaxation:       `2.2_geometry_optimization_with_hydrogen_b/b-Beryllene_h2`
  * space group:      `P1` (`Cm` within 1E-2 Å, no inversion)
  * electronic (SOC): metal, 5 electrons per cell (bands 5/6 half filled)
  * topology (SOC):   Wilson-loop, conditional
    * lowest 4 bands: `Z2 = 0`
    * lowest 6 bands: `Z2 = 0`
  * phonon:           commensurate q stable
    * small imaginary pocket near Γ is a ZA interpolation artefact
    * audit:          `3.5_phonon_dispersion_with_hydrogen_selected/audit_20260927.md`
  * topology folder:  `9.0_topology/b-Beryllene_b`

* β-beryllene double side adhesion
  * material energy:  `E_b_bb = -.13775686E+02`
  * cohesive energy:  `E_b_bb - 2*E_Be - 2*E_H = -11.464067858` (-2.866 eV/atom)
  * adsorption:       `(E_b_bb - E_b_xene - E_H2)/2 = -0.177943200` eV/H
  * lattice constant: `2.5960359691490167`
  * H-Be bond length: `1.44653(0), 1.44642(0), (1:1)`
  * occupied states:  `3/140`
  * thickness:        `0.0948207854203246 + 2*0.70`
  * angle:            `132.726`
  * lattice constant b: `2.596491890616885`
  * relaxation:       `2.2_geometry_optimization_with_hydrogen_b/b-Beryllene_tt`
  * space group:      `P-1` (`C2/m` within 1E-2 Å)
  * electronic (SOC): metal, indirect overlap `-6.15` eV
  * topology (SOC):   conditional `ν = 0` (lowest 6 bands)
  * phonon:           soft mode `-0.646` THz at q = (0.4, 0.6) (old PBE+D3 run)
    * PBE check:      `3.6_phonon_soft_mode_check/b-Beryllene_bb`
  * topology folder:  `9.0_topology/b-Beryllene_bb`

### cubic family

* bulk bcc beryllium
  * material energy:  `E_c_bulk = -.36661415E+01`
  * cohesive energy:  `E_c_bulk - 1*E_Be = -3.627756329`
  * lattice constant: `2.1663874890984025`
  * occupied states:  `1/70`
  * angle:            `90.000`
  * relaxation:       `2.0_geometry_optimization/c0-Beryllium`
  * space group:      `Im-3m` (bcc)

* cubic beryllene trilayer xene
  * material energy:  `E_c3 = -.10008456E+02`
  * cohesive energy:  `E_c3 - 3*E_Be = -9.893300487`
  * lattice constant: `2.18374554`
  * occupied states:  `4/133`
  * thickness         `0.0751622445053417 + 2*1.98`
  * angle:            `90.000`
  * relaxation:       `2.0_geometry_optimization/c3-Beryllene_trilayer`
  * space group:      `P4/mmm`
  * cohesive per atom: `-3.298` eV/atom
  * electronic (SOC): metal, indirect overlap `-3.75` eV
  * topology (SOC):   conditional `ν = 1` (lowest 6 bands)
    * min direct gap: `0.735` meV on Γ–M
    * position:       `+1.19` eV above E_F
  * topology folder:  `9.0_topology/c-Beryllene_trilayer`
