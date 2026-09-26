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
  * Hydrogen atom:    `E_H = -.11174239E+01` (spin polarized, 1 μB, 20 Å box; `1.0_convergence/single_Hydrogen_spin`)
  * Beryllium atom:   `E_Be = -.38385171E-01` (closed-shell 2s², ISPIN = 1 is correct)
  * Hydrogen molecule: `E_H2 = -.67727611E+01` (relaxed H–H 0.750 Å, 20 Å box; `1.0_convergence/H2_molecule`)
  * Correction (2026-09-26): the earlier `E_H = -.14200286E-01` (`1.0_convergence/single_Hydrogen`) was a non-spin-polarized H atom and must not be used; hydrogen adsorption is referenced to H2.

* general information
  * energy cutoff:    `ENCUT = 600`
  * xenes k-points:   `105 105 1`
  * bulks k-points:   `105 105 105`

* cohesive energy formula: `E_tot - m*E_Be - n*E_H` (spin-polarized atoms)
* adsorption energy formula (per H, relative to H2): `(E_tot - E_base - n/2*E_H2)/n`; negative = hydrogenation exothermic

### α-family

* α-Beryllium bulk
  * material energy:  `E_a_bulk = -.75339682E+01`
  * cohesive energy:  `E_a_bulk - 2*E_Be = -7.457197858`
  * lattice constant: `2.2664238526360245`
  * occupied states:  `2/66`
  * angle:            `120.000`

* α-beryllene xene
  * material energy:  `E_a_xene = -.29724647E+01`
  * cohesive energy:  `E_a_xene - 1*E_Be = -2.934079529`
  * lattice constant: `2.1289827591301242`
  * occupied states:  `2/133`
  * thickness:        `2*1.98`
  * angle:            `120.000`

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

### β-family

* β-beryllene xene
  * material energy:  `E_b_xene = -.66470385E+01`
  * cohesive energy:  `E_b_xene - 2*E_Be = -6.570268158`
  * lattice constant: `2.1469504557842742`
  * occupied states:  `3/133`
  * thickness:        `0.0473771641975619 + 2*1.98`
  * angle:            `119.800`

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

* β-beryllene double side adhesion
  * material energy:  `E_b_bb = -.13775686E+02`
  * cohesive energy:  `E_b_bb - 2*E_Be - 2*E_H = -11.464067858` (-2.866 eV/atom)
  * adsorption:       `(E_b_bb - E_b_xene - E_H2)/2 = -0.177943200` eV/H
  * lattice constant: `2.5960359691490167`
  * H-Be bond length: `1.44653(0), 1.44642(0), (1:1)`
  * occupied states:  `3/140`
  * thickness:        `0.0948207854203246 + 2*0.70`
  * angle:            `132.726`

### cubic family

* bulk bcc beryllium
  * material energy:  `E_c_bulk = -.36661415E+01`
  * cohesive energy:  `E_c_bulk - 1*E_Be = -3.627756329`
  * lattice constant: `2.1663874890984025`
  * occupied states:  `1/70`
  * angle:            `90.000`

* cubic beryllene trilayer xene
  * material energy:  `E_c3 = -.10008456E+02`
  * cohesive energy:  `E_c3 - 3*E_Be = -9.893300487`
  * lattice constant: `2.18374554`
  * occupied states:  `4/133`
  * thickness         `0.0751622445053417 + 2*1.98`
  * angle:            `90.000`
