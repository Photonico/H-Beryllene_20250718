#include "symbol.inc"
MODULE hamil_struct_def
  USE prec
  USE pot_struct_def
  ! handle for orbital magnetization and metaGGAs
  TYPE ham_handle
     TYPE(COULOMB_POTENTIAL) :: COULOMB_POT
     COMPLEX(q), POINTER, CONTIGUOUS :: AVTOT(:,:) ! local vector magnetization potential
     RGRID     , POINTER, CONTIGUOUS :: AVEC(:,:)  ! soft part of vector magnetization potential
     COMPLEX(q), POINTER, CONTIGUOUS :: MUTOT(:,:) ! derivative of energy density with respect to kinetic energy density
     RGRID     , POINTER, CONTIGUOUS :: MU(:,:)    ! same as MUTOT, but on GRID instead of GRIDC
  END TYPE ham_handle

END MODULE hamil_struct_def
