# 1 "scala_struct.F"
MODULE scala_struct
  IMPLICIT NONE

!> In most cases DTYPE_ = BLOCK_CYCLIC_2D = 1
  INTEGER,PARAMETER :: DTYPE_=1
!> The BLACS context handle.
!>
!> The BLACS context handle, indicating the BLACS process grid A is
!> distributed over. The context itself is global, but the handle (the
!> integer value) may vary.
  INTEGER,PARAMETER :: CTXT_=2
!> The number of rows in the global array A.
  INTEGER,PARAMETER :: M_=3
!> The number of columns in the global array A.
  INTEGER,PARAMETER :: N_=4
!> The blocking factor used to distribute the rows of the array.
  INTEGER,PARAMETER :: MB_=5
!> The blocking factor used to distribute the columns of the array.
  INTEGER,PARAMETER :: NB_=6
!> The process row over which the first row of the array A is distributed.
  INTEGER,PARAMETER :: RSRC_=7
!> The process column over which the first column of the array A is distributed.
  INTEGER,PARAMETER :: CSRC_ =8
!> The leading dimension of the local array.  LLD_A >= MAX(1,LOCr(M_A)).
  INTEGER,PARAMETER :: LLD_=9
!> dimension of DESCA
  INTEGER,PARAMETER :: DLEN_=9
END MODULE scala_struct
