# 1 "scalapack_wrappers.F"
# 1 "./symbol.inc" 1 

!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define 1                 compile for parallel machine with 1
!------------- end of user part --------------------------------
# 17

# 62

# 91

!
!   charge density: full grid mode
!










# 113

!
!   charge density complex
!






# 133

!
!   wavefunctions: full grid mode
!

# 182

!
!   wavefunctions complex
!








































!
!   common definitions
!








!
!   mpi parallel macros
!







# 268

# 274

# 279

# 286

# 293



# 302







# 317











# 339









!
! OpenMP macros
!
# 354

# 357

# 366









!
! profiling macros
!
# 381




!
! shmem macros
!
# 390

# 393

# 396

!
! quadruple precision
!
# 414











!
! for the 1 interface
!
# 430

!
! CUDA includes
!
# 438

!
! SIMD related definitions
!
# 1 "./simd.inc" 1 
!!#if   defined(__MIC__) || defined(__AVX512F__)
!!#define SIMD512
!!#undef  SIMD256
!!#elif defined(__AVX__) || defined(__AVX2__)
!!#define SIMD256
!!#undef  SIMD512
!!#endif

# 17





# 25





# 35















# 54





# 443 "./symbol.inc" 2 
!
! memalign macros
!
# 456




!
! Macros for PGI/NV HPC compilers version specific code
!
# 466

# 473



# 485










# 497


# 501


!
! OpenACC macros
!
# 527



































# 568












































!
! combined OpenMP and OpenACC macros
!
# 617



!
! routines replaced in LAPACK >=3.6
!
# 625


!
! Macros for HDF5 error check
!


# 634


!
! Macros for GNU version specific code
!
# 641


!
! For machine learning
!




!
! Extra safe initializations (+ overflow protections)
!
# 655




!
! Macros for memory estimation
!




!
! Offloading related macros
!
# 671


# 695

! Line is included when  !OFFLOADING

! Line is a comment when !OFFLOADING



!replace blas and lapack wrapper calls with
!normal calls if no offloading is used:





































# 747









# 759

! Line is a comment when !_OPENACC or   ACC_OFFLOAD



# 769

! Line is a comment when !ACC_OFFLOAD



# 779

! Line is included when  !_OPENACC

! Line is a comment when !_OPENACC



# 789

! Line is a comment when !_OPENMP or   OMP_OFFLOAD



# 799

! Line is a comment when !OMP_OFFLOAD



# 809

! Line is included when  !_OPENMP

! Line is a comment when !_OPENMP


# 2 "scalapack_wrappers.F" 2 

!************************* SUBROUTINE WPDTRTRI *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPDTRTRI(UPLO,DIAG, N,A,IA,JA,DESCA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PDTRTRI_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO, DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL PDTRTRI_OFFLOAD(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)
!!    RETURN
!! ENDIF

      CALL PDTRTRI(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPDTRTRI


!************************* SUBROUTINE PZTRTRI *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPZTRTRI(UPLO,DIAG, N,A,IA,JA,DESCA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PZTRTRI_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO, DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL PZTRTRI_OFFLOAD(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)
!!    RETURN
!! ENDIF

      CALL PZTRTRI(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPZTRTRI

!************************* SUBROUTINE PDPOTRF *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPDPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PDPOTRF_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PDPOTRF_OFFLOAD(UPLO,N,A,IA,JA,DESCA,INFO)
!!    RETURN
!! ENDIF

      CALL PDPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPDPOTRF

!************************* SUBROUTINE PZPOTRF *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPZPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PZPOTRF_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PZPOTRF_OFFLOAD(UPLO,N,A,IA,JA,DESCA,INFO)
!!    RETURN
!! ENDIF

      CALL PZPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPZPOTRF


!************************* SUBROUTINE PSSYEVX *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPSSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                          IFAIL, ICLUSTR, GAP, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PSSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), GAP(*), W(*), WORK(*), Z(*)
      REAL(qs) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PSSYEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
!!                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
!!                         Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
!!                         IFAIL, ICLUSTR, GAP, INFO)
!!    RETURN
!! ENDIF

      CALL PSSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                   IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPSSYEVX


!************************* SUBROUTINE PDSYEVX *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPDSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                          IFAIL, ICLUSTR, GAP, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PDSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), GAP(*), W(*), WORK(*), Z(*)
      REAL(q) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PDSYEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
!!                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
!!                         Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
!!                         IFAIL, ICLUSTR, GAP, INFO)
!!    RETURN
!! ENDIF

      CALL PDSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                   IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPDSYEVX


!************************* SUBROUTINE PCHEEVX *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPCHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                          IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PCHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), WORK(*), Z(*)
      REAL(qs) :: GAP(*), RWORK(*), W(*)
      REAL(qs) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PCHEEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
!!                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
!!                         Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
!!                         IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)
!!    RETURN
!! ENDIF

      CALL PCHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                   IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPCHEEVX


!************************* SUBROUTINE PZHEEVX *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPZHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                          IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PZHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: GAP(*), RWORK(*), W(*)
      REAL(q) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PZHEEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
!!                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
!!                         Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
!!                         IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)
!!    RETURN
!! ENDIF

      CALL PZHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                   IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPZHEEVX


!************************* SUBROUTINE PSSYEVD *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPSSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                           WORK, LWORK, IWORK, LIWORK, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PSSYEVD_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), W(*), WORK(*), Z(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PSSYEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
!!                         WORK, LWORK, IWORK, LIWORK, INFO)
!!    RETURN
!! ENDIF

      CALL PSSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                   WORK, LWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPSSYEVD


!************************* SUBROUTINE PDSYEVD *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPDSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                           WORK, LWORK, IWORK, LIWORK, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PDSYEVD_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), W(*), WORK(*), Z(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PDSYEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
!!                         WORK, LWORK, IWORK, LIWORK, INFO)
!!    RETURN
!! ENDIF

      CALL PDSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                   WORK, LWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPDSYEVD


!************************* SUBROUTINE PCHEEVD *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPCHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PCHEEVD_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), WORK(*), Z(*)
      REAL(qs) :: RWORK(*), W(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PCHEEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
!!                         WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
!!    RETURN
!! ENDIF

      CALL PCHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                   WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPCHEEVD


!************************* SUBROUTINE PZHEEVD *************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPZHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PZHEEVD_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: RWORK(*), W(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PZHEEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
!!                         WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
!!    RETURN
!! ENDIF

      CALL PZHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                   WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPZHEEVD

!************************* SUBROUTINE PDSYEV **************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPDSYEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                         WORK, LWORK, RWORK, LRWORK, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PDSYEV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: RWORK(*), W(*)
      INTEGER :: DESCA(*), DESCZ(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LRWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PDSYEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
!!                        WORK, LWORK, INFO)
!!    RETURN
!! ENDIF

      CALL PDSYEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                  WORK, LWORK, INFO)

      RETURN
      END SUBROUTINE WPDSYEV


!************************* SUBROUTINE PZHEEV **************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPZHEEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                         WORK, LWORK, RWORK, LRWORK, INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PZHEEV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: RWORK(*), W(*)
      INTEGER :: DESCA(*), DESCZ(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LRWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL PZHEEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
!!                        WORK, LWORK, RWORK, LRWORK, INFO)
!!    RETURN
!! ENDIF

      CALL PZHEEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                  WORK, LWORK, RWORK, LRWORK, INFO)

      RETURN
      END SUBROUTINE WPZHEEV


!************************* SUBROUTINE PDTRMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPDTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PDTRMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL PDTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PDTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPDTRMR2D


!************************* SUBROUTINE PSTRMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPSTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PSTRMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL PSTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PSTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPSTRMR2D


!************************* SUBROUTINE PCTRMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPCTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PCTRMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL PCTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PCTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPCTRMR2D


!************************* SUBROUTINE PZTRMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPZTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PZTRMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL PZTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PZTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPZTRMR2D


!************************* SUBROUTINE PSGEMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPSGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PSGEMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

!! IF (OFFLOAD_ON) THEN
!!    CALL PSGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PSGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPSGEMR2D


!************************* SUBROUTINE PDGEMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPDGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PDGEMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

!! IF (OFFLOAD_ON) THEN
!!    CALL PDGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PDGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPDGEMR2D


!************************* SUBROUTINE PCGEMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPCGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PCGEMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

!! IF (OFFLOAD_ON) THEN
!!    CALL PCGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PCGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPCGEMR2D


!************************* SUBROUTINE PZGEMR2D ************************
!
!>
!
!***********************************************************************

      SUBROUTINE WPZGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
!! USE moffload_struct_def
!! USE moffload_scalapack, ONLY : PZGEMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

!! IF (OFFLOAD_ON) THEN
!!    CALL PZGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
!!    RETURN
!! ENDIF

      CALL PZGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPZGEMR2D


