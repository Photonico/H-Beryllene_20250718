# 1 "blas_wrappers.F"
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


# 2 "blas_wrappers.F" 2 

!************************* SUBROUTINE DCOPY ***************************
!
!> Wrapper for BLAS routine DCOPY. (DY := DX)
!> DX and DY must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WDCOPY(N,DX,INCX,DY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DCOPY_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY

!! IF (OFFLOAD_ON) THEN
!!    CALL DCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
!!    RETURN
!! ENDIF

      CALL DCOPY(N,DX,INCX,DY,INCY)

      RETURN
      END SUBROUTINE WDCOPY


!************************* SUBROUTINE SCOPY ***************************
!
!> Wrapper for BLAS routine SCOPY. (DY := DX)
!> DX and DY must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WSCOPY(N,DX,INCX,DY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : SCOPY_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY

!! IF (OFFLOAD_ON) THEN
!!    CALL SCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
!!    RETURN
!! ENDIF

      CALL SCOPY(N,DX,INCX,DY,INCY)

      RETURN
      END SUBROUTINE WSCOPY


!************************* SUBROUTINE CCOPY ***************************
!
!> Wrapper for BLAS routine CCOPY. (ZY := ZX)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WCCOPY(N,ZX,INCX,ZY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : CCOPY_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY

!! IF (OFFLOAD_ON) THEN
!!    CALL CCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
!!    RETURN
!! ENDIF

      CALL CCOPY(N,ZX,INCX,ZY,INCY)

      RETURN
      END SUBROUTINE WCCOPY


!************************* SUBROUTINE ZCOPY ***************************
!
!> Wrapper for BLAS routine ZCOPY. (ZY := ZX)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WZCOPY(N,ZX,INCX,ZY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZCOPY_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY

!! IF (OFFLOAD_ON) THEN
!!    CALL ZCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
!!    RETURN
!! ENDIF

      CALL ZCOPY(N,ZX,INCX,ZY,INCY)

      RETURN
      END SUBROUTINE WZCOPY


!************************* SUBROUTINE DAXPY ***************************
!
!> Wrapper for BLAS routine DAXPY. (DY := DA*DX + DY)
!> DX and DY must be present on the device if offloading is enabled.
!> The scalar DA may be passed from the host.
!
!***********************************************************************

      SUBROUTINE WDAXPY(N,DA,DX,INCX,DY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DAXPY_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*),DY(*)
      REAL(q) :: DA
      INTEGER :: N,INCX,INCY

!! IF (OFFLOAD_ON) THEN
!!    CALL DAXPY_OFFLOAD(N,DA,DX,INCX,DY,INCY)
!!    RETURN
!! ENDIF

      CALL DAXPY(N,DA,DX,INCX,DY,INCY)

      RETURN
      END SUBROUTINE WDAXPY


!************************* SUBROUTINE ZAXPY ***************************
!
!> Wrapper for BLAS routine ZAXPY. (ZY := ZA*ZX + ZY)
!> ZX and ZY must be present on the device if offloading is enabled.
!> The scalar ZA may be passed from the host.
!
!***********************************************************************

      SUBROUTINE WZAXPY(N,ZA,ZX,INCX,ZY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZAXPY_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*),ZY(*)
      COMPLEX(q) :: ZA
      INTEGER :: N,INCX,INCY

!! IF (OFFLOAD_ON) THEN
!!    CALL ZAXPY_OFFLOAD(N,ZA,ZX,INCX,ZY,INCY)
!!    RETURN
!! ENDIF

      CALL ZAXPY(N,ZA,ZX,INCX,ZY,INCY)

      RETURN
      END SUBROUTINE WZAXPY


!************************* SUBROUTINE DSCAL ***************************
!
!> Wrapper for BLAS routine DSCAL. (DX := DA*DX)
!> DX must be present on the device if offloading is enabled.
!> The scalar DA may be passed from the host.
!
!***********************************************************************

      SUBROUTINE WDSCAL(N,DA,DX,INCX)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DSCAL_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*)
      REAL(q) :: DA
      INTEGER :: N,INCX

!! IF (OFFLOAD_ON) THEN
!!    CALL DSCAL_OFFLOAD(N,DA,DX,INCX)
!!    RETURN
!! ENDIF

      CALL DSCAL(N,DA,DX,INCX)

      RETURN
      END SUBROUTINE WDSCAL


!************************* SUBROUTINE ZSCAL ***************************
!
!> Wrapper for BLAS routine ZSCAL. (ZX := ZA*ZX)
!> ZX must be present on the device if offloading is enabled.
!> The scalar ZA may be passed from the host.
!
!***********************************************************************

      SUBROUTINE WZSCAL(N,ZA,ZX,INCX)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZSCAL_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*)
      COMPLEX(q) :: ZA
      INTEGER :: N,INCX

!! IF (OFFLOAD_ON) THEN
!!    CALL ZSCAL_OFFLOAD(N,ZA,ZX,INCX)
!!    RETURN
!! ENDIF

      CALL ZSCAL(N,ZA,ZX,INCX)

      RETURN
      END SUBROUTINE WZSCAL


!************************* SUBROUTINE ZDSCAL **************************
!
!> Wrapper for BLAS routine ZDSCAL. (ZX := DA*ZX)
!> ZX must be present on the device if offloading is enabled.
!> The scalar DA may be passed from the host.
!
!***********************************************************************

      SUBROUTINE WZDSCAL(N,DA,ZX,INCX)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZDSCAL_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*)
      REAL(q) :: DA
      INTEGER :: N,INCX

!! IF (OFFLOAD_ON) THEN
!!    CALL ZDSCAL_OFFLOAD(N,DA,ZX,INCX)
!!    RETURN
!! ENDIF

      CALL ZDSCAL(N,DA,ZX,INCX)

      RETURN
      END SUBROUTINE WZDSCAL


!************************* SUBROUTINE DDOT ****************************
!
!> Wrapper for BLAS function DDOT. (DDOT := DX . DY)
!> DX and DY must be present on the device if offloading is enabled.
!
!***********************************************************************

      REAL(q) FUNCTION WDDOT(N,DX,INCX,DY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DDOT_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY
! local variables
      REAL(q), EXTERNAL :: DDOT

!! IF (OFFLOAD_ON) THEN
!!    WDDOT = DDOT_OFFLOAD(N,DX,INCX,DY,INCY)
!!    RETURN
!! ENDIF

      WDDOT = DDOT(N,DX,INCX,DY,INCY)

      RETURN
      END FUNCTION WDDOT


!************************* SUBROUTINE ZDOTC ***************************
!
!> Wrapper for BLAS function ZDOTC. (ZDOTC := ZX . ZY)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************

      COMPLEX(q) FUNCTION WZDOTC(N,ZX,INCX,ZY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZDOTC_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY
! local variables
      COMPLEX(q), EXTERNAL :: ZDOTC

!! IF (OFFLOAD_ON) THEN
!!    WZDOTC = ZDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
!!    RETURN
!! ENDIF

      WZDOTC = ZDOTC(N,ZX,INCX,ZY,INCY)

      RETURN
      END FUNCTION WZDOTC


!************************* SUBROUTINE CDOTC ***************************
!
!> Wrapper for BLAS function CDOTC. (CDOTC := ZX . ZY)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************

      COMPLEX(qs) FUNCTION WCDOTC(N,ZX,INCX,ZY,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : CDOTC_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY
! local variables
      COMPLEX(qs), EXTERNAL :: CDOTC

!! IF (OFFLOAD_ON) THEN
!!    WCDOTC = CDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
!!    RETURN
!! ENDIF

      WCDOTC = CDOTC(N,ZX,INCX,ZY,INCY)

      RETURN
      END FUNCTION WCDOTC


!************************* SUBROUTINE DGEMV ***************************
!
!> Wrapper for BLAS routine DGEMV.
!> (Y := ALPHA*A*X + BETA*Y, if MODE='N';
!> Y := ALPHA*A'*X + BETA*Y, if MODE='T';
!> Y := ALPHA*conjg(A')*X + BETA*Y, if MODE='C')
!>
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> A and X must be present on the device.
!> Y must be copied to the device if BETA /= 0.
!> Y must be created on the device if BETA = 0.
!> Y must be copied to the host at the end.
!
!***********************************************************************

      SUBROUTINE WDGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DGEMV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),X(*),Y(*)
      REAL(q) :: ALPHA,BETA
      INTEGER :: M,N,LDA,INCX,INCY
      CHARACTER(1) :: MODE

!! IF (OFFLOAD_ON) THEN
!!    CALL DGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!!    RETURN
!! ENDIF

      CALL DGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      RETURN
      END SUBROUTINE WDGEMV


!************************* SUBROUTINE ZGEMV ***************************
!
!> Wrapper for BLAS routine ZGEMV.
!> (Y := ALPHA*A*X + BETA*Y, if MODE='N';
!> Y := ALPHA*A'*X + BETA*Y, if MODE='T';
!> Y := ALPHA*conjg(A')*X + BETA*Y, if MODE='C')
!>
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> A and X must be present on the device.
!> Y must be copied to the device if BETA /= 0.
!> Y must be created on the device if BETA = 0.
!> Y must be copied to the host at the end.
!
!***********************************************************************

      SUBROUTINE WZGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZGEMV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),X(*),Y(*)
      COMPLEX(q) :: ALPHA,BETA
      INTEGER :: M,N,LDA,INCX,INCY
      CHARACTER(1) :: MODE

!! IF (OFFLOAD_ON) THEN
!!    CALL ZGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!!    RETURN
!! ENDIF

      CALL ZGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      RETURN
      END SUBROUTINE WZGEMV


!************************* SUBROUTINE DGEMM ***************************
!
!> Wrapper for BLAS routine DGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WDGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DGEMM_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
      REAL(q) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

!! IF (OFFLOAD_ON) THEN
!!    CALL DGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!    RETURN
!! ENDIF

      CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WDGEMM


!************************* SUBROUTINE SGEMM ***************************
!
!> Wrapper for BLAS routine SGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WSGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : SGEMM_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
      REAL(qs) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

!! IF (OFFLOAD_ON) THEN
!!    CALL SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!    RETURN
!! ENDIF

      CALL SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WSGEMM


!************************* SUBROUTINE ZGEMM ***************************
!
!> Wrapper for BLAS routine ZGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZGEMM_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
      COMPLEX(q) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

!! IF (OFFLOAD_ON) THEN
!!    CALL ZGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!    RETURN
!! ENDIF

      CALL ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WZGEMM


!************************* SUBROUTINE CGEMM ***************************
!
!> Wrapper for BLAS routine CGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WCGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : CGEMM_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
      COMPLEX(qs) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

!! IF (OFFLOAD_ON) THEN
!!    CALL CGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!    RETURN
!! ENDIF

      CALL CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WCGEMM


!************************* SUBROUTINE DTRMM ***************************
!
!> Wrapper for BLAS routine DTRMM.
!>    if SIDE = 'L' or 'l', then B := ALPHA*op(A)*B;
!>    if SIDE = 'R' or 'r', then B := ALPHA*B*op(A);
!>    if TRANSA = 'N' or 'n', then op(A) = A;
!>    if TRANSA = 'T' or 't', then op(A) = A';
!>    if TRANSA = 'C' or 'c', then op(A) = conjg(A');
!> A and B must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WDTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DTRMM_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*)
      REAL(q) :: ALPHA
      INTEGER :: LDA,LDB,M,N
      CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL DTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!!    RETURN
!! ENDIF

      CALL DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      RETURN
      END SUBROUTINE WDTRMM


!************************* SUBROUTINE ZTRMM ***************************
!
!> Wrapper for BLAS routine ZTRMM.
!>    if SIDE = 'L' or 'l', then B := ALPHA*op(A)*B;
!>    if SIDE = 'R' or 'r', then B := ALPHA*B*op(A);
!>    if TRANSA = 'N' or 'n', then op(A) = A;
!>    if TRANSA = 'T' or 't', then op(A) = A';
!>    if TRANSA = 'C' or 'c', then op(A) = conjg(A');
!> A and B must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : ZTRMM_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*)
      COMPLEX(q) :: ALPHA
      INTEGER :: LDA,LDB,M,N
      CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL ZTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!!    RETURN
!! ENDIF

      CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      RETURN
      END SUBROUTINE WZTRMM


!************************* SUBROUTINE DLACPY **************************
!
!> Wrapper for BLAS routine DLACPY. (B := A)
!> Only a full copy via dcopy is implemented for GPU offloading.
!> This means UPLO = F is forced!
!> A and B must be present on the device if offloading is enabled.
!
!***********************************************************************

      SUBROUTINE WDLACPY(UPLO,M,N,A,LDA,B,LDB)

      USE prec
!! USE moffload_struct_def
!! USE moffload_blas, ONLY : DLACPY_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*)
      INTEGER :: LDA,LDB,M,N
      CHARACTER(1) :: UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL DLACPY_OFFLOAD(UPLO,M,N,A,LDA,B,LDB)
!!    RETURN
!! ENDIF

      CALL DLACPY(UPLO,M,N,A,LDA,B,LDB)

      RETURN
      END SUBROUTINE WDLACPY
