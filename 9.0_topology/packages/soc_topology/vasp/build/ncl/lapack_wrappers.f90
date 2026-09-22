# 1 "lapack_wrappers.F"
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


# 2 "lapack_wrappers.F" 2 

!***********************************************************************
! LAPACK
!***********************************************************************

!************************* SUBROUTINE DPOTRF **************************
!
!> Wrapper for the LAPACK routine DPOTRF. (Cholesky factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************

      SUBROUTINE WDPOTRF(UPLO,N,A,LDA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : DPOTRF_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL DPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
!!    RETURN
!! ENDIF

      CALL DPOTRF(UPLO,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WDPOTRF


!************************* SUBROUTINE ZPOTRF **************************
!
!> Wrapper for the LAPACK routine ZPOTRF. (Cholesky factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************

      SUBROUTINE WZPOTRF(UPLO,N,A,LDA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : ZPOTRF_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL ZPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
!!    RETURN
!! ENDIF

      CALL ZPOTRF(UPLO,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WZPOTRF


!************************* SUBROUTINE DGETRF **************************
!
!> Wrapper for the LAPACK routine DGETRF. (LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************

      SUBROUTINE WDGETRF(M,N,A,LDA,IPIV,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : DGETRF_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*)
      INTEGER :: M,N,LDA,INFO,IPIV(*)

!! IF (OFFLOAD_ON) THEN
!!    CALL DGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
!!    RETURN
!! ENDIF

      CALL DGETRF(M,N,A,LDA,IPIV,INFO)

      RETURN
      END SUBROUTINE WDGETRF


!************************* SUBROUTINE ZGETRF **************************
!
!> Wrapper for the LAPACK routine ZGETRF. (LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************

      SUBROUTINE WZGETRF(M,N,A,LDA,IPIV,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : ZGETRF_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*)
      INTEGER :: M,N,LDA,INFO,IPIV(*)

!! IF (OFFLOAD_ON) THEN
!!    CALL ZGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
!!    RETURN
!! ENDIF

      CALL ZGETRF(M,N,A,LDA,IPIV,INFO)

      RETURN
      END SUBROUTINE WZGETRF


!************************* SUBROUTINE DGETRS **************************
!
!> Wrapper for the LAPACK routine DGETRS. (Solves a system of linear
!> equations with a general matrix using the LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************

      SUBROUTINE WDGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : DGETRS_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*)
      INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
      CHARACTER(1) :: TRANS

!! IF (OFFLOAD_ON) THEN
!!    CALL DGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
!!    RETURN
!! ENDIF

      CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      RETURN
      END SUBROUTINE WDGETRS


!************************* SUBROUTINE ZGETRS **************************
!
!> Wrapper for the LAPACK routine ZGETRS. (Solves a system of linear
!> equations with a general matrix using the LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************

      SUBROUTINE WZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : ZGETRS_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*)
      INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
      CHARACTER(1) :: TRANS

!! IF (OFFLOAD_ON) THEN
!!    CALL ZGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
!!    RETURN
!! ENDIF

      CALL ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      RETURN
      END SUBROUTINE WZGETRS


!************************* SUBROUTINE DTRTRI **************************
!
!> Wrapper for the LAPACK routine DTRTRI. (Inverts a triangular matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!
!***********************************************************************

      SUBROUTINE WDTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : DTRTRI_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO,DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL DTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
!!    RETURN
!! ENDIF

      CALL DTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WDTRTRI


!************************* SUBROUTINE ZTRTRI **************************
!
!> Wrapper for the LAPACK routine ZTRTRI. (Inverts a triangular matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!
!***********************************************************************

      SUBROUTINE WZTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : ZTRTRI_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO,DIAG

!! IF (OFFLOAD_ON) THEN
!!    CALL ZTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
!!    RETURN
!! ENDIF

      CALL ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WZTRTRI


!************************* SUBROUTINE DSYGV ***************************
!
!> Wrapper for the LAPACK routine DSYGV. (Solves the generalized
!> eigenvalue problem for a real symmetric-definite matrix pair)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************

      SUBROUTINE WDSYGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : DSYGV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*),W(*),WORK(*)
      INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL DSYGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
!!    RETURN
!! ENDIF

      CALL DSYGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)

      RETURN
      END SUBROUTINE WDSYGV


!************************* SUBROUTINE ZHEGV ***************************
!
!> Wrapper for the LAPACK routine ZHEGV. (Solves the generalized
!> eigenvalue problem for a complex Hermitian-definite matrix pair)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************

      SUBROUTINE WZHEGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : ZHEGV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*),WORK(*)
      REAL(q) :: W(*),RWORK(*)
      INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL ZHEGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
!!    RETURN
!! ENDIF

      CALL ZHEGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)

      RETURN
      END SUBROUTINE WZHEGV


!************************* SUBROUTINE DSYEV ***************************
!
!> Wrapper for the LAPACK routine DSYEV. (Computes all eigenvalues and,
!> optionally, eigenvectors of a real symmetric matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************

      SUBROUTINE WDSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : DSYEV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),W(*),WORK(*)
      INTEGER :: N,LDA,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL DSYEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
!!    RETURN
!! ENDIF

      CALL DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)

      RETURN
      END SUBROUTINE WDSYEV


!************************* SUBROUTINE ZHEEV ***************************
!
!> Wrapper for the LAPACK routine ZHEEV. (Computes all eigenvalues and,
!> optionally, eigenvectors of a complex Hermitian matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************

      SUBROUTINE WZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : ZHEEV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),WORK(*)
      REAL(q) :: W(*),RWORK(*)
      INTEGER :: N,LDA,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL ZHEEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
!!    RETURN
!! ENDIF

      CALL ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)

      RETURN
      END SUBROUTINE WZHEEV


!************************* SUBROUTINE DSYEVX **************************
!
!> Wrapper for the LAPACK routine DSYEVX. (Computes selected eigenvalues
!> and, optionally, eigenvectors of a real symmetric matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************

      SUBROUTINE WDSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : DSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*),VL,VU,ABSTOL
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL DSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
!!    RETURN
!! ENDIF

      CALL DSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WDSYEVX


!************************* SUBROUTINE CHEEVX **************************
!
!>
!
!***********************************************************************

      SUBROUTINE WCHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : CHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(LDA,*),Z(LDZ,*),WORK(*)
      REAL(qs) :: W(*),RWORK(*),VL,VU,ABSTOL
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL CHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
!!    RETURN
!! ENDIF

      CALL CHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WCHEEVX


!************************* SUBROUTINE ZHEEVX **************************
!
!> Wrapper for the LAPACK routine ZHEEVX. (Computes selected eigenvalues
!> and, optionally, eigenvectors of a complex Hermitian matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and Z must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************

      SUBROUTINE WZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : ZHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),Z(LDZ,*),WORK(*)
      REAL(q) :: W(*),RWORK(*),VL,VU,ABSTOL
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL ZHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
!!    RETURN
!! ENDIF

      CALL ZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WZHEEVX


!************************* SUBROUTINE SSYEVX **************************
!
!> Wrapper for the LAPACK routine SSYEVX. (Computes selected eigenvalues
!> and, optionally, eigenvectors of a real symmetric matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************

      SUBROUTINE WSSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      USE prec
!! USE moffload_struct_def
!! USE moffload_lapack, ONLY : SSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*)
      REAL(qs) :: ABSTOL,VL,VU
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

!! IF (OFFLOAD_ON) THEN
!!    CALL SSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
!!    RETURN
!! ENDIF

      CALL SSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WSSYEVX


!***********************************************************************
! 1
!***********************************************************************
