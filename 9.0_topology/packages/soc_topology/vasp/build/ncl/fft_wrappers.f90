# 1 "fft_wrappers.F"
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


# 2 "fft_wrappers.F" 2 

!************************* SUBROUTINE FFTBAS_PLAN **************************
!
!>
!
!***********************************************************************

      SUBROUTINE WFFTBAS(C, GRID, ISN)

      USE prec
      USE mgrid_struct_def
!! USE moffload_struct_def
!! USE moffload_fft, ONLY : FFTBAS_PLAN_OFFLOAD
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISN

!! IF (OFFLOAD_ON) THEN
!! IF (ACC_IS_PRESENT(C) .AND. OFFLOAD_ON) THEN
!!    CALL FFTBAS_PLAN_OFFLOAD(C, GRID, ISN)
!!    RETURN
!! ENDIF

      CALL FFTBAS_PLAN(C, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBAS


!************************* SUBROUTINE FFTBRC_PLAN **************************
!
!>
!
!***********************************************************************

      SUBROUTINE WFFTBRC(C, GRID, ISN)

      USE prec
      USE mgrid_struct_def
!! USE moffload_struct_def
!! USE moffload_fft, ONLY : FFTBRC_PLAN_OFFLOAD

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISN

!! IF (OFFLOAD_ON) THEN
!! IF (ACC_IS_PRESENT(C) .AND. OFFLOAD_ON) THEN
!!    CALL FFTBRC_PLAN_OFFLOAD(C, GRID, ISN)
!!    RETURN
!! ENDIF

      CALL FFTBRC_PLAN(C, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBRC


!************************* SUBROUTINE FFTBAS_PLAN_MU ***********************
!
!>
!
!***********************************************************************

      SUBROUTINE WFFTBAS_MU(N, C, LDC, GRID, ISN)

      USE prec
      USE mgrid_struct_def
!! USE moffload_struct_def
!! USE moffload_fft, ONLY : FFTBAS_PLAN_MU_OFFLOAD

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N, LDC, ISN

!! IF (OFFLOAD_ON) THEN
!! IF (ACC_IS_PRESENT(C) .AND. OFFLOAD_ON) THEN
!!    CALL FFTBAS_PLAN_MU_OFFLOAD(N, C, LDC, GRID, ISN)
!!    RETURN
!! ENDIF

      CALL FFTBAS_PLAN_MU(N, C, LDC, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBAS_MU


!************************* SUBROUTINE FFTBRC_PLAN_MU ***********************
!
!>
!
!***********************************************************************

      SUBROUTINE WFFTBRC_MU(N, C, LDC, GRID, ISN)

      USE prec
      USE mgrid_struct_def
!! USE moffload_struct_def
!! USE moffload_fft, ONLY : FFTBRC_PLAN_MU_OFFLOAD

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N, LDC, ISN

!! IF (OFFLOAD_ON) THEN
!! IF (ACC_IS_PRESENT(C) .AND. OFFLOAD_ON) THEN
!!    CALL FFTBRC_PLAN_MU_OFFLOAD(N, C, LDC, GRID, ISN)
!!    RETURN
!! ENDIF

      CALL FFTBRC_PLAN_MU(N, C, LDC, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBRC_MU


!************************* SUBROUTINE FFTBAS_PLAN_MPI **********************
!
!>
!
!***********************************************************************

      SUBROUTINE WFFTBAS_MPI(A, GRID, ISN)

      USE prec
      USE mgrid_struct_def
!! USE moffload_struct_def
!! USE moffload_fft, ONLY : FFTBAS_PLAN_MPI_OFFLOAD

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: A(*)
      INTEGER :: ISN

!! IF (OFFLOAD_ON) THEN
!! IF (ACC_IS_PRESENT(A) .AND. OFFLOAD_ON) THEN
!!    CALL FFTBAS_PLAN_MPI_OFFLOAD(A, GRID, ISN)
!!    RETURN
!! ENDIF

      CALL FFTBAS_PLAN_MPI(A, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBAS_MPI


!************************* SUBROUTINE FFTBRC_PLAN_MPI **********************
!
!>
!
!***********************************************************************

      SUBROUTINE WFFTBRC_MPI(A, GRID, ISN)

      USE prec
      USE mgrid_struct_def
!! USE moffload_struct_def
!! USE moffload_fft, ONLY : FFTBRC_PLAN_MPI_OFFLOAD

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: A(*)
      INTEGER :: ISN

!! IF (OFFLOAD_ON) THEN
!! IF (ACC_IS_PRESENT(A) .AND. OFFLOAD_ON) THEN
!!    CALL FFTBRC_PLAN_MPI_OFFLOAD(A, GRID, ISN)
!!    RETURN
!! ENDIF

      CALL FFTBRC_PLAN_MPI(A, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBRC_MPI

