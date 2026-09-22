# 1 "extpot.F"
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


# 2 "extpot.F" 2 
!***********************************************************************
!> These are stubs for the Yu/Libisch/Carter PAW DFET embedding routines
!***********************************************************************
module mextpot

    implicit none

contains

   logical function EXTPT_LEXTPOT()
      EXTPT_LEXTPOT = .FALSE.
   end function EXTPT_LEXTPOT

   logical function EXTPT_DEXTPOT()
      EXTPT_DEXTPOT = .FALSE.
   end function EXTPT_DEXTPOT

   subroutine EXTPT_READER(NGXF,NGYF,NGZF,IU0,IU5,IU6)
      integer :: NGXF, NGYF, NGZF, IU0, IU5, IU6
   end subroutine EXTPT_READER

   subroutine EXTPT_EXTERNAL_POT_ADD(GRIDC, LATT_CUR, CVTOT)
      USE prec
      USE lattice
      USE mgrid
      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q) ::      CVTOT(GRIDC%MPLWV)
   end subroutine EXTPT_EXTERNAL_POT_ADD

   subroutine EXTPT_EXTERNAL_POT_ADD_PAW(POTAE, POT, NDIM, NCDIJ, LMMAX, NIP)
      USE prec
      INTEGER             NDIM, NCDIJ, LMMAX, NIP
      REAL(q) POT(:,:,:), POTAE(:,:,:)
   end subroutine EXTPT_EXTERNAL_POT_ADD_PAW

   subroutine EXTPT_CALC_VLM(LATT_CUR, GRIDC, T_INFO, P, LMDIM, WDES, IU0)
      use lattice
      use mgrid
      use poscar
      use pseudo
      use wave
      type (grid_3d)    GRIDC
      type (latt)       LATT_CUR
      type (type_info)  T_INFO
      type (potcar), TARGET :: P(T_INFO%NTYP)
      type (wavedes)    WDES
      integer :: LMDIM, IU0
   end subroutine EXTPT_CALC_VLM

   subroutine EXTPT_DE_DVRN(WDES, GRID_SOFT, GRIDC, SOFT_TO_C, LATT_CUR, &
        P, T_INFO, LMDIM, CRHODE, CHTOT_, CHDEN, IRDMAX )
      use prec
      use mgrid
      use lattice
      use poscar
      use pseudo
      use wave
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      TYPE (grid_3d)     GRID_SOFT
      TYPE (latt)        LATT_CUR
      TYPE (transit)     SOFT_TO_C
      TYPE (wavedes)     WDES
      integer :: LMDIM, IRDMAX
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      COMPLEX(q)   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q), TARGET :: CHTOT_(GRIDC%MPLWV,WDES%NCDIJ)
   end subroutine EXTPT_DE_DVRN

end module mextpot
