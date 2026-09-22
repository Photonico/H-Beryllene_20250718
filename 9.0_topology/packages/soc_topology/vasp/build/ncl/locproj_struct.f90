# 1 "locproj_struct.F"
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


# 2 "locproj_struct.F" 2 
MODULE locproj_struct
  USE prec

!
!> Data structure that stores the definitions of the local functions
!> on which we will project our Bloch orbitals.
!
  TYPE LPRJ_function
! radial part
     INTEGER radial_type    !< 1="optimal" PAW projector, 2=PAW partial wave (PS), 3=Hydrogen-like (=wannier90)
     INTEGER species        !< index of the species at which the local orbital is located
     INTEGER n              !< principal qunatum number
     REAL(q) za             !< real number determining how large the radius of the hydrogenic orbitals should be
     LOGICAL WANNIER90_ORBITAL_DEFINITIONS !< logical governing wether to use wannier90 or LOCPROJ orbital definitions
! spherical part
     INTEGER l              !< l angular quantum number
     INTEGER m              !< m magnetic quantum number
     LOGICAL LROTYLM        !< logical to rotate the radial part of the orbitals
     REAL(q) proj_x(3)      !< x axis to use for the rotated projections
     REAL(q) proj_z(3)      !< z axis to use for the rotated projections
! site part
     INTEGER poscar_site    !< index of the position in the poscar file
     REAL(q) R(3)           !< position of the localized orbital in reduced coordinates
! spin part
     INTEGER spinor         !< spinor index 1 for spin up and 2 for spin down
     REAL(q) spin_qaxis(3)  !< spin quantization axis
! related to PS or PR projectors
     INTEGER iproj          !< index of projector or partial wave to be projected on
     INTEGER ibase          !< first pseudopotential channel with angular moment L=l
     INTEGER nproj          !< number of projectors with same L

!> contains the optimized projectors computed by SPHPRO or SPHPRO_OPTPROJ
!> and stored in the pseudo datastructure
     REAL(q), POINTER :: optproj(:)

!> spline representation of the radial part
     REAL(q), ALLOCATABLE :: SPLINE(:,:,:)
  END TYPE LPRJ_function
END MODULE locproj_struct
