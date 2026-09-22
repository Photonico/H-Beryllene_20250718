# 1 "esf_struct.F"
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


# 2 "esf_struct.F" 2 

MODULE esf_struct_def

  USE prec
  USE poscar_struct_def, ONLY: latt
  USE mkpoints_struct_def, ONLY: skpoints_trans 
  USE wave_struct_def, ONLY: wavedes
  IMPLICIT NONE

! handle for electronic structure factor
  TYPE elstufac_handle
!> number of k-points in x-direction
     INTEGER                 :: NKPX = 0
!> number of k-points in y-direction
     INTEGER                 :: NKPY = 0
!> number of k-points in z-direction
     INTEGER                 :: NKPZ = 0
!> wavefunction descriptor of response function
     TYPE (wavedes), POINTER ::  WGW => NULL()
!> number of irreducible k-points in 1. BZ
     INTEGER                 :: NUMBER_OF_NQ = 0 
!> is symmetry switched on or off
     INTEGER                 :: ISYM = -1
!> kpoint map from IRZ to full BZ
     TYPE( skpoints_trans )  :: KPOINTS_TRANS
!> number of G-vectors
     INTEGER, ALLOCATABLE    :: NG(:)         
!> electronic structure factor
!> similar to diagonal of frequency-integrated responsefunction
!> S(q+G)=\int dw xi(q+G,q+G,w)
!> 2nd and 3rd index correspond to cutoff and lambda value considered
!> electronic stucture factor in q+G domain
     REAL(q),ALLOCATABLE     :: S(:,:)

!> energy cutoff
     REAL(q)                 :: ENCUT=0._q
!> soft energy cutoff
     REAL(q)                 :: ENCUTSOFT =0._q
!> number of G-vectors in 1st direction
     INTEGER                 :: N1 = 0         
!> number of G-vectors in 2nd direction
     INTEGER                 :: N2 = 0    
!> number of G-vectors in 3rd direction
     INTEGER                 :: N3 = 0      
!> cell volume
     REAL(q)                 :: OMEGA=0._q
!> reciprocal Bravais matrix
     REAL(q)                 :: B(3,3)
!> Bravais matrix
     REAL(q)                 :: A(3,3)
!> cutoff
!> splined result
     REAL(q)                 :: ESPLINED=0._q
!> integer error code
     INTEGER                 :: ISTATUS = -1

     TYPE(latt)              :: LATT_CUR       !< stores lattice information

!> singularity correction for exchange potential requires number of electrons
     REAL(q)                 :: NELECT = -1._q 
 
!> auxillary arrays that probably should be removed
     REAL(q),ALLOCATABLE     :: V(:)           !< Coulomb potential in PW basis
     INTEGER,ALLOCATABLE     :: OLD_INDEX(:)   !< stores old index in PW array that corresponds to NE=1 cutoff
  END TYPE elstufac_handle 

END MODULE esf_struct_def 
