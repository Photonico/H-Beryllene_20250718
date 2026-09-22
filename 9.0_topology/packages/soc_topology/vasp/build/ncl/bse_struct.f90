# 1 "bse_struct.F"
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


# 2 "bse_struct.F" 2 
MODULE bse_struct
  USE prec
  IMPLICIT NONE


  INTEGER :: BSE_DESC(9)         ! descriptor for the BSE Hamiltonian
  INTEGER :: VEC_DESC(9)         ! descriptor for non-distributed vector
  INTEGER :: DVEC_DESC(9)        ! descriptor for distributed vector

  TYPE banddesc
    INTEGER :: VBMIN, VBMAX      ! band index for the VB minimum and VB maximum
    INTEGER :: CBMIN, CBMAX      ! band index for the CB minimum and CB maximum
    INTEGER :: CBMIN4, CBMAX4    ! band index for the CB minimum and CB maximum for n4
! this might be reduced for parallelization over bands
    INTEGER :: NGLB              ! block size for conduction bands = total number of conduction bands
    INTEGER :: NGLB4             ! block size for conduction bands for fourth index (N4 usually)
! [ N4 is distributed over nodes when band parallelization is used
!   in this case NGLB4 is equal NGLB/ number of nodes]
  END TYPE banddesc

  TYPE bse_matrix_index
    INTEGER :: NCV                  !< total rank
    INTEGER, ALLOCATABLE ::  N1(:)  !< valence band index
    INTEGER, ALLOCATABLE ::  N3(:)  !< conduction band index
    INTEGER, ALLOCATABLE ::  NK(:)  !< k-point index
    INTEGER, ALLOCATABLE ::  ISP(:) !< spin index
!< index into Hamilton matrix for given valence, conduction, k-point and spin
    INTEGER, ALLOCATABLE ::  INDEX(:,:,:,:)

! used to index full matrix
! column indices K1,K2,ISP1,ISP2,NPOS1, and NPOS3
! index along (1._q,0._q) axis columns
    INTEGER :: COLUMNS
    INTEGER, ALLOCATABLE ::  CK1(:)
    INTEGER, ALLOCATABLE ::  CISP1(:)
    INTEGER, ALLOCATABLE ::  CNPOS1(:)
    INTEGER, ALLOCATABLE ::  CNPOS3(:)

! index in full matrix
! block indices K1,K2,ISP1,ISP2,NPOS1, and NPOS2
    INTEGER :: BLOCKS
    INTEGER, ALLOCATABLE ::  BK1(:)
    INTEGER, ALLOCATABLE ::  BK2(:)
    INTEGER, ALLOCATABLE ::  BISP1(:)
    INTEGER, ALLOCATABLE ::  BISP2(:)
    INTEGER, ALLOCATABLE ::  BNPOS1(:)
    INTEGER, ALLOCATABLE ::  BNPOS2(:)
    INTEGER, ALLOCATABLE ::  BNPOS3(:)
    INTEGER, ALLOCATABLE ::  BNPOS4(:)

! map indices from compact triangular to full ractangular matrix
! this way we can easily find starting and ending position for mask
    INTEGER, ALLOCATABLE ::  BMAP(:)

! mask empty rows/columns
    INTEGER, ALLOCATABLE ::  MASK(:,:)
  END TYPE bse_matrix_index


  TYPE selfenergy_from_bse
     INTEGER :: NOMEGA
     INTEGER :: NBANDS
     INTEGER :: NKPTS
     INTEGER :: ISPIN
     REAL(q), POINTER :: OMEGA(:)  => NULL()
     REAL(q), POINTER :: REAL_PART(:,:,:,:) => NULL() !1st index over frequencies, 2nd index over bands
     REAL(q), POINTER :: IMAG_PART(:,:,:,:) => NULL()!3rd index over kpoints, 4th index over spin
     COMPLEX(q), POINTER :: CELNEW(:,:,:) => NULL()  !bands, kpoints, spin
  END TYPE selfenergy_from_bse

END MODULE bse_struct
