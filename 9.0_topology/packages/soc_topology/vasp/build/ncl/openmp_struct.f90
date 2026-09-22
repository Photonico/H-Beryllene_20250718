# 1 "openmp_struct.F"
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


# 2 "openmp_struct.F" 2 
     MODULE mopenmp_struct_def
! overall threading strategy

!> @ref openmp :
!> number of threads to be used in \"general\":
!> by default this is set to the maximum number of available threads,
!> mopenmp::omp_nthreads=omp_get_max_threads()
      INTEGER :: omp_nthreads = 1

! specific routines

!> @ref openmp :
!> number of threads used in ::m_alltoall_d_omp
!> (experimental and unused for now).
      INTEGER :: omp_nthreads_alltoall = 1

!> @ref openmp :
!> number of threads used in the domain decomposition of the
!> real space projection operators.
!> Default: mopenmp::omp_nthreads_nonlr_rspace=mopenmp::omp_nthreads
      INTEGER :: omp_nthreads_nonlr_rspace = 1

!> @ref openmp :
!> number of threads to be used in connection with OpenACC:
!> in most instances this involves a distribution over orbitals
      INTEGER :: omp_nthreads_acc = 1

!> @ref openmp :
!> if set to true the domain of the real space projection operators
!> is decomposed over (x,y)-planes, in a round-robin fashion, if set
!> to false the projectors are distributed over (z,y)-columns.
      LOGICAL :: omp_nonlr_planewise = .TRUE.

!> @ref openmp :
!> if set to true ::fftmakeplan or ::fftmakeplan_mpi will call
!> dfftw_init_threads and set mopenmp::omp_dfftw_init_threads=.FALSE.
      LOGICAL :: omp_dfftw_init_threads = .TRUE.

# 44

     END MODULE mopenmp_struct_def
