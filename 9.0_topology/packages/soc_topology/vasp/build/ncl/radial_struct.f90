# 1 "radial_struct.F"
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


# 2 "radial_struct.F" 2 
MODULE radial_struct_def
    USE prec
! structure which is used for the logarithmic grid
! the grid points are given by R(i) = RSTART * exp [H (i-1)]
    TYPE rgrid
       REAL(q)  :: RSTART                   ! starting point
       REAL(q)  :: REND                     ! endpoint
       REAL(q)  :: RMAX                     ! radius of augmentation sphere
       REAL(q)  :: D                        ! R(N+1)/R(N) = exp(H)
       REAL(q)  :: H                        !
       REAL(q),POINTER :: R(:)  => NULL()   ! radial grid (r-grid)
       REAL(q),POINTER :: SI(:) => NULL()   ! integration prefactors on r-grid
       INTEGER  :: NMAX                     ! number of grid points
    END TYPE rgrid

! This parameter determines at which magnetization the aspherical contributions
! to the (1._q,0._q) center magnetization are truncated in the non collinear case
!   Without any truncation the aspherical terms for non magnetic atoms
! tend to yield spurious but meaningless contributions to the potential
! so that convergence to the groundstate can not be achieved
! for details see the routines RAD_MAG_DIRECTION and RAD_MAG_DENSITY
    REAL(q), PARAMETER :: MAGMIN=1E-2

! for non collinear calculations, setting
! the parameter USE_AVERAGE_MAGNETISATION  means that the aspherical
! contributions to the (1._q,0._q) center magnetisation are projected onto the
! average magnetization direction in the PAW sphere instead of the
! local moment of the spherical magnetization density at
! each grid-point
! USE_AVERAGE_MAGNETISATION improves the numerical stability significantly
! and must be set
    LOGICAL :: USE_AVERAGE_MAGNETISATION=.TRUE.
END MODULE radial_struct_def
