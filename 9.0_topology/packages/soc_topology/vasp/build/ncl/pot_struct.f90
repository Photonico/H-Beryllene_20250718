# 1 "pot_struct.F"
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


# 2 "pot_struct.F" 2 

module pot_struct_def
    use prec, only: q
    use mgrid_struct_def, only: grid_3d, transit

    implicit none

    type kernel_truncation
        logical :: active                   = .false. !< Decides if the Coulomb truncation method is used
        integer :: dimensionality           = 3       !< Dimensionality of the system.
!< 0 -> atoms and molecules
!< 1 -> nanorods and wires
!< 2 -> surfaces, slabs and 2D materials
!< 3 -> bulk crystals
        integer :: surface_normal_direction = -1      !< Direction in which the surface normal is oriented
        integer :: padding_factor           = huge(1) !< Number of cells that will be padded in any direction
        real(q) :: truncation_factor        = huge(1._q)   !< Truncate the Coulomb kernel with this factor
        real(q) :: truncation_length        = -1._q   !< Truncation length based on the truncation factor
        logical :: coarsen_before_pad       = .false. !< Coarsen the grid before padding it with the zeros in real space.
        logical :: use_simplified_2d_kernel = .false. !< if r_c = L_z/2, use the simplified version of the kernel
        type(grid_3d), pointer :: fine_grid         => NULL() !< The fine grid used in the padding procedure
        type(grid_3d), pointer :: coarse_grid       => NULL() !< The coarse grid used in the padding procedure
        type(transit), pointer :: translation_table => NULL() !< Translation table between fine and coarse grids
    end type

    type coulomb_potential
        type(kernel_truncation) :: kernel_truncate    !< Parameters for truncating the kernel
    end type

end module pot_struct_def
