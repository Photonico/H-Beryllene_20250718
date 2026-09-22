# 1 "nl_struct.F"
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


# 2 "nl_struct.F" 2 
MODULE nonl_struct_def
  USE prec
!
!>  Structure required to support non local projection operators in recip space
!
  TYPE nonl_struct
!only NONL_S
     LOGICAL LRECIP                !< structure set up ?
     INTEGER NTYP                  !< number of types
     INTEGER NIONS                 !< number of ions
     INTEGER NK                    !< kpoint for which CREXP is set up
     INTEGER SELECTED_ION          !< allows to generate a projector for a single ion
     INTEGER, POINTER :: NITYP(:)  !< number of ions for each type
     INTEGER, POINTER :: LMMAX(:)  !< max l-quantum number for each type
     LOGICAL LSPIRAL               !< do we want to calculate spin spirals
!    REAL(q), POINTER, CONTIGUOUS ::QPROJ(:,:,:,:,:) ! projectors in reciprocal space for each k-point
     REAL(q), POINTER ::QPROJ(:,:,:,:,:)
     COMPLEX(q),POINTER, CONTIGUOUS ::CREXP(:,:)  !< phase factor exp (i (G+k) R(ion))
     REAL(q),  POINTER  ::POSION(:,:) !< positions (required for setup)
     REAL(q),  POINTER  ::VKPT_SHIFT(:,:) => NULL() !< k-point shift for each ion
     COMPLEX(q),POINTER ::CQFAK(:,:)  !< i^l

!> complete set of projectors in reciprocal space for a single k-point
!> in over-plane-wave distribution
     COMPLEX(q), POINTER, CONTIGUOUS :: CPROJK(:,:) => NULL()

!> these arrays allow loops over all ions to be restructured from
!> nested "over-types" + "over-ions-of-type" to straightforward
!> loops "over-all-ions" (used in many subroutines).
     INTEGER, POINTER :: ITYP(:)
!> these arrays allow loops over all ions to be restructured from
!> nested "over-types" + "over-ions-of-type" to straightforward
!> loops "over-all-ions" (used in many subroutines).
     INTEGER, POINTER :: LMBASE(:)
  END TYPE nonl_struct
END MODULE nonl_struct_def

MODULE nonlr_struct_def
  USE prec
  USE mpimy
  USE mpi_shmem
!
!>  structures required to support non local projection operators in real space
!

  TYPE nonlr_proj
     REAL(q), POINTER :: PSPRNL(:,:,:)
     INTEGER, POINTER :: LPS(:)
  END TYPE nonlr_proj

!
!> @brief
!> The data structure of non local projection operators in real space
!
!> @details @ref openmp :
!> Under OpenMP the members nonlr_struct::nlimax, nonlr_struct::nli,
!> nonlr_struct::rproj, nonlr_struct::crrexp, and nonlr_struct::nlibase
!> acquire and additional dimension.
!
  TYPE nonlr_struct
!only NONLR_S
     LOGICAL LREAL                     !< structure set up ?
     INTEGER NTYP                      !< number of types
     INTEGER NIONS                     !< number of ions
     INTEGER SELECTED_ION              !< allows to generate a projector for a single ion
     INTEGER IRMAX                     !< maximum number points in sphere
     INTEGER IRALLOC                   !< size for allocation =IRMAX*LMDIM*NIONS
     INTEGER NK                        !< kpoint for which CRREXP is set up
     INTEGER, POINTER :: NITYP(:)      !< number of ions for each type
     INTEGER, POINTER :: ITYP(:)       !< type for each ion
     INTEGER, POINTER :: LMAX(:)       !< max l-quantum number for each type
     INTEGER, POINTER :: LMMAX(:)      !< number lmn-quantum numbers for each type
     INTEGER, POINTER ::CHANNELS(:)    !< number of ln-quantum for each type
     REAL(q), POINTER :: PSRMAX(:)     !< real space cutoff
     REAL(q), POINTER :: RSMOOTH(:)    !< radius for smoothing the projectors around each point
     REAL(q), POINTER :: POSION(:,:)   !< positions (required for setup)
     REAL(q), POINTER :: VKPT_SHIFT(:,:)  !< k-point shift for each ion
     TYPE(nonlr_proj), POINTER :: BETA(:) !< a set of structures containing pointers to
     LOGICAL LSPIRAL               !< do we want to calculate spin spirals

     INTEGER, POINTER :: NLIMAX(:     ) !< maximum index for each ion
     INTEGER, POINTER :: NLI   (:,:   ) !< index for gridpoints
     REAL(q),POINTER, CONTIGUOUS :: RPROJ (:     ) !< projectors on real space grid
     COMPLEX(q),POINTER, CONTIGUOUS::CRREXP(:,:,: ) !< phase factor exp (i k (R(ion)-r(grid)))

!> these arrays allow loops over all ions to be restructured from
!> nested "over-types" + "over-ions-of-type" to straightforward
!> loops "over-all-ions" (used in many subroutines).
     INTEGER, POINTER :: LMBASE(:)
!> these arrays allow loops over all ions to be restructured from
!> nested "over-types" + "over-ions-of-type" to straightforward
!> loops "over-all-ions" (used in many subroutines).
     INTEGER, POINTER :: NLIBASE(: )
# 100

  END TYPE nonlr_struct

  TYPE smoothing_handle
     INTEGER :: N                            !< number of grid points
     REAL(q), POINTER :: WEIGHT(:)           !< weight of each grid point
     REAL(q), POINTER :: X1(:), X2(:), X3(:) !< positions of additional grid points in fractional coordinates
  END TYPE smoothing_handle
END MODULE nonlr_struct_def
