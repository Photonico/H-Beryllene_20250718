# 1 "nonl_high.F"
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


# 2 "nonl_high.F" 2 
!***********************************************************************
!
!> This module contains all high level routines to calculate
!> the non local projection operators, either
!> * on a plane wave grid
!> * on a equally spaced grid
!
!***********************************************************************
MODULE nonl_high
  USE prec
  USE nonlr
  USE nonl
  IMPLICIT NONE

  CONTAINS

!************************** SUBROUTINE PROALL **************************
!
!> Calculates the scalar product of the current wavefunctions
!> stored in W with the projectors
!> ~~~
!>  C_lme ion,n = < b_lme ion | psi_n >
!> ~~~
!> and stores the result in W%CPROJ (wave function character)
!
!***********************************************************************

  SUBROUTINE PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,W)

!! USE moffload

    USE poscar
    USE lattice
    
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
! local
    INTEGER NK

    

# 57


    DO NK=1,W%WDES%NKPTS

       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

       IF (NONLR_S%LREAL) THEN
          CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,W%WDES)
          CALL RPRO(NONLR_S,W%WDES,W,GRID,NK)
       ELSE
          CALL PHASE(W%WDES,NONL_S,NK)
          CALL PROJ(NONL_S,W%WDES,W,NK)
       ENDIF
    ENDDO

# 83


    

    RETURN
  END SUBROUTINE PROALL


!************************** SUBROUTINE W1_PROALL ***********************
!
!> Calculates the scalar product of the current wavefunctions
!> stored in W1(:) with the projectors
!> ~~~
!>  C_lme ion,n = < b_lme ion | psi_n >
!> ~~~
!> and stores the result in W1(:)%CPROJ (wave function character)
!
!***********************************************************************

  SUBROUTINE W1_PROJALL(WDES1, W1, NONLR_S, NONL_S, NMAX)

!! USE moffload_struct_def

    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:)
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    INTEGER, OPTIONAL :: NMAX
! local
    INTEGER NMAX_, NP
!$  INTEGER NSTRIP,NSTRIP_ACT
!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

    

    IF (PRESENT(NMAX)) THEN
       NMAX_=NMAX
    ELSE
       NMAX_=SIZE(W1)
    ENDIF
 
    IF ( NONLR_S%LREAL ) THEN
# 132

       IF (NMAX_ >1 ) THEN
          CALL RPROMU(NONLR_S, WDES1, W1, NMAX_, W1%LDO)
       ELSE
          DO NP=1,NMAX_
             IF (.NOT. W1(NP)%LDO) CYCLE
             CALL RPRO1(NONLR_S, WDES1, W1(NP))
          ENDDO
       ENDIF
    ELSE
# 148

       DO NP=1,NMAX_
          IF (.NOT. W1(NP)%LDO) CYCLE
          CALL PROJ1(NONL_S,WDES1,W1(NP))
       ENDDO
    ENDIF

    

  END SUBROUTINE W1_PROJALL


  SUBROUTINE W1_PROJ(W1, NONLR_S, NONL_S)

!! USE moffload_struct_def

    IMPLICIT NONE
    TYPE (wavefun1) ::  W1
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct)  NONL_S
! local
    IF (NONLR_S%LREAL) THEN
# 174

          CALL RPRO1(NONLR_S,W1%WDES1,W1)
# 178

    ELSE
       CALL PROJ1(NONL_S,W1%WDES1,W1)
    ENDIF
  END SUBROUTINE W1_PROJ

END MODULE nonl_high
