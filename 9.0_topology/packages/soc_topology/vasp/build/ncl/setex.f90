# 1 "setex.F"
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


# 2 "setex.F" 2 
# 4


  MODULE setexm
    USE prec
    USE setexm_struct_def
!
! the exchange stack can be used to save
! the present exchange parameters temporarily
!
    INTEGER, SAVE, PRIVATE :: ISTACK=0,I_XC_DATA_STACK_2=0
    TYPE (xc_info), SAVE, PRIVATE :: XC_DATA_STACK(5)

    CONTAINS

!******************* SUBROUTINE SET_XC_DATA ***********************
!
!> This subroutine determines the information about the functional,
!> like the internal ID, the family (LDA, GGA, or meta-GGA) or
!> the kind (exchange, correlation or exchange-correlation),
!> and store it in XC_DATA.
!
!******************************************************************

    SUBROUTINE SET_XC_DATA(XC)
      USE prec
      USE main_mpi
      USE xc_name
      USE xc_lda_table_name
      USE xc_family_name
# 35

      USE tutor, ONLY : vtutor
      IMPLICIT NONE
      INTEGER :: I,N
      CHARACTER(LEN=40) :: FXCNAME
      TYPE (xc_info) :: XC
# 43


      DO I=1,XC%NXC
         FXCNAME=TRIM(ADJUSTL(XC%NAME(I)))

!no XC
         IF (FXCNAME=='NO') THEN
            XC%ID(I)=ID_XC_NOXC
            XC%FAMILY(I)=XC_FAM_NONE
            XC%KIND(I)=-1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC

!Coulomb hole
         ELSEIF (FXCNAME=='CO') THEN
            XC%ID(I)=ID_XC_COHSM
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC

!LDA
         ELSEIF (FXCNAME=='SL') THEN
            XC%ID(I)=ID_XC_SLATER
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_SLATER
         ELSEIF ((FXCNAME=='PZ_C').OR.(FXCNAME=='CA_C')) THEN
            XC%ID(I)=ID_XC_PZ_C
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PZ_C
         ELSEIF ((FXCNAME=='PZ').OR.(FXCNAME=='CA')) THEN
            XC%ID(I)=ID_XC_PZ
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PZ
         ELSEIF (FXCNAME=='HL') THEN
            XC%ID(I)=ID_XC_HL
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_HL
         ELSEIF (FXCNAME=='VW') THEN
            XC%ID(I)=ID_XC_VWN5
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_VWN5
         ELSEIF (FXCNAME=='WI') THEN
            XC%ID(I)=ID_XC_WI
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_WI
         ELSEIF (FXCNAME=='PW92') THEN
            XC%ID(I)=ID_XC_PW92
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92_RELA
         ELSEIF (FXCNAME=='PW92_C') THEN
            XC%ID(I)=ID_XC_PW92_C
            XC%FAMILY(I)=XC_FAM_LDA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92_C

!GGA
         ELSEIF (FXCNAME=='91') THEN
            XC%ID(I)=ID_XC_PW91
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PZ
         ELSEIF (FXCNAME=='AM') THEN
            XC%ID(I)=ID_XC_AM05
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='BO') THEN
            XC%ID(I)=ID_XC_OPTB88PBE
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='MK') THEN
            XC%ID(I)=ID_XC_OPTB86BPBE
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='ML') THEN
            XC%ID(I)=ID_XC_PW86RPBE
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='CX') THEN
            XC%ID(I)=ID_XC_CX
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='OR') THEN
            XC%ID(I)=ID_XC_OPTPBE
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='PE') THEN
            XC%ID(I)=ID_XC_PBE
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='PBE_X') THEN
            XC%ID(I)=ID_XC_PBE_X
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_X
         ELSEIF (FXCNAME=='PBE_C') THEN
            XC%ID(I)=ID_XC_PBE_C
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92_C
         ELSEIF (FXCNAME=='PS') THEN
            XC%ID(I)=ID_XC_PBESOL
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='RE') THEN
            XC%ID(I)=ID_XC_REVPBE
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='RP') THEN
            XC%ID(I)=ID_XC_RPBE
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_PW92
         ELSEIF (FXCNAME=='B3') THEN
            XC%ID(I)=ID_XC_B3
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_B3
         ELSEIF (FXCNAME=='B5') THEN
            XC%ID(I)=ID_XC_B5
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_B5
         ELSEIF (FXCNAME=='RA') THEN
            XC%ID(I)=ID_XC_ACFDT_RA
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_ACFDT_RA
         ELSEIF (FXCNAME=='PL') THEN
            XC%ID(I)=ID_XC_ACFDT_PL
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_ACFDT_PL
         ELSEIF (FXCNAME=='03') THEN
            XC%ID(I)=ID_XC_ACFDT_03
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_ACFDT_03
         ELSEIF (FXCNAME=='05') THEN
            XC%ID(I)=ID_XC_ACFDT_05
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_ACFDT_05
         ELSEIF (FXCNAME=='10') THEN
            XC%ID(I)=ID_XC_ACFDT_10
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_ACFDT_10
         ELSEIF (FXCNAME=='20') THEN
            XC%ID(I)=ID_XC_ACFDT_20
            XC%FAMILY(I)=XC_FAM_GGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_ACFDT_20
         ELSEIF (FXCNAME=='BF') THEN
# 223

            CALL vtutor%error("VASP needs to be linked against libbeef for Bayesian error estimation &
               &functional support.\nlibbeef sources and binaries can be downloaded from &
               &suncat.stanford.edu")


!MGGA
         ELSEIF (FXCNAME=='M06L') THEN
            XC%ID(I)=ID_XC_M06L
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='M06L_X') THEN
            XC%ID(I)=ID_XC_M06L_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='M06L_C') THEN
            XC%ID(I)=ID_XC_M06L_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='MBJ') THEN
            XC%ID(I)=ID_XC_MBJ
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LLAP=.TRUE.
            IF (XC%NXC>=2) CALL vtutor%error("Error: MBJ can not be combined with another functional.")
         ELSEIF (FXCNAME=='LMBJ') THEN
            XC%ID(I)=ID_XC_LMBJ
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LLAP=.TRUE.
            IF (XC%NXC>=2) CALL vtutor%error("Error: LMBJ can not be combined with another functional.")
         ELSEIF (FXCNAME=='MS0') THEN
            XC%ID(I)=ID_XC_MS0
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.29_q    !MSX_RKAPPA
            XC%PARAM(I,2)=0.28771_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q     !MSX_CFE
            XC%PARAM(I,4)=0.0_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MS0_X') THEN
            XC%ID(I)=ID_XC_MS0_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.29_q    !MSX_RKAPPA
            XC%PARAM(I,2)=0.28771_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q     !MSX_CFE
            XC%PARAM(I,4)=0.0_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MS0_C') THEN
            XC%ID(I)=ID_XC_MS0_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='MS1') THEN
            XC%ID(I)=ID_XC_MS1
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.404_q   !MSX_RKAPPA
            XC%PARAM(I,2)=0.18150_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q     !MSX_CFE
            XC%PARAM(I,4)=0.0_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MS1_X') THEN
            XC%ID(I)=ID_XC_MS1_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.404_q   !MSX_RKAPPA
            XC%PARAM(I,2)=0.18150_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q     !MSX_CFE
            XC%PARAM(I,4)=0.0_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MS1_C') THEN
            XC%ID(I)=ID_XC_MS1_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='MS2') THEN
            XC%ID(I)=ID_XC_MS2
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.504_q   !MSX_RKAPPA
            XC%PARAM(I,2)=0.14601_q !MSX_CFC
            XC%PARAM(I,3)=4.0_q     !MSX_CFE
            XC%PARAM(I,4)=0.0_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MS2_X') THEN
            XC%ID(I)=ID_XC_MS2_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.504_q   !MSX_RKAPPA
            XC%PARAM(I,2)=0.14601_q !MSX_CFC
            XC%PARAM(I,3)=4.0_q     !MSX_CFE
            XC%PARAM(I,4)=0.0_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MS2_C') THEN
            XC%ID(I)=ID_XC_MS2_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='MSPBEL') THEN
            XC%ID(I)=ID_XC_MSPBEL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.804_q     !MSX_RKAPPA
            XC%PARAM(I,2)=0.1036199_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q       !MSX_CFE
            XC%PARAM(I,4)=0.0_q       !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MSRPBEL') THEN
            XC%ID(I)=ID_XC_MSRPBEL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.804_q     !MSX_RKAPPA
            XC%PARAM(I,2)=0.0767086_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q       !MSX_CFE
            XC%PARAM(I,4)=0.0_q       !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='MSB86BL') THEN
            XC%ID(I)=ID_XC_MSB86BL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.804_q      !MSX_RKAPPA
            XC%PARAM(I,2)=0.08809161_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q        !MSX_CFE
            XC%PARAM(I,4)=0.0_q        !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='RMSPBEL') THEN
            XC%ID(I)=ID_XC_RMSPBEL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.804_q     !MSX_RKAPPA
            XC%PARAM(I,2)=0.1036199_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q       !MSX_CFE
            XC%PARAM(I,4)=0.001_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='RMSRPBEL') THEN
            XC%ID(I)=ID_XC_RMSRPBEL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.804_q     !MSX_RKAPPA
            XC%PARAM(I,2)=0.0767086_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q       !MSX_CFE
            XC%PARAM(I,4)=0.001_q     !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='RMSB86BL') THEN
            XC%ID(I)=ID_XC_RMSB86BL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.804_q      !MSX_RKAPPA
            XC%PARAM(I,2)=0.08809161_q !MSX_CFC
            XC%PARAM(I,3)=1.0_q        !MSX_CFE
            XC%PARAM(I,4)=0.001_q      !MSX_TAUREG_X
            XC%NPARAM(I)=4
         ELSEIF (FXCNAME=='R2SCAN') THEN
            XC%ID(I)=ID_XC_R2SCAN
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='R2SCAN_X') THEN
            XC%ID(I)=ID_XC_R2SCAN_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='R2SCAN_C') THEN
            XC%ID(I)=ID_XC_R2SCAN_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='RSCAN') THEN
            XC%ID(I)=ID_XC_RSCAN
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.0001_q !SCAN_TAUREG_X
            XC%PARAM(I,2)=0.0001_q !SCAN_TAUREG_C
            XC%PARAM(I,3)=0.001_q  !SCAN_ALPREG
            XC%NPARAM(I)=3
         ELSEIF (FXCNAME=='RSCAN_X') THEN
            XC%ID(I)=ID_XC_RSCAN_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.0001_q !SCAN_TAUREG_X
            XC%PARAM(I,2)=0.001_q  !SCAN_ALPREG
            XC%NPARAM(I)=2
         ELSEIF (FXCNAME=='RSCAN_C') THEN
            XC%ID(I)=ID_XC_RSCAN_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.0001_q !SCAN_TAUREG_C
            XC%PARAM(I,2)=0.001_q  !SCAN_ALPREG
            XC%NPARAM(I)=2
         ELSEIF (FXCNAME=='RSCANZT') THEN
            XC%ID(I)=ID_XC_RSCANZT
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=1.998902376E-13_q !SCAN_TAUREG_X
            XC%PARAM(I,2)=0.0001_q          !SCAN_TAUREG_C
            XC%PARAM(I,3)=0.001_q           !SCAN_ALPREG
            XC%NPARAM(I)=3
         ELSEIF (FXCNAME=='RSCANZT_X') THEN
            XC%ID(I)=ID_XC_RSCANZT_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=1.998902376E-13_q !SCAN_TAUREG_X
            XC%PARAM(I,2)=0.001_q           !SCAN_ALPREG
            XC%NPARAM(I)=2
         ELSEIF (FXCNAME=='RSCANZT_C') THEN
            XC%ID(I)=ID_XC_RSCANZT_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.0001_q          !SCAN_TAUREG_C
            XC%PARAM(I,2)=0.001_q           !SCAN_ALPREG
            XC%NPARAM(I)=2
         ELSEIF (FXCNAME=='RTPSS') THEN
            XC%ID(I)=ID_XC_RTPSS
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='RTPSS_X') THEN
            XC%ID(I)=ID_XC_RTPSS_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='RTPSS_C') THEN
            XC%ID(I)=ID_XC_RTPSS_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='SCAN') THEN
            XC%ID(I)=ID_XC_SCAN
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.0_q !SCAN_TAUREG_X
            XC%PARAM(I,2)=0.0_q !SCAN_TAUREG_C
            XC%PARAM(I,3)=0.0_q !SCAN_ALPREG
            XC%NPARAM(I)=3
         ELSEIF (FXCNAME=='SCAN_X') THEN
            XC%ID(I)=ID_XC_SCAN_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.0_q !SCAN_TAUREG_X
            XC%PARAM(I,2)=0.0_q !SCAN_ALPREG
            XC%NPARAM(I)=2
         ELSEIF (FXCNAME=='SCAN_C') THEN
            XC%ID(I)=ID_XC_SCAN_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
            XC%PARAM(I,1)=0.0_q !SCAN_TAUREG_C
            XC%PARAM(I,2)=0.0_q !SCAN_ALPREG
            XC%NPARAM(I)=2
         ELSEIF (FXCNAME=='TPSS') THEN
            XC%ID(I)=ID_XC_TPSS
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='TPSS_X') THEN
            XC%ID(I)=ID_XC_TPSS_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='TPSS_C') THEN
            XC%ID(I)=ID_XC_TPSS_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='SCANL') THEN
            XC%ID(I)=ID_XC_SCANL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LLAP=.TRUE.
            XC%LDLAP=.TRUE.
         ELSEIF (FXCNAME=='RSCANL') THEN
            XC%ID(I)=ID_XC_RSCANL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LLAP=.TRUE.
            XC%LDLAP=.TRUE.
         ELSEIF (FXCNAME=='R2SCANL') THEN
            XC%ID(I)=ID_XC_R2SCANL
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LLAP=.TRUE.
            XC%LDLAP=.TRUE.
         ELSEIF (FXCNAME=='OFR2') THEN
            XC%ID(I)=ID_XC_OFR2
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LLAP=.TRUE.
            XC%LDLAP=.TRUE.
         ELSEIF (FXCNAME=='SREGTM1') THEN
            XC%ID(I)=ID_XC_SREGTM1
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='SREGTM2') THEN
            XC%ID(I)=ID_XC_SREGTM2
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='SREGTM3') THEN
            XC%ID(I)=ID_XC_SREGTM3
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='SREGTM2L') THEN
            XC%ID(I)=ID_XC_SREGTM2L
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LLAP=.TRUE.
            XC%LDLAP=.TRUE.
         ELSEIF (FXCNAME=='TASK_X') THEN
            XC%ID(I)=ID_XC_TASK_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='CC_C') THEN
            XC%ID(I)=ID_XC_CC_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='LAK') THEN
            XC%ID(I)=ID_XC_LAK
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='LAK_X') THEN
            XC%ID(I)=ID_XC_LAK_X
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=0
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='LAK_C') THEN
            XC%ID(I)=ID_XC_LAK_C
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=1
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
         ELSEIF (FXCNAME=='PBETEST') THEN
            XC%ID(I)=ID_XC_PBETEST
            XC%FAMILY(I)=XC_FAM_MGGA
            XC%KIND(I)=2
            XC%ID_LDA_TABLE(I)=ID_XC_LDA_TABLE_NOXC
            XC%LTAU=.TRUE.
            XC%LMU=.TRUE.
# 742

         ELSE
            CALL vtutor%error("Error: One of the selected functionals does not exist.")
         ENDIF

         IF (XC%FAMILY(I)>=XC_FAM_LDA) XC%LNOXC=.FALSE.
         IF (XC%FAMILY(I)==XC_FAM_LDA) THEN
            XC%LDOLDA=.TRUE.
         ELSEIF (XC%FAMILY(I)==XC_FAM_GGA) THEN
            XC%LDOGGA=.TRUE.
         ELSEIF (XC%FAMILY(I)==XC_FAM_MGGA) THEN
            XC%LDOMETAGGA=.TRUE.
         ENDIF

      ENDDO

      IF (XC%LUSE_VDW) THEN
         IF (XC%LDOLDA.AND.(.NOT.XC%LDOGGA).AND.(.NOT.XC%LDOMETAGGA)) THEN
            CALL vtutor%error("Error: Nonlocal van der Waals functionals can not be used in combination with a LDA-type functional.")
         ENDIF
      ENDIF

      IF (XC%LUSE_VDW.AND.(XC%IVDW_NL==-1)) THEN
         IF (XC%LDOGGA.AND.(.NOT.XC%LDOMETAGGA)) THEN
            XC%IVDW_NL=1
         ELSEIF ((.NOT.XC%LDOGGA).AND.XC%LDOMETAGGA) THEN
            XC%IVDW_NL=2
         ELSE
            CALL vtutor%error("Error: The type of van der Waals kernel (IVDW_NL) has to be specified.")
         ENDIF
      ENDIF

      IF (XC%LUSE_VDW.AND.XC%LSPIN_VDW.AND.(XC%IVDW_NL/=1).AND.(XC%IVDW_NL/=3).AND.(XC%IVDW_NL/=4)) THEN
         CALL vtutor%error("Error: The spin-polarized vdW DFT is available only for the kernel types IVDW_NL=1, 3 and 4.")
      ENDIF

      RETURN
    END SUBROUTINE SET_XC_DATA

# 882


!******************* SUBROUTINE SETPARAMS_XC ********************
!
!> This subroutine sets the parameters of the functional.
!
!****************************************************************

    SUBROUTINE SETPARAMS_XC(XC)
      USE base, ONLY : in_struct
      USE reader_tags
      USE incar_reader, ONLY: COUNT_ELEMENTS
      USE constant, ONLY : AUTOA
      USE fock_glb, ONLY : AEXX,HFSCREEN
      USE tutor, ONLY : vtutor

      IMPLICIT NONE
      TYPE (in_struct) IO
      TYPE (xc_info) XC
      INTEGER :: I,J,N,NXC,N_ALDAX,N_AGGAX,N_AMGGAX,N_ALDAC,N_AGGAC,N_AMGGAC,IERR
      REAL(q) PARAMTMP
      LOGICAL LOPEN
      CHARACTER (LEN=20) :: ICHAR,JCHAR,CHAR1

      CALL OPEN_INCAR_IF_FOUND(IO%IU5,LOPEN)

      NXC=XC%NXC

      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'XC_C', XC%COEFF, NXC, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
      IF ((N>=1).AND.(SXCTAG/='XC')) THEN
         CALL vtutor%error("SETPARAMS_XC: The XC_C tag can only be used when the functional is specified with the XC tag.")
      ENDIF
      IF ((N/=0).AND.(N/=NXC)) CALL vtutor%error("Error: The number of coefficients XC_C provided is not equal to the number of functional components chosen.")

      N_ALDAX =COUNT_ELEMENTS(INCAR_F,"ALDAX")
      N_AGGAX =COUNT_ELEMENTS(INCAR_F,"AGGAX")
      N_AMGGAX=COUNT_ELEMENTS(INCAR_F,"AMGGAX")
      N_ALDAC =COUNT_ELEMENTS(INCAR_F,"ALDAC")
      N_AGGAC =COUNT_ELEMENTS(INCAR_F,"AGGAC")
      N_AMGGAC=COUNT_ELEMENTS(INCAR_F,"AMGGAC")

      DO I=1,NXC

# 972


         WRITE(ICHAR,*) I
         DO J=1,XC%NPARAM(I)
            WRITE(JCHAR,*) J
# 989

            CHAR1 = 'XC'//TRIM(ADJUSTL(ICHAR))//'_P'//TRIM(ADJUSTL(JCHAR))
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, TRIM(CHAR1), PARAMTMP, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
            IF (N>=1) THEN
               IF (SXCTAG=='XC') THEN
                   XC%PARAM(I,J)=PARAMTMP
               ELSE
                  CALL vtutor%error("SETPARAMS_XC: The XCm_Pn tags can only be used when the functional is specified with the XC tag.")
               ENDIF
            ENDIF
         ENDDO
# 1005

      ENDDO

      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'PARAM1', XC%PARAM1, IERR, WRITEXMLINCAR)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'PARAM2', XC%PARAM2, IERR, WRITEXMLINCAR)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'ZAB_VDW', XC%Zab_VDW, IERR, WRITEXMLINCAR)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'GAMMA_VDW', XC%GAMMA_VDW, IERR, WRITEXMLINCAR)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'ALPHA_VDW', XC%ALPHA_VDW, IERR, WRITEXMLINCAR)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'BPARAM', XC%BPARAM, IERR, WRITEXMLINCAR)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'CPARAM', XC%CPARAM, IERR, WRITEXMLINCAR)

      CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

      RETURN
    END SUBROUTINE SETPARAMS_XC

!******************* SUBROUTINE SETUP_LDA_XC ************************
!
!>  VASP interpolates the XC-energy density from a table (at least
!>  the plane wave part). The required table is generated here.
!
!********************************************************************

    SUBROUTINE SETUP_LDA_XC(XC,ISPIN,IU6,IU0,IDIOT)

!! USE moffload_struct_def

      USE tutor, ONLY: vtutor
      USE ini
      IMPLICIT NONE

      INTEGER  ISPIN            ! spin

! arrays for tutor call
      INTEGER IU6,IU0,IDIOT
      TYPE (xc_info) XC
! temporary
      INTEGER N,NDUMMY
      REAL(q) AMARG
      CHARACTER (1) CSEL
      CHARACTER (2) CEXCH
      IF (XC%NXC==0) THEN
         CALL vtutor%bug("internal ERROR in SETUP_LDA_XC: No functional has been set up", "setex.F", 1047)
      ENDIF
!
! set the exchange correlation type for the internal table
! before the table is used XCTABLE_CHECK should be called
! to ascertain that XC_DATA_TABLE is correct and equivalent to XC
!

      CALL SET_EX_TABLE(XC,NEXCH,EXCTAB%EXCTAB,EXCTAB%NEXCHF,EXCTAB%RHOEXC,IU0)

      AMARG=1E30_q  ! natural boundary conditions required
      CALL SPLCOF(EXCTAB%EXCTAB(1,1,1),EXCTAB%NEXCHF(2),NEXCH,AMARG)

      IF (ISPIN==2) THEN
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,2),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,3),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,4),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,5),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,6),EXCTAB%NEXCHF(2),NEXCH,AMARG)
      ENDIF
      XC_DATA_TABLE=XC

      IF (IU6>=0) WRITE(IU6,7004)TRIM(ADJUSTL(SZXCNAM)),EXCTAB%RHOEXC(1),EXCTAB%NEXCHF(1), &
     &                     EXCTAB%RHOEXC(2),EXCTAB%NEXCHF(2)

 7004 FORMAT(' exchange-correlation table for ',A/ &
     &       '   RHO(1)= ',F8.3,5X,'  N(1)  = ',I8/ &
     &       '   RHO(2)= ',F8.3,5X,'  N(2)  = ',I8)

!$ACC UPDATE DEVICE(EXCTAB) ASYNC(ACC_ASYNC_Q)
!$ACC UPDATE DEVICE(EXCTAB%RHOEXC,EXCTAB%NEXCHF,EXCTAB%EXCTAB) ASYNC(ACC_ASYNC_Q)
      RETURN
    END SUBROUTINE SETUP_LDA_XC

!******************* SUBROUTINE SET_EX_TABLE *********************
!
!> Set up the default xc-table i.e. Ceperly Alder with standard
!> interpolation to spin with relativistic correction.
!
!*****************************************************************

    SUBROUTINE SET_EX_TABLE(XC,N,EXCTAB,NEXCHF,RHOEXC,IU0)

      USE ldalib
      USE string, ONLY: str
      USE xc_lda_table_name
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      INTEGER N, IU0
      REAL(q) :: EXCTAB(N,5,6),RHOEXC(2)
      INTEGER :: NEXCHF(2)
      TYPE (xc_info) XC
! local
      LOGICAL TREL
      CHARACTER (500) CEXCH,CEXCHTMP
      INTEGER J,I,IXC
      REAL(q) :: SLATER, RHOSMA, RHOMAX, RH, EXCP, DEXF, DECF, ALPHA, DRHO, ZETA, FZA, FZB

      TREL=.TRUE.

      CEXCH='('

      NXCLOOP1: DO IXC=1,XC%NXC

      IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_X) THEN
         CEXCHTMP='Slater'
      ELSEIF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_SLATER) THEN
         CEXCHTMP='Slater(with rela. corr.)'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PZ) THEN
         CEXCHTMP='Slater(with rela. corr.)+CA(PZ)'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PZ_C) THEN
         CEXCHTMP='CA(PZ)'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_VWN5) THEN
         CEXCHTMP='Slater(with rela. corr.)+VWN5'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_HL) THEN
         CEXCHTMP='Slater(with rela. corr.)+HL'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_WI) THEN
         CEXCHTMP='Slater(with rela. corr.)+Wigner'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92) THEN
         CEXCHTMP='Slater+PW92'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92_C) THEN
         CEXCHTMP='PW92 correlation'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92_RELA) THEN
         CEXCHTMP='Slater(with rela. corr.)+PW92'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_B3) THEN
         CEXCHTMP='Slater(with rela. corr.)+VWN3'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_B5) THEN
         CEXCHTMP='Slater(with rela. corr.)+VWN5'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_RA) THEN
         CEXCHTMP='RPA of Perdew-Wang'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_03) THEN
         CEXCHTMP='sr-RPA(mu=0.3 A^-1))'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_05) THEN
         CEXCHTMP='sr-RPA(mu=0.5 A^-1)'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_10) THEN
         CEXCHTMP='sr-RPA(mu=1.0 A^-1)'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_20) THEN
         CEXCHTMP='sr-RPA(mu=2.0 A^-1)'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_PL) THEN
         CEXCHTMP='RPA+ of Perdew-Wang'
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_NOXC) THEN
         CEXCHTMP='None'
      ELSE
         CALL vtutor%bug("internal error in SET_EX_TABLE: Wrong exchange-correlation type.","setex.F", 1152)
      ENDIF

      IF (IXC==1) THEN
         CEXCH=TRIM(ADJUSTL(CEXCH))//TRIM(ADJUSTL(CEXCHTMP))
      ELSE IF (IXC/=1) THEN
         CEXCH=TRIM(ADJUSTL(CEXCH))//','//TRIM(ADJUSTL(CEXCHTMP))
      ENDIF

      ENDDO NXCLOOP1

      CEXCH=TRIM(ADJUSTL(CEXCH))//')'

      IF (IU0>=0) THEN
         IF (XC%LFCI==1) THEN 
            WRITE(IU0,*)'LDA part: xc-table for ',TRIM(ADJUSTL(CEXCH)), ', Vosko type interpolation para-ferro'
         ELSE
            WRITE(IU0,*)'LDA part: xc-table for ',TRIM(ADJUSTL(CEXCH)), ', standard interpolation'
         ENDIF
      ENDIF

! Slater parameter
      SLATER=1._q
! standard interpolation   from para- to ferromagnetic corr
      RHOSMA=INT(0.5_q*1000)/1000._q
!#define hugeXCtable
# 1181

      RHOMAX=INT(100.5_q*1000)/1000._q



      RHOEXC(1)=RHOSMA
      RHOEXC(2)=RHOMAX
      NEXCHF(1)=NSMA
      NEXCHF(2)=N

      IF (NSMA/=0) THEN
         RH=RHOSMA/NSMA/2
      ELSE
         RH=RHOMAX/N/100
      ENDIF
      J=1

      CALL EXCHG(XC,RH,EXCP,DEXF,DECF,ALPHA,SLATER,TREL)
        EXCTAB(J,1,1)=RH
        EXCTAB(J,1,2)=RH
        EXCTAB(J,1,3)=RH
        EXCTAB(J,1,4)=RH

        EXCTAB(J,2,1)=EXCP
        EXCTAB(J,2,2)=DEXF
        EXCTAB(J,2,3)=DECF
        EXCTAB(J,2,4)=ALPHA


      IF (NSMA/=0) THEN
         DRHO=RHOSMA/NSMA
         DO 100 I=1,NSMA-1
            J=I+1
            RH=DRHO*I
            CALL EXCHG(XC,RH,EXCP,DEXF,DECF,ALPHA,SLATER,TREL)
          EXCTAB(J,1,1)=RH
          EXCTAB(J,1,2)=RH
          EXCTAB(J,1,3)=RH
          EXCTAB(J,1,4)=RH

          EXCTAB(J,2,1)=EXCP
          EXCTAB(J,2,2)=DEXF
          EXCTAB(J,2,3)=DECF
          EXCTAB(J,2,4)=ALPHA

  100    CONTINUE

         J=NSMA+1
         RH=RHOSMA
         CALL EXCHG(XC,RH,EXCP,DEXF,DECF,ALPHA,SLATER,TREL)
         EXCTAB(J,1,1)=RH
         EXCTAB(J,1,2)=RH
         EXCTAB(J,1,3)=RH
         EXCTAB(J,1,4)=RH

         EXCTAB(J,2,1)=EXCP
         EXCTAB(J,2,2)=DEXF
         EXCTAB(J,2,3)=DECF
         EXCTAB(J,2,4)=ALPHA
      ENDIF

      DRHO=(RHOMAX-RHOSMA)/(N-NSMA)
      DO 200 I=1,N-NSMA-1
         J=I+NSMA+1
         RH=DRHO*I+RHOSMA
         CALL EXCHG(XC,RH,EXCP,DEXF,DECF,ALPHA,SLATER,TREL)
         EXCTAB(J,1,1)=RH
         EXCTAB(J,1,2)=RH
         EXCTAB(J,1,3)=RH
         EXCTAB(J,1,4)=RH

         EXCTAB(J,2,1)=EXCP
         EXCTAB(J,2,2)=DEXF
         EXCTAB(J,2,3)=DECF
         EXCTAB(J,2,4)=ALPHA
  200 CONTINUE

      J=1
      ZETA=0
      FZA =.854960467080682810_q
      FZB =.854960467080682810_q
      EXCTAB(J,1,5)=ZETA
      EXCTAB(J,1,6)=ZETA
      EXCTAB(J,2,5)=FZA
      EXCTAB(J,2,6)=FZB

      DO 800 I=1,N-1
         ZETA=FLOAT(I)/FLOAT(N-1)
         FZA=FZ0(ZETA)/ZETA/ZETA
         FZB=FZ0(ZETA)/ZETA/ZETA
         J=I+1
         EXCTAB(J,1,5)=ZETA
         EXCTAB(J,1,6)=ZETA
         EXCTAB(J,2,5)=FZA
         EXCTAB(J,2,6)=FZB
  800 CONTINUE

      IF (.FALSE.) THEN
      DO 300 I=1,N
         WRITE(97,20) EXCTAB(I,1,1),EXCTAB(I,2,1)
  300 CONTINUE
      DO 400 I=1,N
         WRITE(97,20)  EXCTAB(I,1,2),EXCTAB(I,2,2)
  400 CONTINUE
      DO 500 I=1,N
         WRITE(97,20)  EXCTAB(I,1,3),EXCTAB(I,2,3)
  500 CONTINUE
      DO 600 I=1,N
         WRITE(97,20)  EXCTAB(I,1,4),EXCTAB(I,2,4)
  600 CONTINUE
      DO 700 I=1,N
         WRITE(97,20)  EXCTAB(I,1,5),EXCTAB(I,2,5)
  700 CONTINUE
      DO 710 I=1,N
         WRITE(97,20)  EXCTAB(I,1,6),EXCTAB(I,2,6)
  710 CONTINUE
   20 FORMAT((3(E24.16,2X)))
      ENDIF

      RETURN
   END SUBROUTINE SET_EX_TABLE

!******************* SUBROUTINE EXCHG ****************************************
!
!> EXCHG calculates the LDA part of xc energy density per particle eps_xc(rho)
!>
!> ~~~
!> E_xc = \int eps_xc(rho(r)) rho(r) d^3 r
!> ~~~
!>
!> The subroutines calls various subroutines defined in xclib
!>
!> Before returning, the subroutine divides the results by the density^(1/3)
!> which is the behaviour of the exchange density per particle e_x
!>
!> this routine expects input in Angst, and the output is in eV/Angst
!>
!> Most of the called routines expect atomic units, hence
!> EXCHG calculates the Wigner Seitz radius in a.u. and then calls
!> the various exchange correlation routines (which are expected
!> to return the exchange correlation energy as well as potential
!> in Rydberg)
!
!*****************************************************************************

    SUBROUTINE EXCHG(XC,RHO,EXCP,DEXF,DECF,ALPH,SLATER,TREL)
      USE constant
      USE ldalib
      USE xc_lda_table_name
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      LOGICAL TREL
      REAL(q) :: RHO,EXCP,DEXF,DECF,ALPH,SLATER
      REAL(q) :: RHOP, RHOM, EXP, EXM, DELTA
      TYPE (xc_info) XC
! local
      INTEGER :: IXC
      REAL(q) :: RS, RH, ZETA, FZA, FZB, &
           RHOTHD, SE, ECLDA, ECD1LDA, ECD2LDA, ECLDA_MAG, A0
      REAL(q) :: EXCP_TMP, DEXF_TMP, DECF_TMP, ALPH_TMP, ALDAC

      IF (RHO==0) THEN
         EXCP=0._q
         DEXF=0._q
         DECF=0._q
         ALPH=0._q
         RETURN
      ENDIF

      ALDAC=XC%ALDAC

# 1379


      RHOTHD = RHO**(1/3._q)
      RS = (3._q/(4._q*PI)/RHO)**(1/3._q) /AUTOA

      EXCP=0.0_q
      DEXF=0.0_q
      DECF=0.0_q
      ALPH=0.0_q

      NXCLOOP2: DO IXC=1,XC%NXC

      EXCP_TMP=0.0_q
      DEXF_TMP=0.0_q
      DECF_TMP=0.0_q
      ALPH_TMP=0.0_q

      IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_X) THEN
         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=0._q
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_SLATER) THEN
         EXCP_TMP=SLATER*EX_MOD(XC,RS,1,TREL)
         DEXF_TMP=SLATER*EX_MOD(XC,RS,2,TREL)-EXCP_TMP
         DECF_TMP=0._q
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PZ) THEN
         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECCA(RS,1)*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
         DECF_TMP=ECCA(RS,2)*ALDAC-ECCA(RS,1)*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PZ_C) THEN
         EXCP_TMP=ECCA(RS,1)*ALDAC
         DEXF_TMP=0._q
         DECF_TMP=ECCA(RS,2)*ALDAC-ECCA(RS,1)*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_VWN5) THEN
! Iann Gerber 18/11/04
!        EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECVO(RS,1)
         ZETA=0.0_q
         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+EC_MOD(XC,RS,ZETA,1,TREL)*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
!        DECF_TMP=ECVO(RS,2)*ALDAC-ECVO(RS,1)*ALDAC
! hack: ZETA set to 1.0_q or 0.0_q, respectively
         DECF_TMP=EC_MOD(XC,RS,1._q,2,TREL)*ALDAC-EC_MOD(XC,RS,0._q,1,TREL)*ALDAC
         ALPH_TMP=0._q
! Iann Gerber 18/11/04; updated to scrLSDA by Joachim Paier (28/04/08)
!      ELSE IF (LEXCH_LDA==3) THEN
!         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECGL(RS,1)*ALDAC
!         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
!         DECF_TMP=ECGL(RS,2)*ALDAC-ECGL(RS,1)*ALDAC
!         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_HL) THEN
         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECHL(RS,1)*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
         DECF_TMP=ECHL(RS,2)*ALDAC-ECHL(RS,1)*ALDAC
         ALPH_TMP=0._q
!      ELSE IF (LEXCH_LDA==5) THEN
!         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECBH(RS,1)*ALDAC
!         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
!         DECF_TMP=ECBH(RS,2)*ALDAC-ECBH(RS,1)*ALDAC
!         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_WI) THEN
         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECWI(RS,1)*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
         DECF_TMP=ECWI(RS,2)*ALDAC-ECWI(RS,1)*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)

         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92_C) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)

         EXCP_TMP=ECLDA*ALDAC
         DEXF_TMP=0._q
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92_RELA) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)

         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_B3) THEN
! Joachim Paier
         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECVOIII(RS,1)*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
         DECF_TMP=ECVOIII(RS,2)*ALDAC-ECVOIII(RS,1)*ALDAC
         ALPH_TMP=0._q
! check is VO similar to CA
!        WRITE(77,'(5F14.7)') RS, ECCA(RS,2),ECVO(RS,2)
!        WRITE(78,'(5F14.7)') RS, ECCA(RS,1),ECVO(RS,1)
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_B5) THEN
! Joachim Paier
         EXCP_TMP=EX_MOD(XC,RS,1,TREL)+ECVO(RS,1)*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,TREL)-EX_MOD(XC,RS,1,TREL)
         DECF_TMP=ECVO(RS,2)*ALDAC-ECVO(RS,1)*ALDAC
         ALPH_TMP=0._q
! check is VO similar to CA
!        WRITE(77,'(5F14.7)') RS, ECCA(RS,2),ECVO(RS,2)
!        WRITE(78,'(5F14.7)') RS, ECCA(RS,1),ECVO(RS,1)
! jH-new RPA with Perdew Wang para
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_RA) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_RPA(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA_RPA(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
! Judith Harl
! functionals for range-separated ACFDT (LDA - short range RPA):
! a bit akward at the moment since the range separation parameter
! is hard coded for now. Hopefully this will change ...
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_03) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_15au(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_0_15au(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_05) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_5A(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_0_5A(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_10) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_1_0A(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_1_0A(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_20) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_2_0A(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_2_0A(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q         
! jH-new RPA-plus = ELDAQMC-ELDARPA, exchange term is (0._q,0._q)
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_PL) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP_TMP=EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF_TMP=EX_MOD(XC,RS,2,.FALSE.)-EX_MOD(XC,RS,1,.FALSE.)
         DECF_TMP=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH_TMP=0._q
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_NOXC) THEN
         EXCP_TMP=0._q
         DEXF_TMP=0._q
         DECF_TMP=0._q
         ALPH_TMP=0._q
      ELSE
         CALL vtutor%bug("internal error in EXCHG: No LDA table has been selected.", "setex.F", 1565)
      ENDIF
!
! for Perdew, Burke Ernzerhof we always use the
! recommended  interpolation from nm to magnetic
!
      IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_B3) THEN
         A0=ALPHA0_III(RS)*ALDAC
         ALPH_TMP=DECF_TMP-A0
         DECF_TMP=A0
      ELSE IF ((XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92) .OR. &
     &         (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92_C).OR. &
               (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_PW92_RELA)) THEN
         A0=PBE_ALPHA(RS)*ALDAC
         ALPH_TMP=DECF_TMP-A0
         DECF_TMP=A0
! jH-new new switching for RPA
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_RA) THEN
         A0=RPA_ALPHA(RS)*ALDAC
         ALPH_TMP=DECF_TMP-A0
         DECF_TMP=A0
! jH-new new switching for RPAplus
      ELSE IF (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_ACFDT_PL) THEN          
         A0=(PBE_ALPHA(RS)-RPA_ALPHA(RS))*ALDAC
         ALPH_TMP=DECF_TMP-A0
         DECF_TMP=A0
      ELSE IF ((XC%LFCI==1) .OR. (XC%ID_LDA_TABLE(IXC)==ID_XC_LDA_TABLE_B5)) THEN
         A0=ALPHA0(RS)*ALDAC
         ALPH_TMP=DECF_TMP-A0
         DECF_TMP=A0
      ENDIF

      EXCP=EXCP+XC%COEFF(IXC)*EXCP_TMP
      DEXF=DEXF+XC%COEFF(IXC)*DEXF_TMP
      DECF=DECF+XC%COEFF(IXC)*DECF_TMP
      ALPH=ALPH+XC%COEFF(IXC)*ALPH_TMP

      ENDDO NXCLOOP2

      EXCP=EXCP*RYTOEV/RHOTHD
      DEXF=DEXF*RYTOEV/RHOTHD
      DECF=DECF*RYTOEV/RHOTHD
      ALPH=ALPH*RYTOEV/RHOTHD

      RETURN
    END SUBROUTINE EXCHG

!******************* SUBROUTINE PUSH_XC_TYPE *****************************************
!
!> Subroutine to temporarily use an alternative local exchange-correlation functional.
!
!*************************************************************************************

    SUBROUTINE PUSH_XC_TYPE(XC,XC_DATA_TEMP)

!! USE moffload

      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE fock_glb, ONLY: HFSCREEN, AEXX
      USE wpbe, ONLY : INIT_WPBE
      IMPLICIT NONE
      TYPE (xc_info) :: XC,XC_DATA_TEMP

# 1631


      ISTACK=ISTACK+1
      IF (ISTACK>SIZE(XC_DATA_STACK)) THEN
         CALL vtutor%bug("PUSH_XC_TYPE: stack exhausted " // str(ISTACK), "setex.F", 1635)
      ENDIF

! gK: a comment is in place here
! if LDASCREEN is set to 0, xclib.F falls back to the non-range seperated functional, regardless of LUSE_MODEL_HF
! if this were not the case, (1._q,0._q) would need to set LUSE_MODEL_HF to .FALSE.
! to recover the original PBE functional (see PUSH in paw.F)
      XC_DATA_STACK(ISTACK)=XC
      XC=XC_DATA_TEMP
!gK 18.05.2021 also set some variables in fock_glb; this is way cleaner then the tinkering 1._q before
! this was previously manually set in some places, but in a fairly inconsistent manner
! is was working most of the time, since the calling routine avoided to evaluate the Fock exchange
      AEXX=XC%AEXX
      HFSCREEN=XC%LDASCREEN
      CALL INIT_WPBE(XC)

# 1653

      RETURN
    END SUBROUTINE PUSH_XC_TYPE

!******************* SUBROUTINE POP_XC_TYPE **************
!
!> Subroutine to revert what was 1._q by PUSH_XC_TYPE.
!
!*********************************************************

    SUBROUTINE POP_XC_TYPE(XC)

!! USE moffload

      USE tutor, ONLY: vtutor
      USE fock_glb, ONLY: HFSCREEN, AEXX
      USE wpbe, ONLY : INIT_WPBE
      IMPLICIT NONE
      TYPE (xc_info) :: XC

# 1675


      IF (ISTACK==0) THEN
         CALL vtutor%bug("POP_XC_TYPE: stack is empty, probably PUSH_XC_TYPE was not called.", "setex.F", 1678)
      ENDIF

      XC=XC_DATA_STACK(ISTACK)
!gK 18.05.2021 also set some variables in fock_glb; this is way cleaner then the tinkering 1._q before
      AEXX=XC%AEXX
      HFSCREEN=XC%LDASCREEN
      ISTACK=ISTACK-1
      CALL INIT_WPBE(XC)

# 1690

      RETURN
    END SUBROUTINE POP_XC_TYPE

!******************* SUBROUTINE PUSH_XC_TYPE_METAGGA ****************
!
!> Subroutine to temporarily replace the meta-GGA component(s) of the
!> functional by the functional specified by LEXCH in POTCAR.
!
!********************************************************************

    SUBROUTINE PUSH_XC_TYPE_METAGGA(XC)

!! USE moffload

      USE xc_name
      USE xc_family_name
      USE tutor, ONLY : vtutor
      USE string, ONLY : str

      IMPLICIT NONE
      INTEGER I
      TYPE (xc_info) :: XC

# 1716


      I_XC_DATA_STACK_2=I_XC_DATA_STACK_2+1
      IF (I_XC_DATA_STACK_2>SIZE(XC_DATA_STACK_2)) THEN
         CALL vtutor%bug("PUSH_XC_TYPE_METAGGA: stack exhausted " // str(I_XC_DATA_STACK_2), "setex.F", 1720)
      ENDIF

      XC_DATA_STACK_2(I_XC_DATA_STACK_2)=XC

      DO I=1,XC_DATA_STACK_2(I_XC_DATA_STACK_2)%NXC
         IF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%FAMILY(I)==XC_FAM_MGGA) THEN
            IF (XC_DATA_PP%ID(1)==ID_XC_PZ) THEN
               IF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==0) THEN
                  XC%ID(I)=ID_XC_SLATER
               ELSEIF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==1) THEN
                  XC%ID(I)=ID_XC_PZ_C
               ELSEIF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==2) THEN
                  XC%ID(I)=ID_XC_PZ
               ENDIF
               XC%FAMILY(I)=XC_FAM_LDA
               XC%NPARAM(I)=0
               XC%PARAM(I,:)=0._q
               XC%NAME(I)='PZ'
               XC%LDOLDA=.TRUE.
            ELSEIF (XC_DATA_PP%ID(1)==ID_XC_PW91) THEN
               IF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==0) THEN
!Exchange-only PW91 is not yet implemented
                  CALL vtutor%error("PUSH_XC_TYPE_METAGGA: A POTCAR generated from PW91 can not be used when a exchange-only meta-GGA is specified in INCAR.")
!                  XC%ID(I)=ID_XC_PW91_X
               ELSEIF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==1) THEN
!Correlation-only PW91 is not yet implemented
                  CALL vtutor%error("PUSH_XC_TYPE_METAGGA: A POTCAR generated from PW91 can not be used when a correlation-only meta-GGA is specified in INCAR.")
!                  XC%ID(I)=ID_XC_PW91_C
               ELSEIF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==2) THEN
                  XC%ID(I)=ID_XC_PW91
               ENDIF
               XC%FAMILY(I)=XC_FAM_GGA
               XC%NPARAM(I)=0
               XC%PARAM(I,:)=0._q
               XC%NAME(I)='91'
               XC%LDOGGA=.TRUE.
            ELSEIF (XC_DATA_PP%ID(1)==ID_XC_PBE) THEN
               IF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==0) THEN
                  XC%ID(I)=ID_XC_PBE_X
               ELSEIF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==1) THEN
                  XC%ID(I)=ID_XC_PBE_C
               ELSEIF (XC_DATA_STACK_2(I_XC_DATA_STACK_2)%KIND(I)==2) THEN
                  XC%ID(I)=ID_XC_PBE
               ENDIF
               XC%FAMILY(I)=XC_FAM_GGA
               XC%NPARAM(I)=0
               XC%PARAM(I,:)=0._q
               XC%NAME(I)='PE'
               XC%LDOGGA=.TRUE.
            ENDIF
         ENDIF
      ENDDO
      XC%LDOMETAGGA=.FALSE.

# 1777

      RETURN
    END SUBROUTINE PUSH_XC_TYPE_METAGGA

!******************* SUBROUTINE POP_XC_TYPE_METAGGA ************
!
!> Subroutine to revert what was 1._q by PUSH_XC_TYPE_METAGGA.
!
!***************************************************************

    SUBROUTINE POP_XC_TYPE_METAGGA(XC)

!! USE moffload

      USE tutor, ONLY : vtutor

      IMPLICIT NONE
      TYPE (xc_info) :: XC

# 1798


      IF (I_XC_DATA_STACK_2==0) THEN
         CALL vtutor%bug("POP_XC_TYPE_METAGGA: stack is empty, probably PUSH_XC_TYPE_METAGGA was not called.", "setex.F", 1801)
      ENDIF

      XC=XC_DATA_STACK_2(I_XC_DATA_STACK_2)
      I_XC_DATA_STACK_2=I_XC_DATA_STACK_2-1

# 1809

      RETURN
    END SUBROUTINE POP_XC_TYPE_METAGGA

  END MODULE setexm
