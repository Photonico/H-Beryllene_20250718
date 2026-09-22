# 1 "setex_struct.F"
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


# 2 "setex_struct.F" 2 
      MODULE setexm_struct_def
      USE prec
# 6

      IMPLICIT NONE
!
! the table has two section (1._q,0._q) for low densities
!   (1...NEXCHF(1),0...RHOEXC(1))
! and (1._q,0._q) for high densities
!   (...NEXCHF(2),...RHOEXC(2))
!
      INTEGER, PARAMETER :: NEXCH=4000

      TYPE exctable
        REAL(q) :: EXCTAB(NEXCH,5,6) ! table including spline coeff.
        REAL(q) :: RHOEXC(2)         ! maximal densities
        INTEGER :: NEXCHF(2)         ! number of points
      END TYPE
!
! NSMA is the size of the table EXCTAB
!
      INTEGER, PARAMETER :: NSMA=2000
!
! SXCTAG is the tag (GGA, METAGGA or XC) used in INCAR to specify the functional
! SZXCNAM is the list of the functional components specified in INCAR with the GGA, METAGGA or XC tag
! SZXC_C is the string of the coefficients (XC_C) multiplying the functional components
!
      CHARACTER(LEN=7), SAVE :: SXCTAG
      CHARACTER(LEN=400), SAVE :: SZXCNAM
      CHARACTER(LEN=400), SAVE :: SZXC_C
!
! stores the exchange correlation interpolation table
!
      TYPE (exctable), SAVE :: EXCTAB
!
! for meta-GGA functionals
!
      INTEGER, SAVE :: LMAXTAU=6
      LOGICAL, SAVE :: LMIXTAU=.FALSE.
!
! The maximum number of parameter of a functional
!
      INTEGER, SAVE :: NPARAMMAX=200
!
! xc related variables
!
      TYPE xc_info
         INTEGER :: NXC                            !< Number of components of the functional
         INTEGER :: IVDW_NL                        !< ID number of the nonlocal vdW functional
         INTEGER :: LFCI                           !< Interpolation of correlation from paramagnetic to
!< ferromagnetic case according to
!< Vosko, Wilk and Nusair, CAN. J. PHYS. 58, 1200 (1980).
         INTEGER, ALLOCATABLE :: ID(:)             !< ID number of the components of the functional
         INTEGER, ALLOCATABLE :: FAMILY(:)         !< Family (LDA, GGA or MGGA) of the components of the functional
         INTEGER, ALLOCATABLE :: KIND(:)           !< Kind (exchange, correlation or exchange-correlation) of the components of the functional
         INTEGER, ALLOCATABLE :: ID_LDA_TABLE(:)   !< ID number of the LDA-table functional for each component of the functional
         INTEGER, ALLOCATABLE :: NPARAM(:)         !< Number of parameters for each component of the functional
         REAL(q), ALLOCATABLE :: COEFF(:)          !< Multiplying coefficient for each component of the functional
         REAL(q), ALLOCATABLE :: PARAM(:,:)        !< Parameters for each component of the functional
         REAL(q) :: PARAM1                         !< mu in the exchange GGA optB86b/B86R (GGA=MK) or
!< beta in the exchange GGA optB88 (GGA=BO)
         REAL(q) :: PARAM2                         !< kappa in the exchange GGA optB86b/B86R (GGA=MK) or
!< mu in the exchange GGA optB88 (GGA=BO).
         REAL(q) :: Zab_VDW                        !< Z_ab in the kernel of the nonlocal vdW-DF functional of Dion et al.
         REAL(q) :: GAMMA_VDW                      !< gamma in the kernel of the nonlocal vdW-DF3 functional
         REAL(q) :: ALPHA_VDW                      !< alpha in the kernel of the nonlocal vdW-DF3 functional
         REAL(q) :: BPARAM                         !< b in the kernel of the nonlocal rVV10 functional
         REAL(q) :: CPARAM                         !< C in the kernel of the nonlocal rVV10 functional
         REAL(q) :: ALDAX                          !< Fraction of LDA exchange in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: ALDAC                          !< Fraction of LDA correlation in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AGGAX                          !< Fraction of gradient corrections to the exchange in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AGGAC                          !< Fraction of gradient corrections to the correlation in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AMGGAX                         !< Fraction of meta-GGA exchange in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AMGGAC                         !< Fraction of meta-GGA correlation in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AEXX                           !< Fraction of Hartree-Fock exchange (WARINING: to be used only for the PUSH_XC_TYPE and POP_XC_TYPE subroutines)
         REAL(q) :: LDASCREEN                      !< Screened LDA exchange parameter
         REAL(q) :: LDASCREENC                     !< Screened LDA correlation parameter
         LOGICAL :: LSFBXC                         !< .TRUE. if the sources and drains are removed from the XC B field
         LOGICAL :: LNOXC                          !< .TRUE. if no functional is chosen
         LOGICAL :: LDOLDA                         !< .TRUE. if at least (1._q,0._q) component of the functional is a LDA
         LOGICAL :: LDOGGA                         !< .TRUE. if at least (1._q,0._q) component of the functional is a GGA
         LOGICAL :: LDOMETAGGA                     !< .TRUE. if at least (1._q,0._q) component of the functional is a MGGA
         LOGICAL :: LUSE_LIBXC                     !< .TRUE. if at least (1._q,0._q) component of the functional is from Libxc
         LOGICAL :: LTAU                           !< .TRUE. if at least (1._q,0._q) component of the functional depends on tau
         LOGICAL :: LLAP                           !< .TRUE. if the Laplacian has to be computed (e.g., for MBJ and LL-MGGA)
         LOGICAL :: LMU                            !< .TRUE. if the nonlocal component de_xc/dtau of a MGGA has to be computed
         LOGICAL :: LDLAP                          !< .TRUE. if the component LAP(de_xc/dLAP(rho)) of a LL-MGGA has to be computed
         LOGICAL :: LUSE_VDW                       !< .TRUE. is a nonlocal vdW functional is used
         LOGICAL :: LSPIN_VDW                      !< .TRUE. is the spin-polarized formulation of the vdW functional of Dion et al. is used
         LOGICAL :: LUSE_LONGRANGE_HF              !< Short range LDA exchange interaction only
!< long range contribution is 1._q in HF
!< the default is that HF treats short range, and LDA long range
!< but using this flag the behavior can inverted
!< the flag should be identical to LRHFCALC in fock.F
         LOGICAL :: LRANGE_SEPARATED_CORR          !< .TRUE. for range-separated LDA correlation
!< .FALSE. for complete LDA correlation (default)
         LOGICAL :: LUSE_THOMAS_FERMI              !< Thomas-Fermi screening in local exchange.
!< Should be identical to L_THOMAS_FERMI in fock.F.
         LOGICAL :: LUSE_MODEL_HF                  !< No local exchange in the short range
!< a fraction LDAX of the local exchange in the long range limit
!< should be identical to L_MODEL_HF
         LOGICAL :: LNO_AUG_XC                     !< Compute XC from soft charge density
         CHARACTER(LEN=40), ALLOCATABLE :: NAME(:) !< Name of the components of the functionals
      END TYPE xc_info
      TYPE (xc_info), SAVE :: XC_DATA, &
     &                        XC_DATA_PP, &
     &                        XC_DATA_NOXC, &
     &                        XC_DATA_STACK_2(5), &
     &                        XC_DATA_TABLE, &
     &                        XC_DATA_WPBE_TABLE
!
! libxc related variables
!
      CHARACTER(LEN=35) :: LIBXC_VERSION='     F    Libxc'
# 141


!$ACC DECLARE CREATE(EXCTAB, &
!$ACC&               LMAXTAU, &
!$ACC&               LMIXTAU, &
!$ACC&               EXCTAB%RHOEXC, &
!$ACC&               EXCTAB%NEXCHF, &
!$ACC&               EXCTAB%EXCTAB)

# 160


    CONTAINS

      SUBROUTINE ALLOCATE_XC_DATA(XC_DATA_TEMP,NXC)
      IMPLICIT NONE
      INTEGER :: NXC
      TYPE (xc_info) :: XC_DATA_TEMP

      ALLOCATE(XC_DATA_TEMP%ID(NXC), &
     &         XC_DATA_TEMP%FAMILY(NXC), &
     &         XC_DATA_TEMP%KIND(NXC), &
     &         XC_DATA_TEMP%ID_LDA_TABLE(NXC), &
     &         XC_DATA_TEMP%NPARAM(NXC), &
     &         XC_DATA_TEMP%COEFF(NXC), &
     &         XC_DATA_TEMP%PARAM(NXC,NPARAMMAX), &
     &         XC_DATA_TEMP%NAME(NXC))

      END SUBROUTINE ALLOCATE_XC_DATA

      END MODULE setexm_struct_def

      MODULE xc_name
      ENUM, BIND(C)
         ENUMERATOR :: ID_XC_NOXC,ID_XC_COHSM,ID_XC_SLATER,ID_XC_HL,ID_XC_PZ,ID_XC_VWN5, &
        &   ID_XC_WI,ID_XC_PW91,ID_XC_AM05,ID_XC_OPTB88PBE,ID_XC_OPTB86BPBE,ID_XC_OPTPBE, &
        &   ID_XC_PBE,ID_XC_PBE_X,ID_XC_PBE_C,ID_XC_PZ_C, &
        &   ID_XC_PBESOL,ID_XC_REVPBE,ID_XC_RPBE,ID_XC_B3,ID_XC_B5,ID_XC_ACFDT_RA, &
        &   ID_XC_ACFDT_PL,ID_XC_ACFDT_03,ID_XC_ACFDT_05,ID_XC_ACFDT_10,ID_XC_ACFDT_20, &
        &   ID_XC_BEEF,ID_XC_PW86RPBE,ID_XC_CX,ID_XC_M06L,ID_XC_MBJ,ID_XC_MS0, &
        &   ID_XC_MS1,ID_XC_MS2,ID_XC_R2SCAN,ID_XC_RSCAN,ID_XC_RTPSS,ID_XC_SCAN, &
        &   ID_XC_TPSS,ID_XC_LMBJ,ID_XC_RSCANZT,ID_XC_RSCANZT_X,ID_XC_RSCANZT_C, &
        &   ID_XC_RTPSS_X,ID_XC_RTPSS_C,ID_XC_TPSS_X, &
        &   ID_XC_TPSS_C,ID_XC_M06L_X,ID_XC_M06L_C,ID_XC_MS0_X,ID_XC_MS1_X,ID_XC_MS2_X, &
        &   ID_XC_MS0_C,ID_XC_MS1_C,ID_XC_MS2_C,ID_XC_SCAN_X,ID_XC_SCAN_C,ID_XC_RSCAN_X, &
        &   ID_XC_RSCAN_C,ID_XC_R2SCAN_X,ID_XC_R2SCAN_C,ID_XC_LIBXC,ID_XC_PBETEST, &
        &   ID_XC_SCANL,ID_XC_RSCANL,ID_XC_R2SCANL,ID_XC_OFR2,ID_XC_SREGTM1,ID_XC_SREGTM2, &
        &   ID_XC_SREGTM3,ID_XC_SREGTM2L,ID_XC_TASK_X,ID_XC_CC_C, &
        &   ID_XC_LAK,ID_XC_LAK_X,ID_XC_LAK_C,ID_XC_PW92,ID_XC_PW92_C, &
        &   ID_XC_MSPBEL,ID_XC_MSRPBEL,ID_XC_MSB86BL,ID_XC_RMSPBEL,ID_XC_RMSRPBEL,ID_XC_RMSB86BL
      END ENUM
      END MODULE xc_name

      MODULE xc_lda_table_name
      ENUM, BIND(C)
         ENUMERATOR :: ID_XC_LDA_TABLE_NOXC,ID_XC_LDA_TABLE_SLATER,ID_XC_LDA_TABLE_HL, &
        &   ID_XC_LDA_TABLE_PZ,ID_XC_LDA_TABLE_VWN5,ID_XC_LDA_TABLE_WI,ID_XC_LDA_TABLE_PZ_C, &
        &   ID_XC_LDA_TABLE_PW92,ID_XC_LDA_TABLE_B3,ID_XC_LDA_TABLE_B5,ID_XC_LDA_TABLE_ACFDT_RA, &
        &   ID_XC_LDA_TABLE_ACFDT_03,ID_XC_LDA_TABLE_ACFDT_05,ID_XC_LDA_TABLE_ACFDT_10, &
        &   ID_XC_LDA_TABLE_ACFDT_20,ID_XC_LDA_TABLE_ACFDT_PL,ID_XC_LDA_TABLE_X,ID_XC_LDA_TABLE_PW92_C, &
        &   ID_XC_LDA_TABLE_PW92_RELA
      END ENUM
      END MODULE xc_lda_table_name

      MODULE xc_family_name
      ENUM, BIND(C)
!XC_FAM_XXX is here used because XC_FAMILY_XXX already exist (imported from Libxc).
         ENUMERATOR :: XC_FAM_UNKNOWN,XC_FAM_NONE,XC_FAM_LDA,XC_FAM_GGA,XC_FAM_MGGA
      END ENUM
      END MODULE xc_family_name
