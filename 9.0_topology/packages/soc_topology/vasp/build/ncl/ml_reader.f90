# 1 "ml_reader.F"
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


# 2 "ml_reader.F" 2 

!==============================================================================
!> Contains the reader for machine-learned force fields method.
!==============================================================================
MODULE ml_reader

   USE prec
   IMPLICIT NONE



!===========================================================================
!> Generic procedure to read in a single ML INCAR tag.
!>
!> Each variant searches and assigns given ML tag (+ alternatives) in
!> contents of INCAR file. Actual function body is in ml_reader.inc.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer
!> arrays:
!> ~~~
!>     #define INCLUDE_READ_TAG_BODY
!>     #define TYPE_REAL
!>     #define RANK 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef RANK
!>     #undef TYPE_REAL
!>     #undef INCLUDE_READ_TAG_BODY
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!===========================================================================
   INTERFACE READ_TAG
      PROCEDURE READ_TAG_REAL0
      PROCEDURE READ_TAG_INTEGER0
      PROCEDURE READ_TAG_LOGICAL0
      PROCEDURE READ_TAG_CHARACTER0
      PROCEDURE READ_TAG_REAL1
      PROCEDURE READ_TAG_INTEGER1
   END INTERFACE READ_TAG



   CONTAINS



!===========================================================================
!> Reads INCAR tags for machine-learned force fields.
!>
!> Makes use of ml_reader::READ_TAG generic interface to read in variables
!> of different type.
!===========================================================================
   SUBROUTINE MACHINE_LEARNING_READER (ML_LMLFF,         &
                                       IC,               &
                                       TL,               &
                                       FFM,              &
                                       FF,               &
                                       ML_INPUT_HANDLE,  &
                                       NTYP,             &
                                       IWAVPR,           &
                                       CURRENT_INSTANCE, &
                                       INSTANCES,        &
                                       NSW,              &
                                       IBRION)
      USE tutor, ONLY: vtutor, isError, RandomSeed, newArgument
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE string, ONLY: str
      USE ml_ff_constant, ONLY: EUNIT, FUNIT, SUNIT, AUTOA
      USE ml_ff_taglist
      USE ml_ff_struct, ONLY: FFM_PAR,       &
                              FF_PAR,        &
                              ML_FF_INPUT,   &
                              ML_MB_DEFAULT, &
                              ML_MCONF_DEFAULT
! Input variables
!> Whether ML force fields are enabled (this should always be true
!> if this subroutine is called).
      LOGICAL,           INTENT(IN)    :: ML_LMLFF
!> The INCAR file object potentially containing MLFF related tags.
      TYPE(INCAR_FILE),  INTENT(INOUT) :: IC
!> Stores whether MLFF tags were set to default or user-defined, etc.
      TYPE(TAGLIST),     INTENT(INOUT) :: TL
!> Force field setup.
      TYPE(FFM_PAR),     INTENT(INOUT) :: FFM
!> On-the-fly algorithm setup.
      TYPE(FF_PAR),      INTENT(INOUT) :: FF
!> Type of algorithm for the movement of ions
      INTEGER,           INTENT(IN)    :: IBRION
!> Input structure from VASP.
      TYPE(ML_FF_INPUT), INTENT(INOUT) :: ML_INPUT_HANDLE
! Number of MD steps in VASP
      INTEGER,           INTENT(INOUT) :: NSW
!> Number of atom species.
      INTEGER,           INTENT(IN)    :: NTYP
!> ID of this ML instance.
      INTEGER,           INTENT(IN)    :: CURRENT_INSTANCE
!> Total number of ML instances.
      INTEGER,           INTENT(IN)    :: INSTANCES
      INTEGER,           INTENT(INOUT) :: IWAVPR
! Local variables
      CHARACTER(len=:), ALLOCATABLE :: PF
      CHARACTER(len=:), ALLOCATABLE :: HELPSTRING1
      CHARACTER(len=:), ALLOCATABLE :: HELPSTRING2
      INTEGER :: DTVALUES(8)
      INTEGER :: IDUMLONG
      INTEGER, PARAMETER :: SEED1_MAX = 900000000 
      INTEGER, PARAMETER :: SEED2_MAX = 1000000
      INTEGER :: IWAVPR_HELP
      INTEGER :: N
      INTEGER :: DEFAULT_ISTART
      LOGICAL :: DEFAULT_LFAST
      INTEGER :: DEFAULT_ICRITERIA
      REAL(q) :: DEFAULT_CDOUB
      INTEGER :: DEFAULT_OUTBLOCK
      INTEGER :: DEFAULT_IALGO_LINREG
      REAL(q) :: DEFAULT_EPS_LOW


! Initialize tag list (sets all internal states to TL_UNSET).
      CALL TL%INIT()
! Set the prefix for the reading
      IF (INSTANCES == 1) THEN
         PF = trim("")
      ELSE
         PF = trim("ML_INSTANCE_" // str(CURRENT_INSTANCE) // "/")
      ENDIF

! Allocate arrays necessary to read INCAR file.
      IF (ALLOCATED(ML_INPUT_HANDLE%EATOM_VASP)) THEN
         DEALLOCATE(ML_INPUT_HANDLE%EATOM_VASP)
      ENDIF
      ALLOCATE(ML_INPUT_HANDLE%EATOM_VASP(1:NTYP))
! Reading of the main variable ML_LMLFF is 1._q in READER!!!
! If ML_LMLFF=.TRUE., the parameters related to MLFF are read in the
! INCAR file.
      IF(ML_LMLFF) THEN
! check here if IVAWPR was set
         IF (COUNT_ELEMENTS(IC, PF//'IWAVPR') > 0) THEN
            CALL PROCESS_INCAR(IC, PF//'IWAVPR', IWAVPR_HELP)
         ELSE IF(COUNT_ELEMENTS(IC, 'IWAVPR') > 0) THEN
            CALL PROCESS_INCAR(IC, 'IWAVPR', IWAVPR_HELP)
         ELSE
            IWAVPR = 11
         END IF

! At the beginning ML_MODE is read which sets some defaults
         HELPSTRING2 = "NONE"
         N = READ_TAG(IC, TL, PF, 'ML_MODE', HELPSTRING1, HELPSTRING2)
         FF%MODE=TRIM(ADJUSTL(HELPSTRING1))

         DEFAULT_LFAST        = .FALSE.
         DEFAULT_ICRITERIA    = 1
         DEFAULT_CDOUB        = 2.0_q
         DEFAULT_OUTBLOCK     = 1
         DEFAULT_IALGO_LINREG = 1
         DEFAULT_EPS_LOW      = 1.0E-9_q

         IF (FF%MODE.EQ."train" .OR. FF%MODE.EQ."TRAIN") THEN
            IF (FF%LMLABEXIST) THEN
               DEFAULT_ISTART = 1
            ELSE
               DEFAULT_ISTART = 0
            ENDIF
         ELSE IF (FF%MODE.EQ."refit" .OR. FF%MODE.EQ."REFIT") THEN
            DEFAULT_ISTART       = 4
            NSW                  = 1
            DEFAULT_LFAST        = .TRUE.
            DEFAULT_IALGO_LINREG = 4
            DEFAULT_EPS_LOW      = 1.0E-11_q
         ELSE IF (FF%MODE.EQ."refitbayesian" .OR. &
                  FF%MODE.EQ."REFITBAYESIAN") THEN
            DEFAULT_ISTART  = 4
            NSW             = 1
            DEFAULT_EPS_LOW = 1.0E-11_q
         ELSE IF (FF%MODE.EQ."select" .OR. FF%MODE.EQ."SELECT") THEN
            DEFAULT_ISTART    = 3
            NSW               = 1
            DEFAULT_ICRITERIA = 3
            DEFAULT_CDOUB     = 100.0_q
         ELSE IF (FF%MODE.EQ."run" .OR. FF%MODE.EQ."RUN") THEN
            DEFAULT_ISTART      = 2
            DEFAULT_OUTBLOCK    = 1
         ELSE IF (FF%MODE.EQ."none" .OR. FF%MODE.EQ."NONE") THEN

            DEFAULT_ISTART = 0
# 192

         ELSE
            CALL vtutor%error("ERROR: Option for ML_MODE not recognized. &
                              &Please specify a valid option.")
         ENDIF
! We need to set FF%NSW since it is needed very soon
         FF%NSW=NSW
! Read in whether VASPml C++ library should be used.
! Note (aS): In general, VASPml should be used by default if the mode
! is supported, hence -Dlibvaspml sets the default to .TRUE. here.
! However, for unsupported modes we need to revert to .FALSE. (even
! when we need to override the user INCAR setting). This is 1._q at
! the end of the reader in the POST-READER section.
# 207

         N = READ_TAG(IC, TL, PF, 'ML_LIB', FF%LIB, .FALSE.)

         N = READ_TAG(IC, TL, PF, 'ML_NCSHMEM', FF%NCSHMEM, -1)
! Read FF%CALGO
!    FF%CALGO = 0 : Bayesian error estimation used for error decision
!    FF%CALGO = 1 : Spilling factor used for error decision
         N = READ_TAG(IC, TL, PF, 'ML_CALGO', FF%CALGO, 0)
         IF (FF%CALGO.LT.0 .OR. FF%CALGO.GT.1) THEN
            CALL vtutor%error("ERROR: Option for ML_CALGO not recognized. &
                              &Please specify a valid option.")
         ENDIF
! Read FF%CDOUB.
         N = READ_TAG(IC, TL, PF, 'ML_CDOUB', FF%CDOUB, DEFAULT_CDOUB)

! Read FF%CSIG that is a threshold for variance in the stored
! estimated errors used to determine the criteria when
! FF%ICRITERIA > 0.
         N = READ_TAG(IC, TL, PF, 'ML_CSIG', FF%CSIG, 4E-01_q)

! Read FF%CSLOPE that is a threshold for slope of the stored estimated
! errors used to determine the criteria when FF%ICRITERIA > 0.
         N = READ_TAG(IC, TL, PF, 'ML_CSLOPE', FF%CSLOPE, 2E-01_q)

! Read ML_CTIFOR which controls the selected number of local
! reference configurations
! Depending on ML_CALGO the criterion is either the criterion for
! Bayesian error estimation or spilling factor
         IF (FF%CALGO .EQ. 0) THEN
            N = READ_TAG(IC, TL, PF, 'ML_CTIFOR', FF%CTIFOR, 0.002_q)
         ELSE
            N = READ_TAG(IC, TL, PF, 'ML_CTIFOR', FF%CTIFOR, 0.02_q)
         ENDIF

! Read ML_SCLC_CTIFOR which controls the selected number of local
! reference configurations
         N = READ_TAG(IC, TL, PF, 'ML_SCLC_CTIFOR', FF%SCALE_CTIFOR, 0.6_q)

! Read reference atomic energies, which should be energies of isolated
! atoms.
         N = READ_TAG(IC, TL, PF,                 &
                      'ML_EATOM_REF',             &
                      ML_INPUT_HANDLE%EATOM_VASP, &
                      0.0_q)
         ML_INPUT_HANDLE%EATOM_VASP(1:NTYP) = &
            ML_INPUT_HANDLE%EATOM_VASP(1:NTYP) / EUNIT
! Reading EPS_LOW
         N = READ_TAG(IC, TL, PF, 'ML_EPS_LOW', FF%ML_EPS_LOW, DEFAULT_EPS_LOW)
! Reading EPS_REG
         N = READ_TAG(IC, TL, PF, 'ML_EPS_REG', FF%ML_EPS_REG, 1E-15_q)
! Read FF%IMAT_SPARS, which specifies the type of cutoff function.
!    FF%IMAT_SPARS = 1 : Covariance matrix.
!    FF%IMAT_SPARS = 2 : SOAP matrix.
         N = READ_TAG(IC, TL, PF, 'ML_IMAT_SPARS', FF%IMAT_SPARS, 1)
! Read FF%ISTART.
!    ML_ISTART = 0 : Start the calculation without any previous data.
!    ML_ISTART = 1 : Restart the calculation after reading the
!                    previous database file.
!    ML_ISTART = 2 : Only force field is used without learning.
!    ML_ISTART = 3 : New force-field created only from ML_AB file.
         N = READ_TAG(IC, TL, PF, 'ML_ISTART', FF%ISTART, DEFAULT_ISTART)

# 275


         IF (FF%ISTART == 3) THEN
            CALL vtutor%alert( &
               "ML_MODE = select was selected. Structure in POSCAR is ignored, &
               &hence no prediction of error, forces and stress will be &
               &available on that structure at the end of this calculation.")
         ENDIF

! If ML_FF file is missing for ML_ISTART=2 we exit
         IF (FF%ISTART == 2) THEN
            IF (.NOT.FF%LMLFFEXIST) THEN
               CALL vtutor%error("ERROR: ML_FF is missing, exiting...")
            ENDIF
         ENDIF
! If ML_AB file is missing for ML_ISTART=1,3,4 we exit
         IF (FF%ISTART == 1 .OR. FF%ISTART == 3 .OR. FF%ISTART == 4) THEN
            IF (.NOT.FF%LMLABEXIST) THEN
               CALL vtutor%error("ERROR: ML_AB is missing, exiting...")
            ENDIF
         ENDIF
! Set LMLONLY_FF to .FALSE. except for ML_ISTART=2
         FF%LMLONLY = .FALSE.
         IF (FF%ISTART == 2) THEN
            FF%LMLONLY = .TRUE.
         ENDIF
! Read FF%ISVD.
         N = READ_TAG(IC, TL, PF, 'ML_ISVD', FF%ISVD, 1)
! Read FF%IUPDATE_CRITERIA.
! This variable should be always 1, do not change it !!!
            N = READ_TAG(IC, TL, PF,            &
                         'ML_IUPDATE_CRITERIA', &
                         FF%IUPDATE_CRITERIA,   &
                         1)
! Read FF%LBASIS_DISCARD.
!   If FF%LBASIS_DISCARD = .TRUE., the basis sets are thrown away if
!   NB > FFM%MB.
         N = READ_TAG(IC, TL, PF,          &
                      'ML_LBASIS_DISCARD', &
                      FF%LBASIS_DISCARD,   &
                      .TRUE.)
! Read LCCONF_DISCARD_FF.
!   If FF%LCONF_DISCARD = .TRUE., the training sets are thrown away if
!   NCONF > FF%MCONF.
         N = READ_TAG(IC, TL, PF, 'ML_LCONF_DISCARD', FF%LCONF_DISCARD, .FALSE.)
! Read LUSE_NAMES
         N = READ_TAG(IC, TL, PF, 'ML_LUSE_NAMES', FF%LUSE_NAMES,.FALSE.)
! Read FF%ICRITERIA.
         IF (FF%CALGO.EQ.0) THEN
            N = READ_TAG(IC, TL, PF,     &
                         'ML_ICRITERIA', &
                         FF%ICRITERIA,   &
                         DEFAULT_ICRITERIA)
         ELSE IF (FF%CALGO.EQ.1) THEN
            N = READ_TAG(IC, TL, PF,     &
                         'ML_ICRITERIA', &
                         FF%ICRITERIA,   &
                         0)
         ENDIF
         IF (FF%CALGO.EQ.1 .AND. FF%ICRITERIA.GT.0) THEN
            CALL vtutor%error("ERROR: For ML_CALGO=1, &
                              &ML_ICRITERIA=0 must be set.")
         ENDIF
! Read FF%LMLMB.
!   FF%LMLMB = .FALSE. : No many-body interactions. But currently,
!   this program does not work when FF%LMLMB is set as .FALSE.
!   FF%LMLMB = .TRUE.  : Execute force field generation with the
!   many-body interactions.
         N = READ_TAG(IC, TL, PF, 'ML_LMLMB', FF%LMLMB, .TRUE.)
! Read ML_INPUT_HANDLE%NDIM_SCALAPACK_FF
         N = READ_TAG(IC, TL, PF, 'ML_NDIM_SCALAPACK',                &
                       ML_INPUT_HANDLE%NDIM_SCALAPACK_FF, &
                       2)
! Read ML_INPUT_HANDLE%NWRITE_FF.
         N = READ_TAG(IC, TL, PF, 'ML_FF_NWRITE', ML_INPUT_HANDLE%NWRITE_FF, 2)
! Read FF%IWEIGHT
         N = READ_TAG(IC, TL, PF, 'ML_IWEIGHT', FF%IWEIGHT, 3)
! Tag to control wether ML_NMDINT is multiplied with a random
! number between 0 and 1
! If .TRUE. then NMDINT is varying between 1 and NMDINT
         N = READ_TAG(IC, TL, PF,          &
                      'ML_LNMDINT_RANDOM', &
                      FF%LNMDINT_RANDOM,   &
                      .FALSE.)
! Read FF%LTEST.
!   FF%LTEST = .TRUE.  : Execute ab initio calculations to get test
!                        data.
!   FF%LTEST = .FALSE. : Do not execute ab initio calculations to get
!                        test data.
         N = READ_TAG(IC, TL, PF, 'ML_LTEST', FF%LTEST, .FALSE.)
! Read FF%LTRJ.
!   FF%LTRJ = .TRUE.  : MD trajectory is always obtained using ab
!                       initio forces and stress tensors.
!   FF%LTRJ = .FALSE. : MD trajectory is obtained using generating
!                       force field when the estimated error is small.
!                       Otherwise, ab initio forces and stress tensors
!                       are used.
         N = READ_TAG(IC, TL, PF, 'ML_LTRJ', FF%LTRJ, .FALSE.)
! READ FF%LFAST
!   FF%LFAST = .TRUE.  : Execute fast ISTART=2, no error estimate
!   FF%LFAST = .FALSE. : Execute slow ISTART=2, withe error estimate
         N = READ_TAG(IC, TL, PF, 'ML_LFAST', FF%LFAST, DEFAULT_LFAST)
         IF (FF%LFAST) THEN
            IF (FF%ISTART == 1 .AND. FF%NSW > 1) THEN
               CALL vtutor%error("ERROR: ML_LFAST only available for &
                                 &ML_ISTART=1 and NSW=0,1")
            ENDIF
            IF (FF%ISTART == 0) THEN
               CALL vtutor%error("ERROR: ML_LFAST not available for &
                                 &ML_ISTART=0")
            ENDIF
            IF (FF%ISTART == 3) THEN
               CALL vtutor%error("ERROR: ML_LFAST not available for &
                                 &ML_ISTART=3")
            ENDIF
         ENDIF
! READ FF%LFORCESLOW
!   FF%LFORCESLOW = .TRUE.  : Slow algorithm is enforced for force field
!    trained with ML_LFAST. No error estimation available.
         N = READ_TAG(IC, TL, PF, 'ML_LFORCESLOW', FF%LFORCESLOW, .FALSE.)
! Choose which algorithm is used for linear regression
!    ML_IALGO_LINREG = 1 : Solve L2-norm (ridge regression)
!    ML_IALGO_LINREG = 2 : Solve L1-norm with QR factorization
!    ML_IALGO_LINREG = 3 : Solve L1-norm with truncated SVD
!    ML_IALGO_LINREG = 4 : SVD + Tikhonov regularization
         N = READ_TAG(IC, TL, PF,        &
                      'ML_IALGO_LINREG', &
                      FF%IALGO_LINREG,   &
                      DEFAULT_IALGO_LINREG)
         IF (FF%LFAST .AND. FF%IALGO_LINREG < 3. .AND. FF%ISTART.NE.2) THEN
            CALL vtutor%error("ERROR: ML_LFAST only available for &
                              &ML_IALGO_LINREG>2")
         ENDIF
! READ FF%LTOTEN_SYSTEM
!    FF%LTOTEN_SYSTEM = .TRUE. : Total energy of system instead
!                                energy/atom is learned
!    FF%LTOTEN_SYSTEM = .FALSE. : Energy/atom is learned
         N = READ_TAG(IC, TL, PF, 'ML_LTOTEN_SYSTEM', FF%LTOTEN_SYSTEM, .FALSE.)
! Read FF%MCONF that is the maximum number of configurations used as
! training data.
         N = READ_TAG(IC, TL, PF, 'ML_MCONF', FF%MCONF, ML_MCONF_DEFAULT)
! Read FF%MCONF_NEW that is the maximum number of new configurations
! stored as candidates of training data temporally.
         N = READ_TAG(IC, TL, PF, 'ML_MCONF_NEW', FF%MCONF_NEW, 5)
! Read FF%MHIS that is the number of estimated errors stored in memory
! to judge the criteria.
         N = READ_TAG(IC, TL, PF, 'ML_MHIS', FF%MHIS, 10)
! Read FF%NMDINT that is the number of interval where machine do not
! collect any samples.
         IF (FF%ISTART == 3) THEN
            N = READ_TAG(IC, TL, PF, 'ML_NMDINT', FF%NMDINT, 1)
         ELSE
            N = READ_TAG(IC, TL, PF, 'ML_NMDINT', FF%NMDINT, 10)
         ENDIF
! Read FF%NTEST.
! If FF%LTEST = .TRUE., test data will be obtained at every NTEST MD
! step.
         IF(FF%LTEST) THEN
            N = READ_TAG(IC, TL, PF, 'ML_NTEST', FF%NTEST, 10)
         ENDIF
! Read FF%OUTPUT_MODE
         IF (FF%ISTART == 2) THEN
            N = READ_TAG(IC, TL, PF, 'ML_OUTPUT_MODE', FF%OUTPUT_MODE, 0)
         ELSE
            N = READ_TAG(IC, TL, PF, 'ML_OUTPUT_MODE', FF%OUTPUT_MODE, 1)
         ENDIF
! Read FF%OUTBLOCK
         N = READ_TAG(IC, TL, PF, 'ML_OUTBLOCK', FF%OUTBLOCK, DEFAULT_OUTBLOCK)
         IF (FF%OUTBLOCK < 1) THEN
            CALL vtutor%error("ML_INTERFACE: ERROR: ML_OUTBLOCK must be &
                              &larger than or equal to 1.")
         ENDIF
! ML_OUTBLOCK should not be larger than MAX(NSW,1)
         IF (FF%OUTBLOCK>MAX(NSW,1)) THEN
            CALL vtutor%alert("ML_OUTBLOCK > max(NSW,1), and will be reset to &
                              &ML_OUTBLOCK = "//str(MAX(NSW,1))//".")
            FF%OUTBLOCK = MAX(NSW,1)
         ENDIF
! ML_OUTBLOCK/=1 is only supported for ML_ISTART=2 and IBRION=0
         IF (FF%ISTART /= 2 .AND. FF%OUTBLOCK /= 1) THEN
            CALL vtutor%alert("ML_OUTBLOCK > 1 is only supported for ML_MODE = run, &
                              &defaulting back to ML_OUTBLOCK = 1 ...")
            FF%OUTBLOCK = 1
         ENDIF
         IF (FF%OUTBLOCK>1 .AND. IBRION/=0) THEN
            CALL vtutor%alert("ML_OUTBLOCK > 1 is only supported for IBRION = 0, &
                              &defaulting back to ML_OUTBLOCK = 1 ...")
            FF%OUTBLOCK = 1
         ENDIF

! Read FF%IERR.
         IF (FF%ISTART == 2) THEN
            N = READ_TAG(IC, TL, PF, 'ML_IERR', FF%IERR, 0)
         ELSE
            N = READ_TAG(IC, TL, PF, 'ML_IERR', FF%IERR, 1)
         ENDIF

         IF (FF%NSW .LT. FF%IERR) THEN
            CALL vtutor%alert("The number of ionic steps (NSW) &
                 &is smaller than ML_IERR. Hence no error &
                 &estimation will be calculated &
                 &and written out.")
         ENDIF
         IF (FF%ISTART /= 2 .AND. FF%IERR/=1) THEN
            CALL vtutor%error("ERROR: ML_IERR=1 must be set for &
                              &ML_MODE/=run.")
         ENDIF
! Set here FF%IERR correctly
         IF (FF%ISTART == 2) THEN
            IF (FF%IERR .GT. 0) THEN
               IF (FF%IERR .LT. FF%OUTBLOCK) THEN
                  CALL vtutor%error("ML_IERR is smaller than &
                     &ML_OUTBLOCK. ML_IERR needs to be set as an &
                     &integer multiple of ML_OUTBLOCK.")  
               ELSE
                  IF (MOD(FF%IERR,FF%OUTBLOCK) .NE. 0) THEN
                     CALL vtutor%error("ERROR: ML_IERR must be &
                              &integer divisible by ML_OUTBLOCK. &
                              &ML_IERR needs to be set as an &
                              &integer multiple of ML_OUTBLOCK.")
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (FF%LFORCESLOW .AND. FF%IERR>0) THEN
            CALL vtutor%error("ERROR: ML_IERR=0 must be set for &
                              &ML_LFORCESLOW = .TRUE.")
         ENDIF

! Read FF%WTIFOR, FF%WTOTEN and FF%WTSIF
! in (eV Angst^-1), (eV atom^-1) and (kbar), respectively.
         IF(FF%IWEIGHT == 1) THEN
            N = READ_TAG(IC, TL, PF, 'ML_WTIFOR', FF%WTIFOR, 5E-02_q)
            N = READ_TAG(IC, TL, PF, 'ML_WTOTEN', FF%WTOTEN, 5E-03_q)
            N = READ_TAG(IC, TL, PF, 'ML_WTSIF', FF%WTSIF, 5.0_q)
         ELSE
            N = READ_TAG(IC, TL, PF, 'ML_WTIFOR', FF%WTIFOR, 1.0_q)
            N = READ_TAG(IC, TL, PF, 'ML_WTOTEN', FF%WTOTEN, 1.0_q)
            N = READ_TAG(IC, TL, PF, 'ML_WTSIF', FF%WTSIF, 1.0_q)
         ENDIF
! Read FF%CX (parameter x for criteria update).
         N = READ_TAG(IC, TL, PF, 'ML_CX', FF%CX, 0.0_q)
! If FF%LMLMB = .TRUE., read variables necessary for the many-body
! fitting.
         IF (FF%LMLMB) THEN
            IF (FF%LFAST) THEN
               FFM%LSIC=.FALSE.
            ELSE
               FFM%LSIC=.TRUE.
            ENDIF
! Read FFM%W1 that
            N = READ_TAG(IC, TL, PF, 'ML_W1', FFM%W1, 0.1_q)
! Read FFM%W2 that
! We disable the reading of ML_W2 for the moment because
! everything is determined by ML_W1
            FFM%W2 = 1.0_q - FFM%W1
!N = READ_TAG(IC, TL, PF, 'ML_W2', FFM%W2, 1.0_q - FFM%W1)
! Read ML_DESC_TYPE
! ML_DESC_TYPE = 0: Standard descriptor
! ML_DESC_TYPE = 1: Linear-scaling with element types
! ML_DESC_TYPE = 2: dual: no-element + default
! ML_DESC_TYPE = 3: Linear-scaling with MRB
! ML_DESC_TYPE = 4: dual: default + default
            N = READ_TAG(IC, TL, PF, 'ML_DESC_TYPE', FFM%DESC_TYPE, 0)
! Linear-scaling descriptor only allowed for LSIC=.FALSE.
            IF (FFM%DESC_TYPE>0) THEN
               FFM%LSIC=.FALSE.
!               IF (FF%MODE.EQ."refitbayesian" .OR. FF%MODE.EQ."REFITBAYESIAN") THEN
!                  CALL vtutor%error("ML_INTERFACE: ERROR: ML_MODE = refitbayesian &
!                                   &selected together with ML_DESC_TYPE > 0. Only &
!                                   &ML_DESC_TYPE = 0 is available for this ML_MODE. &
!                                   &ML_DESC_TYPE > 0 is available for ML_MODE = refit.")
!               ENDIF
!               IF (FF%MODE.EQ."train" .OR. FF%MODE.EQ."TRAIN" .OR. &
!                 FF%MODE.EQ."select" .OR. FF%MODE.EQ."SELECT") THEN
!                  CALL vtutor%error("ML_INTERFACE: ML_DESC_TYPE > 0 is &
!                  &not allowed together with ML_MODE = train or &
!                  &ML_MODE = select.")
!               ENDIF
            ENDIF
            IF (FFM%DESC_TYPE>0 .AND. FFM%W2.LE.0.0_q) THEN
               CALL vtutor%error("ML_INTERFACE: ERROR: ML_DESC_TYPE must be 0 &
                                    &for ML_W2=0.0.")
            ENDIF
            IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4) THEN
               IF (FF%MODE.EQ."train" .OR. FF%MODE.EQ."TRAIN") THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: ML_DESC_TYPE=2,4 not &
                                    &available for ML_MODE=TRAIN.")
               ENDIF
               IF (FF%MODE.EQ."SELECT" .OR. FF%MODE.EQ."select") THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: ML_DESC_TYPE=2,4 not &
                                    &available for ML_MODE=SELECT.")
               ENDIF
               IF (FF%ISTART.EQ.0) THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: ML_DESC_TYPE=2,4 not &
                                    &available for ML_ISTART=0.")
               ENDIF
               IF (FF%ISTART.EQ.1) THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: ML_DESC_TYPE=2,4 not &
                                    &available for ML_ISTART=1.")
               ENDIF
               IF (FF%ISTART.EQ.3) THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: ML_DESC_TYPE=2,4 not &
                                    &available for ML_ISTART=3.")
               ENDIF
               IF (FF%MODE.EQ."refit" .OR. FF%MODE.EQ."REFIT") THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: ML_DESC_TYPE=2,4 not &
                                    &available for ML_MODE=REFIT. Please select &
                                    &ML_MODE=REFITBAYESIAN.")
               ENDIF
            ENDIF
! READ FFM%RV%DESC_RATIO
            N = READ_TAG(IC, TL, PF, 'ML_DESC_RATIO_DUAL', FFM%RV%RATIO_DESC, 0.5_q)
            FFM%RV%RATIO_DESC_LOW_RES = 1.0_q - FFM%RV%RATIO_DESC
! READ
            N = READ_TAG(IC, TL, PF, 'ML_DESC_FACTORE_TESTE', FFM%RV%FACTOR_DESC, 1.0_q)
! Read FFM%LAFILT2.
            N = READ_TAG(IC, TL, PF, 'ML_LAFILT2', FFM%LAFILT2, .TRUE.)
            IF (FFM%LAFILT2) THEN
               N = READ_TAG(IC, TL, PF, 'ML_IAFILT2', FFM%IAFILT2, 2)
               IF (FFM%IAFILT2 == 2) THEN
                  N = READ_TAG(IC, TL, PF, 'ML_AFILT2', FFM%AFILT2, 2E-03_q)
               ENDIF
            ENDIF
! Read ML_IBROADX (X=1 or 2).
!  ! this variable determines the way to braoden the atomic distributions.
! ML_IBROADX=1 : No broadening.
! ML_IBROADX=2 : Element-independent user specified broadening.
! ML_IBROADX=3 : Element-dependent broadening
! Probably only ML_IBROADX=2 works!!!
            N = READ_TAG(IC, TL, PF, 'ML_IBROAD1', FFM%IBROAD1, 2)
            N = READ_TAG(IC, TL, PF, 'ML_IBROAD2', FFM%IBROAD2, 2)
! Set basis function type
! ML_BASIS_TYPEX=0 : Modified spherical Bessel functions
! ML_BASIS_TYPEX=1 : Gauss functions
            N = READ_TAG(IC, TL, PF, 'ML_BASIS_TYPE2', FFM%BASIS_TYPE2, 0)
            N = READ_TAG(IC, TL, PF, 'ML_BASIS_TYPE2_RED', FFM%RV%BASIS_TYPE2_RED, 4)
! Read ML_ICUTX (X=1 or 2), which specifies the type of
! cutoff function.
!   ML_ICUTX = 1 : Behler-Parinello type cutoff function.
!   ML_ICUTX = 2 : Miwa-Ohno type cutoff function.
            N = READ_TAG(IC, TL, PF, 'ML_ICUT1', FFM%ICUT1, 1)
            N = READ_TAG(IC, TL, PF, 'ML_ICUT2', FFM%ICUT2, 1)
! this variable determines the method to calculate inverse of SOAP
! matrix.
!   FFM%INVERSE_SOAP = 1 : Inverse matrix is calculated by using
!                          eigen-values and eigen-vectors.
!   FFM%INVERSE_SOAP = 2 : Inverse matrix is calculated by LU
!                          factorization.
            N = READ_TAG(IC, TL, PF, 'ML_INVERSE_SOAP', FFM%INVERSE_SOAP, 2)
! Read FFM%IREG
            N = READ_TAG(IC, TL, PF, 'ML_IREG', FFM%IREG, 2)
! Read FFM%ISCALE_TOTEN
            N = READ_TAG(IC, TL, PF, 'ML_ISCALE_TOTEN', FFM%ISCALE_TOTEN, 2)
! Read FFM%LSUPERVEC
            N = READ_TAG(IC, TL, PF, 'ML_LSUPERVEC', FFM%LSUPERVEC, .TRUE.)
!IF (FFM%DESC_TYPE.EQ.2 .AND. FFM%LSUPERVEC) THEN
!   CALL vtutor%error("ML_INTERFACE: ERROR: ML_LSUPERVEC must be .FALSE. &
!                        &for ML_DESC_TYPE=2.")
!ENDIF
! Read FFM%LEATOM
            N = READ_TAG(IC, TL, PF, 'ML_LEATOM', FFM%LEATOM, .FALSE.)
! Read FFM%LHEAT
            N = READ_TAG(IC, TL, PF, 'ML_LHEAT', FFM%LHEAT, .FALSE.)
! Read ML_LMETRICX.
            N = READ_TAG(IC, TL, PF, 'ML_LMETRIC1', FFM%LMETRIC1, .FALSE.)
            N = READ_TAG(IC, TL, PF, 'ML_LMETRIC2', FFM%LMETRIC2, .FALSE.)
! Read FFM%LSPARSDES and related tags.
            N = READ_TAG(IC, TL, PF, 'ML_LSPARSDES', FFM%LSPARSDES, .FALSE.)
            IF (FFM%LSPARSDES) THEN
               IF (FF%MODE.EQ."train" .OR. FF%MODE.EQ."TRAIN" .OR. & 
                 FF%MODE.EQ."select" .OR. FF%MODE.EQ."SELECT") THEN
                  CALL vtutor%error("ML_INTERFACE: ML_LSPARSDES = .TRUE. is &
                  &not allowed together with ML_MODE = train or &
                  &ML_MODE = select.")
               ENDIF
            ENDIF
            N = READ_TAG(IC, TL, PF, 'ML_NRANK_SPARSDES', FFM%NRANK_SPARSDES, 5)
            N = READ_TAG(IC, TL, PF,         &
                         'ML_RDES_SPARSDES', &
                         FFM%RDES_SPARSDES,  &
                         0.5_q)
! Sparsification of descriptors is only allowed in ISTART = 1 mode,
! or with ISTART = 4 ( refit ) mode.
! Hence overwrite user setting if set incorrectly. Related tags can
! be left unchanged, they will only be used if LSPARSDES = .TRUE.
            IF ( (FF%ISTART /= 1 .AND. FF%ISTART /= 4) .AND. &
                TL%GET_STATE('ML_LSPARSDES') /= TL_DEFAULT) THEN
               FFM%LSPARSDES = .FALSE.
               CALL TL%SET_STATE('ML_LSPARSDES', TL_OVERRIDE) 
            END IF
! Read ML_LVARTRANX.
            N = READ_TAG(IC, TL, PF, 'ML_LVARTRAN1', FFM%LVARTRAN1, .FALSE.)
            N = READ_TAG(IC, TL, PF, 'ML_LVARTRAN2', FFM%LVARTRAN2, .FALSE.)
! Read ML_NMETRICX.
            IF (FFM%LMETRIC1) THEN
               N = READ_TAG(IC, TL, PF, 'ML_NMETRIC1', FFM%NMETRIC1, 6)
            ENDIF
            IF (FFM%LMETRIC2) THEN
               N = READ_TAG(IC, TL, PF, 'ML_NMETRIC2', FFM%NMETRIC2, 6)
            ENDIF
! Read ML_NVARTRANX.
            IF (FFM%LVARTRAN1) THEN
               N = READ_TAG(IC, TL, PF, 'ML_NVARTRAN1', FFM%NVARTRAN1, 6)
            ENDIF
            IF (FFM%LVARTRAN2) THEN
               N = READ_TAG(IC, TL, PF, 'ML_NVARTRAN2', FFM%NVARTRAN2, 6)
            ENDIF
! Read LWINDOWX_FFY (X=1, 2 or nothing, and Y=M or Q)
!   ML_LWINDOWX_Y = .TRUE. : window function is used to smoothen
!                               the expanded atomic distribution.
!   ML_LWINDOWX_Y = .FLASE.: window function is not used.
            N = READ_TAG(IC, TL, PF, 'ML_LWINDOW1', FFM%LWINDOW1, .FALSE.)
            N = READ_TAG(IC, TL, PF, 'ML_LWINDOW2', FFM%LWINDOW2, .FALSE.)
! Read IWINDOWX_FFY (X=1, 2 or nothing, Y=M or Q) if LWINDOWX_FFY =
! .TRUE.
            IF (FFM%LWINDOW1) THEN
               N = READ_TAG(IC, TL, PF, 'ML_IWINDOW1', FFM%IWINDOW1, 6)
            ENDIF
            IF (FFM%LWINDOW2) THEN
               N = READ_TAG(IC, TL, PF, 'ML_IWINDOW2', FFM%IWINDOW2, 6)
            ENDIF
! Read FFM%LMAX2 that defines the maximum angular momentum number
! of spherical harmonics used to expand atomic distributions.
            IF (FFM%LAFILT2) THEN
               N = READ_TAG(IC, TL, PF, 'ML_LMAX2', FFM%LMAX2, 3)
            ELSE
               N = READ_TAG(IC, TL, PF, 'ML_LMAX2', FFM%LMAX2, 6)
            ENDIF
            N = READ_TAG(IC, TL, PF, 'ML_LMAX2_RED', FFM%RV%LMAX2_RED, 3)
! Read FFM%LNORM1 that specifies the execution of normalization.
            N = READ_TAG(IC, TL, PF, 'ML_LNORM1', FFM%LNORM1, .TRUE.)
! Read FFM%LNORM2 that specifies the execution of normalization.
            N = READ_TAG(IC, TL, PF, 'ML_LNORM2', FFM%LNORM2, .TRUE.)
! Read FFM%MB that is the maximum number of basis sets used for
! describing many-body interactions.
            N = READ_TAG(IC, TL, PF, 'ML_MB', FFM%MB, ML_MB_DEFAULT)
! Read FFM%MB_MIN which corresponds to the minimum number of local
! reference configurations required for the learning algorithm.
            N = READ_TAG(IC, TL, PF, 'ML_MB_MIN', FFM%MB_MIN, 3)
            IF (FFM%MB_MIN < 2) THEN
               CALL vtutor%error("ML_INTERFACE: The minimum number of local &
                  &reference configurations (ML_MB_MIN) must be at least 2.")
            ENDIF
! Read ML_MSPLX (X=1 or 2) that is the number of radial grids
! for spline-interpolations.
            N = READ_TAG(IC, TL, PF, 'ML_MSPL1', FFM%MSPL1, 1000)
            N = READ_TAG(IC, TL, PF, 'ML_MSPL2', FFM%MSPL2, FFM%MSPL1)
! Read NHYP that corresponds to the polynomial parameters of SOAP
! 1-point kernel.
            IF (FFM%LSUPERVEC) THEN
               N = READ_TAG(IC, TL, PF, 'ML_NHYP', FFM%NHYP1, 4)
            ELSE
               N = READ_TAG(IC, TL, PF, 'ML_NHYP', FFM%NHYP1, 1)
            ENDIF
! Read NHYP that corresponds to the polynomial parameters of SOAP
! 2-point kernel.
            N = READ_TAG(IC, TL, PF, 'ML_NHYP2', FFM%NHYP2, 4)
            IF (FFM%LSUPERVEC) THEN
               FFM%NHYP2 = FFM%NHYP1
            ENDIF
! Read ML_NRX (X=1 or 2) that is the number of radial grids
! used to execute radial integrations to compute the many-body
! descriptors.
            N = READ_TAG(IC, TL, PF, 'ML_NR1', FFM%NR1, 1000)
            N = READ_TAG(IC, TL, PF, 'ML_NR2', FFM%NR2, FFM%NR1)
! Read ML_RCUTX (X=1 or 2) that is the cutoff radius used for
! calculating many-body descriptors.
            N = READ_TAG(IC, TL, PF, 'ML_RCUT1', FFM%RCUT1, 8.0_q)
            N = READ_TAG(IC, TL, PF, 'ML_RCUT2', FFM%RCUT2, 5.0_q)
            IF (FFM%W1 == 0.0_q) THEN
               FFM%RCUT1 = FFM%RCUT2
            ENDIF
            IF (FFM%W2 == 0.0_q) THEN
               FFM%RCUT2 = FFM%RCUT1
            ENDIF
! Read ML_RMETRICX (X=1 or 2).
            IF (FFM%LMETRIC1) THEN
               N = READ_TAG(IC, TL, PF, 'ML_RMETRIC1', FFM%RMETRIC1, 1.0_q)
            ENDIF
            IF (FFM%LMETRIC2) THEN
               N = READ_TAG(IC, TL, PF, 'ML_RMETRIC2', FFM%RMETRIC2, 1.0_q)
            ENDIF
! Read FFM%SIGV0
            N = READ_TAG(IC, TL, PF, 'ML_SIGV0', FFM%SIGV0, 1.0_q)
! Read FFM%SIGW0
            IF (FF%IALGO_LINREG == 4) THEN
               N = READ_TAG(IC, TL, PF, 'ML_SIGW0', FFM%SIGW0, 1E-7_q)
            ELSE
               N = READ_TAG(IC, TL, PF, 'ML_SIGW0', FFM%SIGW0, 1.0_q)
            ENDIF
! Read ML_SIONX (X=1 or 2) that specifies the width of
! Gaussian functions used for broadening the atomic distributions.
            IF (FFM%IBROAD1 /= 1) THEN
               N = READ_TAG(IC, TL, PF, 'ML_SION1', FFM%SION1, 5.0E-01_q)
            ENDIF
            IF (FFM%IBROAD2 /= 1) THEN
               N = READ_TAG(IC, TL, PF, 'ML_SION2', FFM%SION2, FFM%SION1)
            ENDIF
! Read FFM%LCOUPLE that specifies whether the coupling parameter is
! introduced for the calculation of the chemical potential.
            N = READ_TAG(IC, TL, PF, 'ML_LCOUPLE', FFM%LCOUPLE, .FALSE.)
            IF (FFM%LCOUPLE) THEN
               N = READ_TAG(IC, TL, PF,         &
                            'ML_NATOM_COUPLED', &
                            FFM%NATOM_COUPLED,  &
                            0)
               IF (ALLOCATED(FFM%ICOUPLE)) THEN
                  DEALLOCATE(FFM%ICOUPLE)
               ENDIF
               ALLOCATE(FFM%ICOUPLE(1:FFM%NATOM_COUPLED))
               N = READ_TAG(IC, TL, PF, 'ML_ICOUPLE', FFM%ICOUPLE, 0)
               IF (N > 0 .AND. N /= FFM%NATOM_COUPLED) THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: Number of atoms for &
                                    &coupling .NE. ML_NATOM_COUPLED.")
               ENDIF
               N = READ_TAG(IC, TL, PF, 'ML_RCOUPLE', FFM%RCOUPLE, 1.0_q)
            ELSE
               FFM%NATOM_COUPLED = 0
               CALL TL%SET_STATE('ML_NATOM_COUPLED', TL_DEFAULT)
               IF(ALLOCATED(FFM%ICOUPLE)) THEN
                  DEALLOCATE(FFM%ICOUPLE)
               ENDIF
               ALLOCATE(FFM%ICOUPLE(1:FFM%NATOM_COUPLED))
               FFM%ICOUPLE = 0
               CALL TL%SET_STATE('ML_ICOUPLE', TL_DEFAULT)
               FFM%RCOUPLE = 1.0_q
               CALL TL%SET_STATE('ML_RCOUPLE', TL_DEFAULT)
            ENDIF
! Read ML_MRBX (X=1 or 2) that is the number of radial basis
! sets used to expand the atomic distributions.
! Automatically scale values of FFM%MRB2 if no value is supplied
            N = READ_TAG(IC, TL, PF, 'ML_MRB1', FFM%MRB1, 12)
            N = READ_TAG(IC, TL, PF, 'ML_MRB2', FFM%MRB2,  8)
            N = READ_TAG(IC, TL, PF, 'ML_MRB2_MAX', FFM%MRB2_MAX,  MAX(8,FFM%MRB2))
            IF (FFM%MRB2_MAX.LT.FFM%MRB2) THEN
               CALL vtutor%error("ML_INTERFACE: ERROR: ML_MRB2_MAX is smaller than ML_MRB2.  &
                                  &Please set ML_MRB2_MAX>=ML_MRB2.")
            ENDIF
            N = READ_TAG(IC, TL, PF, 'ML_MRB2_RED', FFM%RV%MRB2_RED,  2)
            N = READ_TAG(IC, TL, PF, 'ML_MRB2_MAX_RED', FFM%RV%MRB2_MAX_RED,  MAX(8,FFM%RV%MRB2_RED))
            IF (FFM%RV%MRB2_MAX_RED.LT.FFM%RV%MRB2_RED) THEN
               CALL vtutor%error("ML_INTERFACE: ERROR: ML_MRB2_MAX_RED is smaller than ML_MRB2_RED.  &
                                  &Please set ML_MRB2_MAX_RED>=ML_MRB2_RED.")
            ENDIF
! Use random seed if available
            N = READ_TAG(IC, TL, PF, 'ML_RANDOM_SEED', FF%SEED, 0)
            IF (N == 0) THEN
               CALL DATE_AND_TIME(VALUES=DTVALUES)
! with this choice we'r sure not to exceed
! SEED1_MAX (maximal value we can get here is 893581000)
               IDUMLONG = DTVALUES(3) * 24 * 60 * 60 * 300 &
                        + DTVALUES(5) * 60 * 60 * 1000 &
                        + DTVALUES(6) * 60 * 1000 &
                        + DTVALUES(7) * 1000 &
                        + DTVALUES(8)
               FF%SEED(1) = MOD(IDUMLONG, SEED1_MAX)
            ELSE IF ((FF%SEED(1) < 0) .OR. (FF%SEED(1) > SEED1_MAX) .OR. &
                     (FF%SEED(2) < 0) .OR. (FF%SEED(2) > SEED2_MAX) .OR. &
                     (FF%SEED(3) < 0)) THEN
               CALL vtutor%write(isError, RandomSeed, newArgument(ival = FF%SEED))     
            ENDIF
! Decide whether configurations not giving basis functions will be
! thrown away
            IF (FF%ISTART == 3) THEN
               N = READ_TAG(IC, TL, PF,                                &
                            'ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS', &
                            FF%LDISCARD_STRUCTURES_NOT_GIVING_BASIS,   &
                            .FALSE.)
            ENDIF

!======================================================================================
!
! POST-READER SANITY CHECKS AND WARNINGS
!
!======================================================================================
            IF (FF%ISTART == 3 .AND. FF%ICRITERIA == 3 .AND. ( &
                TL%GET_STATE("ML_CTIFOR") == TL_INCAR     .OR. &
                TL%GET_STATE("ML_CTIFOR") == TL_INCAR_ALT .OR. &
                TL%GET_STATE("ML_CX")     == TL_INCAR     .OR. &
                TL%GET_STATE("ML_CX")     == TL_INCAR_ALT .OR. &
                TL%GET_STATE("ML_CSLOPE") == TL_INCAR     .OR. &
                TL%GET_STATE("ML_CSLOPE") == TL_INCAR_ALT .OR. &
                TL%GET_STATE("ML_CSIG")   == TL_INCAR     .OR. &
                TL%GET_STATE("ML_CSIG")   == TL_INCAR_ALT .OR. &
                TL%GET_STATE("ML_MHIS")   == TL_INCAR     .OR. &
                TL%GET_STATE("ML_MHIS")   == TL_INCAR_ALT)) THEN
               CALL vtutor%alert( &
                  "An INCAR tag related to the Bayesian error &
                  &threshold was set (ML_CTIFOR, ML_CX, ML_CSLOPE, ML_CSIG, &
                  &ML_MHIS). However, these tags will not be used because &
                  &ML_ICRITERIA=3 which means that the Bayesian error &
                  &threshold will be read directly from the ML_AB file.")
            ENDIF
            IF (FF%ISTART /= 3 .AND. FF%ICRITERIA == 3) THEN
               CALL vtutor%error( &
                  "Reading Bayesian error threshold from ML_AB file &
                  &(ML_ICRITERIA=3) is only possible for ML_MODE=select. &
                  &Please select a different Bayesian error update method (&
                  &ML_ICRITERIA tag).")
            ENDIF

! Disable VASPml for unsupported modes.
! =====================================
! Currently supported:
! --------------------
! * ML_MODE = run, only with FF trained with ML_LFAST = .TRUE.
!
            IF (FF%LIB) THEN
! Prediction mode we can test here, the use of fast mode is only
! determined when reading the ML_FF file. This is 1._q later
! inside VASPml which must fail with error when a ML_FF with
! ML_LFAST = .FALSE. is read in.
               IF (.NOT. FF%ISTART == 2) THEN
! If user requested VASPml explicitly, issue warning
! and override INCAR settings.
                  IF (TL%GET_STATE("ML_LIB") == TL_INCAR) THEN
                     CALL vtutor%alert( &
                        "The selected mode of operation is not available &
                        &in the exerimental VASPml library, yet ML_LIB=.TRUE. &
                        &was set in the INCAR file. Resetting now &
                        &ML_LIB to .FALSE. and continuing without VASPml lib.")
                     CALL TL%SET_STATE('ML_LIB', TL_OVERRIDE)
                  ENDIF
                  FF%LIB = .FALSE.
               ENDIF
            ENDIF


! Set some important stuff
! Set some important variables
            FF%NHIS = 0
            FFM%LMAX1 = 0
            FFM%MB_INPUT = FFM%MB
            FFM%MB_ALLOCATE = FFM%MB
            FFM%NSPL1 = FFM%MSPL1
            FFM%NSPL2 = FFM%MSPL2

! Unit conversions
            IF (FF%IWEIGHT == 1) THEN
               FF%WTIFOR = FF%WTIFOR / FUNIT
               FF%WTOTEN = FF%WTOTEN / EUNIT
               FF%WTSIF = FF%WTSIF / SUNIT
            ENDIF
            IF (FF%WTSIF .LT. 1.E-13_q) THEN
               FF%WTSIF = 1.E-13_q
            ENDIF
            FF%CTIFOR = FF%CTIFOR / FUNIT
            FFM%RCUT1 = FFM%RCUT1 / AUTOA
            FFM%RCUT2 = FFM%RCUT2 / AUTOA
            FFM%SION1 = FFM%SION1 / AUTOA
            FFM%SION2 = FFM%SION2 / AUTOA
            IF (FFM%LMETRIC1) THEN
               FFM%RMETRIC1 = FFM%RMETRIC1 / AUTOA
            ENDIF
            IF (FFM%LMETRIC2) THEN
               FFM%RMETRIC2 = FFM%RMETRIC2 / AUTOA
            ENDIF
         ELSE
            FFM%IAFILT2 = 0
            FFM%INVERSE_SOAP = 0
            FFM%IREG = 0
            FFM%LMAX1 = 0
            FFM%LMAX2 = 0
            FFM%MB = 0
            FFM%MB_INPUT = 0
            FFM%MB_ALLOCATE = 0
            FFM%MSPL1 = 0
            FFM%MSPL2 = 0
            FFM%NHYP1 = 0
            FFM%NHYP2 = 0
            FFM%NR1 = 0
            FFM%NR2 = 0
            FFM%NSPL1 = 0
            FFM%NSPL2 = 0
            FFM%MRB1 = 0
            FFM%MRB2 = 0
            FFM%NMETRIC1 = 0
            FFM%NMETRIC2 = 0
            FFM%NVARTRAN1 = 0
            FFM%NVARTRAN2 = 0
            FFM%RCUT1 = 0.0_q
            FFM%RCUT2 = 0.0_q
            FFM%RMETRIC1 = 0.0_q
            FFM%RMETRIC2 = 0.0_q
            FFM%SIGV0 = 0.0_q
            FFM%SIGW0 = 0.0_q
            FFM%SION1 = 0.0_q
            FFM%SION2 = 0.0_q
            FFM%W1 = 0.0_q
            FFM%W2 = 0.0_q
         ENDIF
      ENDIF


! defined(1) || defined(MPI_CHAIN)

   END SUBROUTINE MACHINE_LEARNING_READER



! Overloaded READ_TAG functions generated with "template-like" inc file.




!> Read INCAR tag with (1._q,0._q) real value.
   INTEGER FUNCTION READ_TAG_REAL0 &
# 1 "./ml_reader.inc" 1 
! Fortran "template" workaround from incar_reader.inc.


# 10


# 16


# 20



!==============================================================================
!> Search and assign given ML tag (+ alternatives) in contents of INCAR file.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer arrays:
!> ~~~
!>     #define 
!>     #define 
!>     #define 0 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef 0
!>     #undef 
!>     #undef 
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!==============================================================================
!  INTEGER FUNCTION READ_TAG_...       &
                            (INCAR,    &
                             TAG_LIST, &
                             PREFIX,   &
                             TAG,      &
                             VARIABLE, &
                             DEFAULT_VALUE) RESULT(NARGS)
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE ml_ff_taglist
      IMPLICIT NONE

!> String representation of INCAR file.
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
!> List with all tag states.
      TYPE(TAGLIST),    INTENT(INOUT) :: TAG_LIST
!> INCAR tag prefix (in case of multiple ML instances).
      CHARACTER(len=*), INTENT(IN)    :: PREFIX
!> INCAR tag to be processed.
      CHARACTER(len=*), INTENT(IN)    :: TAG
!> Variable where to store default or INCAR tag argument.
      REAL(q),          INTENT(INOUT) :: VARIABLE 
!> Default value for variable.
      REAL(q),          INTENT(IN)    :: DEFAULT_VALUE

! Alternative tags
      CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
      INTEGER                         :: I

! Set default value.
      VARIABLE = DEFAULT_VALUE
      CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

! Collect alternative tag list.
      ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

      NARGS = 0
! Loop over alternative tags (includes current "preferred" tag as
! first entry).
      DO I = 1, SIZE(ALT_TAGS)
! Check for per-instance INCAR tag first.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR,                                &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF

! Now comes the same check without prefix.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))), VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF
      END DO

! If found anywhere, set tag list state accordingly.
      IF (NARGS > 0) THEN
         IF (I > 1) THEN
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
         ELSE
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
         END IF
      END IF

!  END FUNCTION READ_TAG_...



# 985 "ml_reader.F" 2 
   END FUNCTION READ_TAG_REAL0





!> Read INCAR tag with (1._q,0._q) integer value.
   INTEGER FUNCTION READ_TAG_INTEGER0 &
# 1 "./ml_reader.inc" 1 
! Fortran "template" workaround from incar_reader.inc.


# 10




# 20



!==============================================================================
!> Search and assign given ML tag (+ alternatives) in contents of INCAR file.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer arrays:
!> ~~~
!>     #define 
!>     #define TYPE_REAL
!>     #define 0 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef 0
!>     #undef TYPE_REAL
!>     #undef 
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!==============================================================================
!  INTEGER FUNCTION READ_TAG_...       &
                            (INCAR,    &
                             TAG_LIST, &
                             PREFIX,   &
                             TAG,      &
                             VARIABLE, &
                             DEFAULT_VALUE) RESULT(NARGS)
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE ml_ff_taglist
      IMPLICIT NONE

!> String representation of INCAR file.
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
!> List with all tag states.
      TYPE(TAGLIST),    INTENT(INOUT) :: TAG_LIST
!> INCAR tag prefix (in case of multiple ML instances).
      CHARACTER(len=*), INTENT(IN)    :: PREFIX
!> INCAR tag to be processed.
      CHARACTER(len=*), INTENT(IN)    :: TAG
!> Variable where to store default or INCAR tag argument.
      INTEGER,          INTENT(INOUT) :: VARIABLE 
!> Default value for variable.
      INTEGER,          INTENT(IN)    :: DEFAULT_VALUE

! Alternative tags
      CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
      INTEGER                         :: I

! Set default value.
      VARIABLE = DEFAULT_VALUE
      CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

! Collect alternative tag list.
      ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

      NARGS = 0
! Loop over alternative tags (includes current "preferred" tag as
! first entry).
      DO I = 1, SIZE(ALT_TAGS)
! Check for per-instance INCAR tag first.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR,                                &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF

! Now comes the same check without prefix.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))), VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF
      END DO

! If found anywhere, set tag list state accordingly.
      IF (NARGS > 0) THEN
         IF (I > 1) THEN
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
         ELSE
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
         END IF
      END IF

!  END FUNCTION READ_TAG_...



# 994 "ml_reader.F" 2 
   END FUNCTION READ_TAG_INTEGER0





!> Read INCAR tag with (1._q,0._q) logical value.
   INTEGER FUNCTION READ_TAG_LOGICAL0 &
# 1 "./ml_reader.inc" 1 
! Fortran "template" workaround from incar_reader.inc.


# 10


# 14


# 20



!==============================================================================
!> Search and assign given ML tag (+ alternatives) in contents of INCAR file.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer arrays:
!> ~~~
!>     #define 
!>     #define TYPE_REAL
!>     #define 0 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef 0
!>     #undef TYPE_REAL
!>     #undef 
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!==============================================================================
!  INTEGER FUNCTION READ_TAG_...       &
                            (INCAR,    &
                             TAG_LIST, &
                             PREFIX,   &
                             TAG,      &
                             VARIABLE, &
                             DEFAULT_VALUE) RESULT(NARGS)
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE ml_ff_taglist
      IMPLICIT NONE

!> String representation of INCAR file.
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
!> List with all tag states.
      TYPE(TAGLIST),    INTENT(INOUT) :: TAG_LIST
!> INCAR tag prefix (in case of multiple ML instances).
      CHARACTER(len=*), INTENT(IN)    :: PREFIX
!> INCAR tag to be processed.
      CHARACTER(len=*), INTENT(IN)    :: TAG
!> Variable where to store default or INCAR tag argument.
      LOGICAL,          INTENT(INOUT) :: VARIABLE 
!> Default value for variable.
      LOGICAL,          INTENT(IN)    :: DEFAULT_VALUE

! Alternative tags
      CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
      INTEGER                         :: I

! Set default value.
      VARIABLE = DEFAULT_VALUE
      CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

! Collect alternative tag list.
      ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

      NARGS = 0
! Loop over alternative tags (includes current "preferred" tag as
! first entry).
      DO I = 1, SIZE(ALT_TAGS)
! Check for per-instance INCAR tag first.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR,                                &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF

! Now comes the same check without prefix.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))), VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF
      END DO

! If found anywhere, set tag list state accordingly.
      IF (NARGS > 0) THEN
         IF (I > 1) THEN
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
         ELSE
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
         END IF
      END IF

!  END FUNCTION READ_TAG_...



# 1003 "ml_reader.F" 2 
   END FUNCTION READ_TAG_LOGICAL0





!> Read INCAR tag with (1._q,0._q) string value.
   INTEGER FUNCTION READ_TAG_CHARACTER0 &
# 1 "./ml_reader.inc" 1 
! Fortran "template" workaround from incar_reader.inc.


# 10


# 18





!==============================================================================
!> Search and assign given ML tag (+ alternatives) in contents of INCAR file.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer arrays:
!> ~~~
!>     #define 
!>     #define TYPE_REAL
!>     #define 0 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef 0
!>     #undef TYPE_REAL
!>     #undef 
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!==============================================================================
!  INTEGER FUNCTION READ_TAG_...       &
                            (INCAR,    &
                             TAG_LIST, &
                             PREFIX,   &
                             TAG,      &
                             VARIABLE, &
                             DEFAULT_VALUE) RESULT(NARGS)
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE ml_ff_taglist
      IMPLICIT NONE

!> String representation of INCAR file.
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
!> List with all tag states.
      TYPE(TAGLIST),    INTENT(INOUT) :: TAG_LIST
!> INCAR tag prefix (in case of multiple ML instances).
      CHARACTER(len=*), INTENT(IN)    :: PREFIX
!> INCAR tag to be processed.
      CHARACTER(len=*), INTENT(IN)    :: TAG
!> Variable where to store default or INCAR tag argument.
      CHARACTER (len=:), ALLOCATABLE,          INTENT(INOUT) :: VARIABLE 
!> Default value for variable.
      CHARACTER (len=:), ALLOCATABLE,          INTENT(IN)    :: DEFAULT_VALUE

! Alternative tags
      CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
      INTEGER                         :: I

! Set default value.
      VARIABLE = DEFAULT_VALUE
      CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

! Collect alternative tag list.
      ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

      NARGS = 0
! Loop over alternative tags (includes current "preferred" tag as
! first entry).
      DO I = 1, SIZE(ALT_TAGS)
! Check for per-instance INCAR tag first.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR,                                &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF

! Now comes the same check without prefix.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))), VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF
      END DO

! If found anywhere, set tag list state accordingly.
      IF (NARGS > 0) THEN
         IF (I > 1) THEN
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
         ELSE
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
         END IF
      END IF

!  END FUNCTION READ_TAG_...



# 1012 "ml_reader.F" 2 
   END FUNCTION READ_TAG_CHARACTER0





!> Read INCAR tag with 1-dimensional real array.
   INTEGER FUNCTION READ_TAG_REAL1 &
# 1 "./ml_reader.inc" 1 
! Fortran "template" workaround from incar_reader.inc.
# 4


# 10


# 16


# 20



!==============================================================================
!> Search and assign given ML tag (+ alternatives) in contents of INCAR file.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer arrays:
!> ~~~
!>     #define 
!>     #define 
!>     #define 1 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef 1
!>     #undef 
!>     #undef 
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!==============================================================================
!  INTEGER FUNCTION READ_TAG_...       &
                            (INCAR,    &
                             TAG_LIST, &
                             PREFIX,   &
                             TAG,      &
                             VARIABLE, &
                             DEFAULT_VALUE) RESULT(NARGS)
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE ml_ff_taglist
      IMPLICIT NONE

!> String representation of INCAR file.
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
!> List with all tag states.
      TYPE(TAGLIST),    INTENT(INOUT) :: TAG_LIST
!> INCAR tag prefix (in case of multiple ML instances).
      CHARACTER(len=*), INTENT(IN)    :: PREFIX
!> INCAR tag to be processed.
      CHARACTER(len=*), INTENT(IN)    :: TAG
!> Variable where to store default or INCAR tag argument.
      REAL(q),          INTENT(INOUT) :: VARIABLE (:)
!> Default value for variable.
      REAL(q),          INTENT(IN)    :: DEFAULT_VALUE

! Alternative tags
      CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
      INTEGER                         :: I

! Set default value.
      VARIABLE = DEFAULT_VALUE
      CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

! Collect alternative tag list.
      ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

      NARGS = 0
! Loop over alternative tags (includes current "preferred" tag as
! first entry).
      DO I = 1, SIZE(ALT_TAGS)
! Check for per-instance INCAR tag first.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR,                                &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF

! Now comes the same check without prefix.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))), VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF
      END DO

! If found anywhere, set tag list state accordingly.
      IF (NARGS > 0) THEN
         IF (I > 1) THEN
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
         ELSE
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
         END IF
      END IF

!  END FUNCTION READ_TAG_...



# 1021 "ml_reader.F" 2 
   END FUNCTION READ_TAG_REAL1





!> Read INCAR tag with 1-dimensional integer array.
   INTEGER FUNCTION READ_TAG_INTEGER1 &
# 1 "./ml_reader.inc" 1 
! Fortran "template" workaround from incar_reader.inc.
# 4


# 10




# 20



!==============================================================================
!> Search and assign given ML tag (+ alternatives) in contents of INCAR file.
!>
!> Example usage to get a READ_TAG function for 1-dimensional integer arrays:
!> ~~~
!>     #define 
!>     #define TYPE_REAL
!>     #define 1 1
!>        INTEGER FUNCTION READ_TAG_REAL1 &
!>     #include "ml_interface.inc"
!>        END FUNCTION READ_TAG_REAL1
!>     #undef 1
!>     #undef TYPE_REAL
!>     #undef 
!> ~~~
!> @attention Don't forget the '&' at the end of the function definition!
!>
!> @returns NARGS ... Number of elements found for the given tag.
!==============================================================================
!  INTEGER FUNCTION READ_TAG_...       &
                            (INCAR,    &
                             TAG_LIST, &
                             PREFIX,   &
                             TAG,      &
                             VARIABLE, &
                             DEFAULT_VALUE) RESULT(NARGS)
      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
      USE ml_ff_taglist
      IMPLICIT NONE

!> String representation of INCAR file.
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
!> List with all tag states.
      TYPE(TAGLIST),    INTENT(INOUT) :: TAG_LIST
!> INCAR tag prefix (in case of multiple ML instances).
      CHARACTER(len=*), INTENT(IN)    :: PREFIX
!> INCAR tag to be processed.
      CHARACTER(len=*), INTENT(IN)    :: TAG
!> Variable where to store default or INCAR tag argument.
      INTEGER,          INTENT(INOUT) :: VARIABLE (:)
!> Default value for variable.
      INTEGER,          INTENT(IN)    :: DEFAULT_VALUE

! Alternative tags
      CHARACTER(len=:), ALLOCATABLE   :: ALT_TAGS(:)
      INTEGER                         :: I

! Set default value.
      VARIABLE = DEFAULT_VALUE
      CALL TAG_LIST%SET_STATE(TAG, TL_DEFAULT)

! Collect alternative tag list.
      ALT_TAGS = TAGLIST_GET_ALT_TAGS(TAG)

      NARGS = 0
! Loop over alternative tags (includes current "preferred" tag as
! first entry).
      DO I = 1, SIZE(ALT_TAGS)
! Check for per-instance INCAR tag first.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR,                                &
                               ADJUSTL(TRIM(PREFIX // ALT_TAGS(I))), &
                               VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF

! Now comes the same check without prefix.
         NARGS = COUNT_ELEMENTS(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))))
         IF (NARGS > 0) THEN
            CALL PROCESS_INCAR(INCAR, ADJUSTL(TRIM(ALT_TAGS(I))), VARIABLE)
            EXIT ! Skip alternative tags if (1._q,0._q) is already found.
         END IF
      END DO

! If found anywhere, set tag list state accordingly.
      IF (NARGS > 0) THEN
         IF (I > 1) THEN
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR_ALT, I)
         ELSE
            CALL TAG_LIST%SET_STATE(TAG, TL_INCAR, I)
         END IF
      END IF

!  END FUNCTION READ_TAG_...



# 1030 "ml_reader.F" 2 
   END FUNCTION READ_TAG_INTEGER1






! defined(1) || defined(MPI_CHAIN)

! 

END MODULE ml_reader
