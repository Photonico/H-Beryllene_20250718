# 1 "ml_ff_constant.F"
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


# 2 "ml_ff_constant.F" 2 
!************************************************************************
!
!  this module contains some control data structures for
!  machine-learning force field.
!
!***********************************************************************
      MODULE ML_FF_CONSTANT
      USE ML_FF_PREC
! Constants used for machine learning force field

!  Some important Parameters, to convert to a.u.
!  - AUTOA  = 1. a.u. in Angstroem
!  - RYTOEV = 1 Ry in Ev
!  - EVTOJ  = 1 eV in Joule
!  - AMTOKG = 1 atomic mass unit ("proton mass") in kg

      REAL(q), PARAMETER :: AUTOA=0.529177249_q,RYTOEV=13.605826_q
      REAL(q), PARAMETER :: EVTOJ=1.60217733E-19_q,AMTOKG=1.6605402E-27_q
      REAL(q), PARAMETER :: PI =3.141592653589793238_q

! EUNIT = 1 Hartree in ev units
! MUNIT = Electron mass in atomic units
! RUNIT = Bohr Radius in Angstrom1
! FUNIT = Unit of force in Hartree/Bohr
! SUNIT = Unit of stress tensor
! TUNIT = Unit of time in atomic unit

      REAL(q), PARAMETER :: EUNIT = 2.0_q*RYTOEV, &
                            MUNIT = 1.054571726E-034_q**2/2.0_q/(RYTOEV*EVTOJ*AUTOA**2*1.0E-20_q)/AMTOKG, &
                            FUNIT = 2.0_q*RYTOEV/AUTOA, &
                            SUNIT = 2.0_q*RYTOEV/(AUTOA*1.0E-10_q)**3*EVTOJ/1.0E8_q, &
                            TUNIT = 1.054571726E-034_q/(2.0_q*RYTOEV*EVTOJ)*1.0E15_q

! Parameters for Normalization (to avoid it if norm becomes (0._q,0._q))

      REAL(q), PARAMETER :: TOLERANCE_NORM = 1.0E-08_q

      INTEGER, PARAMETER :: NBLOCK_SCALAPACK = 32
!> Length in bytes of the ML_FF ASCII header.
      INTEGER, PARAMETER :: ML_FF_HEADER_SIZE = 4096
!======================================================================
!> ML_LOGFILE version number.
!>
!> Version history:
!> * 0.1.0: Initial version.
!> * 0.2.0: Minor updates, changes include:
!>   - Additional STATUS category "predfast".
!>   - New tags output:
!>     ML_MODE, ML_LFAST, ML_LSIC, ML_EPS_REG, ML_OUTPUT_MODE, ML_OUTBLOCK,
!>     ML_LSUPERVEC, ML_LBASIS_DISCARD.
!>   - New info lines NDESC and NDESC_SIC.
!> * 0.2.1: Minor updates, changes:
!>   - New tag output: ML_DESC_TYPE.
!>   - Maximum buffer size below tag ML_MB added.
!> * 0.2.2: Minor updates, changes:
!>   - New tag output: ML_MB_MIN.
!>   - New info lines: MSG.
!>   - Skipping per-step log lines in fast prediction mode.
!> * 0.2.3: Minor updates, changes:
!>   - New log line: SF(F).
!>   - New tag output: ML_IERR, ML_CALGO, ML_LIB
!======================================================================
      INTEGER, PARAMETER :: ML_LOGFILE_VERSION(3) = [0, 2, 3]
!======================================================================
!> ML_FF version number.
!>
!> Version history:
!> * 0.0.1: Original ASCII file.
!> * 0.1.0: Binary version of previous version.
!> * 0.2.0: Added ASCII header to binary file.
!> * 0.2.1: Added ML_DESC_TYPE=1 output.
!> * 0.2.2: Added ML_DESC_TYPE=2,3 output.
!> * 0.2.3: Added inverse kernel matrix SINV.
!> * 0.2.4: Added LFAST field in ML_FF file.
!======================================================================
      INTEGER, PARAMETER :: ML_FF_VERSION_WRITE(3) = [0, 2, 4]

      END MODULE
