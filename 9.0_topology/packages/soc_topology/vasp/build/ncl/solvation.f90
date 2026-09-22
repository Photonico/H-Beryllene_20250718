# 1 "solvation.F"
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


# 2 "solvation.F" 2 
!**********************************************************************************
!**********************************************************************************
! This file is only a placeholder for the actual solvation.F file:
! its public procedures are only stubs and the public variables are set so
! as not to interfere with the workings of VASP, but allows you to compile
! VASP with the hooks to the solvation model code VASPsol in place.
!
! The solvation.F file with the actual VASPsol code can be obtained at:
!
!  http://vaspsol.mse.cornell.edu/
!
!**********************************************************************************
!**********************************************************************************


!******************** MODULE SOLVATION ********************************************
!
!
!> Interfaces the solvation engine with the rest of vasp.
!
!
!**********************************************************************************
MODULE solvation

  USE prec

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SOL_READER, SOL_WRITER, XML_WRITE_SOL, SOL_VCORRECTION

  LOGICAL, SAVE :: LSOL=.FALSE.

  REAL(q), PUBLIC, SAVE :: Ediel_SOL=0._q
  REAL(q), PUBLIC, ALLOCATABLE, SAVE :: EIFOR_SOL(:,:)

CONTAINS

!******************** SUBROUTINE SOL_READER ***************************************
!
!
!> Reads in the solvation model parameters.
!
!
!**********************************************************************************
  SUBROUTINE SOL_READER(NIONS,EDIFF,IO)
    USE base
    IMPLICIT NONE

    TYPE (in_struct), INTENT(in) :: IO
    REAL(q), INTENT(in) :: EDIFF
    INTEGER, INTENT(in) :: NIONS
   
! this has to be 1._q ALWAYS
    IF (ALLOCATED(EIFOR_SOL)) DEALLOCATE(EIFOR_SOL)
    ALLOCATE(EIFOR_SOL(3,NIONS))

    EIFOR_SOL=0._q; Ediel_SOL=0._q
     
    RETURN
  END SUBROUTINE SOL_READER


!******************** SUBROUTINE SOL_WRITER ***************************************
!
!
!> Writes the solvation model parameters to the OUTCAR file.
!
!
!**********************************************************************************
  SUBROUTINE SOL_WRITER(IO)
    USE base
    TYPE (in_struct), INTENT(in) :: IO
    RETURN
  END SUBROUTINE SOL_WRITER


!******************** SUBROUTINE XML_WRITE_SOL ************************************
!
!
!> Writes the solvation model parameters to vasprun.xml.
!
!
!**********************************************************************************
  SUBROUTINE XML_WRITE_SOL
    RETURN
  END SUBROUTINE XML_WRITE_SOL


!******************** SUBROUTINE SOL_VCORRECTION *********************************
!
!
!> Computes the potential, energy and force corrections due to solvation.
!
!
!********************************************************************************
  SUBROUTINE SOL_VCORRECTION(INFO, T_INFO, LATT_CUR, P, WDES, GRIDC, CHTOT, CVTOT)
    USE base
    USE poscar
    USE lattice
    USE pseudo
    USE mgrid
    USE wave
    USE mdipol

    TYPE (info_struct), INTENT(in) :: INFO
    TYPE (type_info), INTENT(in) :: T_INFO
    TYPE (latt), INTENT(IN) :: LATT_CUR
    TYPE (potcar), INTENT(IN) :: P(T_INFO%NTYP)
    TYPE (wavedes), INTENT(IN) :: WDES
    TYPE (grid_3d), INTENT(IN) :: GRIDC
    
    COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ)
    COMPLEX(q) CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
    
    RETURN
  END SUBROUTINE SOL_VCORRECTION

END MODULE solvation
