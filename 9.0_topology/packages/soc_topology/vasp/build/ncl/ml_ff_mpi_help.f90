# 1 "ml_ff_mpi_help.F"

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


# 3 "ml_ff_mpi_help.F" 2 
!****************************************************************************************************
! Module for some mpi routines that are needed very early
!****************************************************************************************************

      MODULE MPI_HELP

        USE ML_FF_CONSTANT
        USE ML_FF_PREC
        IMPLICIT NONE

        INCLUDE "mpif.h"
        TYPE ML_MPI_PAR 
            INTEGER :: MPI_COMM = 0 !< MPI_Communicator
            INTEGER :: NCPU     = 0 !< total number of proc in this communicator
            INTEGER :: NODE_ME  = 0 !< node id starting from 1 ... NCPU
            INTEGER :: COLOR    = 0 !< color of process within (sub)communicator
        END TYPE ML_MPI_PAR


        INTERFACE M_BCAST
           MODULE PROCEDURE M_BCAST_I_SCALAR
           MODULE PROCEDURE M_BCAST_R_SCALAR
           MODULE PROCEDURE M_BCAST_L_SCALAR
           MODULE PROCEDURE M_BCAST_I_1D
           MODULE PROCEDURE M_BCAST_R_1D
           MODULE PROCEDURE M_BCAST_L_1D
           MODULE PROCEDURE M_BCAST_I_2D
           MODULE PROCEDURE M_BCAST_R_2D
           MODULE PROCEDURE M_BCAST_L_2D
           MODULE PROCEDURE M_BCAST_I_3D
           MODULE PROCEDURE M_BCAST_R_3D
           MODULE PROCEDURE M_BCAST_L_3D
        END INTERFACE M_BCAST

        CONTAINS

!******************************************************************************************
! 1 Broadcast routines
!******************************************************************************************
        SUBROUTINE M_BCAST_I_SCALAR(COMM,SCAL,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: SCAL
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 54

           CALL MPI_BCAST(SCAL,1,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
# 58


        END SUBROUTINE

        SUBROUTINE M_BCAST_R_SCALAR(COMM,SCAL,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: SCAL
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 74

           CALL MPI_BCAST(SCAL,1,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
# 78


        END SUBROUTINE

        SUBROUTINE M_BCAST_L_SCALAR(COMM,SCAL,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: SCAL
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 95

           CALL MPI_BCAST(SCAL,1,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
# 99


        END SUBROUTINE

        SUBROUTINE M_BCAST_I_1D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: VEC(:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 116

           CALL MPI_BCAST(VEC(1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
# 120


        END SUBROUTINE

        SUBROUTINE M_BCAST_R_1D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: VEC(:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 137

           CALL MPI_BCAST(VEC(1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
# 141


        END SUBROUTINE

        SUBROUTINE M_BCAST_L_1D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: VEC(:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 158

           CALL MPI_BCAST(VEC(1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
# 162


        END SUBROUTINE

        SUBROUTINE M_BCAST_I_2D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: VEC(:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 179

           CALL MPI_BCAST(VEC(1,1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
# 183


        END SUBROUTINE

        SUBROUTINE M_BCAST_R_2D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: VEC(:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 200

           CALL MPI_BCAST(VEC(1,1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
# 204


        END SUBROUTINE

        SUBROUTINE M_BCAST_L_2D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: VEC(:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 221

           CALL MPI_BCAST(VEC(1,1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
# 225


        END SUBROUTINE

        SUBROUTINE M_BCAST_I_3D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: VEC(:,:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 242

           CALL MPI_BCAST(VEC(1,1,1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
# 246


        END SUBROUTINE

        SUBROUTINE M_BCAST_R_3D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: VEC(:,:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 263

           CALL MPI_BCAST(VEC(1,1,1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
# 267


        END SUBROUTINE

        SUBROUTINE M_BCAST_L_3D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: VEC(:,:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
!Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
# 284

           CALL MPI_BCAST(VEC(1,1,1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
# 288


        END SUBROUTINE

      END MODULE MPI_HELP

