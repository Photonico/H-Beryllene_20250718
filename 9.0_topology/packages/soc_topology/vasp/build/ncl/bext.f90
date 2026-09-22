# 1 "bext.F"
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


# 2 "bext.F" 2 
      MODULE bexternal
      USE prec
      IMPLICIT NONE

      PUBLIC BEXT_READER,LBEXTERNAL,BEXT,BEXT_ADDV

      PRIVATE

      REAL(q), SAVE :: BEXT(3)=0

      LOGICAL, SAVE :: LBEXT=.FALSE.

      CONTAINS

!***********************************************************************
!******************** PUBLIC PROCEDURES ********************************
!***********************************************************************

!******************** SUBROUTINE BEXT_READER ***************************
!
!> Reads BEXT from the INCAR file
!
!***********************************************************************

      SUBROUTINE BEXT_READER(IU0,IU5,ISPIN,LNONCOLLINEAR)
      USE base
      USE vaspxml
      USE reader_tags
      USE tutor, ONLY: vtutor
      USE string, ONLY: str

      INTEGER       :: IU5,IU6,IU0
! local variables
      INTEGER       :: IDUM, N, IERR
      REAL(q)       :: RDUM
      COMPLEX(q)    :: CDUM
      LOGICAL       :: LOPEN,LDUM
      CHARACTER (1) :: CHARAC

      REAL(q)       :: BNORM2
      INTEGER       :: ISPIN,NREQ
      LOGICAL       :: LSORBIT,LNONCOLLINEAR

      NREQ=0
      IF (LNONCOLLINEAR) THEN
         NREQ=3
      ELSEIF (ISPIN==2) THEN
         NREQ=1
      ENDIF

      IF (NREQ/=0) THEN
! read in BEXT tag
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BEXT', BEXT(1:NREQ), NREQ, IERR, WRITEXMLINCAR, LCONTINUE=.TRUE., FOUNDNUMBER=N)
! error handling
         IF (N==NREQ) THEN
           LBEXT=.TRUE.

           BNORM2=BEXT(1)*BEXT(1)+BEXT(2)*BEXT(2)+BEXT(3)*BEXT(3)
           IF (BNORM2<1.E-8_q) THEN
             CALL vtutor%advice("BEXT is tiny. "//&
               "Are you sure it will make a difference?")
           ENDIF
         ELSEIF (N/=NREQ.AND.N/=0) THEN
           CALL vtutor%error("ERROR: You have set "//str(N)//&
             " value(s) for BEXT; however "//&
             "for LNONCOLLINEAR="//str(LNONCOLLINEAR)//" and ISPIN="//&
             str(ISPIN)//" BEXT takes "//str(NREQ)//" value(s)." )
         ENDIF
      ENDIF

      CALL CLOSE_INCAR_IF_FOUND(IU5)
      RETURN
      END SUBROUTINE BEXT_READER


!******************** FUNCTION LBEXTERNAL ******************************
!
!***********************************************************************

      FUNCTION LBEXTERNAL()
      LOGICAL LBEXTERNAL
      LBEXTERNAL=LBEXT
      END FUNCTION LBEXTERNAL


!******************** SUBROUTINE BEXT_ADDV *****************************
!
!***********************************************************************

      SUBROUTINE BEXT_ADDV(CVTOT,GRIDC,NCDIJ)
      USE mgrid
      
      TYPE (grid_3d) GRIDC

      COMPLEX(q)  :: CVTOT(GRIDC%MPLWV,NCDIJ)
      INTEGER     :: NCDIJ
! local variables
      INTEGER     :: I

      DO I=2,NCDIJ
         CALL ADD2VG0(CVTOT(1,I),GRIDC,CMPLX(BEXT(I-1),0._q,KIND=q))
      ENDDO

      RETURN
      END SUBROUTINE BEXT_ADDV

      END MODULE bexternal


!******************** SUBROUTINE ADD2VG0 *******************************
!
!***********************************************************************

      SUBROUTINE ADD2VG0(CVTOT,GRIDC,C)

!! USE moffload_struct_def

      USE mgrid
      TYPE (grid_3d) GRIDC
      COMPLEX(q)  :: CVTOT(GRIDC%RC%NROW,GRIDC%RC%NCOL)
      COMPLEX(q)  :: C
! local variables
      INTEGER N1,N2,N3,NC

      N1=1; N2=1; N3=1

!$ACC PARALLEL LOOP PRESENT(CVTOT) 
      DO NC=1,GRIDC%RC%NCOL
         IF (GRIDC%RC%I2(NC)==N2 .AND. GRIDC%RC%I3(NC)==N3) THEN
            CVTOT(NC,N1)=CVTOT(NC,N1)+C
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE ADD2VG0
