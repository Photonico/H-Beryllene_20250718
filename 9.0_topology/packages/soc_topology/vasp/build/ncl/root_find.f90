# 1 "root_find.F"
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


# 2 "root_find.F" 2 

!> @brief General routines for finding the zeros of a function
!
!> All the arguments necessary to the function are passed using an array of C pointers.
!> A wrapper function must defined to transform the C pointers into Fortran ones,
!> call the function and return. An arbitrary large number of arguments is possible,
!> and they can be of any type.
!> An example of how these routines can be used is present in lcao_bare.F
MODULE root_find
    USE prec
    USE tutor, ONLY: vtutor
    USE iso_c_binding
    IMPLICIT NONE

    abstract interface
        function func_template(x,nargs,args) RESULT(FX)
        use prec
        use iso_c_binding
        REAL(q),INTENT(IN) :: X
        INTEGER,INTENT(IN) :: NARGS
        TYPE(C_PTR),VALUE,INTENT(IN) :: ARGS
        REAL(q) :: FX
        end function func_template
    end interface


    PUBLIC :: ZBRAC ! Bracket outwards
    PUBLIC :: ZBRAK ! Bracket inwards
    PUBLIC :: BISECT ! Find root using Bissection method
    PUBLIC :: ZBRENT ! Find root using Brent's method
   
    CONTAINS

!>@brief find the bracket to search the root
!> Given a function func and an initial guessed range x1 to x2,
!> the routine expands the range geometrically until a root is bracketed
!> by the returned values x1 and x2 (in which case succes returns as .true.)
!> or until the range becomes unacceptably large (in which case succes returns as .false.).
!> Taken from "Numerical recipes in Fortran90"
  SUBROUTINE ZBRAC(FUNC,NARGS,ARGS,X1,X2,SUCCESS)
    PROCEDURE(func_template) :: FUNC    !< function as input
    INTEGER, INTENT(IN) :: NARGS !< number of arguments
    TYPE(C_PTR), VALUE, INTENT(IN) :: ARGS !< array of c pointers containing the arguments to F
    REAL(q), INTENT(INOUT)  :: X1,X2 !< limits in which to find zeros
    LOGICAL, INTENT(OUT) :: SUCCESS !< whether a (0._q,0._q) was found
! local
    REAL(q),PARAMETER :: FACTOR = 1.01 !< Factor at which the interval is expanded
    INTEGER,PARAMETER :: NTRY=20  !< The maximum expansion of the interval is FACTOR**NTRY*ABS(X2-X1)
    INTEGER :: I
    REAL(q) :: F1,F2

    F1=FUNC(X1,NARGS,ARGS)
    F2=FUNC(X2,NARGS,ARGS)
    SUCCESS = .TRUE.
    DO I=1,NTRY
       IF(F1*F2<0) RETURN
       X1=X1-FACTOR*(X2-X1)
       F1=FUNC(X1,NARGS,ARGS)
       IF(F1*F2<0) RETURN
       X2=X2+FACTOR*(X2-X1)
       F2=FUNC(X2,NARGS,ARGS)
    ENDDO
    SUCCESS=.FALSE.
  END SUBROUTINE ZBRAC

!>@brief find the bracket to search the root
!> Given a function fx defined on the interval from x1-x2 subdivide the
!> interval into n equally spaced segments, and search for (0._q,0._q) crossings
!> of the function. nb is input as the maximum number of roots sought, and
!> is reset to the number of bracketing pairs xb1(1:nb),xb2(1:nb) that
!> are found.
!> Taken from "Numerical recipes in Fortran90"
  SUBROUTINE ZBRAK(FUNC,NARGS,ARGS,X1,X2,N,XB1,XB2,NB)
    PROCEDURE(func_template) :: FUNC    !< function as input
    INTEGER, INTENT(IN) :: NARGS !< number of arguments
    TYPE(C_PTR), VALUE, INTENT(IN) :: ARGS !< array of c pointers containing the arguments to F
    REAL(q), INTENT(IN)   :: X1,X2 !< limits in which to find zeros
    INTEGER, INTENT(IN)   :: N !< Number of divisions between X1 and X2
    INTEGER, INTENT(INOUT) :: NB !< Maximum number of roots on entrance, number of roots on exit
    REAL(q), INTENT(OUT)  :: XB1(NB),XB2(NB) !< Lower and upper bounds for for intervals bracketing zeros
! local
    INTEGER :: I,NBB
    REAL(q) :: DX,FC,FP,X

    NBB=0
    X=X1
    DX=(X2-X1)/N
    FP=FUNC(X,NARGS,ARGS)

    DO I=1,N
        X=X+DX
        FC=FUNC(X,NARGS,ARGS)
        IF(FC*FP<0) THEN
           NBB=NBB+1
           XB1(NBB)=X-DX
           XB2(NBB)=X
           IF(NBB==NB) EXIT
        ENDIF
        FP=FC
    ENDDO
    NB=NBB
  END SUBROUTINE ZBRAK

!> @brief Bissection method to find roots of a function in interval A,B
  SUBROUTINE BISECT(FUNC,NARGS,ARGS,A,B,TOL,C,FC)
    PROCEDURE(func_template) :: FUNC    !< function as input
    INTEGER, INTENT(IN) :: NARGS !< number of arguments
    TYPE(C_PTR), VALUE, INTENT(IN) :: ARGS !< array of c pointers containing the arguments to F
    REAL(q), INTENT(INOUT) :: A,B !< limits in which to find zeros
    REAL(q), INTENT(IN) :: TOL !< Tolerance to find the (0._q,0._q)
    REAL(q), INTENT(OUT) :: C !< Solution C
    REAL(q), INTENT(OUT) :: FC !< Value of FUNC(C)
! local variables
    INTEGER, PARAMETER :: MAX_ITER = 500
    REAL(q) :: FA,FB
    INTEGER :: ITER
    FA = FUNC(A,NARGS,ARGS)
    FB = FUNC(B,NARGS,ARGS)
! Start bissection method to find the energy that yields the correct logaritmic derivative
    DO ITER=1,MAX_ITER
       C = 0.5_q*(A+B)
!C=(A*FB-B*FA)/(FB-FA)
       FC = FUNC(C,NARGS,ARGS)
!WRITE(*,'(I3,7F18.8)') ITER,A,C,B,FA,FC,FB,0.5_q*(B-A)

       IF(FC==0.OR.0.5_q*(B-A)<TOL) EXIT

       IF(SIGN(1.0_q,FC)==SIGN(1.0_q,FA)) THEN
          A  = C
          FA = FC
       ELSE
          B  = C
          FB = FC
       ENDIF
    ENDDO
    IF (ITER>MAX_ITER) CALL vtutor%error('BISECT: exceeded maximum iterations')
  END SUBROUTINE BISECT

!> @brief Brent's (0._q,0._q) finding alorithm.
!
!> Based on the routine already present in pade_fit.F
!> Taken from "Numerical Recipes in Fortran90"
  SUBROUTINE ZBRENT(FUNC,NARGS,ARGS,X1,X2,TOL,C,FC)
    PROCEDURE(func_template) :: FUNC    !< function as input
    INTEGER, INTENT(IN) :: NARGS !< number of arguments
    TYPE(C_PTR), VALUE, INTENT(IN) :: ARGS !< array of c pointers containing the arguments to F
    REAL(q), INTENT(IN)   :: X1,X2,TOL
! local
    INTEGER, PARAMETER :: MAX_ITER = 500 !< Maximum number of iterations
    REAL(q), PARAMETER :: EPS=1E-14 !< Machine precision
    INTEGER(q) :: ITER
    REAL(q)    :: A,B,C,D,E,FA,FB,FC
    REAL(q)    :: P,Q1,R,S,TOL1,XM

    A=X1
    B=X2
    FA=FUNC(A,NARGS,ARGS)
    FB=FUNC(B,NARGS,ARGS)
    IF ((FA > 0.0 .AND. FB > 0.0) .OR. (FA < 0.0 .AND. FB < 0.0)) THEN
        CALL vtutor%error('ZBRENT: root must be bracketed for ZBRENT')
    ENDIF
    C=B
    FC=FB
    DO ITER=1,MAX_ITER
       IF ((FB > 0.0 .AND. FC > 0.0) .OR. (FB < 0.0 .AND. FC < 0.0)) THEN
          C=A
          FC=FA
          D=B-A
          E=D
       END IF
       IF (ABS(FC) < ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
       END IF
       TOL1=2.0_q*EPS*ABS(B)+0.5_q*TOL
       XM=0.5_q*(C-B)
       IF (ABS(XM) <= TOL1 .OR. FB == 0.0) THEN
          C=B
          FC=FB
          RETURN
       END IF
       IF (ABS(E) >= TOL1 .AND. ABS(FA) > ABS(FB)) THEN
          S=FB/FA
          IF (A == C) THEN
             P=2.0_q*XM*S
             Q1=1.0_q-S
          ELSE
             Q1=FA/FC
             R=FB/FC
             P=S*(2.0_q*XM*Q1*(Q1-R)-(B-A)*(R-1.0_q))
             Q1=(Q1-1.0_q)*(R-1.0_q)*(S-1.0_q)
          END IF
          IF (P > 0.0) Q1=-Q1
          P=ABS(P)
          IF (2.0_q*P < MIN(3.0_q*XM*Q1-ABS(TOL1*Q1),ABS(E*Q1))) THEN
             E=D
             D=P/Q1
          ELSE
             D=XM
             E=D
          END IF
       ELSE
          D=XM
          E=D
       END IF
       A=B
       FA=FB
       IF(ABS(D)>TOL1) THEN
          B=B+D
       ELSE
          B=B+SIGN(TOL1,XM)
       ENDIF
       FB=FUNC(B,NARGS,ARGS)
    END DO
    CALL vtutor%error('ZBRENT: exceeded maximum iterations')
  END SUBROUTINE ZBRENT

END MODULE root_find
