# 1 "lattice.F"
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


# 2 "lattice.F" 2 
      MODULE LATTICE
      USE prec
      USE poscar_struct_def, ONLY: latt

      CONTAINS

!**************** SUBROUTINE LATTIC  ***********************************
! RCS:  $Id: lattice.F,v 1.2 2001/04/05 10:34:10 kresse Exp $
!
!>  subroutine for calculating the reciprocal lattice from the direct
!>  lattice
!>
!>  in addition the norm of the lattice-vectors and the volume of
!>  the basis-cell is calculated
!***********************************************************************

      SUBROUTINE LATTIC(Mylatt)
      USE prec

      IMPLICIT NONE

      TYPE(LATT) Mylatt
      REAL(q) Omega
      INTEGER I,J
      INTRINSIC SUM

      CALL EXPRO(Mylatt%B(1:3,1),Mylatt%A(1:3,2),Mylatt%A(1:3,3))
      CALL EXPRO(Mylatt%B(1:3,2),Mylatt%A(1:3,3),Mylatt%A(1:3,1))
      CALL EXPRO(Mylatt%B(1:3,3),Mylatt%A(1:3,1),Mylatt%A(1:3,2))

      Omega =Mylatt%B(1,1)*Mylatt%A(1,1)+Mylatt%B(2,1)*Mylatt%A(2,1) &
     &      +Mylatt%B(3,1)*Mylatt%A(3,1)

      DO I=1,3
      DO J=1,3
        Mylatt%B(I,J)=Mylatt%B(I,J)/Omega
      ENDDO
      ENDDO

      DO I=1,3
        Mylatt%ANORM(I)=SQRT(SUM(Mylatt%A(:,I)*Mylatt%A(:,I)))
        Mylatt%BNORM(I)=SQRT(SUM(Mylatt%B(:,I)*Mylatt%B(:,I)))
      ENDDO
      Mylatt%Omega=Omega
      RETURN
      END SUBROUTINE

      SUBROUTINE LATOLD(OMEGA,A,B)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION A(3,3),B(3,3)

      CALL EXPRO(B(1,1),A(1,2),A(1,3))
      CALL EXPRO(B(1,2),A(1,3),A(1,1))
      CALL EXPRO(B(1,3),A(1,1),A(1,2))

      OMEGA =B(1,1)*A(1,1)+B(2,1)*A(2,1)+B(3,1)*A(3,1)

      DO I=1,3
      DO J=1,3
        B(I,J)=B(I,J)/OMEGA
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE


!**************** FUNCTION LATTCHK *************************************
!
!> checks whether the crystalline and reciprocal lattices are
!> consistently classified by LATTYP
!
!***********************************************************************

      FUNCTION LATTCHK(LATT_CUR,LCONTINUE) RESULT(RES)
      USE prec
      USE tutor, ONLY : vtutor
      IMPLICIT NONE
      TYPE (latt) :: LATT_CUR
      LOGICAL, OPTIONAL :: LCONTINUE
      INTEGER :: RES(3)
! local
      INTEGER :: IBRAVA,IBRAVB,IBRAVA2B
      REAL(q) :: VTMP(3,3),CELDIM(6)
      LOGICAL :: LCONT
      CHARACTER(LEN=43) :: STRIBAV

      LCONT=.FALSE. ; IF (PRESENT(LCONTINUE)) LCONT=LCONTINUE

      VTMP=LATT_CUR%A; CALL LATTYP(VTMP(1,1),VTMP(1,2),VTMP(1,3),IBRAVA,CELDIM,-1)
      VTMP=LATT_CUR%B; CALL LATTYP(VTMP(1,1),VTMP(1,2),VTMP(1,3),IBRAVB,CELDIM,-1)

      IF (IBRAVB==IBRAVA2B(IBRAVA) .OR. LCONT) THEN
         RES(1)=IBRAVA ; RES(2)=IBRAVB ; RES(3)=IBRAVA2B(IBRAVA)
      ELSE
         CALL vtutor%alert( &
            'Inconsistent Bravais lattice types found for crystalline and reciprocal lattice:\n\n&
            &   Crystalline: ' // TRIM(STRIBAV(IBRAVA)) // '\n' // &
            '   Reciprocal : ' // TRIM(STRIBAV(IBRAVB)) // '\n' // &
            '               (instead of '// TRIM(STRIBAV(IBRAVA2B(IBRAVA))) // ')\n\n&
            &In most cases this is due to inaccuracies in the specification of the crytalline lattice vectors.\n\n&
            &Suggested SOLUTIONS:\n&
            & ) Refine the lattice parameters of your structure,\n&
            & ) and/or try changing SYMPREC.' )
      ENDIF

      RETURN
      END FUNCTION LATTCHK


!**************** SUBROUTINE EXPRO   ***********************************
! EXPRO
!> calculates the x-product of two vectors
!
!***********************************************************************

      SUBROUTINE EXPRO(H,U1,U2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION H(3),U1(3),U2(3)

      H(1)=U1(2)*U2(3)-U1(3)*U2(2)
      H(2)=U1(3)*U2(1)-U1(1)*U2(3)
      H(3)=U1(1)*U2(2)-U1(2)*U2(1)

      RETURN
      END SUBROUTINE


!**************** SUBROUTINE KARDIR ************************************
!> transform a set of vectors from cartesian coordinates to
!> ) direct lattice      (BASIS must be equal to B reciprocal lattice)
!> ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************

      SUBROUTINE KARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION V(3,NMAX),BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE


!**************** SUBROUTINE DIRKAR ************************************
!> transform a set of vectors from
!> ) direct lattice      (BASIS must be equal to A direct lattice)
!> ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
!> to cartesian coordinates
!***********************************************************************

      SUBROUTINE DIRKAR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION V(3,NMAX),BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(1,2)+V(3,N)*BASIS(1,3)
        V2=V(1,N)*BASIS(2,1)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(2,3)
        V3=V(1,N)*BASIS(3,1)+V(2,N)*BASIS(3,2)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE

!**************** SUBROUTINE TOPRIM ************************************
!> bring all ions into the primitive cell
!***********************************************************************

      SUBROUTINE TOPRIM(NIONS,POSION)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION POSION(3,NIONS)

      DO I=1,NIONS
      POSION(1,I)=MOD(POSION(1,I)+60,1._q)
      POSION(2,I)=MOD(POSION(2,I)+60,1._q)
      POSION(3,I)=MOD(POSION(3,I)+60,1._q)
      ENDDO
      END SUBROUTINE

      END MODULE


!**************** SUBROUTINE KARDIR ************************************
!> transform a set of vectors from cartesian coordinates to
!> ) direct lattice      (BASIS must be equal to B reciprocal lattice)
!> ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************

      SUBROUTINE CKARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(3,NMAX),V1,V2,V3
      DIMENSION  BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE


!**************** SUBROUTINE DIRKAR ************************************
!> transform a set of vectors from
!> ) direct lattice      (BASIS must be equal to A direct lattice)
!> ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
!> to cartesian coordinates
!***********************************************************************

      SUBROUTINE CDIRKAR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(3,NMAX),V1,V2,V3
      DIMENSION  BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(1,2)+V(3,N)*BASIS(1,3)
        V2=V(1,N)*BASIS(2,1)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(2,3)
        V3=V(1,N)*BASIS(3,1)+V(2,N)*BASIS(3,2)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE

