# 1 "stepver.F"
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


# 2 "stepver.F" 2 
!*********************************************************************
! RCS:  $Id: stepver.F,v 1.2 2002/04/16 07:28:52 kresse Exp $
!
!>  Subroutine STEP performs a verlet step on the ions.
!>
!>  This part of the program uses slightly different reduced units:
!>
!>  particules coordinates are in direct lattice   range [0,1] \n
!>  unit time     UT  =  POTIM * 1E-15 \n
!>    in reduced units (1._q,0._q) timestep for the ion movents takes 1 UT \n
!>  unit energy   UE  = 1 eV
!>
!>  These parameters must be supplied by the main porgram \n
!>  NIONS    number of IONS \n
!>  A        direct lattice \n
!>  TS       soll-temperature \n
!>  POTIM    length of timestep \n
!>  X        on entry: current positions  (time step n) \n
!>           on exit:  updated positions  (time step n+1) \n
!>  V        on entry: velocities of ions in r.u.
!>           these are the velocities at the time step n-1/2
!>           (on first Acall initial velocities) \n
!>           on exit: updated velocities at time step  n+1/2 \n
!>  D2C      D2C must contain half the accelations of the ions in r.u.
!>           evaluated at the current positions (time step n) \n
!>  SMASS    mass Parameter for nose dynamic, it has the dimension
!>           of a Energy*time**2 = mass*length**2 and is supplied in
!>           atomic mass * ULX**2 (this makes it  easy to define
!>           the parameter) \n
!>         <0     microcanonical ensemble
!>                if ISCALE is 1 scale velocities to TS \n
!>         -2     dont change velocities at all
!>
!>  Except for initialisation only D2C must be recalculated by
!>  the main program.
!>
!>  Following varibles are used by step, and should not be changed \n
!>  S        nose-parameter \n
!>  D2       used as temporary work array to calculate consistent velocity. \n
!>  D3C unused
!>
!>  Following Quantities are returned \n
!>  XC       position are equal to those at the previous timestep n \n
!>  X        updated position of ions (see above) \n
!>  EKIN     kinetic energy of system evaluated after C-step \n
!>  EPS      kinetic energy of fictivious Nose-parameter after C-step \n
!>  ES       potential energy of fictivious Nose-parameter after C-step \n
!>  DISMAX   maximum distance the ions have moved in r.u.
!>
!> The standard verlet is given by : x(n+1)=x(n)+(x(n)-x(n-1))+a(n) \n
!> This is equivalent to the leap-frog verlet \n
!>           v(n+1/2)=v(n-1/2)+a(n)
!>           x(n+1)=x(n)+v(n+1/2)
!
!*********************************************************************

      SUBROUTINE STEP(INIT,ISCALE,NIONS,A,ANORM,D2C,SMASS,X,XC, &
     &      POTIM,POMASS,NTYP,ITYP,TS,V,D2,D3,S, &
     &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IU6)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION  X(3,NIONS),XC(3,NIONS), &
     &           V(3,NIONS),D2(3,NIONS),D3(3,NIONS)
      DIMENSION  ITYP(NIONS),POMASS(NTYP)
      DIMENSION  S(4)
      DIMENSION  A(3,3),ANORM(3)
!=======================================================================
!   Some important Parameters, to convert to UL,UT,UE
!   BOLK is the Bolzmannk.,BOLKEV is the Bolzmannk. in EV
!=======================================================================
      PARAMETER (EVTOJ=1.60217733E-19_q,AMTOKG=1.6605402E-27_q, &
     &           BOLK=1.380658E-23_q,BOLKEV=8.6173857E-5_q)

!=======================================================================
!   D2C ..are the calculated vectors for the accelerations
!=======================================================================
      DIMENSION  D2C(3,NIONS)

!=======================================================================
!   set unit length,unit time ...
!=======================================================================
      UL =1E-10_q*ANORM(1)
      UT =POTIM*1E-15_q
!=======================================================================
!   convert SMASS to more usefull units
!   SMASS is in         (atomic mass* UL**2)
!   what we want now is (eV*UT**2)
!
!=======================================================================
      FACT=(AMTOKG/EVTOJ)*(UL/UT)**2

      IF (SMASS>0) THEN
        SQQ=SMASS *FACT
      ELSE
        SQQ=1.0_q
        IF (SMASS==-2) THEN
! SMASS=-2 set D2C and D2 and D3 to 0 (simulate no accelerations at all)
        DO N=1,3
        DO  I=1,NIONS
          D2C(N,I)=0
          D2(N,I) =0
          D3(N,NIONS)=0
        ENDDO
        ENDDO
      ENDIF
      ENDIF

!=======================================================================
!   INIT=0 initialize variables and perform first step
!   copy D2C to     D2
!        X   to     XC
!   set  D3  to     0
!   also initialize S
!=======================================================================
      IF (INIT==0) THEN
        INIT=1
!-------use kinetic Energy to calculate initial values for S
        CALL EKINC(EKIN,NIONS,NTYP,ITYP,POMASS,POTIM,A,V)
        S(1)=1._q
        S(2)=0._q
        S(3)=0._q
        IF (SMASS>0) S(3)=(EKIN-NDEGREES_OF_FREEDOM*BOLKEV*TS/2)*S(1)/SQQ
        S(4)=0._q
      ENDIF
!=======================================================================
!   calculate velocitites kinetic energy and temperature at
!   current timestep
!   D2= v(current timestep) = v(n)+a(v)/2
!   this requires an iterative refinement, because a(v) depends on the
!   v(current timestep) in the Nose formulation
!=======================================================================
      DO N=1,3
      DO I=1,NIONS
         D2(N,I)=V(N,I)
      ENDDO
      ENDDO
      SVEL=S(2)

      DO 210 ITER=0,10
      IF (SMASS>0) THEN
        SC=SVEL/S(1)
      ELSE
        SC=0
      ENDIF
      CALL EKINC(EKIN,NIONS,NTYP,ITYP,POMASS,POTIM,A,D2)

      ERROR=0
      DO N=1,3
      DO I=1,NIONS
         ERROR=ERROR+(V(N,I)+D2C(N,I)-SC*0.5_q*D2(N,I)-D2(N,I))**2
         D2(N,I)=V(N,I)+D2C(N,I)-SC*0.5_q*D2(N,I)
      ENDDO
      ENDDO

      IF (SMASS>0) THEN
      ERROR=ERROR+ &
     &  (S(2)+S(1)*((EKIN-NDEGREES_OF_FREEDOM*BOLKEV*TS/2)/SQQ+0.5_q*SC**2)-SVEL)**2
      SVEL=S(2)+S(1)*((EKIN-NDEGREES_OF_FREEDOM*BOLKEV*TS/2)/SQQ+0.5_q*SC**2)
      ENDIF
      IF (SQRT(ERROR)<1E-10_q) GOTO 220
  210 CONTINUE
  220 CONTINUE
!=======================================================================
!   scale velocities if requested
!=======================================================================

      IF (ISCALE==1) THEN
      SCALE=SQRT(TS/(2*EKIN/(BOLKEV*NDEGREES_OF_FREEDOM)))
      IF (IU6>=0) &
      WRITE(IU6,11) SCALE
   11 FORMAT('scaling velocities: factor=',F10.4)

      DO I =1,3
      DO NI=1,NIONS
         V(I,NI) =V(I,NI)*SCALE
         D2(I,NI)=D2(I,NI)*SCALE
      ENDDO
      ENDDO

      CALL EKINC(EKIN,NIONS,NTYP,ITYP,POMASS,POTIM,A,D2)

      ENDIF
!*********************************************************************
!   final verlet step
!*********************************************************************
      IF (SMASS>0) THEN
        SC=SVEL/S(1)
        EPS=0.5_q*SQQ*(SVEL/S(1))**2
        ES =NDEGREES_OF_FREEDOM*BOLKEV*TS*LOG(S(1))
      ELSE
        SC=0
        EPS=0
        ES=0
      ENDIF

      DO N=1,3
      DO I=1,NIONS
         XC(N,I)=X(N,I)
         V(N,I) =V(N,I)+2*D2C(N,I)-SC*D2(N,I)
         X(N,I) =X(N,I)+V(N,I)
         IF(X(N,I)>=1._q) X(N,I)=X(N,I)-1._q
         IF(X(N,I)<0._q) X(N,I)=X(N,I)+1._q
      ENDDO
      ENDDO
      IF (SMASS>0) THEN
      S(2)=S(2)+S(1)*((2*EKIN-BOLKEV*TS*NDEGREES_OF_FREEDOM)/SQQ+SC**2)
      S(1)=S(1)+S(2)
      ENDIF
!=======================================================================
!   calculate the maximum distance moved by the ions
!=======================================================================

      DISMAX=0
      DO 120 I=1,NIONS
         DISMAX=MAX(DISMAX, &
     &        V(1,I)**2+ &
     &        V(2,I)**2+ &
     &        V(3,I)**2)
  120 CONTINUE

      DISMAX=SQRT(DISMAX)

!=======================================================================
!   now the routine drops back to MAIN
!   MAIN evaluates the accelerations D2C at the predictor
!   coordinates X
!   (these accelerations will serve as corrector-derivatives of second
!   order).
!=======================================================================

      RETURN
      END
