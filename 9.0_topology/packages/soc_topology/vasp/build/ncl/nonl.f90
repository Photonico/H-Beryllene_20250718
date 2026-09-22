# 1 "nonl.F"
!#define debug
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


# 3 "nonl.F" 2 
!***********************************************************************
!
!> This module contains all the routines required to support
!> reciprocal space presentation of the non local projection operators
!> on a plane wave grid
!>
!> @note
!> The wavefunction character is strictly defined as
!> ~~~
!>  C_n = e (i k R) \sum_G p_n,G+k e (i G R) C_G
!> ~~~
!> VASP however does not include the term <TT>e (i k R)</TT>, since
!> it does not change the occupancy matrix or energy. Hence
!> ~~~
!>  W%CPROJ(n) = \sum_G p_n,G+k e (i G R) C_G
!> ~~~
!> the derivative of C_n with respect to R can either include the
!> term from e (i k R)
!> ~~~
!>    i k \sum_G  p_n,G+k e (i G R) C_G = i k C_n
!> ~~~
!> or not.
!> This contribution does not contribute to the (1._q,0._q) center occupancy
!> matrix or the total energy either, since it is "complex"
!> (expectation values involving C_m* D_mn C_n are complex)
!> the Berry phase seems to require <TT>-Dshift_der_k</TT> however
!> (nonlr.F uses a similar convention)
!***********************************************************************


MODULE nonl

!! USE moffload

  USE prec
  USE wave
  USE nonl_struct_def

  INTERFACE
     SUBROUTINE VNLAC0(NONL_S,WDES1,CPROJ_LOC,CACC)
       USE nonl_struct_def
       USE wave_struct_def
       TYPE (nonl_struct) NONL_S
       TYPE (wavedes1)    WDES1
       COMPLEX(q)    CPROJ_LOC
       COMPLEX(q) CACC
     END SUBROUTINE VNLAC0
  END INTERFACE

  INTERFACE
     SUBROUTINE VNLAC0_DER(NONL_S,WDES1,CPROJ_LOC,CACC, LATT_CUR, ION, IDIR)
       USE nonl_struct_def
       USE wave_struct_def
       USE lattice
       TYPE (nonl_struct) NONL_S
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CACC
       COMPLEX(q)         CPROJ_LOC
       TYPE (latt)  LATT_CUR
       INTEGER      ION
       INTEGER      IDIR
     END SUBROUTINE VNLAC0_DER
  END INTERFACE

CONTAINS

!****************** SUBROUTINE NONL_ALLOC ******************************
! RCS:  $Id: nonl.F,v 1.2 2002/08/14 13:59:41 kresse Exp $
!
!> Allocates required arrays based on T_INFO and P structure,
!>
!> The number of ions and types is taken from T_INFO and
!> LMDIM is taken from the P structure
!
!***********************************************************************

  SUBROUTINE NONL_ALLOC(NONL_S,T_INFO,P,WDES, LREAL)
    USE prec
    USE pseudo
    USE poscar
    USE ini
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    LOGICAL LREAL

! local var
    INTEGER NIONS,NTYPD,LMDIM,NRPLWV,NKPTS,NT,NIS,NI

    NIONS =T_INFO%NIONS
    NTYPD =T_INFO%NTYPD
    LMDIM =P(1)%LMDIM
    NRPLWV=WDES%NGDIM
    NKPTS =WDES%NKPTS   ! number of k-points in the IBZ
! (required for HF type calculations)
    NONL_S%LRECIP=.NOT. LREAL
    NONL_S%NTYP  =T_INFO%NTYP
    NONL_S%NK    =0
    NONL_S%NIONS =T_INFO%NIONS
    NONL_S%NITYP =>T_INFO%NITYP
    NONL_S%ITYP  =>T_INFO%ITYP
    NONL_S%POSION=>T_INFO%POSION
    NONL_S%LSPIRAL=WDES%LSPIRAL
    NULLIFY(NONL_S%VKPT_SHIFT)

    ALLOCATE(NONL_S%LMMAX(NONL_S%NTYP),NONL_S%LMBASE(NONL_S%NIONS+1))

    NONL_S%LMBASE(1)=0; NIS=1
    DO NT=1,NONL_S%NTYP
       NONL_S%LMMAX(NT)=P(NT)%LMMAX
       DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
          NONL_S%LMBASE(NI+1)=NONL_S%LMBASE(NI)+NONL_S%LMMAX(NT)
       ENDDO
       NIS=NIS+NONL_S%NITYP(NT)
    ENDDO

    NULLIFY(NONL_S%CREXP,NONL_S%QPROJ,NONL_S%CQFAK)

    IF ((.NOT.LREAL).AND.(.NOT.NONL_S%LSPIRAL)) THEN
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
         NONL_S%CQFAK(LMDIM,NTYPD))
         CALL REGISTER_ALLOCATE(16._q*SIZE(NONL_S%CREXP)+8._q*SIZE(NONL_S%QPROJ), "nonl-proj")
    ENDIF
    IF ((.NOT.LREAL).AND. NONL_S%LSPIRAL) THEN
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
         NONL_S%CQFAK(LMDIM,NTYPD))
         CALL REGISTER_ALLOCATE(16._q*SIZE(NONL_S%CREXP)+8._q*SIZE(NONL_S%QPROJ), "nonl-proj")
    ENDIF

    RETURN
  END SUBROUTINE NONL_ALLOC


  SUBROUTINE NONL_DEALLOC(NONL_S)
    USE ini
    IMPLICIT NONE
    TYPE (nonl_struct) NONL_S
    LOGICAL LREAL

    LREAL=.NOT.NONL_S%LRECIP

    IF (.NOT.LREAL) THEN
       CALL DEREGISTER_ALLOCATE(16._q*SIZE(NONL_S%CREXP)+8._q*SIZE(NONL_S%QPROJ), "nonl-proj")
       DEALLOCATE(NONL_S%CREXP,NONL_S%QPROJ,NONL_S%CQFAK)
       DEALLOCATE(NONL_S%LMMAX)
       DEALLOCATE(NONL_S%LMBASE)
   ENDIF
       

    RETURN
  END SUBROUTINE NONL_DEALLOC

!****************** SUBROUTINE NONL_ALLOC_SPHPRO ***********************
!
!> allocates required arrays to describe projectors for (1._q,0._q) ion
!> and (1._q,0._q) k-point
!
!***********************************************************************

  SUBROUTINE NONL_ALLOC_SPHPRO(NONL_S,P,WDES)
    USE pseudo
    USE poscar
    IMPLICIT NONE


    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P
    TYPE (wavedes)     WDES
! local var
    INTEGER NIONS,NTYPD,LMDIM,NRPLWV,NKPTS,NT,NIS,NI

    NIONS =1
    NTYPD =1
    LMDIM =P%LMDIM
    NRPLWV=WDES%NGDIM
    NKPTS =WDES%NKPTS

    NONL_S%LRECIP=.TRUE.
    NONL_S%NTYP  =NTYPD
    NONL_S%NK    =0
    NONL_S%NIONS =NIONS
    NONL_S%LSPIRAL=WDES%LSPIRAL
    ALLOCATE(NONL_S%NITYP(NTYPD)); NONL_S%NITYP(1)=1
    ALLOCATE(NONL_S%ITYP(NTYPD)) ; NONL_S%ITYP(1)=1
    ALLOCATE(NONL_S%LMMAX(NTYPD)); NONL_S%LMMAX(1)=P%LMMAX

    ALLOCATE(NONL_S%LMBASE(NIONS+1))
    NONL_S%LMBASE(1)=0; NIS=1
    DO NT=1,NONL_S%NTYP
       DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
          NONL_S%LMBASE(NI+1)=NONL_S%LMBASE(NI)+NONL_S%LMMAX(NT)
       ENDDO
       NIS=NIS+NONL_S%NITYP(NT)
    ENDDO

    IF (.NOT.NONL_S%LSPIRAL) &
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
         NONL_S%CQFAK(LMDIM,NTYPD))
    IF (NONL_S%LSPIRAL) &
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
         NONL_S%CQFAK(LMDIM,NTYPD))
    RETURN
  END SUBROUTINE NONL_ALLOC_SPHPRO

!****************** SUBROUTINE NONL_DEALLOC_SPHPRO *********************
!
!> Deallocate required arrays to describe projectors for (1._q,0._q) ion
!> and (1._q,0._q) k-point
!
!***********************************************************************

  SUBROUTINE NONL_DEALLOC_SPHPRO(NONL_S)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    DEALLOCATE(NONL_S%CREXP,NONL_S%QPROJ,NONL_S%CQFAK)
    DEALLOCATE(NONL_S%NITYP,NONL_S%ITYP)
    DEALLOCATE(NONL_S%LMMAX,NONL_S%LMBASE)

    RETURN
  END SUBROUTINE NONL_DEALLOC_SPHPRO


!****************** SUBROUTINE SPHER *********************************
!
!>  Calculates the sperical harmonics multiplied
!>  by the pseudopotential QPROJ on the compressed reciprocal
!>  lattice grid
!>
!>  The full non local pseudopotential is given by
!> ~~~
!>  <G|V(K)|GP> = SUM(L,LP,site,R)
!>                     L                                    LP
!>    QPROJ(G,L site) i * DIJ(L,LP site) QPROJ(GP,LP site) i *
!>                        EXP(-i(G-GP) R)
!> ~~~
!>  i^L is stored in NONL_S\%CQFAK (nonl_struct::CQFAK)
!>
!> @param IZERO = -1   set NONL_S to the negative projectors\n
!>              =  1   set NONL_S to the positive projectors\n
!>              =  0   add to NONL_S the positive projectors
!
!*********************************************************************

  SUBROUTINE SPHER(GRID, NONL_S, P, WDES, LATT_CUR, IZERO, NK_SELECT, DK)

    USE pseudo
    USE mpimy
    USE mgrid
    USE lattice
    USE constant
    USE asa
    USE string, ONLY: str
    USE tutor, ONLY: vtutor

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    INTEGER  IZERO
    INTEGER, OPTIONAL :: NK_SELECT
    REAL(q), OPTIONAL :: DK(3)

! work arrays
    REAL(q),ALLOCATABLE :: GLEN(:),XS(:),YS(:),ZS(:),VPS(:),FAKTX(:),YLM(:,:),GLENP(:)

    

!! 

    LYDIM=MAXL(NONL_S%NTYP,P)
    LMYDIM=(LYDIM+1)**2          ! number of lm pairs

    NA=WDES%NGDIM

    ALLOCATE(GLEN(NA),XS(NA),YS(NA),ZS(NA),VPS(NA),FAKTX(NA),YLM(NA,LMYDIM))

    IF (PRESENT(DK)) ALLOCATE(GLENP(NA))

# 297

    IF (WDES%LGAMMA) THEN
       CALL vtutor%bug("internal error in SPHER: WDES%LGAMMA is incorrect " // str(WDES%LGAMMA), "nonl.F", 299)
    ENDIF


!-----------------------------------------------------------------------
! spin spiral propagation vector in cartesian coordinates
! is simply (0._q,0._q) when LSPIRAL=.FALSE.
!-----------------------------------------------------------------------
    QX= (WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3))/2
    QY= (WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3))/2
    QZ= (WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3))/2

!=======================================================================
! Loop over NSPINORS: here only in case of spin spirals NRSPINOR=2
!=======================================================================
    IF (WDES%LSPIRAL) THEN
       NSPINORS=2
    ELSE
       NSPINORS=1
    ENDIF
    spinor: DO ISPINOR=1,NSPINORS
!=======================================================================
! main loop over all special points
!=======================================================================
       kpoint: DO NK=1,WDES%NKPTS
          IF (PRESENT(NK_SELECT)) THEN
             IF ( NK_SELECT/=NK) CYCLE
          ENDIF
!         unfortunately this presently breaks the sc-GW code
!          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!=======================================================================
! now calculate the necessary tables:
! containing the length, phasefactor exp(i m phi) and
! sin(theta),cos(theta)
!=======================================================================
          DO IND=1,WDES%NGVECTOR(NK)

             N1=MOD(WDES%IGX(IND,NK)+GRID%NGX,GRID%NGX)+1
             N2=MOD(WDES%IGY(IND,NK)+GRID%NGY,GRID%NGY)+1
             N3=MOD(WDES%IGZ(IND,NK)+GRID%NGZ,GRID%NGZ)+1

             G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
             G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
             G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))

             IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
                G1= G1 +NONL_S%VKPT_SHIFT(1,NK)
                G2= G2 +NONL_S%VKPT_SHIFT(2,NK)
                G3= G3 +NONL_S%VKPT_SHIFT(3,NK)
             ENDIF

             FACTM=1._q
             IF (WDES%LGAMMA .AND. (N1/=1 .OR. N2/=1 .OR. N3/=1)) FACTM=SQRT(2._q)

             GX= (G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)-QX) *TPI
             GY= (G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)-QY) *TPI
             GZ= (G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)-QZ) *TPI


             GLEN(IND)=MAX(SQRT(GX*GX+GY*GY+GZ*GZ),1E-10_q)
             FAKTX(IND)=FACTM
             XS(IND)  =GX/GLEN(IND)
             YS(IND)  =GY/GLEN(IND)
             ZS(IND)  =GZ/GLEN(IND)

             IF (PRESENT(DK)) THEN
                G1P= G1 +DK(1)
                G2P= G2 +DK(2)
                G3P= G3 +DK(3)
                GXP= (G1P*LATT_CUR%B(1,1)+G2P*LATT_CUR%B(1,2)+G3P*LATT_CUR%B(1,3)-QX) *TPI
                GYP= (G1P*LATT_CUR%B(2,1)+G2P*LATT_CUR%B(2,2)+G3P*LATT_CUR%B(2,3)-QY) *TPI
                GZP= (G1P*LATT_CUR%B(3,1)+G2P*LATT_CUR%B(3,2)+G3P*LATT_CUR%B(3,3)-QZ) *TPI
                GLENP(IND)=MAX(SQRT(GXP*GXP+GYP*GYP+GZP*GZP),1E-10_q)
             ENDIF

          ENDDO

          INDMAX=IND-1
!=======================================================================
! now calculate the tables containing the spherical harmonics
! multiply by the pseudopotential and 1/(OMEGA)^(1/2)
!=======================================================================
          CALL SETYLM(LYDIM,INDMAX,YLM,XS,YS,ZS)

          typ: DO NT=1,NONL_S%NTYP

             LMIND=1
             l_loop: DO L=1,P(NT)%LMAX
!=======================================================================
! first interpolate the non-local pseudopotentials
! and multiply by 1/(OMEGA)^(1/2)
!=======================================================================
                FAKT= 1/SQRT(LATT_CUR%OMEGA)
                ARGSC=NPSNL/P(NT)%PSMAXN
!$OMP SIMD PRIVATE(ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3)
                DO IND=1,INDMAX
                   ARG=(GLEN(IND)*ARGSC)+1
                   NADDR=INT(ARG)

                   VPS(IND)=0._q

                   IF (ASSOCIATED(P(NT)%PSPNL_SPLINE)) THEN
                      IF (NADDR<NPSNL) THEN
                         NADDR  =MIN(INT(GLEN(IND)*ARGSC)+1,NPSNL-1)
                         REM=GLEN(IND)-P(NT)%PSPNL_SPLINE(NADDR,1,L)
                         VPS(IND)=(P(NT)%PSPNL_SPLINE(NADDR,2,L)+REM*(P(NT)%PSPNL_SPLINE(NADDR,3,L)+ &
                           &         REM*(P(NT)%PSPNL_SPLINE(NADDR,4,L)+REM*P(NT)%PSPNL_SPLINE(NADDR,5,L))))*FAKT*FAKTX(IND)
                      ENDIF
                   ELSE IF (NADDR<NPSNL-2) THEN
                      REM=MOD(ARG,1.0_q)
                      V1=P(NT)%PSPNL(NADDR-1,L)
                      V2=P(NT)%PSPNL(NADDR,L  )
                      V3=P(NT)%PSPNL(NADDR+1,L)
                      V4=P(NT)%PSPNL(NADDR+2,L)
                      T0=V2
                      T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
                      T2=(V1+V3-(2*V2))/2._q
                      T3=(V4-V1+(3*(V2-V3)))/6._q
                      VPS(IND)=(T0+REM*(T1+REM*(T2+REM*T3)))*FAKT*FAKTX(IND)
                   ENDIF

                   IF (VPS(IND)/=0._q.AND.PRESENT(DK)) THEN
                      ARG=(GLENP(IND)*ARGSC)+1
                      NADDR=INT(ARG)
                      IF (NADDR>NPSNL-1.OR.(.NOT.ASSOCIATED(P(NT)%PSPNL_SPLINE).AND.NADDR>NPSNL-3)) VPS(IND)=0._q
                   ENDIF
                ENDDO
!=======================================================================
! initialize to 0
!=======================================================================
                LL=P(NT)%LPS(L)
                MMAX=2*LL

                IF (IZERO==1 .OR. IZERO==-1) THEN
                   DO LM=0,MMAX
!DIR$ IVDEP
!OCL NOVREC
                      DO IND=1,INDMAX
                         NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)=0
                      ENDDO
                   ENDDO
                ENDIF
                IF (IZERO==-1) THEN
                   DO IND=1,INDMAX
                      VPS(IND)=-VPS(IND)
                   ENDDO
                ENDIF

!=======================================================================
! now multiply with the spherical harmonics
! and set "phase factor" CSET
!=======================================================================
                IF (LL==0) THEN
                   CSET=1.0_q
                ELSE IF (LL==1) THEN
                   CSET=(0.0_q,1.0_q)
                ELSE IF (LL==2) THEN
                   CSET=-1.0_q
                ELSE IF (LL==3) THEN
                   CSET=(0.0_q,-1.0_q)
                ENDIF

                DO LM=0,MMAX
                   NONL_S%CQFAK(LMIND+LM,NT)=CSET
                ENDDO
                LMBASE=LL**2+1

                DO LM=0,MMAX
                   DO IND=1,INDMAX
                      NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)= NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)+ &
                           VPS(IND)*YLM(IND,LM+LMBASE)
                   ENDDO
                ENDDO

                LMIND=LMIND+MMAX+1

             ENDDO l_loop
             IF (LMIND-1/=P(NT)%LMMAX) THEN
                CALL vtutor%bug("internal ERROR: SPHER: LMMAX is wrong " // str(P(NT)%LMMAX) // " " // &
                   str(LMIND-1), "nonl.F", 478)
             ENDIF
          ENDDO typ

!=======================================================================
! and of loop over special-points
!=======================================================================
       ENDDO kpoint

! conjugate phase alteration for spin down: -q/2 -> q/2
       QX=-QX
       QY=-QY
       QZ=-QZ
    ENDDO spinor

    DEALLOCATE(GLEN,XS,YS,ZS,VPS,FAKTX,YLM)

    IF (PRESENT(DK)) DEALLOCATE(GLENP)

!! 
!$ACC UPDATE DEVICE(NONL_S%QPROJ,NONL_S%CQFAK) IF_PRESENT 

    

    RETURN
  END SUBROUTINE SPHER


!****************** SUBROUTINE SPHER_DER *****************************
!
!> Calculates the derivative of the non local projector
!> with respect to k for the cartesian
!> component IDIR times -i = -i nabla_k <p_k|
!>
!> This is equivalent to the projector times r_i  (i=IDIR)
!> ~~~
!> <p r_i| G+k> =             \int p(r) r_i Y_lm(r) e^ir(G+k) dr
!>              = 1/i  d/dk_i \int p(r)     Y_lm(r) e^ir(G+k) dr
!>              = -i d/dk_i <p | G+k>
!> ~~~
!> For simplicity finite differences are used
!
!*********************************************************************

  SUBROUTINE SPHER_DER(GRID,NONL_S,P,WDES,LATT_CUR,  IDIR)
    USE pseudo
    USE mpimy
    USE mgrid
    USE lattice
    USE constant
    USE asa

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    INTEGER IDIR
! local
    REAL(q),POINTER  :: VKPT_ORIG(:,:)
    REAL(q) :: V(3)
! displacement for finite differences
    REAL(q), PARAMETER :: DIS=fd_displacement
    INTEGER IDIS,NK

    

! reset all entries
    NONL_S%QPROJ=0

! first generate a pointer to the original k-point array
    VKPT_ORIG=>WDES%VKPT
    NULLIFY(WDES%VKPT)
    ALLOCATE(WDES%VKPT(3,WDES%NKPTS))

! k-displacement
    V=0
    V(IDIR)=DIS/TPI
    CALL KARDIR(1,V(1),LATT_CUR%A)

! negative displacement
    DO NK=1,WDES%NKPTS
       WDES%VKPT(:,NK)=VKPT_ORIG(:,NK)-V
    ENDDO
!   CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,-1)
    CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,-1,DK= 2*V)

! positive displacement
    DO NK=1,WDES%NKPTS
       WDES%VKPT(:,NK)=VKPT_ORIG(:,NK)+V
    ENDDO
!   CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR, 0)
    CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR, 0,DK=-2*V)

    NONL_S%QPROJ=NONL_S%QPROJ/DIS/2

    NONL_S%CQFAK=NONL_S%CQFAK*(0.0_q,-1.0_q)


    DEALLOCATE(WDES%VKPT)
    WDES%VKPT=>VKPT_ORIG

    

  END SUBROUTINE SPHER_DER


  SUBROUTINE NONL_ALLOC_DER(NONL_S, NONL_D)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S, NONL_D
! local var
    INTEGER NTYPD, NRPLWV, LMDIM, NKPTS

    NRPLWV=SIZE(NONL_S%QPROJ,1)
    LMDIM =SIZE(NONL_S%QPROJ,2)
    NTYPD =SIZE(NONL_S%QPROJ,3)
    NKPTS =SIZE(NONL_S%QPROJ,4)

    NONL_D=NONL_S

    NULLIFY(NONL_D%QPROJ,NONL_D%CQFAK)
    IF (NONL_S%LSPIRAL) THEN
       ALLOCATE(NONL_D%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
            NONL_D%CQFAK(LMDIM,NTYPD))
    ELSE
       ALLOCATE(NONL_D%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
            NONL_D%CQFAK(LMDIM,NTYPD))
    ENDIF
! NONL_D%QPROJ = 0.0_q; NONL_D%CQFAK = (0.0_q, 0.0_q) ! Breaks ???, operations on entire uninitialized arrays in SPHER_DER.

    RETURN
  END SUBROUTINE NONL_ALLOC_DER

  SUBROUTINE NONL_DEALLOC_DER(NONL_D)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_D
    DEALLOCATE(NONL_D%QPROJ,NONL_D%CQFAK)

    RETURN
  END SUBROUTINE NONL_DEALLOC_DER


!************************ SUBROUTINE PHASE *****************************
!
!> Calculates the phasefactor CREXP (<TT>exp(ig.r)</TT>) for (1._q,0._q) k-point,
!> on the compressed grid of lattice vectors
!>
!> @note CREXP is only calculated if NK changes, if the ionc positions
!> change PHASE must be called with NK=0 to force the routine to recalculate
!> the phase-factor in the next call
!
!***********************************************************************

  SUBROUTINE PHASE(WDES,NONL_S,NK)
    USE constant

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavedes)     WDES
    TYPE (nonl_struct) NONL_S

    

!=======================================================================
! check if special point changed
!=======================================================================
    IF (NK==0 .OR.NK==NONL_S%NK) THEN
       NONL_S%NK=NK
!$ACC UPDATE DEVICE(NONL_S%NK) 
       
       RETURN
    ENDIF
    NONL_S%NK=NK
!$ACC UPDATE DEVICE(NONL_S%NK) 

! number of G-vectors
    NPL= WDES%NGVECTOR(NK)

! set the phase-factor
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(NONL_S,WDES) &
!$ACC PRIVATE(GXDX,GYDY,GZDZ,CGDR) 
    ion:  DO NI=1,NONL_S%NIONS
  GXDX=NONL_S%POSION(1,NI)
  GYDY=NONL_S%POSION(2,NI)
  GZDZ=NONL_S%POSION(3,NI)
 !$OMP SIMD PRIVATE(CGDR)
       DO M=1,NPL
!!     GXDX=NONL_S%POSION(1,NI)
!!     GYDY=NONL_S%POSION(2,NI)
!!     GZDZ=NONL_S%POSION(3,NI)
          CGDR=CITPI*(WDES%IGX(M,NK)*GXDX+WDES%IGY(M,NK)*GYDY+WDES%IGZ(M,NK)*GZDZ)
          NONL_S%CREXP(M,NI)=EXP(CGDR)
       ENDDO
    ENDDO ion

!$ACC UPDATE SELF(NONL_S%CREXP) 

    

    RETURN
  END SUBROUTINE PHASE


!************************ SUBROUTINE VNLACC ****************************
!
!> Calculates the non local contribution of the Hamiltonian,
!> using the reciprocal space projection scheme
!>
!> The projection of the wavefunction onto the projection operators
!> must be supplied in W1\%CPROJ
!
!***********************************************************************

  SUBROUTINE VNLACC(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    REAL(q) :: EVALUE
    COMPLEX(q) :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
! local
    COMPLEX(q) :: CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) 
!$ACC KERNELS PRESENT(CACC) 
    CACC=0
!$ACC END KERNELS
    CALL OVERL1(W1%WDES1, SIZE(CDIJ,1),CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))
!$ACC EXIT DATA DELETE(CRESUL) 
    RETURN
  END SUBROUTINE VNLACC


  SUBROUTINE VNLACC_ADD(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    REAL(q) :: EVALUE
    COMPLEX(q) :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
! work array
    COMPLEX(q) :: CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) 
    CALL OVERL1(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))
!$ACC EXIT DATA DELETE(CRESUL) 
    RETURN
  END SUBROUTINE VNLACC_ADD


  SUBROUTINE VNLACC_C(NONL_S, W1, CDIJ, CQIJ, ISP, CEVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER :: LMDIM
    COMPLEX(q) :: CEVALUE
    COMPLEX(q) :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
! work array
    COMPLEX(q) :: CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) 
!$ACC KERNELS PRESENT(CACC) 
    CACC=0
!$ACC END KERNELS
    CALL OVERL1_C(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), CEVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))
!$ACC EXIT DATA DELETE(CRESUL) 

    RETURN
  END SUBROUTINE VNLACC_C


  SUBROUTINE VNLACC_ADD_C(NONL_S,W1, CDIJ, CQIJ, ISP, CEVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER :: LMDIM
    COMPLEX(q) :: CEVALUE
    COMPLEX(q) :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
! work array
    COMPLEX(q) :: CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) 
    CALL OVERL1_C(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), CEVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))
!$ACC EXIT DATA DELETE(CRESUL) 

    RETURN
  END SUBROUTINE VNLACC_ADD_C


! here CDIJ and QCIJ are always complex
  SUBROUTINE VNLACC_ADD_CCDIJ(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    REAL(q) :: EVALUE
    COMPLEX(q) :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
! work array
    COMPLEX(q) :: CRESUL(W1%WDES1%NPRO)

    CALL OVERL1_CCDIJ(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))

    RETURN
  END SUBROUTINE VNLACC_ADD_CCDIJ


!************************ SUBROUTINE PROJ ******************************
!
!> Calculates the projection of all bands of (1._q,0._q) specific
!> k-point onto the reciprocal space projection operators
!
!***********************************************************************

  SUBROUTINE PROJ(NONL_S,WDES,W,NK)
    USE mopenmp_struct_def
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    INTEGER NK
! local
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    INTEGER ISP, N

    

!$ACC ENTER DATA CREATE(WDES1,W1) 
    CALL SETWDES(WDES,WDES1,NK)

    DO ISP=1,WDES%ISPIN
       DO N=1,WDES%NBANDS
          CALL SETWAV(W,W1,WDES1,N,ISP)
          CALL PROJ1(NONL_S,WDES1,W1)
       ENDDO
    ENDDO
# 841

    

  END SUBROUTINE PROJ


!************************ SUBROUTINE PROJ1 *****************************
!
!> Calculates the scalar product of (1._q,0._q) wavefunction with
!> all projector functions in reciprocal space
!> (thesis gK Eq. 10.34)
!>
!> @details @ref openmp :
!> under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is contracted to a single loop over \"all ions\", and the
!> latter is distributed over all available threads.
!
!***********************************************************************

  SUBROUTINE PROJ1(NONL_S,WDES1,W1)
    IMPLICIT NONE

    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes1)    :: WDES1
    TYPE (wavefun1), TARGET :: W1

# 872

    CALL PROJ1_(NONL_S,WDES1,W1)

    RETURN
  END SUBROUTINE PROJ1

  SUBROUTINE PROJ1_(NONL_S,WDES1,W1)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1), TARGET :: W1

! local
    INTEGER :: NPL, LMBASE, ISPIRAL, ISPINOR, NT, LMMAXC, NI, K, LM
    REAL(q) :: SUMR, SUMI
    REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
    COMPLEX(q)    :: CPROJ(WDES1%NPRO_TOT)

    

    IF (WDES1%NK /= NONL_S%NK) THEN
       CALL vtutor%bug("internal error in PROJ1: PHASE not properly set up " // str(WDES1%NK) // " " &
          // str(NONL_S%NK), "nonl.F", 895)
    ENDIF

    NPL=WDES1%NGVECTOR

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
!=======================================================================
! performe loop over ions
!=======================================================================
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NI,NT,LMMAXC,LMBASE,WORK,TMP,LM,SUMR,SUMI)
       ion: DO NI=1,NONL_S%NIONS
          NT=NONL_S%ITYP(NI)
          LMMAXC=NONL_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion
          LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)
!=======================================================================
! multiply with phasefactor and divide into real and imaginary part
!=======================================================================
          CALL CREXP_MUL_WAVE( NPL, NONL_S%CREXP(1,NI), W1%CPTWFP(NPL*ISPINOR+1), WORK(1))
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
# 930

          
          CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, WORK(1) , 1 , 0._q ,  TMP(1,1), 1)
          CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, WORK(1+NPL) , 1 , 0._q ,  TMP(1,2), 1)
          

          l_loop: DO LM=1,LMMAXC
             SUMR=TMP(LM,1)
             SUMI=TMP(LM,2)
             CPROJ(LM+LMBASE)=(CMPLX( SUMR , SUMI ,KIND=q) *NONL_S%CQFAK(LM,NT))
          ENDDO l_loop

       ENDDO ion
!$OMP END PARALLEL DO
    ENDDO spinor

! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

    

    RETURN
  END SUBROUTINE PROJ1_


!************************ SUBROUTINE PROJ1_ACC *************************
!
!***********************************************************************
# 1046


!************************ SUBROUTINE PROJMU ****************************
!
!***********************************************************************

  SUBROUTINE PROJMU(NONL_S,WDES1,W1)
    IMPLICIT NONE
    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes1)    :: WDES1
    TYPE (wavefun1)    :: W1(:)

! local
    INTEGER :: NSIM,NP

    

    NSIM=SIZE(W1)

# 1070

    DO NP=1,NSIM
       IF (W1(NP)%LDO) CALL PROJ1(NONL_S,WDES1,W1(NP))
    ENDDO

    

    RETURN
  END SUBROUTINE PROJMU


!************************ SUBROUTINE PROJMU_ACC ************************
!
!***********************************************************************
# 1181


!************************ SUBROUTINE PROJXYZ ***************************
!
!> Calculates the first order change of the wave function character
!> upon moving the ions for (1._q,0._q) selected k-point and spin component
!>
!> The results are stored in CPROJXYZ
!>
!> A factor 2 stemming from variations of the bra or kat is included
!>
!> @param LADDK (optional) LADDK=.TRUE.
!> forces the treatment usually 1._q only for <TT>#define </TT>.
!> This allows for correct results for the Born effective charges
!> even if <TT></TT> is not defined
!>
!> @details @ref openmp :
!> loops over the basis vectors of the orbitals are distributed
!> over all available threads.
!
!***********************************************************************

  SUBROUTINE PROJXYZ(NONL_S, WDES, W, LATT_CUR, ISP, NK, CPROJXYZ, LADDK)
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes)     :: WDES
    TYPE (wavespin)    :: W
    TYPE (latt)        :: LATT_CUR

    COMPLEX(q)    :: CPROJXYZ(WDES%NPROD, WDES%NBANDS, 3)
    INTEGER :: ISP, NK

    LOGICAL, OPTIONAL :: LADDK

! local variable
    TYPE (wavedes1) :: WDES1  ! descriptor for (1._q,0._q) k-point

    COMPLEX(q) :: CMUL
    INTEGER    :: NPL,N,LMBASE,LMMAXC,ISPINOR,ISPIRAL,K,KK,NI,NT,LM

    COMPLEX(q) :: CX_,CY_,CZ_,CVAL
# 1228

    COMPLEX(q) :: CX(WDES%NPRO_TOT),CY(WDES%NPRO_TOT),CZ(WDES%NPRO_TOT)


    LOGICAL :: LVKPT

    

    NPL=WDES%NGVECTOR(NK)

    LVKPT=.FALSE. ; IF (PRESENT(LADDK)) LVKPT=LADDK

!$ACC ENTER DATA CREATE(CX,CY,CZ) 

!$ACC ENTER DATA CREATE(WDES1) 
    CALL SETWDES(WDES,WDES1,NK)
    CALL PHASE(WDES,NONL_S,NK)

# 1250


 !!!$OMP PARALLEL DO PRIVATE(N,CX,CY,CZ,LMBASE,ISPIRAL,ISPINOR,NI,NT,LMMAXC,LM,CMUL,CX_,CY_,CZ_,K,KK,CVAL) &
 !!!$OMP DEFAULT(SHARED) SCHEDULE(static)
!$ACC PARALLEL LOOP GANG COLLAPSE(3) PRESENT(NONL_S,W,WDES,LATT_CUR,CX,CY,CZ) &
!$ACC PRIVATE(ISPIRAL,NT,LMMAXC,LMBASE) 
    band: DO N=1,WDES%NBANDS

  CX=0; CY=0; CZ=0

       spinor: DO ISPINOR=0,WDES%NRSPINORS-1
          ion: DO NI=1,NONL_S%NIONS

             ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1

             NT=NONL_S%ITYP(NI)
             LMMAXC=NONL_S%LMMAX(NT)
             IF (LMMAXC==0) CYCLE ion
             LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)

!$ACC LOOP PRIVATE(CMUL,CX_,CY_,CZ_)
             l_loop: DO LM=1,LMMAXC
                CMUL=NONL_S%CQFAK(LM,NT)*CITPI
                CX_=0; CY_=0; CZ_=0

                IF (LVKPT) THEN
!$ACC LOOP VECTOR PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
 !$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(static) &
 !$OMP PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CVAL=2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*W%CPTWFP(KK,N,NK,ISP)*NONL_S%CREXP(K,NI)*CMUL
                      CX_=CX_+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                      CY_=CY_+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                      CZ_=CZ_+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
                   ENDDO
 !$OMP END PARALLEL DO SIMD
                ELSE
!$ACC LOOP VECTOR PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
 !$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(static) &
 !$OMP PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CVAL=2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*W%CPTWFP(KK,N,NK,ISP)*NONL_S%CREXP(K,NI)*CMUL

                      CX_=CX_+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                      CY_=CY_+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                      CZ_=CZ_+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
# 1302

                   ENDDO
 !$OMP END PARALLEL DO SIMD
                ENDIF

# 1311

                CX(LM+LMBASE)=CX_*LATT_CUR%B(1,1)+CY_*LATT_CUR%B(1,2)+CZ_*LATT_CUR%B(1,3)
                CY(LM+LMBASE)=CX_*LATT_CUR%B(2,1)+CY_*LATT_CUR%B(2,2)+CZ_*LATT_CUR%B(2,3)
                CZ(LM+LMBASE)=CX_*LATT_CUR%B(3,1)+CY_*LATT_CUR%B(3,2)+CZ_*LATT_CUR%B(3,3)


             ENDDO l_loop
          ENDDO ion
       ENDDO spinor


       CALL DIS_PROJ(WDES1,CX,CPROJXYZ(1,N,1))
       CALL DIS_PROJ(WDES1,CY,CPROJXYZ(1,N,2))
       CALL DIS_PROJ(WDES1,CZ,CPROJXYZ(1,N,3))


    ENDDO band
 !!!$OMP END PARALLEL DO

# 1339


    

    RETURN
  END SUBROUTINE PROJXYZ


!************************ SUBROUTINE PROJXYZ_WA ************************
!
!> Same as #PROJXYZ, however, is accepts a single array for the
!> orbitals' PW coefficients instead  of a wavefunction array W
!>
!> This subroutine is used by the module greens_real_space
!>
!> @details @ref openmp :
!> loops over the basis vectors of the orbitals are distributed
!> over all available threads.
!
!***********************************************************************

  SUBROUTINE PROJXYZ_WA(NONL_S, WDES, CPTWFP, LATT_CUR, NK, CPROJXYZ, LADDK)
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    COMPLEX(q),CONTIGUOUS :: CPTWFP(:,:)
    TYPE (latt)        LATT_CUR
    INTEGER            ISP, NK
    COMPLEX(q), CONTIGUOUS :: CPROJXYZ(:,:,:)
    LOGICAL, OPTIONAL :: LADDK
! local variable
    TYPE (wavedes1)    WDES1         ! descriptor for (1._q,0._q) k-point
    COMPLEX(q)       :: CX(WDES%NPRO_TOT) ,CY(WDES%NPRO_TOT) ,CZ(WDES%NPRO_TOT)
    COMPLEX(q)       :: CX_, CY_, CZ_
    INTEGER    :: NPRO, NPL, N, LMBASE, LMMAXC, ISPINOR, ISPIRAL, K, KK, NI, NIS, NT, LM
    COMPLEX(q) :: CMUL
    COMPLEX(q)       :: CVAL

    

    NPRO=WDES%NPRO_TOT
    NPL= WDES%NGVECTOR(NK)

    CALL SETWDES(WDES,WDES1,NK)
    CALL PHASE(WDES,NONL_S,NK)

!!!$OMP PARALLEL DO PRIVATE(N,CX,CY,CZ,LMBASE,ISPIRAL,ISPINOR,NIS,NT,LMMAXC,NI,LM,CMUL,CX_,CY_,CZ_,K,KK,CVAL) &
!!!$OMP DEFAULT(SHARED) SCHEDULE(static)
    band: DO N=1,WDES%NBANDS

       CX=0
       CY=0
       CZ=0

       LMBASE= 0

       ISPIRAL=1
       spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
          NIS=1
          typ:  DO NT=1,NONL_S%NTYP
             LMMAXC=NONL_S%LMMAX(NT)
             IF (LMMAXC==0) GOTO 100

             ions:   DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
                l_loop: DO LM=1,LMMAXC
                   CMUL=NONL_S%CQFAK(LM,NT)*CITPI
                   CX_=0
                   CY_=0
                   CZ_=0
                   IF (PRESENT (LADDK)) THEN
!$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(static) &
!$OMP PRIVATE(KK,CVAL)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CVAL =  2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*CPTWFP(KK,N)*NONL_S%CREXP(K,NI)*CMUL
                      CX_=CX_+CVAL*(WDES%IGX(K,NK)+WDES%VKPT(1,NK))
                      CY_=CY_+CVAL*(WDES%IGY(K,NK)+WDES%VKPT(2,NK))
                      CZ_=CZ_+CVAL*(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))
                   ENDDO
!$OMP END PARALLEL DO SIMD
                   ELSE
!$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(static) &
!$OMP PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
! factor 2 because either the bra or kat can vary
                      CVAL =  2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*CPTWFP(KK,N)*NONL_S%CREXP(K,NI)*CMUL

                      CX_=CX_+CVAL*(WDES%IGX(K,NK)+WDES%VKPT(1,NK))
                      CY_=CY_+CVAL*(WDES%IGY(K,NK)+WDES%VKPT(2,NK))
                      CZ_=CZ_+CVAL*(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))
# 1438

                   ENDDO
!$OMP END PARALLEL DO SIMD
                   ENDIF

                   CX(LM+LMBASE)=CX_*LATT_CUR%B(1,1)+CY_*LATT_CUR%B(1,2)+CZ_*LATT_CUR%B(1,3)
                   CY(LM+LMBASE)=CX_*LATT_CUR%B(2,1)+CY_*LATT_CUR%B(2,2)+CZ_*LATT_CUR%B(2,3)
                   CZ(LM+LMBASE)=CX_*LATT_CUR%B(3,1)+CY_*LATT_CUR%B(3,2)+CZ_*LATT_CUR%B(3,3)

                ENDDO l_loop
                LMBASE= LMMAXC+LMBASE
             ENDDO ions

100          NIS = NIS+NONL_S%NITYP(NT)
          ENDDO typ

          IF (NONL_S%LSPIRAL) ISPIRAL=2
       ENDDO spinor

       CALL DIS_PROJ(WDES1,CX,CPROJXYZ(1,N,1))
       CALL DIS_PROJ(WDES1,CY,CPROJXYZ(1,N,2))
       CALL DIS_PROJ(WDES1,CZ,CPROJXYZ(1,N,3))

    ENDDO band
!!!$OMP END PARALLEL DO

    

    RETURN
  END SUBROUTINE PROJXYZ_WA


!************************ SUBROUTINE PROJLAT_DER ***********************
!
!> Calculates the first order change of the wave function character
!> upon changing the lattice
!>
!> The results are stored in CPROJXYZ
!
!***********************************************************************

  SUBROUTINE PROJLAT_DER(P, NONL_S, WDES, W, LATT_CUR, ISP, NK, CPROJXYZ)
    USE lattice
    USE constant
    USE pseudo
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc

    IMPLICIT NONE

    TYPE (potcar)      P(:)
    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (latt)        LATT_CUR
    INTEGER            ISP, NK
    COMPLEX(q), TARGET :: CPROJXYZ(WDES%NPROD, WDES%NBANDS, 6)
! local variable
    TYPE (wavedes1)    WDES1         ! descriptor for (1._q,0._q) k-point
    TYPE (wavefun1)    W1
    INTEGER    :: N
    INTEGER :: IDIR, JDIR, IJDIR
    TYPE (latt)        LATT_FIN1,LATT_FIN2
    REAL(q) :: DIS=fd_displacement

    

!$ACC ENTER DATA CREATE(WDES1) 
    CALL SETWDES(WDES,WDES1,NK)
    CALL PHASE(WDES,NONL_S,NK)

    IJDIR=0
    DO IDIR=1,3
       DO JDIR=1,IDIR
          IJDIR=IJDIR+1
          LATT_FIN1=LATT_CUR
          LATT_FIN1%A(IDIR,:)=LATT_CUR%A(IDIR,:)+DIS*LATT_CUR%A(JDIR,:)

          LATT_FIN2=LATT_CUR
          LATT_FIN2%A(IDIR,:)=LATT_CUR%A(IDIR,:)-DIS*LATT_CUR%A(JDIR,:)

          CALL LATTIC(LATT_FIN1)
          CALL LATTIC(LATT_FIN2)

          CALL SPHER(WDES%GRID,NONL_S,P,WDES,LATT_FIN2, IZERO=-1, NK_SELECT=NK)
          CALL SPHER(WDES%GRID,NONL_S,P,WDES,LATT_FIN1, IZERO=0 , NK_SELECT=NK)

!!!! !$OMP PARALLEL DO PRIVATE(W1) NUM_THREADS(omp_nthreads_acc)
          band: DO N=1,WDES%NBANDS

!!             

             CALL SETWAV(W,W1,WDES1,N,ISP); W1%CPROJ => CPROJXYZ(:,N,IJDIR)
# 1532

             CALL PROJ1(NONL_S,WDES1,W1) ! calculate W1%CPROJ (linked to CPROJXYZ)

!$ACC KERNELS PRESENT(CPROJXYZ) 
             CPROJXYZ(:,N,IJDIR)=CPROJXYZ(:,N,IJDIR)/DIS
!$ACC END KERNELS
# 1540

          ENDDO band
!!!! !$OMP END PARALLEL DO
       ENDDO
    ENDDO

    CALL SPHER(WDES%GRID, NONL_S, P, WDES, LATT_CUR,  1, NK_SELECT=NK)

# 1550


    

    RETURN
  END SUBROUTINE PROJLAT_DER


!************************ SUBROUTINE PROJ_DER **************************
!
!> Calculates the first order change of the wave function character
!> upon moving (1._q,0._q) ion ION in the direction IDIR
!>
!> A term  (i k)  is included as well to "center" the derivate at the
!> ion
!
!***********************************************************************

  SUBROUTINE PROJ_DER(NONL_S, WDES1, W1, LATT_CUR, ION, IDIR)

!! USE moffload_struct_def

    USE lattice
    USE constant
    USE iso_c_binding

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    TYPE (latt)        LATT_CUR
    INTEGER            ION, IDIR
! local variable
    COMPLEX(q)       :: CVAL
    INTEGER    :: NPL, N, LMBASE, LMMAXC, ISPINOR, ISPIRAL, K, KK, NT, LM
    COMPLEX(q) :: CMUL
    COMPLEX(q),TARGET :: CTMP(WDES1%NRPLWV)
    REAL(q),CONTIGUOUS,POINTER :: CTMP_REAL(:)
    COMPLEX(q),TARGET       :: C(WDES1%NPRO_TOT)
    REAL(q) :: SUMR, SUMI
    REAL(q) :: TMP(101,2)

    

    NPL= WDES1%NGVECTOR
    CALL C_F_POINTER(C_LOC(CTMP),CTMP_REAL,[2*WDES1%NRPLWV])

!$ACC ENTER DATA CREATE(C,CTMP) 
!$ACC KERNELS 
    C=0
!$ACC END KERNELS

    ISPIRAL=1

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       NT=NONL_S%ITYP(ION)
       LMMAXC=NONL_S%LMMAX(NT)
       IF (LMMAXC==0) CYCLE
       LMBASE=NONL_S%LMBASE(ION)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)

       IF (IDIR==0) THEN
!$ACC PARALLEL LOOP PRIVATE(KK) PRESENT(W1,W1%CPTWFP,NONL_S,CTMP) 
!DIR$ IVDEP
!OCL NOVREC
          DO K=1,NPL
             KK=K+NPL*ISPINOR
             CTMP(K)=W1%CPTWFP(KK)*NONL_S%CREXP(K,ION)
          ENDDO
       ELSE
!$ACC PARALLEL LOOP PRIVATE(KK) PRESENT(W1,W1%CPTWFP,NONL_S,WDES1,LATT_CUR,CTMP) 
!DIR$ IVDEP
!OCL NOVREC
          DO K=1,NPL
             KK=K+NPL*ISPINOR
             CTMP(K)=W1%CPTWFP(KK)*NONL_S%CREXP(K,ION)* &

                 ((WDES1%IGX(K)+WDES1%VKPT(1))*LATT_CUR%B(IDIR,1)+ &
                  (WDES1%IGY(K)+WDES1%VKPT(2))*LATT_CUR%B(IDIR,2)+ &
                  (WDES1%IGZ(K)+WDES1%VKPT(3))*LATT_CUR%B(IDIR,3))
# 1634

          ENDDO
       ENDIF

# 1660

       CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, CTMP_REAL(1), 2, 0._q,  TMP(1,1), 1)
       CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, CTMP_REAL(2), 2, 0._q,  TMP(1,2), 1)
       IF (IDIR==0) THEN
           DO LM=1,LMMAXC
               SUMR = TMP(LM,1)
               SUMI = TMP(LM,2)
               C(LM+LMBASE)=(CMPLX( SUMR, SUMI, KIND=q)*NONL_S%CQFAK(LM,NT))
           ENDDO
       ELSE
           DO LM=1,LMMAXC
               SUMR = TMP(LM,1)
               SUMI = TMP(LM,2)
               C(LM+LMBASE)=(CMPLX( SUMR, SUMI, KIND=q)*NONL_S%CQFAK(LM,NT)*CITPI)
           ENDDO
       ENDIF

       IF (NONL_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor

    CALL DIS_PROJ(WDES1,C,W1%CPROJ(1))
!$ACC EXIT DATA DELETE(C,CTMP) 

    

    RETURN
  END SUBROUTINE PROJ_DER


!************************ SUBROUTINE VNLACC_DER ************************
!
!> Calculates the non local contribution of the first derivative
!> of the Hamiltonian upon moving (1._q,0._q) ion ION in the direction IDIR
!>
!> More specifically the contribution
!> ~~~
!>    | d p_i/ d R > (D_ij - epsilon Q_ij) c_j
!> ~~~
!> For IDIR=0 the derivative is replaced by | p_i >
!
!***********************************************************************

  SUBROUTINE VNLACC_DER(NONL_S, W1, &
       &     CDIJ,CQIJ,EVALUE,  CACC, LATT_CUR, ION, IDIR)
    USE lattice

    IMPLICIT NONE
    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER LMDIM
    REAL(q) EVALUE
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)   CACC(:)
    TYPE (latt)  LATT_CUR
    INTEGER      ION, IDIR
! work arrays
    COMPLEX(q)         CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) 
    CALL OVERL1(W1%WDES1,SIZE(CDIJ,1),CDIJ(1,1,1,1),CQIJ(1,1,1,1), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0_DER(NONL_S,W1%WDES1,CRESUL(1),CACC(1), LATT_CUR, ION, IDIR )
!$ACC EXIT DATA DELETE(CRESUL) 

    RETURN
  END SUBROUTINE VNLACC_DER

  SUBROUTINE VNLACC_DER_C(NONL_S, W1, &
       &     CDIJ, CQIJ,EVALUE,  CACC, LATT_CUR, ION, IDIR)
    USE lattice

    IMPLICIT NONE
    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER LMDIM
    COMPLEX(q)   EVALUE
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)   CACC(:)
    TYPE (latt)  LATT_CUR
    INTEGER      ION, IDIR
! work arrays
    COMPLEX(q)         CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) 
    CALL OVERL1_C(W1%WDES1,SIZE(CDIJ,1),CDIJ(1,1,1,1),CQIJ(1,1,1,1), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0_DER(NONL_S,W1%WDES1,CRESUL(1),CACC(1), LATT_CUR, ION, IDIR )
!$ACC EXIT DATA DELETE(CRESUL) 

    RETURN
  END SUBROUTINE VNLACC_DER_C


!************************ SUBROUTINE FORNL *****************************
!
!> Calculates the forces related to the non local pseudopotential
!>
!> The projection of the wavefunction onto the projector functions
!> must be stored in W\%CPROJ on entry
!>
!> Algorithm:
!> ~~~
!> ACC(G') = SUM(L,L',R)
!>         SUM(G)   QPROJ(G ,L)  EXP( iG  R) C(G)  iG *
!>         D(L,L') + Evalue  Q(L,L')
!>         SUM(GP)  QPROJ(GP,LP) EXP(-iGP R) C(GP)
!> ~~~
!
!***********************************************************************

  SUBROUTINE FORNL(NONL_S,WDES,W,LATT_CUR, LMDIM,CDIJ,CQIJ,EINL)

!! USE moffload

    USE lattice
    USE constant

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (latt)        LATT_CUR
    TYPE (wavedes1)    WDES1         ! descriptor for (1._q,0._q) k-point

    DIMENSION EINL(3,NONL_S%NIONS)
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! work arrays
    REAL(q) ::  ENL(NONL_S%NIONS)
    COMPLEX(q) :: CM
    COMPLEX(q),ALLOCATABLE :: CX(:) ,CY(:) ,CZ(:), C(:)
    COMPLEX(q),ALLOCATABLE :: CXL(:),CYL(:),CZL(:),CL(:)

    

    N =WDES%NPRO_TOT
    NL=WDES%NPRO

    ALLOCATE(CX(N),CY(N),CZ(N),C(N),CXL(NL),CYL(NL),CZL(NL),CL(NL))
!$ACC ENTER DATA CREATE(WDES1,CX,CY,CZ,C,CXL,CYL,CZL,CL,EINL,ENL) 

!$ACC KERNELS PRESENT(EINL,ENL) 
    EINL=0._q
    ENL =0._q
!$ACC END KERNELS
!=======================================================================
! loop over special points, and bands
!=======================================================================

    kpoint: DO NK=1,WDES%NKPTS

       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

       NPL= WDES%NGVECTOR(NK)
       CALL SETWDES(WDES,WDES1,NK)

       CALL PHASE(WDES,NONL_S,NK)

       spin:   DO ISP=1,WDES%ISPIN
          band: DO N=1,WDES%NBANDS

             EVALUE=W%CELEN(N,NK,ISP)
             WEIGHT=W%FERWE(N,NK,ISP)*WDES%WTKPT(NK)*WDES%RSPIN
!=======================================================================
! first build up CX, CY, CZ for tables
!=======================================================================
!$ACC KERNELS PRESENT(CX,CY,CZ,C) 
             CX=0
             CY=0
             CZ=0
             C=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,NONL_S,W,CX,CY,CZ) &
!$ACC& PRIVATE(NT,LMMAXC,LMBASE,CMUL,CTMP1,CTMP2,CTMP3,CTMP) 
             spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
           ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
                ions: DO NI=1,NONL_S%NIONS
!!              ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
                   NT=NONL_S%ITYP(NI)
                   LMMAXC=NONL_S%LMMAX(NT)
                   IF (LMMAXC==0) CYCLE
                   LMBASE =NONL_S%LMBASE(NI)+ISPINOR *NONL_S%LMBASE(NONL_S%NIONS+1)

                   l_loop: DO LM=1,LMMAXC
                      CMUL=NONL_S%CQFAK(LM,NT)*CITPI

                      CTMP=0; CTMP1=0; CTMP2=0; CTMP3=0

!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP VECTOR REDUCTION(+:CTMP1,CTMP2,CTMP3,CTMP) PRIVATE(KK,CVAL)
                      DO K=1,NPL
                         KK=K+NPL*ISPINOR
                         CVAL=(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W%CPTWFP(KK,N,NK,ISP))*CMUL
! new version calculate the projected wave function character on the fly
                         CTMP=CTMP+(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W%CPTWFP(KK,N,NK,ISP))*NONL_S%CQFAK(LM,NT)

! new version add i k which contributes nothing to the energy derivative
                         CTMP1=CTMP1+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                         CTMP2=CTMP2+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                         CTMP3=CTMP3+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
# 1868

                      ENDDO
                      C (LM+LMBASE)=C (LM+LMBASE)+CTMP
                      CX(LM+LMBASE)=CX(LM+LMBASE)+CTMP1
                      CY(LM+LMBASE)=CY(LM+LMBASE)+CTMP2
                      CZ(LM+LMBASE)=CZ(LM+LMBASE)+CTMP3
                   ENDDO l_loop
                ENDDO ions
             ENDDO spinor

             CALL DIS_PROJ(WDES1,C ,CL )
             CALL DIS_PROJ(WDES1,CX,CXL)
             CALL DIS_PROJ(WDES1,CY,CYL)
             CALL DIS_PROJ(WDES1,CZ,CZL)
!=======================================================================
! sum up local contributions
! calculate SUM_LP  D(LP,L)-E Q(LP,L) * C(LP)
!=======================================================================
!$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(WDES,NONL_S,W,CX,CY,CZ) &
!$ACC& PRIVATE(NT,LMMAXC,LMBASE,LMBASE_,NIP,CTMP,CTMP1,CTMP2,CTMP3) 
             spinor2: DO ISPINOR=0,WDES%NRSPINORS-1
                DO ISPINOR_=0,WDES%NRSPINORS-1

                   ions2: DO NI=1,WDES%NIONS
                      NT=WDES%ITYP(NI)
                      LMMAXC=WDES%LMMAX(NT)
                      IF (LMMAXC==0) CYCLE
                      LMBASE =WDES%LMBASE(NI)+ISPINOR *WDES%NPRO/2
                      LMBASE_=WDES%LMBASE(NI)+ISPINOR_*WDES%NPRO/2

                      NIP=NI_GLOBAL(NI, WDES%COMM_INB) !  local storage index

                      CTMP=0; CTMP1=0; CTMP2=0; CTMP3=0

!$ACC LOOP COLLAPSE(2) REDUCTION(+:CTMP) REDUCTION(-:CTMP1,CTMP2,CTMP3) PRIVATE(CM)
                      l_loop2: DO LM=1,LMMAXC
                         DO LMP=1,LMMAXC
                            CM=(CDIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR))* &
! new version use the calculated projected wave function character
                               CONJG(CL(LM+LMBASE))
!                                 CONJG(W%CPROJ(LM+LMBASE,N,NK,ISP))
                            CTMP=CTMP + (W%CPROJ(LMP+LMBASE_,N,NK,ISP)*CM)*WEIGHT
                            CTMP1=CTMP1-(2*WEIGHT)*(CXL(LMP+LMBASE_)*CM)
                            CTMP2=CTMP2-(2*WEIGHT)*(CYL(LMP+LMBASE_)*CM)
                            CTMP3=CTMP3-(2*WEIGHT)*(CZL(LMP+LMBASE_)*CM)
                         ENDDO
                      ENDDO l_loop2
!$ACC ATOMIC UPDATE
                      ENL (NIP)  =ENL (NIP)  +REAL(CTMP, KIND=q)
!$ACC ATOMIC UPDATE
                      EINL(1,NIP)=EINL(1,NIP)+REAL(CTMP1,KIND=q)
!$ACC ATOMIC UPDATE
                      EINL(2,NIP)=EINL(2,NIP)+REAL(CTMP2,KIND=q)
!$ACC ATOMIC UPDATE
                      EINL(3,NIP)=EINL(3,NIP)+REAL(CTMP3,KIND=q)
                   ENDDO ions2
                ENDDO
             ENDDO spinor2
!=======================================================================
          ENDDO band
       ENDDO spin
    ENDDO kpoint
# 1933

!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
    CALL M_sum_d(WDES%COMM, EINL(1,1),NONL_S%NIONS*3)
    CALL M_sum_d(WDES%COMM, ENL (1)  ,NONL_S%NIONS)
    CALL  DIRKAR(NONL_S%NIONS,EINL,LATT_CUR%B)
    DEALLOCATE(CX,CY,CZ,C,CXL,CYL,CZL,CL)

    

    RETURN
  END SUBROUTINE FORNL

!************************ SUBROUTINE STRENL *****************************
!
!> Calculates non-local contributions to stress
!>
!> Easiest to implement and definitly quite fast, is an approach based
!> on finite differences
!
!***********************************************************************

  SUBROUTINE STRENL(GRID,NONL_S,P,W,WDES,LATT_CUR, &
       LMDIM,CDIJ,CQIJ, ISIF,FNLSIF)

!! USE moffload

    USE pseudo
    USE mpimy
    USE mgrid
    USE lattice
    USE constant

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (wavedes1)    WDES1
    TYPE (wavespin)    W
    TYPE (wavefun1)    W1
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR,LATT_FIN

    DIMENSION FNLSIF(3,3)
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! work arrays
    COMPLEX(q),ALLOCATABLE,TARGET ::  CWORK(:)
    REAL(q) FNLSIFTMP

    

    DIS=fd_displacement
!TEST which step should be used in the finite differences
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST
    ALLOCATE(CWORK(WDES%NPROD))

    FNLSIF=0
!$ACC ENTER DATA CREATE(CWORK,WDES1) 

    DO IDIR=1,3
       DO JDIR=1,3
          FNLSIFTMP=0
!=======================================================================
! use central differences to calculate the stress
! set up QPROJ so that central differences can be  calculated
!=======================================================================
          LATT_FIN=LATT_CUR
          IF (ISIF==1) THEN
!  only isotrop pressure
             DO I=1,3; DO J=1,3
                LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
             ENDDO; ENDDO
          ELSE
!  all directions
             DO I=1,3
                LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
             ENDDO
          ENDIF
          CALL LATTIC(LATT_FIN)

          CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,  IZERO=-1)

          LATT_FIN=LATT_CUR
          IF (ISIF==1) THEN
!  only isotrop pressure
             DO I=1,3; DO J=1,3
                LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1-DIS/3)
             ENDDO; ENDDO
          ELSE
!  all directions
             DO I=1,3
                LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)-DIS*LATT_CUR%A(JDIR,I)
             ENDDO
          ENDIF
          CALL LATTIC(LATT_FIN)

          CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,  IZERO=0)

!=======================================================================
! loop over all k-points spin and bands
!=======================================================================
          kpoint: DO NK=1,WDES%NKPTS

             IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

             CALL PHASE(WDES,NONL_S,NK)
             CALL SETWDES(WDES,WDES1,NK)

             spin: DO ISP=1,WDES%ISPIN
                band: DO N=1,WDES%NBANDS
                   CALL SETWAV(W,W1,WDES1,N,ISP); W1%CPROJ => CWORK
# 2053

                   CALL PROJ1(NONL_S,WDES1,W1)  ! calculate W1%CPROJ (linked to CWORK)

                   WEIGHT=WDES%WTKPT(NK)*WDES%RSPIN*W1%FERWE
                   EVALUE=W1%CELEN

!$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(WDES,CWORK,W,W%CPROJ,CDIJ,CQIJ) &
!$ACC PRIVATE(LBASE,LBASE_,NT,LMMAXC) REDUCTION(+:FNLSIFTMP) 
                   spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
                      DO ISPINOR_=0,WDES%NRSPINORS-1

                         DO NI=1,WDES%NIONS
                            NT=WDES%ITYP(NI)
                            LMMAXC=WDES%LMMAX(NT)
                            IF (LMMAXC==0) CYCLE
                            LBASE =WDES%LMBASE(NI)+ISPINOR *WDES%NPRO/2
                            LBASE_=WDES%LMBASE(NI)+ISPINOR_*WDES%NPRO/2

!$ACC LOOP COLLAPSE(2) REDUCTION(+:FNLSIFTMP)
                            DO L=1 ,LMMAXC
 !$OMP SIMD REDUCTION(+:FNLSIFTMP)
                               DO LP=1,LMMAXC
                                  FNLSIFTMP = FNLSIFTMP+WEIGHT*CWORK(LBASE_+LP)* &
                                       CONJG(W%CPROJ(LBASE+L,N,NK,ISP))* &
                                       (CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO spinor
# 2085

                ENDDO band
             ENDDO spin
          ENDDO kpoint
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
          FNLSIF(IDIR,JDIR)=FNLSIFTMP
!
!  only isotrop pressure finish now
!
          IF (ISIF==1) THEN
             FNLSIF(2,2)= FNLSIF(1,1)
             FNLSIF(3,3)= FNLSIF(1,1)
             GOTO 310  ! terminate (not very clean but who cares)
          ENDIF
!=======================================================================
! next direction
!=======================================================================
       ENDDO
    ENDDO

310 CONTINUE
# 2109

    CALL M_sum_d(WDES%COMM, FNLSIF, 9)

    FNLSIF=FNLSIF/DIS

# 2118

!=======================================================================
! recalculate the projection operators
! (the array was used as a workspace)
!=======================================================================
    IZERO=1
    CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,  IZERO)

    DEALLOCATE(CWORK)

    

    RETURN
  END SUBROUTINE STRENL


END MODULE nonl

!************************ SUBROUTINE VNLAC0 ****************************
!
!> Calculates a linear combination of projection operators in reciprocal
!> space
!>
!> The result is *added* to CACC
!
!***********************************************************************

  SUBROUTINE VNLAC0(NONL_S,WDES1,CPROJ_LOC,CACC)

!! USE moffload

    USE nonl_struct_def
    USE wave_struct_def

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    COMPLEX(q)    CPROJ_LOC(WDES1%NPRO)
    COMPLEX(q) CACC(WDES1%NRPLWV)

# 2164

    CALL VNLAC0_(NONL_S,WDES1,CPROJ_LOC,CACC)
    RETURN
  END SUBROUTINE VNLAC0

  SUBROUTINE VNLAC0_(NONL_S,WDES1,CPROJ_LOC,CACC)
    USE nonl_struct_def
    USE wave_struct_def

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    COMPLEX(q)    CPROJ_LOC(WDES1%NPRO)
    COMPLEX(q) CACC(WDES1%NRPLWV)
! local
    INTEGER :: NPL, LMBASE, ISPIRAL, ISPINOR, NT, LMMAXC, NI, LM, K
    COMPLEX(q) :: CTMP

    REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
    COMPLEX(q)    :: CPROJ(WDES1%NPRO_TOT)
# 2190

    

    IF (WDES1%NK /= NONL_S%NK) THEN
       WRITE(*,*) 'internal error in VNLAC0: PHASE not properly set up',WDES1%NK, NONL_S%NK
    ENDIF

    NPL=WDES1%NGVECTOR

! merge projected wavefunctions from all nodes
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
!=======================================================================
! performe loop over ions
!=======================================================================
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NI,NT,LMMAXC,LMBASE,WORK,CTMP,TMP,LM) REDUCTION(+:CACC)
       ion: DO NI=1,NONL_S%NIONS
          NT=NONL_S%ITYP(NI)
          LMMAXC=NONL_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion
          LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)
          DO LM=1,LMMAXC
             CTMP= CPROJ(LMBASE+LM)*CONJG(NONL_S%CQFAK(LM,NT))
             TMP(LM,1)= REAL( CTMP ,KIND=q)
             TMP(LM,2)= AIMAG(CTMP)
          ENDDO
# 2235

          

          CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, TMP(1,1) , 1 , 0._q , WORK(1), 1)
          CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, TMP(1,2) , 1 , 0._q , WORK(1+NPL), 1)

          

!=======================================================================
! add acceleration from this ion to CACC
!=======================================================================
          CALL WORK_MUL_CREXP( NPL, WORK(1), NONL_S%CREXP(1,NI), CACC(1+NPL*ISPINOR))

       ENDDO ion
    ENDDO spinor

# 2255


    

    RETURN
  END SUBROUTINE VNLAC0_

# 2321


!************************ SUBROUTINE VNLAC0_DER *************************
!
!> Calculates the non local contribution of the first derivative
!> of the Hamiltonian upon moving (1._q,0._q) ion ION in the direction IDIR
!>
!> More specifically the contribution
!> ~~~
!>    | d p_i/ d R > c_i
!> ~~~
!> For IDIR=0 the derivative is replaced by | p_i >
!>
!> A term (i k) is included as well to "center" the derivate at the
!> ion
!
!***********************************************************************

  SUBROUTINE VNLAC0_DER(NONL_S,WDES1,CPROJ_LOC,CACC, LATT_CUR, ION, IDIR)

!! USE moffload_struct_def

    USE nonl_struct_def
    USE wave_struct_def
    USE lattice
    USE constant
    USE string, ONLY: str
    USE tutor, ONLY: vtutor

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    COMPLEX(q) CACC(WDES1%NRPLWV)
    COMPLEX(q)       CPROJ_LOC(WDES1%NPRO)
    TYPE (latt)  LATT_CUR
    INTEGER      ION, IDIR

! work arrays
    COMPLEX(q)    :: CPROJ(WDES1%NPRO_TOT)
    INTEGER :: NPL, LMBASE, LMMAXC, ISPINOR,  ISPIRAL, K, KK, NI, NIS, NT, LM, L
    COMPLEX(q) CMUL
    COMPLEX(q) CACC_ADD(WDES1%NRPLWV)

    

    IF (WDES1%NK /= NONL_S%NK) &
       CALL vtutor%bug("VNLAC0_DER: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONL_S%NK),"nonl.F",2369)

    NPL=WDES1%NGVECTOR

!$ACC ENTER DATA CREATE(CPROJ,CACC_ADD) 
! merge projected wavefunctions from all nodes
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))
!=======================================================================
! performe loop over ions
!=======================================================================
    LMBASE= 0
    ISPIRAL = 1

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       NIS=1
       typ: DO NT=1,NONL_S%NTYP
          LMMAXC=NONL_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ions: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
!note(ar): the following IF is cuts down all the preceeding loops to a single pass (no loops at all)
             IF (NI==ION) THEN
!$ACC KERNELS PRESENT(CACC_ADD) 
                CACC_ADD=0
!$ACC END KERNELS
!$ACC PARALLEL LOOP SEQ PRIVATE(CMUL) PRESENT(CACC_ADD,NONL_S,CPROJ) 
                l_loop: DO LM=1,LMMAXC
                   CMUL=CPROJ(LMBASE+LM)*CONJG(NONL_S%CQFAK(LM,NT))
!$ACC LOOP VECTOR GANG
!DIR$ IVDEP
!OCL NOVREC
                   DO K=1,NPL
                      CACC_ADD(K)=CACC_ADD(K)+NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)*CMUL* &
                           (CONJG(NONL_S%CREXP(K,NI)))
                   ENDDO
                ENDDO l_loop
!=======================================================================
! add acceleration from this ion
!=======================================================================

                IF (IDIR==0) THEN
!$ACC PARALLEL LOOP PRESENT(CACC,CACC_ADD) PRIVATE(KK) 
!DIR$ IVDEP
!OCL NOVREC
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CACC(KK)=CACC(KK)+CACC_ADD(K)
                   ENDDO
                ELSE
!$ACC PARALLEL LOOP PRESENT(CACC,CACC_ADD,WDES1,LATT_CUR) PRIVATE(KK) 
!DIR$ IVDEP
!OCL NOVREC
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CACC(KK)= CACC(KK)+CACC_ADD(K)*CONJG(CITPI)* &

                          ((WDES1%IGX(K)+WDES1%VKPT(1))*LATT_CUR%B(IDIR,1)+ &
                           (WDES1%IGY(K)+WDES1%VKPT(2))*LATT_CUR%B(IDIR,2)+ &
                           (WDES1%IGZ(K)+WDES1%VKPT(3))*LATT_CUR%B(IDIR,3))
# 2433

                   ENDDO
                ENDIF

             ENDIF
             LMBASE= LMMAXC+LMBASE
          ENDDO ions

600       NIS = NIS+NONL_S%NITYP(NT)
       ENDDO typ
       IF (NONL_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$ACC EXIT DATA DELETE(CPROJ,CACC_ADD) 

    

    RETURN
  END SUBROUTINE VNLAC0_DER


!***********************************************************************
!
!> small f77 helper routine to multiply with phase factor and divide
!> into real and imaginary part
!
!***********************************************************************

  SUBROUTINE CREXP_MUL_WAVE( NPL, CREXP, CPTWFP, WORK)
!$ACC ROUTINE VECTOR
    USE prec
    IMPLICIT NONE
    INTEGER NPL
    COMPLEX(q) :: CREXP(NPL), CPTWFP(NPL)
    REAL(q)    :: WORK(2*NPL)
    COMPLEX(q) :: CTMP
! local
    INTEGER K

 

!$ACC LOOP VECTOR
 !$OMP PARALLEL DO SIMD SCHEDULE(static) DEFAULT(SHARED) PRIVATE(CTMP)
    DO K=1,NPL
       CTMP=    CREXP(K) *CPTWFP(K)
       WORK(K)    = REAL( CTMP ,KIND=q)
       WORK(K+NPL)= AIMAG(CTMP)
    ENDDO
 !$OMP END PARALLEL DO SIMD

 

  END SUBROUTINE CREXP_MUL_WAVE


  SUBROUTINE WORK_MUL_CREXP( NPL, WORK, CREXP, CPTWFP)
!$ACC ROUTINE VECTOR
    USE prec

!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

    IMPLICIT NONE

    INTEGER NPL
    COMPLEX(q) :: CREXP(NPL), CPTWFP(NPL)
    REAL(q)    :: WORK(2*NPL)
    COMPLEX(q) :: CTMP
! local
    INTEGER K

 

!$ACC LOOP VECTOR PRIVATE(CTMP)
 !$OMP PARALLEL DO SIMD SCHEDULE(static) DEFAULT(SHARED) PRIVATE(CTMP)
    DO K=1,NPL

       CPTWFP(K)=CPTWFP(K)+ CMPLX( WORK(K) , WORK(K+NPL) ,KIND=q) *CONJG(CREXP(K))
# 2512

    ENDDO
 !$OMP END PARALLEL DO SIMD

 

  END SUBROUTINE WORK_MUL_CREXP
