#include "symbol.inc"
!******************** SUBROUTINE FEXCG ****************************
!
!> This subroutine calculates derivatives of the electron density
!> and calls the subroutines that calculate xc quantities
!> (GGA, meta-GGA, nonlocal van der Waals, dDsC) for the energy,
!> potential, forces, stress tensor, etc.
!
!******************************************************************

      SUBROUTINE FEXCG(&
     &   XC,NCDIJ,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
     &   KINEDEN,CHTOT,DENCOR,CWORK,DMUWORK,XDM,IVDW,NPDDSC,NPMU)

GPACC USE moffload

      USE prec
      USE base, only : CAST_COMPLEX_TO_RGRID
      USE lattice
      USE mpimy
      USE mgrid
      USE constant
      USE xc_driver
      USE mbj, ONLY : CMBJ, GET_CMBJ_PW, LscMBJ, SUM_GRHO_OVER_RHO_PW, UPDATE_CMBJ_AUX
      USE tau_mu, ONLY : tau_handle
      USE setexm_struct_def, ONLY : xc_info
      USE xc_name
      USE vdw_ll, ONLY : VDW_NONLOC, VDW_NONLOC_SPIN
      USE VDW_rVV10, ONLY : VDW_NONLOC_RVV10
      USE ggalib, ONLY : ggaapprox

      IMPLICIT NONE

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      TYPE (tau_handle)  KINEDEN
      TYPE (xc_info)     XC

      INTEGER :: IVDW,NPDDSC,NCDIJ,I,IDIR,JDIR,IDUMP,IONODE,ISP,NPMU,IDH
      INTEGER :: N1,N2,N3,NC,NODE_ME
      INTEGER :: NP1,NP4,NPVDW

      REAL(q) XCSIF(3,3),STRESS(3,3),XDM(NPDDSC)
      REAL(q) MAG_NORM,MAG_NORM_I,NABMAG(3),TAU_NORM,TAU_DOT_MDIR,LAP_NORM,LAP_DOT_MDIR
      REAL(q) LAPLUP,LAPLDW,XCENC,EXC,ABSNAB,ABSNABUP,ABSNABDW
      REAL(q) ANAB1U,ANAB2U,ANAB3U,ANAB1D,ANAB2D,ANAB3D
      REAL(q) DEXC1,DEXC2,DVC_,DVXC1,DVXC2,EVTOH,EXCL,XCENCC,VXC,VXC1,VXC2
      REAL(q) GG,GX,GY,GZ,GI,GJ
      REAL(q) MUUP,MUDW,EDCMU
      REAL(q) DLAPU,DLAPD
      REAL(q) CVMBJ
      REAL(q) tXDM
      REAL(q) RHO,RHO1,RHO2,RINPL
      REAL(q) SIF11,SIF12,SIF22,SIF23,SIF31,SIF33
      REAL(q) SIFLAP11,SIFLAP12,SIFLAP22,SIFLAP23,SIFLAP31,SIFLAP33
      REAL(q) TAUD,TAUU
      REAL(q) :: CHGMIN=1E-10_q, NABMIN=1E-10_q, TAUMIN=1E-10_q, LAPMIN=1E-10_q
      REAL(q), ALLOCATABLE :: DWORKG(:,:),DWORK1(:,:),DWORK2(:,:),DWORK3(:,:),DVC(:),LAPLACIAN(:,:)
      REAL(q), ALLOCATABLE :: DWORK4(:),DWORK4S(:),DWORK6(:),DWORK6S(:),DWORK7(:),DWORK7S(:)
      REAL(q), ALLOCATABLE :: DWORKL(:,:),DENSHESS(:,:,:)
      REAL(q) DVXCFAC

      COMPLEX(q) CVZERO
      COMPLEX(q), ALLOCATABLE, TARGET :: CWGRAD(:,:),CKINEDEN(:,:)
      COMPLEX(q), TARGET :: CHTOT(GRIDC%MPLWV,NCDIJ),CWORK(GRIDC%MPLWV,NCDIJ)

      RGRID DENCOR(GRIDC%RL%NP)
      RGRID DMUWORK(DIMREAL(NPMU),NCDIJ)
      RGRID, ALLOCATABLE :: DWORK5(:),DWORK5S(:),TAUC(:)
      RGRID, POINTER :: DWGRAD(:,:) => NULL()
      RGRID, POINTER :: DHTOT(:,:) => NULL()
      RGRID, POINTER :: DWORK(:,:) => NULL()
      RGRID, POINTER :: DKINEDEN(:,:) => NULL()

      LOGICAL LCHGMINSMALL

!$ACC ROUTINE (GGA_DRIVER) SEQ
!$ACC ROUTINE (GGASPIN_DRIVER) SEQ
!$ACC ROUTINE (GGAAPPROX) SEQ
!$ACC ROUTINE (METAGGA_DRIVER) SEQ
!$ACC ROUTINE (LscMBJ) SEQ
!$ACC ROUTINE (GET_CMBJ_PW) SEQ

      IF (XC%LDOMETAGGA) THEN
         DVXCFAC=1.0_q
      ELSE
         DVXCFAC=0.5_q
      ENDIF

      NP1=GRIDC%RL%NP
      NP4=1
      NPVDW=1
      IF ((NCDIJ==1).OR.XC%LUSE_VDW) NP4=NP1
      IF (XC%LUSE_VDW) NPVDW=NP1

      ALLOCATE(CWGRAD(GRIDC%MPLWV,NCDIJ),DWORKG(NP1,MIN(NCDIJ,2)), &
        DWORK1(NP1,NCDIJ),DWORK2(NP1,NCDIJ),DWORK3(NP1,NCDIJ),DVC(NP1), &
        DWORK4(NP4),DWORK5(NPVDW),DWORK6(NPVDW),DWORK7(NPVDW))
!$ACC ENTER DATA CREATE(CWGRAD,DWORKG,DWORK1,DWORK2,DWORK3,DVC,DWORK4,DWORK5,DWORK6,DWORK7) __IF_ASYNC__

      IF (XC%LSPIN_VDW) THEN
         ALLOCATE(DWORK4S(NP4),DWORK5S(NPVDW),DWORK6S(NPVDW),DWORK7S(NPVDW))
!$ACC ENTER DATA CREATE(DWORK4S,DWORK5S,DWORK6S,DWORK7S) __IF_ASYNC__
      ENDIF

      IF (XC%LTAU) THEN
         ALLOCATE(CKINEDEN(GRIDC%MPLWV,NCDIJ),TAUC(NP1))
!$ACC ENTER DATA CREATE(KINEDEN,CKINEDEN,TAUC) __IF_ASYNC__
!$ACC KERNELS PRESENT(CKINEDEN,TAUC) __IF_ASYNC__
         CKINEDEN=KINEDEN%TAU
         TAUC=KINEDEN%TAUC
!$ACC END KERNELS
      ENDIF

      IF (XC%LLAP) THEN
         ALLOCATE(LAPLACIAN(NP1,NCDIJ))
!$ACC ENTER DATA CREATE(LAPLACIAN) __IF_ASYNC__
!$ACC KERNELS PRESENT(LAPLACIAN) __IF_ASYNC__
         LAPLACIAN=0._q
!$ACC END KERNELS
      ENDIF

      IF (XC%LDLAP) THEN
         ALLOCATE(DWORKL(GRIDC%RL%NP,NCDIJ),DENSHESS(GRIDC%RL%NP,6,NCDIJ))
!$ACC ENTER DATA CREATE(DWORKL,DENSHESS) __IF_ASYNC__
!$ACC KERNELS PRESENT(DWORKL,DENSHESS) __IF_ASYNC__
         DWORKL=0
         DENSHESS=0
!$ACC END KERNELS
      ENDIF

      CALL CAST_COMPLEX_TO_RGRID(CWGRAD,DWGRAD)
      CALL CAST_COMPLEX_TO_RGRID(CHTOT,DHTOT)
      CALL CAST_COMPLEX_TO_RGRID(CWORK,DWORK)
      IF (XC%LTAU) CALL CAST_COMPLEX_TO_RGRID(CKINEDEN,DKINEDEN)
!$ACC ENTER DATA CREATE(DWORK) __IF_ASYNC__

      ! Flip CHTOT from (rho,mag) to (rho_up,rho_down)
      IF ((NCDIJ==1).OR.((NCDIJ==2).AND.XC%LDOMETAGGA)) CALL RL_FLIP(CHTOT,GRIDC,NCDIJ,.TRUE.)

      ! Non-collinear magnetism: gradient corrections in the
      ! noncollinear case are calculated a bit differently than in the collinear case
      ! CHTOT and KINEDEN should be in (rho,mag) representation

!     TAUMIN=0 for MBJ is ok since there is no division by TAU
      IF ((XC%ID(1)==ID_XC_MBJ) .OR. (XC%ID(1)==ID_XC_LMBJ)) TAUMIN=0._q

! set to one for error-dumps
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
      IDUMP=0
#ifdef debug
      IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
      NODE_ME=0
      IONODE =0
      IDUMP=0
#ifdef debug
      IDUMP=1
#endif
#endif

#ifdef vector
      IF ((XC%NXC==1).AND.(XC%ID(1)==ID_XC_PW91)) THEN
        CALL vtutor%error("ERROR: only PW91 is supported if -Dvector is specified in the makefile \n &
           &remove -Dvector from the makefile")
      ENDIF
#endif
!
! only for the analytical PBE functional we can use a low charge density threshold
! since numerical differentiation leads to errors for the other functionals
!
      LCHGMINSMALL=.TRUE.
      DO I=1,XC%NXC
         IF ((XC%ID(I)/=ID_XC_PBE) .AND. &
           & (XC%ID(I)/=ID_XC_PBE_X) .AND. &
           & (XC%ID(I)/=ID_XC_PBE_C) .AND. &
           & (XC%ID(I)/=ID_XC_PBETEST) .AND. &
           & (XC%ID(I)/=ID_XC_RPBE) .AND. &
           & (XC%ID(I)/=ID_XC_PBESOL) .AND. &
           & (XC%ID(I)/=ID_XC_BEEF) .AND. &
           & (XC%ID(I)/=ID_XC_REVPBE) .AND. &
           & (XC%ID(I)/=ID_XC_OPTPBE) .AND. &
           & (XC%ID(I)/=ID_XC_OPTB88PBE) .AND. &
           & (XC%ID(I)/=ID_XC_OPTB86BPBE) .AND. &
           & (XC%ID(I)/=ID_XC_PW86RPBE) .AND. &
           & (XC%ID(I)/=ID_XC_CX)) LCHGMINSMALL=.FALSE.
      ENDDO
      IF (LCHGMINSMALL) THEN
         CHGMIN=1E-10
      ELSE
         CHGMIN=1E-7
      ENDIF

      RINPL=1._q/GRIDC%NPLWV                    ! Scaling of Energy
      EVTOH=1._q/(2.*HSQDTM)*AUTOA5             ! KinEDens eV to Hartree

!=======================================================================
! First phase: Transform DENCOR (core charge) and
!  CHTOT (pseudo chargedensity) to real space
!=======================================================================
      spin: DO ISP=1,NCDIJ
      IF ((NCDIJ==1).OR.(NCDIJ==2)) THEN
! set CWORK to total real charge in reciprocal space
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK,DENCOR,DHTOT,RINPL,LATT_CUR,NCDIJ,ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DENCOR,DHTOT,LATT_CUR) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            DWORK(I,ISP)=(DENCOR(I)/NCDIJ+DHTOT(I,ISP))*RINPL/LATT_CUR%OMEGA
         ENDDO
NOACC !$OMP END PARALLEL DO
      ELSEIF (NCDIJ==4) THEN
         IF (ISP==1) THEN
! Add core charge density to pseudo charge density
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK,DENCOR,DHTOT,RINPL,LATT_CUR,ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DENCOR,DHTOT,LATT_CUR) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
               DWORK(I,ISP)=(DENCOR(I)+DHTOT(I,ISP))*RINPL/LATT_CUR%OMEGA
            ENDDO
NOACC !$OMP END PARALLEL DO
         ELSE
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK,DHTOT,RINPL,LATT_CUR,ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DHTOT,LATT_CUR) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
               DWORK(I,ISP)=DHTOT(I,ISP)*RINPL/LATT_CUR%OMEGA
            ENDDO
NOACC !$OMP END PARALLEL DO
         ENDIF
      ENDIF
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3D(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

!=======================================================================
! now calculate the gradient of the charge density
!=======================================================================
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,CWGRAD,CWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(CWGRAD,CWORK) __IF_ASYNC__
      DO  I=1,GRIDC%RC%NP
        CWGRAD(I,ISP)=CWORK(I,ISP)
      ENDDO
NOACC !$OMP END PARALLEL DO
! x-component:
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GX) SHARED(GRIDC,LATT_CUR,CWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX) __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         CWORK(I,ISP)=CWORK(I,ISP)*GX*CITPI
       ENDDO
NOACC !$OMP END PARALLEL DO
! grad_x in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3D(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK1,DWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK1,DWORK) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
        DWORK1(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
NOACC !$OMP END PARALLEL DO

! y-component:
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,CWORK,CWGRAD,ISP)
!$ACC PARALLEL LOOP PRESENT(CWGRAD,CWORK) __IF_ASYNC__
      DO  I=1,GRIDC%RC%NP
        CWORK(I,ISP)=CWGRAD(I,ISP)
      ENDDO
NOACC !$OMP END PARALLEL DO
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GY) SHARED(GRIDC,LATT_CUR,CWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,LATT_CUR) PRIVATE(N1,N2,N3,NC,GY) __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         CWORK(I,ISP)=CWORK(I,ISP)*GY*CITPI
      ENDDO
NOACC !$OMP END PARALLEL DO
! grad_y in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3D(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK2,DWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK2,DWORK) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
        DWORK2(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
NOACC !$OMP END PARALLEL DO

! z-component:
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,CWORK,CWGRAD,ISP)
!$ACC PARALLEL LOOP PRESENT(CWGRAD,CWORK) __IF_ASYNC__
      DO  I=1,GRIDC%RC%NP
        CWORK(I,ISP)=CWGRAD(I,ISP)
      ENDDO
NOACC !$OMP END PARALLEL DO
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GZ) SHARED(GRIDC,LATT_CUR,CWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,LATT_CUR) PRIVATE(N1,N2,N3,NC,GZ) __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
         CWORK(I,ISP)=CWORK(I,ISP)*GZ*CITPI
      ENDDO
NOACC !$OMP END PARALLEL DO
! grad_z in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3D(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK3,DWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK3,DWORK) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
        DWORK3(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
NOACC !$OMP END PARALLEL DO

      IF (XC%LLAP) THEN
        IF (XC%LDLAP) THEN

          ! and if laplacian-dependent functional is called, also compute density hessian
          ! indices for density hessian: 1 = xx, 2 = xy, 3 = xz, 4 = yy, 5 = yz, 6 = zz
          IDH = 0
          DO IDIR=1,3
            DO JDIR = IDIR,3
              IDH = IDH + 1
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GI, GJ) SHARED(GRIDC,LATT_CUR,CWORK,CWGRAD,ISP,IDIR,JDIR)
!$ACC PARALLEL LOOP PRESENT(GRIDC,LATT_CUR,CWGRAD,CWORK) PRIVATE(N1,N2,N3,NC,GI,GJ) __IF_ASYNC__
              DO I=1,GRIDC%RC%NP
                N1= MOD((I-1),GRIDC%RC%NROW) +1
                NC= (I-1)/GRIDC%RC%NROW+1
                N2= GRIDC%RC%I2(NC)
                N3= GRIDC%RC%I3(NC)

                GI = GRIDC%LPCTX(N1)*LATT_CUR%B(IDIR,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(IDIR,2)&
                &      +GRIDC%LPCTZ(N3)*LATT_CUR%B(IDIR,3)
                GJ = GRIDC%LPCTX(N1)*LATT_CUR%B(JDIR,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(JDIR,2)&
                &      +GRIDC%LPCTZ(N3)*LATT_CUR%B(JDIR,3)
                CWORK(I,ISP) = CWGRAD(I,ISP) * GI*GJ*CITPI*CITPI
              ENDDO
NOACC !$OMP END PARALLEL DO
              CALL SETUNB(CWORK(1,ISP),GRIDC)
              CALL FFT3D(CWORK(1,ISP),GRIDC,1)
              CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DENSHESS,DWORK,ISP,IDH)
!$ACC PARALLEL LOOP PRESENT(DENSHESS,DWORK) __IF_ASYNC__
              DO I=1,GRIDC%RL%NP
                 DENSHESS(I,IDH,ISP) = REAL( DWORK(I,ISP) ,KIND=q)
              ENDDO
NOACC !$OMP END PARALLEL DO
            ENDDO
          ENDDO
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,LAPLACIAN,DENSHESS,ISP)
!$ACC PARALLEL LOOP PRESENT(LAPLACIAN,DENSHESS) __IF_ASYNC__
          DO I=1,GRIDC%RL%NP
            LAPLACIAN(I,ISP) = DENSHESS(I,1,ISP) + DENSHESS(I,4,ISP) + DENSHESS(I,6,ISP)
          ENDDO
NOACC !$OMP END PARALLEL DO
        ELSE

          DO IDIR=1,3
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GG) SHARED(GRIDC,LATT_CUR,CWORK,CWGRAD,ISP,IDIR)
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,LATT_CUR,CWGRAD) PRIVATE(N1,N2,N3,NC,GG) __IF_ASYNC__
            DO I=1,GRIDC%RC%NP  ! loop over all grid points NP in the reciprocal (RC) grid
               ! index of grid point along three reciprocal lattice vectors (N1, N2, N3)
               N1= MOD((I-1),GRIDC%RC%NROW) +1
               NC= (I-1)/GRIDC%RC%NROW+1
               N2= GRIDC%RC%I2(NC)
               N3= GRIDC%RC%I3(NC)
               ! convert to lattice vector component in direction x, y or z (corresponding to index J)
               GG=(GRIDC%LPCTX(N1)*LATT_CUR%B(IDIR,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(IDIR,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(IDIR,3))
               CWORK(I,ISP)=CWGRAD(I,ISP)*GG*GG*CITPI*CITPI
               !write(*,*)'xcgrad:cwgrad ', CWORK4(I), CWGRAD(I), G1, G2, CITPI
            ENDDO
NOACC !$OMP END PARALLEL DO
            CALL SETUNB(CWORK(1,ISP),GRIDC)
            CALL FFT3D(CWORK(1,ISP),GRIDC,1)
            CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,LAPLACIAN,DWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(LAPLACIAN,DWORK) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
               LAPLACIAN(I,ISP)=LAPLACIAN(I,ISP)+REAL( DWORK(I,ISP) ,KIND=q)
            ENDDO
NOACC !$OMP END PARALLEL DO
          ENDDO

        ENDIF
      ENDIF

      ENDDO spin


!=======================================================================
!  grad rho    d    f_xc
! ---------- * ------------      (Phys.Rev.B 50,7 (1994) 4954)
! |grad rho|   d |grad rho|
!
!  MIND: the factor OMEGA is difficult to understand:
!   1/N sum_r energy_density * rho *OMEGA = Energy
!   1/N sum_r energy_density * \bar rho   = Energy (\bar rho=rho*LATT_CUR%OMEGA)
!=======================================================================

!$ACC KERNELS PRESENT(DWORK,DVC,DWORKG) __IF_ASYNC__
      DWORK = 0._q ; DVC = 0._q ; DWORKG = 0._q
!$ACC END KERNELS

      EXC=0._q; EDCMU=0._q; CVMBJ=CMBJ

NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP SHARED(GRIDC,NCDIJ,CHGMIN,NABMIN,LAPMIN,TAUMIN,IVDW,XDM,DMUWORK,DVXCFAC,DHTOT,DENCOR,LATT_CUR,XC,EVTOH,DWORK1,DWORK2,DWORK3,DWORK4,DWORK4S,DWORK7,DWORK7S,TAUC,DKINEDEN,DWORK,DWORKG,DVC,LAPLACIAN,DWORKL) &
NOACC !$OMP PRIVATE(I,MAG_NORM,MAG_NORM_I,RHO1,RHO2,NABMAG,ABSNABUP,ABSNABDW,ABSNAB,TAU_NORM,TAU_DOT_MDIR,TAUU,TAUD,LAP_NORM,LAP_DOT_MDIR,LAPLUP,LAPLDW,EXCL,DEXC1,DEXC2,DVXC1,DVXC2,DVC_,DLAPU,DLAPD,MUUP,MUDW,RHO,TXDM) &
NOACC !$OMP FIRSTPRIVATE(CVMBJ) &
NOACC !$OMP REDUCTION(+:EXC,EDCMU)
!$ACC PARALLEL LOOP PRESENT(DHTOT,DENCOR,LATT_CUR,XC,DWORK1,DWORK2,DWORK3,DWORK4,DWORK4S,DWORK7,DWORK7S,TAUC,DKINEDEN,DWORK,DWORKG,DVC,LAPLACIAN,DWORKL) &
!$ACC PRIVATE(MAG_NORM,MAG_NORM_I,RHO1,RHO2,NABMAG,ABSNABUP,ABSNABDW,ABSNAB,TAU_NORM,TAU_DOT_MDIR,TAUU,TAUD, &
!$ACC LAP_NORM,LAP_DOT_MDIR,LAPLUP,LAPLDW,EXCL,DEXC1,DEXC2,DVXC1,DVXC2,DVC_,DLAPU,DLAPD,MUUP,MUDW,RHO,TXDM) &
!$ACC FIRSTPRIVATE(CVMBJ) REDUCTION(+:EXC,EDCMU) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP

         IF (NCDIJ==1) THEN
         ! nonspin-polarized case
            RHO1=MAX(REAL((DHTOT(I,1)+DENCOR(I))*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
            RHO2=RHO1

            ABSNAB=SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1) &
                 +DWORK3(I,1)*DWORK3(I,1))

            ABSNABUP= 0.5_q*ABSNAB
            ABSNABDW= 0.5_q*ABSNAB

!           not needed for MBJ
            IF ((XC%ID(1)/=ID_XC_MBJ) .AND. (XC%ID(1)/=ID_XC_LMBJ)) THEN
               ABSNABUP=MAX(ABSNABUP,NABMIN)
               ABSNABDW=MAX(ABSNABDW,NABMIN)
               ABSNAB=MAX(ABSNAB,2*NABMIN)
            ENDIF

            IF (XC%LLAP) THEN
               LAPLUP=0.5_q*REAL(LAPLACIAN(I,1),KIND=q)
            ELSE
               LAPLUP=0._q
            ENDIF
            LAPLDW=LAPLUP

            IF (XC%LTAU) THEN
               TAUU=MAX(0.5_q*REAL(DKINEDEN(I,1)+TAUC(I),KIND=q), TAUMIN)
            ELSE
               TAUU=0._q
            ENDIF
            TAUD=TAUU

            DWORK4(I)=2.0_q*RHO1

            IF (XC%LUSE_VDW) DWORK7(I)=ABSNAB

            ! DDsC dispersion correction
            IF (IVDW==4) CALL GGAapprox(ABS(RHO1)*AUTOA3,ABSNABUP*AUTOA4,XDM(I))

         ELSEIF (NCDIJ==2) THEN
         ! Spin-polarized case
            RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)*0.5_q)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
            RHO2= MAX(REAL((DHTOT(I,2)+DENCOR(I)*0.5_q)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)

            ABSNABUP=SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1) &
                          +DWORK3(I,1)*DWORK3(I,1))

            ABSNABDW=SQRT(DWORK1(I,2)*DWORK1(I,2)+DWORK2(I,2)*DWORK2(I,2) &
                          +DWORK3(I,2)*DWORK3(I,2))

            ABSNAB= (DWORK1(I,1)+DWORK1(I,2))*(DWORK1(I,1)+DWORK1(I,2))+ &
                    (DWORK2(I,1)+DWORK2(I,2))*(DWORK2(I,1)+DWORK2(I,2))+ &
                    (DWORK3(I,1)+DWORK3(I,2))*(DWORK3(I,1)+DWORK3(I,2))
            ABSNAB=SQRT(ABSNAB)

!           not needed for MBJ
            IF ((XC%ID(1)/=ID_XC_MBJ) .AND. (XC%ID(1)/=ID_XC_LMBJ)) THEN
               ABSNABUP=MAX(ABSNABUP,NABMIN)
               ABSNABDW=MAX(ABSNABDW,NABMIN)
               ABSNAB=MAX(ABSNAB,2*NABMIN)
            ENDIF
!
! presently VASP converges very badly if correlation_ABS_DRHOUP_ABS_DRHOD
! is not set
! it is probably the exchange potential, which is partly canceled by
! the correlation energy if ABSNAB=ABSNABUP+ABSNABDW
! so let us stick to ABSNAB=ABSNABUP+ABSNABDW for the time beeing
!
!#define  correlation_ABS_DRHOUP_ABS_DRHOD
#ifdef correlation_ABS_DRHOUP_ABS_DRHOD
            ABSNAB=ABSNABUP+ABSNABDW
#endif
            IF (XC%LLAP) THEN
               LAPLUP=REAL(LAPLACIAN(I,1),KIND=q)
               LAPLDW=REAL(LAPLACIAN(I,2),KIND=q)
            ELSE
               LAPLUP=0._q
               LAPLDW=0._q
            ENDIF

            IF (XC%LTAU) THEN
               TAUU=MAX(REAL(DKINEDEN(I,1)+TAUC(I)*0.5_q,KIND=q), TAUMIN)
               TAUD=MAX(REAL(DKINEDEN(I,2)+TAUC(I)*0.5_q,KIND=q), TAUMIN)
            ELSE
               TAUU=0._q
               TAUD=0._q
            ENDIF

            IF (XC%LUSE_VDW) THEN
               IF (XC%LSPIN_VDW) THEN
                  DWORK4(I)=RHO1
                  DWORK4S(I)=RHO2
                  DWORK7(I)=ABSNABUP
                  DWORK7S(I)=ABSNABDW
               ELSE
                  DWORK4(I)=RHO1+RHO2
                  DWORK7(I)=ABSNAB
               ENDIF
            ENDIF

            ! DDsC dispersion correction
            IF (IVDW.eq.4) then
               CALL GGAapprox(ABS(RHO1)*AUTOA3,ABSNABUP*AUTOA4,XDM(I))
               CALL GGAapprox(ABS(RHO2)*AUTOA3,ABSNABDW*AUTOA4,tXDM)
               XDM(I)=(XDM(I)+tXDM)*0.5_q
            ENDIF

         ELSEIF (NCDIJ==4) THEN
         ! non-collinear case
!
! | m |
!
            MAG_NORM=SQRT(ABS( &
           &    DHTOT(I,2)*DHTOT(I,2)+DHTOT(I,3)*DHTOT(I,3)+DHTOT(I,4)*DHTOT(I,4)))
            MAG_NORM_I=1._q/MAX(MAG_NORM,CHGMIN)
!
! RHO1: \rho_up   = ( \rho_tot + | m | )/2
! RHO2: \rho_down = ( \rho_tot - | m | )/2
!
            RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)+MAG_NORM)*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
            RHO2= MAX(REAL((DHTOT(I,1)+DENCOR(I)-MAG_NORM)*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
!
! \nabla | m |
!
            NABMAG(1)=(DWORK1(I,2)*DHTOT(I,2)+DWORK1(I,3)*DHTOT(I,3)+DWORK1(I,4)*DHTOT(I,4))*MAG_NORM_I
            NABMAG(2)=(DWORK2(I,2)*DHTOT(I,2)+DWORK2(I,3)*DHTOT(I,3)+DWORK2(I,4)*DHTOT(I,4))*MAG_NORM_I
            NABMAG(3)=(DWORK3(I,2)*DHTOT(I,2)+DWORK3(I,3)*DHTOT(I,3)+DWORK3(I,4)*DHTOT(I,4))*MAG_NORM_I
!
! | ( \nabla \rho + \nabla | m | )/2 |
!
            ABSNABUP=SQRT( &
           &            (DWORK1(I,1)+NABMAG(1))*(DWORK1(I,1)+NABMAG(1)) + &
           &             (DWORK2(I,1)+NABMAG(2))*(DWORK2(I,1)+NABMAG(2)) + &
           &              (DWORK3(I,1)+NABMAG(3))*(DWORK3(I,1)+NABMAG(3)) ) * 0.5_q
!
! | ( \nabla \rho - \nabla | m | )/2 |
!
            ABSNABDW=SQRT( &
           &            (DWORK1(I,1)-NABMAG(1))*(DWORK1(I,1)-NABMAG(1)) + &
           &             (DWORK2(I,1)-NABMAG(2))*(DWORK2(I,1)-NABMAG(2)) + &
           &              (DWORK3(I,1)-NABMAG(3))*(DWORK3(I,1)-NABMAG(3)) ) *0.5_q
!
! | \nabla \rho |
!
            ABSNAB=SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1)+DWORK3(I,1)*DWORK3(I,1))
!
! presently VASP converges very badly if correlation_ABS_DRHOUP_ABS_DRHOD
! is not set
! it is probably the exchange potential, which is partly canceled by
! the correlation energy if ABSNAB=ABSNABUP+ABSNABDW
! so let us stick to ABSNAB=ABSNABUP+ABSNABDW for the time beeing
!
!#define  correlation_ABS_DRHOUP_ABS_DRHOD
#ifdef correlation_ABS_DRHOUP_ABS_DRHOD
            ABSNAB=ABSNABUP+ABSNABDW
#endif
!
! Refill DWORK[1..3](:,2) with ( \nabla \rho - \nabla | m | )/2
!
            DWORK1(I,2)=(DWORK1(I,1)-NABMAG(1))*0.5_q
            DWORK2(I,2)=(DWORK2(I,1)-NABMAG(2))*0.5_q
            DWORK3(I,2)=(DWORK3(I,1)-NABMAG(3))*0.5_q
!
! Refill DWORK[1..3](:,1) with ( \nabla \rho + \nabla | m | )/2
!
            DWORK1(I,1)=(DWORK1(I,1)+NABMAG(1))*0.5_q
            DWORK2(I,1)=(DWORK2(I,1)+NABMAG(2))*0.5_q
            DWORK3(I,1)=(DWORK3(I,1)+NABMAG(3))*0.5_q

            IF (XC%LLAP) THEN
               !
               ! | m |
               !
!               LAP_NORM=MAX(SQRT(ABS( &
!              &    LAPLACIAN(I,2)*LAPLACIAN(I,2)+LAPLACIAN(I,3)*LAPLACIAN(I,3)+LAPLACIAN(I,4)*LAPLACIAN(I,4))), LAPMIN)
               LAP_DOT_MDIR=(DHTOT(I,2)*LAPLACIAN(I,2)+DHTOT(I,3)*LAPLACIAN(I,3)+DHTOT(I,4)*LAPLACIAN(I,4))*MAG_NORM_I
               !
               ! Laplacian of the density
!              LAPLUP=REAL((LAPLACIAN(I,1)+LAP_NORM)*0.5_q,KIND=q)
!              LAPLDW=REAL((LAPLACIAN(I,1)-LAP_NORM)*0.5_q,KIND=q)
               LAPLUP=REAL((LAPLACIAN(I,1)+LAP_DOT_MDIR)*0.5_q,KIND=q)
               LAPLDW=REAL((LAPLACIAN(I,1)-LAP_DOT_MDIR)*0.5_q,KIND=q)
            ELSE
               LAPLUP=0._q
               LAPLDW=0._q
            ENDIF

            IF (XC%LTAU) THEN
               !
               ! | m |
               !
!               TAU_NORM=MAX(SQRT(ABS( &
!              &    DKINEDEN(I,2)*DKINEDEN(I,2)+DKINEDEN(I,3)*DKINEDEN(I,3)+DKINEDEN(I,4)*DKINEDEN(I,4))), TAUMIN)
               TAU_DOT_MDIR=(DHTOT(I,2)*DKINEDEN(I,2)+DHTOT(I,3)*DKINEDEN(I,3)+DHTOT(I,4)*DKINEDEN(I,4))*MAG_NORM_I
               !
               ! kinetic energy density
               ! TAUU: \tau_up   = ( \tau_tot + | m | )/2
               ! TAUD: \tau_down = ( \tau_tot - | m | )/2
               !
!              TAUU=MAX(REAL((DKINEDEN(I,1)+TAU_NORM)*0.5_q,KIND=q), TAUMIN)
!              TAUD=MAX(REAL((DKINEDEN(I,1)-TAU_NORM)*0.5_q,KIND=q), TAUMIN)
               TAUU=MAX(REAL((DKINEDEN(I,1)+TAU_DOT_MDIR)*0.5_q+TAUC(I)*0.5_q,KIND=q), TAUMIN)
               TAUD=MAX(REAL((DKINEDEN(I,1)-TAU_DOT_MDIR)*0.5_q+TAUC(I)*0.5_q,KIND=q), TAUMIN)
            ELSE
               TAUU=0._q
               TAUD=0._q
            ENDIF
!
! store total charge and total gradient for non-local vdW-DF
            IF (XC%LUSE_VDW) THEN
               IF (XC%LSPIN_VDW) THEN
                  DWORK4(I)=RHO1
                  DWORK4S(I)=RHO2
                  DWORK7(I)=ABSNABUP
                  DWORK7S(I)=ABSNABDW
               ELSE
                  DWORK4(I)=RHO1+RHO2
                  DWORK7(I)=ABSNAB
               ENDIF
            ENDIF

            ! DDsC dispersion correction
            IF (IVDW.eq.4) then
               CALL GGAapprox(ABS(RHO1)*AUTOA3,ABSNABUP*AUTOA4,XDM(I))
               CALL GGAapprox(ABS(RHO2)*AUTOA3,ABSNABDW*AUTOA4,tXDM)
               XDM(I)=(XDM(I)+tXDM)*0.5_q
            ENDIF

         ENDIF

         IF (XC%LDOMETAGGA) THEN

            ! The CMBJ parameter may depend on the spatial position,
            ! if so, it is set here to the right entry in CMBJ_AUX
            CALL GET_CMBJ_PW(XC,I,CVMBJ)

            CALL METAGGA_DRIVER(&
           &   XC,RHO1*AUTOA3,RHO2*AUTOA3,ABSNABUP*AUTOA4,ABSNABDW*AUTOA4,ABSNAB*AUTOA4, &
           &   LAPLUP*AUTOA5,LAPLDW*AUTOA5,TAUU*EVTOH,TAUD*EVTOH, &
           &   EXCL,DEXC1,DEXC2,DVXC1,DVXC2,DVC_,MUUP,MUDW,DLAPU,DLAPD,CVMBJ)

            RHO=RHO1+RHO2

            IF (LscMBJ()) CALL SUM_GRHO_OVER_RHO_PW(RHO*AUTOA3,ABSNAB*AUTOA4)

            EXC=EXC+EXCL*RHO*RYTOEV*LATT_CUR%OMEGA
            DVXC1=DVXC1*RYTOEV*AUTOA
            DVXC2=DVXC2*RYTOEV*AUTOA
            DVC_ =DVC_ *RYTOEV*AUTOA
            IF (XC%LDLAP) THEN
              DLAPU=DLAPU*RYTOEV*AUTOA2
              DLAPD=DLAPD*RYTOEV*AUTOA2
            ENDIF

            IF (XC%LMU) THEN
               ! Double counting contribution stemming from dExc/d\mu
               EDCMU=EDCMU+0.5_q*(MUUP*TAUU+MUDW*TAUD-(MUUP+MUDW)*TAUC(I)*0.5_q)*LATT_CUR%OMEGA
               MUUP=0.5_q*MUUP*RYTOEV*AUTOA2 ! HSQDTM = (plancks CONSTANT/(2*PI))**2/(2*ELECTRON MASS)
               MUDW=0.5_q*MUDW*RYTOEV*AUTOA2
            ENDIF

#ifdef correlation_ABS_DRHOUP_ABS_DRHOD
            DVXC1=DVXC1+DVC_
            DVXC2=DVXC2+DVC_
            DVC_=0
#endif

            ! Store d f_xc / d \tau_up,down in DMUWORK if required
            IF (XC%LMU) THEN
               DMUWORK(I,1)=MUUP
               IF ((NCDIJ==2).OR.(NCDIJ==4)) DMUWORK(I,2)=MUDW
            ENDIF

            ! Store d f_x/ d (|\nabla \rho_up,down| ) / |\nabla \rho_up,down| in DWORK
            ! and  d f_xc/ d \rho_up,down  in DWORKG
            DWORK(I,1)  = DVXC1 / MAX(ABSNABUP, NABMIN)
            DWORKG(I,1) = DEXC1*RYTOEV
            IF (XC%LDLAP) DWORKL(I,1) = DLAPU
            IF ((NCDIJ==2).OR.(NCDIJ==4)) THEN
               DWORK(I,2)  = DVXC2 / MAX(ABSNABDW, NABMIN)
               DWORKG(I,2) = DEXC2*RYTOEV
               IF (XC%LDLAP) DWORKL(I,2) = DLAPD
            ENDIF
            ! Store d f_c/ d(|\nabla \rho|) / |\nabla \rho| in DVC
            DVC(I)=DVC_/MAX(ABSNAB, NABMIN)

         ENDIF

         IF (NCDIJ==1) THEN
            RHO=RHO1+RHO2
            CALL GGA_DRIVER(XC,RHO*AUTOA3,ABSNAB*AUTOA4,EXCL,DEXC1,DVXC1,.FALSE.)
            EXC=EXC+EXCL*RHO*RYTOEV*LATT_CUR%OMEGA
            DWORK(I,1)  = DWORK(I,1) + DVXCFAC*DVXC1*RYTOEV*AUTOA/MAX(ABSNABUP,NABMIN)
            DWORKG(I,1) = DWORKG(I,1) + DEXC1*RYTOEV

         ELSEIF ((NCDIJ==2).OR.(NCDIJ==4)) THEN
            CALL GGASPIN_DRIVER(XC,RHO1*AUTOA3,RHO2*AUTOA3, &
        &      ABSNABUP*AUTOA4,ABSNABDW*AUTOA4,ABSNAB*AUTOA4, &
        &      EXCL,DEXC1,DEXC2,DVXC1,DVXC2,DVC_,.FALSE.)

            RHO=RHO1+RHO2

            EXC=EXC+EXCL*RHO*RYTOEV*LATT_CUR%OMEGA
            DVXC1=DVXC1*RYTOEV*AUTOA
            DVXC2=DVXC2*RYTOEV*AUTOA
            DVC_ =DVC_ *RYTOEV*AUTOA
#ifdef correlation_ABS_DRHOUP_ABS_DRHOD
            DVXC1=DVXC1+DVC_
            DVXC2=DVXC2+DVC_
            DVC_=0
#endif

!
!   store d f/ d (|d rho| ) / |d rho|  in DWORK
!
            DWORK(I,1)  = DWORK(I,1) + DVXC1 / MAX(ABSNABUP,NABMIN)
            DWORK(I,2)  = DWORK(I,2) + DVXC2 / MAX(ABSNABDW,NABMIN)
            DVC(I)      = DVC(I) + DVC_  / MAX(ABSNAB,NABMIN)
!
!   store d f/ d rho  in DWORKG
!
            DWORKG(I,1) = DWORKG(I,1) + DEXC1*RYTOEV
            DWORKG(I,2) = DWORKG(I,2) + DEXC2*RYTOEV

         ENDIF

      ENDDO
NOACC !$OMP END PARALLEL DO

      IF (XC%LLAP) THEN
!$ACC EXIT DATA DELETE(LAPLACIAN) __IF_ASYNC__
         DEALLOCATE(LAPLACIAN)
      ENDIF
GPACC CALL VTUTOR_OFFLOAD_CHECK()

!=======================================================================
! non-local vdW-DF contributions
!=======================================================================
      IF (XC%LUSE_VDW) THEN
         IF (NCDIJ==1) THEN
            IF ((XC%IVDW_NL==1).OR.(XC%IVDW_NL==3).OR.(XC%IVDW_NL==4)) THEN
!$ACC KERNELS PRESENT(DWORK) __IF_ASYNC__
               DWORK = DWORK/2._q/DVXCFAC
!$ACC END KERNELS
               CALL VDW_NONLOC(DWORK4,DWORK(:,1),DWORKG(:,1),DWORK7,EXC,GRIDC,LATT_CUR,STRESS,1,XC%IVDW_NL,XC%ZAB_VDW,XC%GAMMA_VDW,XC%ALPHA_VDW)
!$ACC KERNELS PRESENT(DWORK) __IF_ASYNC__
               DWORK = DWORK*2._q*DVXCFAC
!$ACC END KERNELS
            ELSEIF (XC%IVDW_NL==2) THEN
!$ACC KERNELS PRESENT(DWORK) __IF_ASYNC__
               DWORK = DWORK/2._q/DVXCFAC
!$ACC END KERNELS
               CALL VDW_NONLOC_RVV10(DWORK4,DWORK(:,1),DWORKG(:,1),DWORK7,EXC,GRIDC,LATT_CUR,STRESS,1,XC%BPARAM,XC%CPARAM)
!$ACC KERNELS PRESENT(DWORK) __IF_ASYNC__
               DWORK = DWORK*2._q*DVXCFAC
!$ACC END KERNELS
            ENDIF
         ELSEIF ((NCDIJ==2).OR.(NCDIJ==4)) THEN
!$ACC KERNELS PRESENT(DWORK5,DWORK6,DWORK5S,DWORK6S) __IF_ASYNC__
            DWORK5=0._q
            DWORK6=0._q
            IF (XC%LSPIN_VDW) THEN
               DWORK5S=0._q
               DWORK6S=0._q
            ENDIF
!$ACC END KERNELS
            IF ((XC%IVDW_NL==1).OR.(XC%IVDW_NL==3).OR.(XC%IVDW_NL==4)) THEN
               IF (XC%LSPIN_VDW) THEN
                  CALL VDW_NONLOC_spin( DWORK4, DWORK4S, DWORK5, DWORK5S, &
                     & DWORK6, DWORK6S, DWORK7, DWORK7S, EXC,GRIDC,LATT_CUR,stress,XC%IVDW_NL,XC%ZAB_VDW,XC%GAMMA_VDW,XC%ALPHA_VDW)
               ELSE
                  CALL VDW_NONLOC(DWORK4,DWORK5,DWORK6,DWORK7,EXC,GRIDC,LATT_CUR,STRESS,2,XC%IVDW_NL,XC%ZAB_VDW,XC%GAMMA_VDW,XC%ALPHA_VDW)
               ENDIF
            ELSEIF (XC%IVDW_NL==2) THEN
               CALL VDW_NONLOC_RVV10(DWORK4,DWORK5,DWORK6,DWORK7,EXC,GRIDC,LATT_CUR,STRESS,2,XC%BPARAM,XC%CPARAM)
            ENDIF
            IF (XC%LSPIN_VDW) THEN
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK,DWORK5,DWORK5S,DWORK1,DWORK2,DWORK3)
!$ACC PARALLEL PRESENT(DWORK,DWORK5,DWORK5S,DWORK1,DWORK2,DWORK3,DWORKG,DWORK6,DWORK6S) __IF_ASYNC__
!$ACC LOOP GANG VECTOR
               DO I=1,GRIDC%RL%NP
                  DWORK(I,1)=DWORK(I,1)+DWORK5(I)/ MAX(SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1) &
                            +DWORK3(I,1)*DWORK3(I,1)),1.E-10_q)
                  DWORK(I,2)=DWORK(I,2)+DWORK5S(I)/ MAX(SQRT(DWORK1(I,2)*DWORK1(I,2)+DWORK2(I,2)*DWORK2(I,2) &
                            +DWORK3(I,2)*DWORK3(I,2)),1.E-10_q)
               ENDDO
NOACC !$OMP END PARALLEL DO
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORKG,DWORK6,DWORK6S)
!$ACC LOOP GANG VECTOR
               DO I=1,GRIDC%RL%NP
                  DWORKG(I,1)=DWORKG(I,1)+DWORK6(I)
                  DWORKG(I,2)=DWORKG(I,2)+DWORK6S(I)
               ENDDO
NOACC !$OMP END PARALLEL DO
!$ACC END PARALLEL
            ELSE
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK,DWORK5,DWORK1,DWORK2,DWORK3)
!$ACC PARALLEL PRESENT(DWORK,DWORK5,DWORK1,DWORK2,DWORK3,DWORKG,DWORK6) __IF_ASYNC__
!$ACC LOOP GANG VECTOR
               DO I=1,GRIDC%RL%NP
                  DWORK(I,1)=DWORK(I,1)+DWORK5(I)/ MAX(SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1) &
              &             +DWORK3(I,1)*DWORK3(I,1)),1.E-10_q)
                  DWORK(I,2)=DWORK(I,2)+DWORK5(I)/ MAX(SQRT(DWORK1(I,2)*DWORK1(I,2)+DWORK2(I,2)*DWORK2(I,2) &
              &             +DWORK3(I,2)*DWORK3(I,2)),1.E-10_q)
               ENDDO
NOACC !$OMP END PARALLEL DO
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORKG,DWORK6)
!$ACC LOOP GANG VECTOR
               DO I=1,GRIDC%RL%NP
                  DWORKG(I,1)=DWORKG(I,1)+DWORK6(I)
                  DWORKG(I,2)=DWORKG(I,2)+DWORK6(I)
               ENDDO
NOACC !$OMP END PARALLEL DO
!$ACC END PARALLEL
            ENDIF
         ENDIF
      ENDIF

!**********************************************************************
! TR: update the values for local SC CMBJ
!**********************************************************************

      IF ((XC%ID(1)==ID_XC_MBJ).OR.(XC%ID(1)==ID_XC_LMBJ)) &
     &   CALL UPDATE_CMBJ_AUX(GRIDC,LATT_CUR,XC,DHTOT,DENCOR,DWORK1,DWORK2,DWORK3,NCDIJ)

!=======================================================================
! gradient terms in stress tensor
!          d    f_xc     grad rho  x grad rho
! sum_r   ------------   --------------------- * LATT_CUR%OMEGA
!         d |grad rho|        |grad rho|
!=======================================================================
      SIF11=0
      SIF22=0
      SIF33=0
      SIF12=0
      SIF23=0
      SIF31=0
      IF ((XC%ID(1)/=ID_XC_MBJ) .AND. (XC%ID(1)/=ID_XC_LMBJ)) THEN
      IF (NCDIJ==1) THEN
         IF (XC%LDOMETAGGA) THEN
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK1, DWORK2, DWORK3, DWORK, DVC) &
NOACC !$OMP REDUCTION(+: SIF11, SIF22, SIF33, SIF12, SIF23, SIF31)
!$ACC PARALLEL LOOP PRESENT(DWORK1,DWORK2,DWORK3,DWORK,DVC) &
!$ACC REDUCTION(+:SIF11,SIF22,SIF33,SIF12,SIF23,SIF31) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
               SIF11=SIF11+DWORK1(I,1)*DWORK1(I,1)*DWORK(I,1)*0.5_q
               SIF22=SIF22+DWORK2(I,1)*DWORK2(I,1)*DWORK(I,1)*0.5_q
               SIF33=SIF33+DWORK3(I,1)*DWORK3(I,1)*DWORK(I,1)*0.5_q
               SIF12=SIF12+DWORK1(I,1)*DWORK2(I,1)*DWORK(I,1)*0.5_q
               SIF23=SIF23+DWORK2(I,1)*DWORK3(I,1)*DWORK(I,1)*0.5_q
               SIF31=SIF31+DWORK3(I,1)*DWORK1(I,1)*DWORK(I,1)*0.5_q

               SIF11=SIF11+DWORK1(I,1)*(DWORK1(I,1))*DVC(I)
               SIF22=SIF22+DWORK2(I,1)*(DWORK2(I,1))*DVC(I)
               SIF33=SIF33+DWORK3(I,1)*(DWORK3(I,1))*DVC(I)
               SIF12=SIF12+DWORK1(I,1)*(DWORK2(I,1))*DVC(I)
               SIF23=SIF23+DWORK2(I,1)*(DWORK3(I,1))*DVC(I)
               SIF31=SIF31+DWORK3(I,1)*(DWORK1(I,1))*DVC(I)
            ENDDO
NOACC !$OMP END PARALLEL DO
         ELSE
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK1, DWORK2, DWORK3, DWORK) &
NOACC !$OMP REDUCTION(+: SIF11, SIF22, SIF33, SIF12, SIF23, SIF31)
!$ACC PARALLEL LOOP PRESENT(DWORK1,DWORK2,DWORK3,DWORK) &
!$ACC& REDUCTION(+:SIF11,SIF22,SIF33,SIF12,SIF23,SIF31) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
               SIF11=SIF11+DWORK1(I,1)*DWORK1(I,1)*DWORK(I,1)
               SIF22=SIF22+DWORK2(I,1)*DWORK2(I,1)*DWORK(I,1)
               SIF33=SIF33+DWORK3(I,1)*DWORK3(I,1)*DWORK(I,1)
               SIF12=SIF12+DWORK1(I,1)*DWORK2(I,1)*DWORK(I,1)
               SIF23=SIF23+DWORK2(I,1)*DWORK3(I,1)*DWORK(I,1)
               SIF31=SIF31+DWORK3(I,1)*DWORK1(I,1)*DWORK(I,1)
            ENDDO
NOACC !$OMP END PARALLEL DO
         ENDIF
      ELSEIF ((NCDIJ==2).OR.(NCDIJ==4)) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(DWORK1,DWORK2,DWORK3,DWORK,DVC) &
!$ACC REDUCTION(+:SIF11,SIF22,SIF33,SIF12,SIF23,SIF31) __IF_ASYNC__
         DO ISP=1,MIN(NCDIJ,2)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK1, DWORK2, DWORK3, DWORK, DVC, ISP) &
NOACC !$OMP REDUCTION(+: SIF11, SIF22, SIF33, SIF12, SIF23, SIF31)
         DO I=1,GRIDC%RL%NP
            SIF11=SIF11+DWORK1(I,ISP)*DWORK1(I,ISP)*DWORK(I,ISP)
            SIF22=SIF22+DWORK2(I,ISP)*DWORK2(I,ISP)*DWORK(I,ISP)
            SIF33=SIF33+DWORK3(I,ISP)*DWORK3(I,ISP)*DWORK(I,ISP)
            SIF12=SIF12+DWORK1(I,ISP)*DWORK2(I,ISP)*DWORK(I,ISP)
            SIF23=SIF23+DWORK2(I,ISP)*DWORK3(I,ISP)*DWORK(I,ISP)
            SIF31=SIF31+DWORK3(I,ISP)*DWORK1(I,ISP)*DWORK(I,ISP)

            SIF11=SIF11+DWORK1(I,ISP)*(DWORK1(I,1)+DWORK1(I,2))*DVC(I)
            SIF22=SIF22+DWORK2(I,ISP)*(DWORK2(I,1)+DWORK2(I,2))*DVC(I)
            SIF33=SIF33+DWORK3(I,ISP)*(DWORK3(I,1)+DWORK3(I,2))*DVC(I)
            SIF12=SIF12+DWORK1(I,ISP)*(DWORK2(I,1)+DWORK2(I,2))*DVC(I)
            SIF23=SIF23+DWORK2(I,ISP)*(DWORK3(I,1)+DWORK3(I,2))*DVC(I)
            SIF31=SIF31+DWORK3(I,ISP)*(DWORK1(I,1)+DWORK1(I,2))*DVC(I)
         ENDDO
NOACC !$OMP END PARALLEL DO
         ENDDO
      ENDIF

! waith here for device offloading since the addition is done on the host
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      SIF11=SIF11*RINPL*LATT_CUR%OMEGA
      SIF22=SIF22*RINPL*LATT_CUR%OMEGA
      SIF33=SIF33*RINPL*LATT_CUR%OMEGA
      SIF12=SIF12*RINPL*LATT_CUR%OMEGA
      SIF23=SIF23*RINPL*LATT_CUR%OMEGA
      SIF31=SIF31*RINPL*LATT_CUR%OMEGA
      ENDIF   ! XC%ID(1)/=ID_XC_MBJ and XC%ID(1)/=ID_XC_LMBJ

!=======================================================================
! calculate
!              d    f_xc     grad rho
!        div  (------------  --------  )
!              d |grad rho| |grad rho|
!
! in reciprocal space
!=======================================================================

      IF (NCDIJ==1) THEN
         IF (XC%LDOMETAGGA) THEN
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I,ANAB1U,ANAB2U,ANAB3U,ANAB1D,ANAB2D,ANAB3D) SHARED(GRIDC,DWORK1,DWORK2,DWORK3,DWORK,DVC)
!$ACC PARALLEL LOOP PRESENT(DWORK1,DWORK2,DWORK3,DWORK,DVC) &
!$ACC& PRIVATE(ANAB1U,ANAB2U,ANAB3U,ANAB1D,ANAB2D,ANAB3D) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
               ANAB1U= DWORK1(I,1)*0.5_q
               ANAB2U= DWORK2(I,1)*0.5_q
               ANAB3U= DWORK3(I,1)*0.5_q
               ANAB1D= ANAB1U
               ANAB2D= ANAB2U
               ANAB3D= ANAB3U
               DWORK1(I,1) = ANAB1U* DWORK(I,1) + (ANAB1U+ANAB1D) * DVC(I)
               DWORK2(I,1) = ANAB2U* DWORK(I,1) + (ANAB2U+ANAB2D) * DVC(I)
               DWORK3(I,1) = ANAB3U* DWORK(I,1) + (ANAB3U+ANAB3D) * DVC(I)
            ENDDO
NOACC !$OMP END PARALLEL DO
         ELSE
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC,DWORK1,DWORK2,DWORK3,DWORK)
!$ACC PARALLEL LOOP PRESENT(DWORK1,DWORK2,DWORK3,DWORK) __IF_ASYNC__
            DO I=1,GRIDC%RL%NP
               DWORK1(I,1) = DWORK1(I,1)* REAL( DWORK(I,1) ,KIND=q)
               DWORK2(I,1) = DWORK2(I,1)* REAL( DWORK(I,1) ,KIND=q)
               DWORK3(I,1) = DWORK3(I,1)* REAL( DWORK(I,1) ,KIND=q)
            ENDDO
NOACC !$OMP END PARALLEL DO
         ENDIF
      ELSEIF ((NCDIJ==2).OR.(NCDIJ==4)) THEN
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I,ANAB1U,ANAB2U,ANAB3U,ANAB1D,ANAB2D,ANAB3D) SHARED(GRIDC,DWORK1,DWORK2,DWORK3,DWORK,DVC)
!$ACC PARALLEL LOOP PRESENT(DWORK1,DWORK2,DWORK3,DWORK,DVC) &
!$ACC& PRIVATE(ANAB1U,ANAB2U,ANAB3U,ANAB1D,ANAB2D,ANAB3D) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            ANAB1U= DWORK1(I,1)
            ANAB2U= DWORK2(I,1)
            ANAB3U= DWORK3(I,1)
            ANAB1D= DWORK1(I,2)
            ANAB2D= DWORK2(I,2)
            ANAB3D= DWORK3(I,2)

            DWORK1(I,1) = ANAB1U* DWORK(I,1) + (ANAB1U+ANAB1D) * DVC(I)
            DWORK2(I,1) = ANAB2U* DWORK(I,1) + (ANAB2U+ANAB2D) * DVC(I)
            DWORK3(I,1) = ANAB3U* DWORK(I,1) + (ANAB3U+ANAB3D) * DVC(I)

            DWORK1(I,2) = ANAB1D* DWORK(I,2) + (ANAB1U+ANAB1D) * DVC(I)
            DWORK2(I,2) = ANAB2D* DWORK(I,2) + (ANAB2U+ANAB2D) * DVC(I)
            DWORK3(I,2) = ANAB3D* DWORK(I,2) + (ANAB3U+ANAB3D) * DVC(I)
         ENDDO
NOACC !$OMP END PARALLEL DO
      ENDIF

      spin2: DO ISP=1,MIN(NCDIJ,2)
! x-component:
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK, DWORK1, ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DWORK1) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK1(I,ISP)
      ENDDO
NOACC !$OMP END PARALLEL DO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3D(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GX) SHARED(GRIDC, LATT_CUR, CWGRAD, CWORK, ISP)
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CWGRAD,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX) __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         CWGRAD(I,ISP)=CWORK(I,ISP)*GX*CITPI
      ENDDO
NOACC !$OMP END PARALLEL DO

! y-component:
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK, DWORK2, ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DWORK2) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK2(I,ISP)
      ENDDO
NOACC !$OMP END PARALLEL DO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3D(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GY) SHARED(GRIDC, LATT_CUR, CWGRAD, CWORK, ISP)
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CWGRAD,LATT_CUR) PRIVATE(N1,N2,N3,NC,GY) __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         CWGRAD(I,ISP)=CWGRAD(I,ISP)+CWORK(I,ISP)*GY*CITPI
      ENDDO
NOACC !$OMP END PARALLEL DO

! z-component:
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK, DWORK3, ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DWORK3) __IF_ASYNC__
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK3(I,ISP)
      ENDDO
NOACC !$OMP END PARALLEL DO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3D(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GZ) SHARED(GRIDC, LATT_CUR, CWGRAD, CWORK, ISP)
!$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CWGRAD,LATT_CUR) PRIVATE(N1,N2,N3,NC,GZ) __IF_ASYNC__
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
         CWGRAD(I,ISP)=CWGRAD(I,ISP)+CWORK(I,ISP)*GZ*CITPI
      ENDDO
NOACC !$OMP END PARALLEL DO

      CALL SETUNB(CWGRAD(1,ISP),GRIDC)
      CALL FFT3D(CWGRAD(1,ISP),GRIDC,1)
      CALL OPSYNC(CWGRAD(1,ISP),DWGRAD(1,ISP),GRIDC%NPLWV)

      ENDDO spin2


!=======================================================================
! laplacian terms in stress tensor
!                     d  f_xc         d^2 n_sigma
! sum_(r,sigma)    ------------      ------------- * LATT_CUR%OMEGA
!                  d lap rho_sigma     d r_i d r_j
!=======================================================================

! if ispin == 1, then denshess(:,:,1) is the Hessian of the total density
! dworkl is still d f_xc / d lap rho_sigma, so the sum will still yield
! twice the sum over a single spin channel. No additional factor needed

! if ispin == 2, then denshess(:,:,sigma) is the Hessian for each spin channel sigma,
! no correcting factor needed

! this has to appear here rather than with the gradient contribution to
! the stress tensor because DWORK is required for the gradient contributions
! to the XC potential

      IF (XC%LDLAP) THEN
        SIFLAP11=0
        SIFLAP22=0
        SIFLAP33=0
        SIFLAP12=0
        SIFLAP23=0
        SIFLAP31=0
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(DWORKL,DENSHESS) &
!$ACC REDUCTION(+:SIFLAP11,SIFLAP22,SIFLAP33,SIFLAP12,SIFLAP23,SIFLAP31) __IF_ASYNC__
        SLST: DO ISP=1,MIN(NCDIJ,2)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORKL, DENSHESS, ISP) &
NOACC !$OMP REDUCTION(+: SIFLAP11, SIFLAP22, SIFLAP33, SIFLAP12, SIFLAP23, SIFLAP31)
          DO I=1,GRIDC%RL%NP
            SIFLAP11=SIFLAP11 + 2._q*DWORKL(I,ISP)*DENSHESS(I,1,ISP)
            SIFLAP22=SIFLAP22 + 2._q*DWORKL(I,ISP)*DENSHESS(I,4,ISP)
            SIFLAP33=SIFLAP33 + 2._q*DWORKL(I,ISP)*DENSHESS(I,6,ISP)
            SIFLAP12=SIFLAP12 + 2._q*DWORKL(I,ISP)*DENSHESS(I,2,ISP)
            SIFLAP23=SIFLAP23 + 2._q*DWORKL(I,ISP)*DENSHESS(I,5,ISP)
            SIFLAP31=SIFLAP31 + 2._q*DWORKL(I,ISP)*DENSHESS(I,3,ISP)
          ENDDO
NOACC !$OMP END PARALLEL DO
        ENDDO SLST

! waith here for device offloading since the addition is done on the host
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
        SIF11=SIF11 + SIFLAP11*RINPL*LATT_CUR%OMEGA
        SIF22=SIF22 + SIFLAP22*RINPL*LATT_CUR%OMEGA
        SIF33=SIF33 + SIFLAP33*RINPL*LATT_CUR%OMEGA
        SIF12=SIF12 + SIFLAP12*RINPL*LATT_CUR%OMEGA
        SIF23=SIF23 + SIFLAP23*RINPL*LATT_CUR%OMEGA
        SIF31=SIF31 + SIFLAP31*RINPL*LATT_CUR%OMEGA

!=======================================================================
! calculate
!               d    f_xc
!        lap ( ------------ )
!               d lap rho
!
! in reciprocal space
!=======================================================================

        SPINLAP: DO ISP=1,MIN(NCDIJ,2)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK, DWORKL, ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DWORKL) __IF_ASYNC__
          DO I=1,GRIDC%RL%NP
            DWORK(I,ISP) = DWORKL(I,ISP)
          ENDDO
NOACC !$OMP END PARALLEL DO

          CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
          CALL FFT3D(CWORK(1,ISP),GRIDC,-1)
          CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, N1, NC, N2, N3, GX, GY, GZ, GG) SHARED(GRIDC,LATT_CUR,CWORK,ISP)
!$ACC PARALLEL LOOP PRESENT(GRIDC,LATT_CUR,CWORK) PRIVATE(N1,N2,N3,NC,GX,GY,GZ,GG) __IF_ASYNC__
          DO I=1,GRIDC%RC%NP
             N1= MOD((I-1),GRIDC%RC%NROW) +1
             NC= (I-1)/GRIDC%RC%NROW+1
             N2= GRIDC%RC%I2(NC)
             N3= GRIDC%RC%I3(NC)

             GX = GRIDC%LPCTX(N1)*LATT_CUR%B(1,1) &
            &    + GRIDC%LPCTY(N2)*LATT_CUR%B(1,2) &
            &    + GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
             GY = GRIDC%LPCTX(N1)*LATT_CUR%B(2,1) &
            &    + GRIDC%LPCTY(N2)*LATT_CUR%B(2,2) &
            &    + GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
             GZ = GRIDC%LPCTX(N1)*LATT_CUR%B(3,1) &
            &    + GRIDC%LPCTY(N2)*LATT_CUR%B(3,2) &
            &    + GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)
             GG = GX**2 + GY**2 + GZ**2

             CWORK(I,ISP) = GG*CITPI*CITPI*CWORK(I,ISP)
          ENDDO
NOACC !$OMP END PARALLEL DO

          CALL SETUNB(CWORK(1,ISP),GRIDC)
          CALL FFT3D(CWORK(1,ISP),GRIDC,1)
          CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I) SHARED(GRIDC, DWORK, DWORKL, ISP)
!$ACC PARALLEL LOOP PRESENT(DWORK,DWORKL) __IF_ASYNC__
          DO I=1,GRIDC%RL%NP
            DWORKL(I,ISP) = REAL(DWORK(I,ISP),KIND=q)
          ENDDO
NOACC !$OMP END PARALLEL DO

        ENDDO SPINLAP
      ENDIF

!=======================================================================
! Now prepare the rest:
! (store rho in DWORK3 and quantity of above in DWORK1)
!=======================================================================
      XCENC=0._q
      CVZERO=0._q
      XCENCC=0._q

      IF (NCDIJ==1) THEN
      ! nonspin-polarized case
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, RHO, VXC) SHARED(GRIDC, DWORK4, DWORKG, DWGRAD, RINPL, DWORK, DWORKL, LATT_CUR, DHTOT, XC) &
NOACC !$OMP REDUCTION(-:XCENCC, XCENC, CVZERO)
!$ACC PARALLEL LOOP PRESENT(DWORK4,DWORKG,DWGRAD,DWORK,LATT_CUR,XC,DHTOT,DWORKL) &
!$ACC& REDUCTION(+:CVZERO,XCENCC,XCENC) PRIVATE(RHO,VXC) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            RHO=REAL( DWORK4(I) ,KIND=q)
            VXC=DWORKG(I,1)- REAL( DWGRAD(I,1) ,KIND=q) *RINPL
            IF (XC%LDLAP) THEN
               VXC=VXC + DWORKL(I,1)*RINPL
            ENDIF
            DWORK(I,1)=VXC
            CVZERO=CVZERO+VXC
            XCENCC=XCENCC-VXC*RHO*LATT_CUR%OMEGA
            XCENC=XCENC  -VXC* REAL( DHTOT(I,1) ,KIND=q)
         ENDDO
NOACC !$OMP END PARALLEL DO
      ELSEIF (NCDIJ==2) THEN
      ! Spin-polarized case
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, RHO1, RHO2, VXC1, VXC2, VXC) SHARED(GRIDC, DENCOR, NCDIJ, DWORKG, DWGRAD, RINPL, DWORK, DWORKL, LATT_CUR, DHTOT, XC, CHGMIN) &
NOACC !$OMP REDUCTION(-:XCENCC, XCENC, CVZERO)
!$ACC PARALLEL LOOP PRESENT(DHTOT,DENCOR,LATT_CUR,XC,DWORKG,DWGRAD,DWORK,DWORKL) &
!$ACC& REDUCTION(+:CVZERO,XCENCC,XCENC) PRIVATE(RHO1,RHO2,VXC1,VXC2,VXC) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
            RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)/NCDIJ)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
            RHO2= MAX(REAL((DHTOT(I,2)+DENCOR(I)/NCDIJ)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)

            VXC1=DWORKG(I,1)- REAL( DWGRAD(I,1) ,KIND=q) *RINPL
            VXC2=DWORKG(I,2)- REAL( DWGRAD(I,2) ,KIND=q) *RINPL
            IF (XC%LDLAP) THEN
               VXC1=VXC1 + DWORKL(I,1)*RINPL
               VXC2=VXC2 + DWORKL(I,2)*RINPL
            ENDIF
            DWORK(I,1)=VXC1
            DWORK(I,2)=VXC2
            VXC = 0.5_q*(VXC1+VXC2)
            CVZERO=CVZERO+VXC
            XCENCC=XCENCC-VXC1*RHO1*LATT_CUR%OMEGA-VXC2*RHO2*LATT_CUR%OMEGA
            XCENC=XCENC  -VXC1* REAL( DHTOT(I,1) ,KIND=q) -VXC2* REAL( DHTOT(I,2) ,KIND=q)
         ENDDO
NOACC !$OMP END PARALLEL DO
      ELSEIF (NCDIJ==4) THEN
      ! non-collinear case
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(I, MAG_NORM, RHO1, RHO2, VXC1, VXC2, VXC) SHARED(GRIDC, DENCOR, DWORKG, DWGRAD, RINPL, DWORK, DWORKL, LATT_CUR, DHTOT, XC, CHGMIN) &
NOACC !$OMP REDUCTION(-:XCENCC, XCENC, CVZERO)
!$ACC PARALLEL LOOP PRESENT(DHTOT,DENCOR,LATT_CUR,XC,DWORKG,DWGRAD,DWORK,DWORKL) &
!$ACC& REDUCTION(+:CVZERO,XCENCC,XCENC) PRIVATE(MAG_NORM,RHO1,RHO2,VXC1,VXC2,VXC) __IF_ASYNC__
         DO I=1,GRIDC%RL%NP
!
! | m |
!
            MAG_NORM=SQRT(ABS(DHTOT(I,2)*DHTOT(I,2)+ DHTOT(I,3)*DHTOT(I,3) + DHTOT(I,4)*DHTOT(I,4)))
!
! RHO1: \rho_up   = ( \rho_tot + | m | )/2
! RHO2: \rho_down = ( \rho_tot - | m | )/2
!
            RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)+MAG_NORM)*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
            RHO2= MAX(REAL((DHTOT(I,1)+DENCOR(I)-MAG_NORM)*0.5_q/LATT_CUR%OMEGA ,KIND=q), CHGMIN)

            VXC1=DWORKG(I,1)- REAL( DWGRAD(I,1) ,KIND=q) *RINPL
            VXC2=DWORKG(I,2)- REAL( DWGRAD(I,2) ,KIND=q) *RINPL
            IF (XC%LDLAP) THEN
               VXC1=VXC1 + DWORKL(I,1)*RINPL
               VXC2=VXC2 + DWORKL(I,2)*RINPL
            ENDIF
            DWORK(I,1)=VXC1
            DWORK(I,2)=VXC2
            VXC = 0.5_q*(VXC1+VXC2)
            CVZERO=CVZERO+VXC
            XCENCC=XCENCC-VXC1*RHO1*LATT_CUR%OMEGA-VXC2*RHO2*LATT_CUR%OMEGA
            XCENC=XCENC  -VXC1* REAL( (DHTOT(I,1)+MAG_NORM)*0.5_q ,KIND=q) &
           &             -VXC2* REAL( (DHTOT(I,1)-MAG_NORM)*0.5_q ,KIND=q)
         ENDDO
NOACC !$OMP END PARALLEL DO
      ENDIF
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(DWORKL,DENSHESS) __IF_ASYNC__

      CVZERO=CVZERO*RINPL
      XCENC =(XCENC+EXC-EDCMU)*RINPL
      XCENCC=(XCENCC+EXC)*RINPL
      EXC   =EXC*RINPL
      EDCMU =EDCMU*RINPL

      SIF11=SIF11-XCENCC
      SIF22=SIF22-XCENCC
      SIF33=SIF33-XCENCC
      XCSIF(1,1)=SIF11
      XCSIF(2,2)=SIF22
      XCSIF(3,3)=SIF33
      XCSIF(1,2)=SIF12
      XCSIF(2,1)=SIF12
      XCSIF(2,3)=SIF23
      XCSIF(3,2)=SIF23
      XCSIF(3,1)=SIF31
      XCSIF(1,3)=SIF31

      IF (XC%LUSE_VDW) XCSIF = XCSIF - STRESS

      IF (.NOT.XC%LMU) THEN
         CALLMPI( M_sum_2(GRIDC%COMM, EXC, XCENC))
      ELSE
         CALLMPI( M_sum_3(GRIDC%COMM, EXC, XCENC, EDCMU))
      ENDIF
      CALLMPI( M_sum_d(GRIDC%COMM, XCSIF, 9))
      CALLMPI( M_sum_z(GRIDC%COMM,CVZERO,1))

      ! Flip CHTOT back to (rho,mag)
      IF ((NCDIJ==1).OR.((NCDIJ==2).AND.XC%LDOMETAGGA)) CALL RL_FLIP(CHTOT,GRIDC,NCDIJ,.FALSE.)

!$ACC EXIT DATA DELETE(CWGRAD,DWORK,DWORKG,DWORK1,DWORK2,DWORK3,DVC,DWORK4,DWORK5,DWORK6,DWORK7) __IF_ASYNC__
      DEALLOCATE(CWGRAD,DWORKG,DWORK1,DWORK2,DWORK3,DVC,DWORK4,DWORK5,DWORK6,DWORK7)

      IF (XC%LSPIN_VDW) THEN
!$ACC EXIT DATA DELETE(DWORK4S,DWORK5S,DWORK6S,DWORK7S) __IF_ASYNC__
         DEALLOCATE(DWORK4S,DWORK5S,DWORK6S,DWORK7S)
      ENDIF

      IF (XC%LTAU) THEN
!$ACC EXIT DATA DELETE(KINEDEN,CKINEDEN,TAUC) __IF_ASYNC__
         DEALLOCATE(CKINEDEN,TAUC)
      ENDIF

      RETURN
      END SUBROUTINE FEXCG
