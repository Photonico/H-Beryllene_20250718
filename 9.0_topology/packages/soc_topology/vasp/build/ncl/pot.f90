# 1 "pot.F"
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


# 2 "pot.F" 2 
      MODULE pot
      USE prec
      USE base, ONLY: write_pot
      USE charge
      USE lattice, ONLY: latt
      USE pot_struct_def
      USE pot_electrostat
      USE coulomb_cutoff

      COMPLEX(q), ALLOCATABLE, PRIVATE :: EXTERNAL_POTENTIAL(:)
      PRIVATE ADD_EXTERNAL_POTENTIAL

      INTERFACE str
         MODULE PROCEDURE WRT_POTENTIAL_TO_STRING
      END INTERFACE str

      CONTAINS
!************************ SUBROUTINE POTLOK ****************************
!
!> This subroutine calculates the total local potential CVTOT
!> which is the sum of the hartree potential, the exchange-correlation
!> potential and the ionic local potential. The routine also calculates
!> the total local potential SV on the small grid.
!>
!> On entry: \n
!>  CHTOT(:,1)    density \n
!>  CHTOT(:,2)    respectively CHTOT(:,2:4) contain the magnetization \n
!> On return (LNONCOLLINEAR=.FALSE.): \n
!>  CVTOT(:,1)    potential for up \n
!>  CVTOT(:,2)    potential for down \n
!>
!> On return (LNONCOLLINEAR=.TRUE.): \n
!>  CVTOT(:,1:4)  spinor representation of potential
!>
!> @details @ref openacc :
!> On entering this subroutine OpenACC execution is switched on.
!> On exit the OpenACC execution mode reverts to its original status,
!> and SV, CHTOT, CVTOT are copied back to the host.
!>
!
! RCS:  $Id: pot.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!***********************************************************************
    SUBROUTINE POTLOK(KINEDEN,GRID,GRIDC,GRID_SOFT,COMM_INTER,WDES,  &
                  INFO,XC,P,T_INFO,E,LATT_CUR,  &
                  CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,MUTOT,MU,SOFT_TO_C,XCSIF, &
                  COULOMB_POT)

!! USE moffload

      USE prec
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, only: fexcf, fexcp
      USE tau_mu, ONLY : tau_handle
      USE base
      USE wave
      USE mdipol
      USE Constrained_M_modular
      USE main_mpi, ONLY: COMM
      USE scpc, ONLY: SCPC_APPLY
      USE plugins, ONLY: PLUGINS_LOCAL_POTENTIAL
! solvation__
      USE solvation
! solvation__
! bexternal__
      USE bexternal
! bexternal__
! embedding__
      USE mextpot
! embedding__
      IMPLICIT NONE

      TYPE (tau_handle)  KINEDEN
      TYPE (xc_info)     XC
      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT
      TYPE (wavedes)     WDES
      TYPE (transit)     SOFT_TO_C
      TYPE (info_struct) INFO
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (energy)      E
      TYPE (latt)        LATT_CUR
      TYPE (communic)    COMM_INTER
      TYPE (coulomb_potential), OPTIONAL :: COULOMB_POT

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP), CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      COMPLEX(q), INTENT(OUT) :: SV(GRID%MPLWV, WDES%NCDIJ)
      COMPLEX(q), INTENT(OUT) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      COMPLEX(q), INTENT(IN), TARGET :: CHTOT(GRIDC%MPLWV, WDES%NCDIJ)
      COMPLEX(q)      DENCOR(GRIDC%RL%NP)
      REAL(q)    XCSIF(3,3),TMPSIF(3,3),XDMTMP(1)
! work arrays (allocated after call to FEXCG)
      COMPLEX(q), ALLOCATABLE::  CWORK1(:),CWORK(:,:)
      REAL(q) ELECTROSTATIC,EXC,EXCG,RINPL,XCENCG,FACTM
      LOGICAL, EXTERNAL :: L_NO_LSDA_GLOBAL
      INTEGER I,ISP,MWORK1,N1,N2,N3,NC,NG, ASYNC_PRE
      COMPLEX(q) CVZERG
      COMPLEX(q) , POINTER :: MU(:,:)
      COMPLEX(q), POINTER :: MUTOT(:,:)

      COMPLEX(q) , POINTER :: MU_(:,:)
      COMPLEX(q), POINTER :: MUTOT_(:,:)
      TYPE(COULOMB_POTENTIAL) :: COULOMB_POT_PERIODIC

# 118


# 127


      COMPLEX(q), POINTER :: CTMP(:,:)

      

# 150


      IF (XC%LNO_AUG_XC) THEN
         ALLOCATE(CTMP(GRIDC%MPLWV,WDES%NCDIJ))
!$ACC ENTER DATA COPYIN(CHDEN) 
      ELSE
         CTMP => CHTOT
      ENDIF
!$ACC ENTER DATA CREATE(CTMP) 

      IF (XC%LMU) THEN
         IF (.NOT.ASSOCIATED(MU) .OR. .NOT.ASSOCIATED(MUTOT)) &
            CALL vtutor%bug("POTLOK: MU and/or MUTOT not associated.", "pot.F", 162)
         MU_ => MU ; MUTOT_ => MUTOT
!$ACC ENTER DATA CREATE(MU_,MUTOT_) 
!$ACC KERNELS PRESENT(MU_,MUTOT_) 
         MU_=0 ; MUTOT_=0
!$ACC END KERNELS
      ELSE
         ALLOCATE(MU_(1,WDES%NCDIJ),MUTOT_(1,WDES%NCDIJ))
!$ACC ENTER DATA CREATE(MU_,MUTOT_) 
      ENDIF

      MWORK1=MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV)
      ALLOCATE(CWORK1(MWORK1),CWORK(GRIDC%MPLWV,WDES%NCDIJ))
!$ACC ENTER DATA CREATE(CWORK,CWORK1) 

!-----------------------------------------------------------------------
!
!  calculate the exchange correlation potential and the dc. correction
!
!-----------------------------------------------------------------------
      EXC     =0
      E%XCENC =0
      E%EXCG  =0
      E%CVZERO=0
      XCSIF   =0

!$ACC KERNELS PRESENT(CVTOT) 
      CVTOT   =0
!$ACC END KERNELS

# 202


      IF (.NOT.XC%LNOXC) THEN
! transform the charge density to real space
         EXCG  =0
         XCENCG=0
         CVZERG=0
         TMPSIF=0

         IF (XC%LNO_AUG_XC) THEN
! We will use the soft charge density instead of the augmented (1._q,0._q)
!$ACC KERNELS PRESENT(CTMP) 
            CTMP=0
!$ACC END KERNELS
            DO ISP=1,WDES%NCDIJ
               CALL ADD_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CHDEN(1,ISP),CTMP(1,ISP))
               CALL SETUNB_COMPAT(CTMP(1,ISP),GRIDC)
               CALL FFT3D(CTMP(1,ISP),GRIDC,1)
            ENDDO
         ENDIF

! Bring the augmented charge density to real space
         DO ISP=1,WDES%NCDIJ
            CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
         ENDDO

# 238


         IF (WDES%ISPIN==2) THEN

! get the charge and the total magnetization
            CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)
! do LDA+U instead of LSDA+U
            IF (L_NO_LSDA_GLOBAL()) THEN
!$ACC KERNELS PRESENT(CWORK) 
               CWORK(:,2)=0
!$ACC END KERNELS
            ENDIF

            IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
! gradient corrections to LDA
! unfortunately FEXCGS requires (up,down) density
! instead of (rho,mag)
               IF (.NOT.XC%LDOMETAGGA) CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
! GGA potential
               CALL FEXCG(XC,2,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                    KINEDEN,CWORK,DENCOR,CVTOT,MUTOT_,XDMTMP,0,1,SIZE(MUTOT_,1))
               IF (.NOT.XC%LDOMETAGGA) CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
            ENDIF

! add LDA part of potential
            CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
               CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
               E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
!gk COH
! add Coulomb hole
            CALL COHSM1_RGRID(XC, 2, CWORK(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
!gK COHend
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

! get the proper direction vx = v0 + hat m delta v
            CALL MAG_DIRECTION(CTMP(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)

         ELSEIF (WDES%LNONCOLLINEAR) THEN

            IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
               CALL FEXCG(XC,4,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                    KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDMTMP,0,1,SIZE(MUTOT_,1))
            ENDIF

            IF (XC%LDOMETAGGA) THEN
               IF (XC%LMU) THEN
! MUTOT comes out as two component (up,down) quantities, where
! "up" and "down" are taken w.r.t. the local magnetization direction,
! this is now cast into (rho, mag) form.
                  CALL MAG_DIRECTION(CTMP(1,1), MUTOT_(1,1), GRIDC, WDES%NCDIJ)
! and rearranged from (rho,mag) to spinor representation
                  CALL POT_FLIP_RL(MUTOT_, GRIDC, WDES%NCDIJ)
               ENDIF
            ENDIF

! FEXCF requires (rho,mag) density instead of (up,down)
            CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)
! quick hack to do LDA+U instead of LSDA+U
            IF (L_NO_LSDA_GLOBAL()) THEN
!$ACC KERNELS PRESENT(CWORK) 
               CWORK(:,2)=0
!$ACC END KERNELS
            ENDIF
! end of hack
! add LDA part of potential
            CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
               CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
               E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
!gk COH
! add Coulomb hole
            CALL COHSM1_RGRID(XC, 2, CWORK(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
!gK COHend
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)
! get the proper direction vx = v0 + hat m delta v

            CALL MAG_DIRECTION(CTMP(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)

         ELSE

            IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
! gradient corrections to LDA
               CALL FEXCG(XC,1,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                    KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDMTMP,0,1,SIZE(MUTOT_,1))
            ENDIF

! LDA part of potential
            CALL FEXCP(GRIDC,XC,LATT_CUR%OMEGA, &
                 CTMP,DENCOR,CVTOT,CWORK,E%CVZERO,EXC,E%XCENC,XCSIF,.TRUE.)
!gk COH
! add Coulomb hole
            CALL COHSM1_RGRID(XC, 1, CTMP(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
!gK COHend

         ENDIF

# 379


         XCSIF=XCSIF+TMPSIF
         E%EXCG=EXC+EXCG
         E%XCENC=E%XCENC+XCENCG
         E%CVZERO=E%CVZERO+CVZERG

      ELSE
         DO ISP=1,WDES%NCDIJ
            CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
         ENDDO
      ENDIF

      IF (XC%LNO_AUG_XC) THEN
! Store CVTOT in CTMP
!$ACC KERNELS PRESENT(CTMP,CVTOT) 
         CTMP=CVTOT
!$ACC END KERNELS
      ENDIF

      IF (XC%LSFBXC) CALL REMOVE_SOURCES_FROM_BXC(E, CVTOT, CHTOT, GRIDC, WDES, LATT_CUR)

!-MM- changes to accomodate constrained moments
!-----------------------------------------------------------------------
! add constraining potential
!-----------------------------------------------------------------------


      IF (M_CONSTRAINED()) THEN
! NB. at this point both CHTOT and CVTOT must be given
! in (charge,magnetization) convention in real space
!$ACC UPDATE SELF(CHTOT,CVTOT) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
         CALL M_INT(CHTOT,GRIDC,WDES)
         CALL ADD_CONSTRAINING_POT(CVTOT,GRIDC,WDES)
!$ACC UPDATE DEVICE(CVTOT) 
      ENDIF


!-MM- end of addition

!-----------------------------------------------------------------------
! calculate the total potential
!-----------------------------------------------------------------------
! add external electrostatic potential
      DIP%ECORR=0
      DIP%E_ION_EXTERN=0

      IF (DIP%LCOR_DIP) THEN
! get the total charge and store it in CWORK
         IF  ( WDES%NCDIJ > 1) THEN
            CALL MAG_DENSITY(CHTOT,CWORK, GRIDC, WDES%NCDIJ)
         ELSE
            CALL RL_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
         ENDIF

!$ACC UPDATE SELF(CWORK,CVTOT) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
           CALL CDIPOL(GRIDC, LATT_CUR,P,T_INFO, &
             CWORK,CSTRF,CVTOT(1,1), WDES%NCDIJ, INFO%NELECT )
!$ACC UPDATE DEVICE(CVTOT) 

         CALL EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT(1,1))
      ELSE
         CALL EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT(1,1))
      ENDIF
! embedding__
      IF (EXTPT_LEXTPOT()) CALL EXTPT_EXTERNAL_POT_ADD(GRIDC, LATT_CUR, CVTOT(1,1))
! embedding__

      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO
!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3D(CVTOT(1,ISP),GRIDC,-1)
      ENDDO

      IF (XC%LMU) THEN
!-----------------------------------------------------------------------
! The same procedure to calculate MU from MUTOT
!-----------------------------------------------------------------------
         RINPL=1._q/GRIDC%NPLWV
         DO  ISP=1,WDES%NCDIJ
            CALL RL_ADD(MUTOT_(1,ISP),RINPL,MUTOT_(1,ISP),0._q,MUTOT_(1,ISP),GRIDC)
            CALL FFT3D(MUTOT_(1,ISP),GRIDC,-1)
         ENDDO
!-----------------------------------------------------------------------
! copy MUTOT_ to MU_ and set contribution of unbalanced lattice-vectors
! to (0._q,0._q), then FFT of MU_ and MUTOT_ to real space
!-----------------------------------------------------------------------
         DO ISP=1,WDES%NCDIJ
            CALL SETUNB_COMPAT(MUTOT_(1,ISP),GRIDC)
            CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,MUTOT_(1,ISP),CWORK1)
            CALL SETUNB(CWORK1,GRID_SOFT)
            CALL FFT3D(CWORK1,GRID_SOFT, 1)
            CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,MU_(1,ISP),GRID_SOFT)

! final result is only correct for first in-band-group
! (i.e. proc with nodeid 1 in COMM_INTER)
! copy to other in-band-groups using COMM_INTER
! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
# 486

            CALL M_bcast_z(COMM_INTER, MU_(1,ISP), GRID%RL%NP)

            CALL FFT3D(MUTOT_(1,ISP),GRIDC,1)
         ENDDO
      ENDIF
!-----------------------------------------------------------------------
! add the hartree potential and the double counting corrections
!-----------------------------------------------------------------------
       IF (PRESENT(COULOMB_POT)) THEN
          CALL POTHAR(GRIDC, LATT_CUR, COULOMB_POT, CHTOT, CWORK,E%DENC)
          CALL GET_RHO0(GRIDC,CWORK,E%POTHARZERO)
       ELSE
          CALL POTHAR(GRIDC, LATT_CUR, COULOMB_POT_PERIODIC, CHTOT, CWORK,E%DENC)
       END IF
# 507

!$ACC PARALLEL LOOP PRESENT(CVTOT,CWORK) 
      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
!-----------------------------------------------------------------------
! add external potential in reciprocal space
!-----------------------------------------------------------------------
      CALL ADD_EXTERNAL_POTENTIAL(GRIDC, CVTOT)
! solvation__
!-----------------------------------------------------------------------
! add the dielectric corrections to CVTOT and the energy
!-----------------------------------------------------------------------
      CALL SOL_Vcorrection(INFO,T_INFO,LATT_CUR,P,WDES,GRIDC,CHTOT,CVTOT)
! solvation__
!-----------------------------------------------------------------------
!  add local pseudopotential potential
!-----------------------------------------------------------------------
      IF(INFO%TURBO==0)THEN
         IF (PRESENT(COULOMB_POT)) THEN
            CALL POTION(GRIDC,P,LATT_CUR,T_INFO,COULOMB_POT,CWORK,CSTRF,E%PSCENC)
            CALL GET_RHO0(GRIDC,CWORK,E%POTIONZERO)
         ELSE
            CALL POTION(GRIDC,P,LATT_CUR,T_INFO,COULOMB_POT_PERIODIC,CWORK,CSTRF,E%PSCENC)
         END IF
      ELSE
!!    
         CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,COULOMB_POT,CWORK(:,1),E%PSCENC,E%TEWEN)
!!    
!$ACC UPDATE DEVICE(CWORK) 
      ENDIF

# 545


      ELECTROSTATIC=0
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRIDC,CWORK,CHTOT) &
!$ACC& REDUCTION(+:ELECTROSTATIC) PRIVATE(N2,N3,NG,FACTM) 
      col: DO NC=1,GRIDC%RC%NCOL
 N2= GRIDC%RC%I2(NC)
 N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
!!   N2= GRIDC%RC%I2(NC)
!!   N3= GRIDC%RC%I3(NC)
        NG = N1 + (NC-1)*GRIDC%RC%NROW
        FACTM=1
        

        ELECTROSTATIC=ELECTROSTATIC+  CWORK(NG,1)*CONJG(CHTOT(NG,1))
      ENDDO row
      ENDDO col
      ELECTROSTATIC=ELECTROSTATIC+E%PSCENC-E%DENC+E%TEWEN

      E%PSCENC=E%PSCENC + DIP%ECORR + DIP%E_ION_EXTERN

! apply self consistent potential correction
      CALL SCPC_APPLY(GRIDC, LATT_CUR, COULOMB_POT_PERIODIC, CHTOT, P, T_INFO, CSTRF, CVTOT, E%ESCPC)

!$ACC PARALLEL LOOP PRESENT(CVTOT,CWORK) 
      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
! bexternal__
      IF (LBEXTERNAL()) CALL BEXT_ADDV(CVTOT,GRIDC,SIZE(CVTOT,2))
! bexternal__

# 610


!=======================================================================
! if overlap is used :
! copy CVTOT to SV and set contribution of unbalanced lattice-vectors
! to (0._q,0._q),  then  FFT of SV and CVTOT to real space
!=======================================================================

      CALL POT_FLIP(CVTOT, GRIDC,WDES%NCDIJ )

      DO ISP=1,WDES%NCDIJ
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
         CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CVTOT(1,ISP),CWORK1)
         CALL SETUNB(CWORK1,GRID_SOFT)
         CALL FFT3D(CWORK1,GRID_SOFT, 1)
         CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,SV(1,ISP),GRID_SOFT)

!  final result is only correct for first in-band-group
! (i.e. proc with nodeid 1 in COMM_INTER)
!  copy to other in-band-groups using COMM_INTER
! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
# 633

         CALL M_bcast_z(COMM_INTER, SV(1,ISP), GRID%RL%NP)

         CALL FFT3D(CVTOT(1,ISP),GRIDC,1)
      ENDDO

      IF (XC%LNO_AUG_XC) THEN
         CALL POT_FLIP(CTMP, GRIDC,WDES%NCDIJ )
! Subtract CTMP from CVTOT
         DO ISP=1,WDES%NCDIJ
            CALL RL_ADD(CVTOT(1,ISP),1.0_q,CTMP(1,ISP),-1.0_q,CVTOT(1,ISP),GRIDC)
         ENDDO
      ENDIF

!$ACC EXIT DATA DELETE(CTMP,CWORK1,CWORK) 
      IF (XC%LNO_AUG_XC) THEN
!$ACC EXIT DATA DELETE(CHDEN) 
         DEALLOCATE(CTMP)
      ENDIF

! Ensure that MU and MUTOT are up-to-date on the host
!$ACC UPDATE SELF(MU_,MUTOT_) IF(OFFLOAD_ON.AND.XC%LMU) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(MU_,MUTOT_) IF(OFFLOAD_ON.AND.XC%LMU) ASYNC(ACC_ASYNC_Q)
      IF (.NOT.XC%LMU) THEN
!$ACC EXIT DATA DELETE(MU_,MUTOT_) 
         DEALLOCATE(MU_,MUTOT_)
      ENDIF
      NULLIFY(MU_,MUTOT_)

      DEALLOCATE(CWORK1,CWORK)

# 680


      

      RETURN
    END SUBROUTINE POTLOK


!***********************************************************************
!
!> Small helper routine to set the local potential on the coarse
!> plane wave grid from the full dense (augmentation) grid.
!>
!> CVTOT must be supplied in reciprocal space (not usually the case)
!> and SV is returned in real space.
!
!***********************************************************************

    SUBROUTINE SET_SV( GRID, GRIDC, GRID_SOFT, COMM_INTER, SOFT_TO_C, NCDIJ, SV, CVTOT)

      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT NONE

      INTEGER NCDIJ
      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT
      TYPE (transit)     SOFT_TO_C

      COMPLEX(q)   SV(GRID%MPLWV, NCDIJ)
      COMPLEX(q) CVTOT(GRIDC%MPLWV, NCDIJ)
      TYPE (communic)    COMM_INTER
! work arrays
      COMPLEX(q) ::  CWORK1(GRID_SOFT%MPLWV)
      INTEGER ISP


      DO ISP=1,NCDIJ
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
         CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CVTOT(1,ISP),CWORK1)
         CALL SETUNB(CWORK1,GRID_SOFT)
! transform to real-space representation
         CALL FFT3D(CWORK1,GRID_SOFT, 1)
         CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,SV(1,ISP),GRID_SOFT)

!  final result is only correct for first in-band-group
! (i.e. proc with nodeid 1 in COMM_INTER)
!  copy to other in-band-groups using COMM_INTER
! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
# 731

         CALL M_bcast_z(COMM_INTER, SV(1,ISP), GRID%RL%NP)

         CALL FFT3D(CVTOT(1,ISP),GRIDC,1)
      ENDDO
    END SUBROUTINE SET_SV

    SUBROUTINE ADD_EXTERNAL_POTENTIAL(GRIDC, CVTOT)
! adds the external potential stored as module variable
      USE mgrid, ONLY: grid_3d
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(grid_3d), INTENT(IN) :: GRIDC
      COMPLEX(q), INTENT(INOUT) :: CVTOT(:,:)
!
      IF (.NOT.ALLOCATED(EXTERNAL_POTENTIAL)) RETURN
!
! sanity checks
      IF (SIZE(EXTERNAL_POTENTIAL, 1) < GRIDC%RC%NP) &
         CALL vtutor%bug("Size of external potential incompatible with grid.", "pot.F", 750)
      IF (SIZE(CVTOT, 2) > 1) &
         CALL vtutor%bug("Implementation only tested without spin polarization.", "pot.F", 752)
!
      CALL ZAXPY(GRIDC%RC%NP, (1._q,0._q), EXTERNAL_POTENTIAL, 1, CVTOT, 1)
!
    END SUBROUTINE ADD_EXTERNAL_POTENTIAL

    SUBROUTINE ADD_TO_EXTERNAL_POTENTIAL(POTENTIAL)
! add another contribution to existing external potential
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      COMPLEX(q), INTENT(IN) :: POTENTIAL(:)
!
      IF (.NOT.ALLOCATED(EXTERNAL_POTENTIAL)) THEN
         ALLOCATE(EXTERNAL_POTENTIAL(SIZE(POTENTIAL)))
         EXTERNAL_POTENTIAL = POTENTIAL
      ELSE
         IF (SIZE(POTENTIAL) /= SIZE(EXTERNAL_POTENTIAL)) &
            CALL vtutor%bug("Size of potentials inconsistent.", "pot.F", 769)
         CALL ZAXPY(SIZE(POTENTIAL), (1._q,0._q), POTENTIAL, 1, EXTERNAL_POTENTIAL, 1)
      END IF
!
    END SUBROUTINE ADD_TO_EXTERNAL_POTENTIAL

    SUBROUTINE RESET_EXTERNAL_POTENTIAL()
! clear the external potential
      IF (ALLOCATED(EXTERNAL_POTENTIAL)) DEALLOCATE(EXTERNAL_POTENTIAL)
    END SUBROUTINE RESET_EXTERNAL_POTENTIAL

!************************ SUBROUTINE POTXC  ****************************
!
!> This subroutine to calculate the XC-potential including gradient
!> corrections.
!>
!> This routine is required to calculate the partial core
!> corrections to the forces.
!>
!> On entry: \n
!>  CHTOT(:,1)    density \n
!>  CHTOT(:,2)    respectively CHTOT(:,2:4) contain the magnetization \n
!> On return (LNONCOLLINEAR=.FALSE.): \n
!>  CVTOT(:,1)    average potential \n
!>  CVTOT(:,2:4)  magnetic field
!
!***********************************************************************

      SUBROUTINE POTXC(KINEDEN,GRIDC,GRID_SOFT,SOFT_TO_C,INFO,XC,WDES, &
     &   LATT_CUR,CHDEN,CHTOT,DENCOR,CVTOT,MUTOT,IVDW,XDM,NPDDSC)

!! USE moffload_struct_def

      USE prec
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, only: fexcf, fexcp
      USE tau_mu, ONLY : tau_handle
      USE mpimy
      USE mgrid
      USE lattice
      USE base
      USE wave

      IMPLICIT NONE

      TYPE (tau_handle)  KINEDEN
      TYPE (xc_info)     XC
      TYPE (grid_3d)     GRIDC,GRID_SOFT
      TYPE (wavedes)     WDES
      TYPE (transit)     SOFT_TO_C
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR
      INTEGER IVDW,NPDDSC,ISP
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      COMPLEX(q) CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      COMPLEX(q), TARGET :: CHTOT(GRIDC%MPLWV,WDES%NCDIJ)
      COMPLEX(q)      DENCOR(GRIDC%RL%NP)
      COMPLEX(q), POINTER :: MUTOT(:,:)
! work arrays
      REAL(q)    XCSIF(3,3),XDM(NPDDSC),EXC,EXCG,RINPL,XCENC,XCENCG
      COMPLEX(q), ALLOCATABLE:: CWORK(:,:)
      COMPLEX(q) CVZERG,CVZERO
      COMPLEX(q), POINTER :: MUTOT_(:,:)

      COMPLEX(q), POINTER :: CTMP(:,:)

      

!$ACC ENTER DATA COPYIN(CHTOT) 
      IF (XC%LNO_AUG_XC) THEN
         ALLOCATE(CTMP(GRIDC%MPLWV,WDES%NCDIJ))
!$ACC ENTER DATA COPYIN(CHDEN) 
      ELSE
         CTMP => CHTOT
      ENDIF
!$ACC ENTER DATA CREATE(CTMP) 

      IF (XC%LMU) THEN
         IF (.NOT.ASSOCIATED(MUTOT)) &
            CALL vtutor%bug("POTXC: MUTOT not associated.", "pot.F", 848)
         MUTOT_ => MUTOT
!$ACC ENTER DATA CREATE(MUTOT_) 
!$ACC KERNELS PRESENT(MUTOT_) 
         MUTOT_=0
!$ACC END KERNELS
      ELSE
         ALLOCATE(MUTOT_(1,WDES%NCDIJ))
!$ACC ENTER DATA CREATE(MUTOT_) 
      ENDIF

!$ACC ENTER DATA CREATE(CVTOT) 
!$ACC KERNELS PRESENT(CVTOT) 
      CVTOT = 0
!$ACC END KERNELS

      ALLOCATE(CWORK(GRIDC%MPLWV,WDES%NCDIJ))
!$ACC ENTER DATA CREATE(CWORK) 

      IF (XC%LNO_AUG_XC) THEN
! We will use the soft charge density instead of the augmented (1._q,0._q)
!$ACC KERNELS PRESENT(CTMP) 
         CTMP=0
!$ACC END KERNELS
         DO ISP=1,WDES%NCDIJ
            CALL ADD_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CHDEN(1,ISP),CTMP(1,ISP))
            CALL SETUNB_COMPAT(CTMP(1,ISP),GRIDC)
            CALL FFT3D(CTMP(1,ISP),GRIDC,1)
         ENDDO
      ENDIF

! Bring the charge density to real space
      DO ISP=1,WDES%NCDIJ
         CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
      ENDDO

      IF (WDES%ISPIN==2) THEN

! get the charge and the total magnetization
         CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)

         IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
! gradient corrections to LDA
! unfortunately FEXCGS requires (up,down) density
! instead of (rho,mag)
            IF (.NOT.XC%LDOMETAGGA) CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
            CALL FEXCG(XC,2,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                 KINEDEN,CWORK,DENCOR,CVTOT,MUTOT_,XDM,IVDW,NPDDSC,SIZE(MUTOT_,1))
            IF (.NOT.XC%LDOMETAGGA) CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
         ENDIF

! add LDA part of potential
         CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
            CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
            CVZERO,EXC,XCENC,XCSIF, .TRUE.)
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

! get the proper direction vx = v0 + hat m delta v
         CALL MAG_DIRECTION(CTMP(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)

      ELSEIF (WDES%LNONCOLLINEAR) THEN

!-MM- gradient corrections in the noncollinear case are calculated
!     a bit differently than in the collinear case
         IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
            CALL FEXCG(XC,4,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                 KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDM,IVDW,NPDDSC,SIZE(MUTOT_,1))
         ENDIF

! FEXCF requires (rho,mag) density instead of (up,down)
         CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)
! add LDA part of potential
         CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
            CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
            CVZERO,EXC,XCENC,XCSIF, .TRUE.)
! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)
! get the proper direction vx = v0 + hat m delta v
         IF (.NOT.XC%LDOMETAGGA) CALL MAG_DIRECTION(CTMP(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
!-MM- end of changes to calculation of gga in noncollinear case

         IF (XC%LDOMETAGGA) THEN
            IF (XC%LMU) THEN
! mM: just until all is 1._q
!              WRITE(*,*) 'non-collinear mode for revTPSS not yet fully implemented, sorry ...'
!              CALL M_stop('VASP aborting ...'); stop
! mM
! MUTOT and CWORK come out as two component (up,down) quantities, where
! "up" and "down" are taken w.r.t. the local magnetization direction,
! this is now cast into (rho, mag) form.
               CALL MAG_DIRECTION_KINDENS(CTMP(1,1), KINEDEN%TAU(1,1), MUTOT_(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ, LATT_CUR)
! and rearranged from (rho,mag) to spinor representation
               CALL POT_FLIP_RL(MUTOT_(1,1), GRIDC, WDES%NCDIJ)
            ELSE
! CWORK(potential) comes out with (spin up and down) representation
! It needs to flip to (total, mag) representation
               CALL MAG_DIRECTION(CTMP(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
            ENDIF
! rearrange the potential in (spinor) representation
            CALL POT_FLIP_RL(CVTOT(1,1), GRIDC, WDES%NCDIJ)
         ENDIF

      ELSE

         IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
! gradient corrections to LDA
            CALL FEXCG(XC,1,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                 KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDM,IVDW,NPDDSC,SIZE(MUTOT_,1))
         ENDIF

! LDA part of potential
         CALL FEXCP(GRIDC,XC,LATT_CUR%OMEGA, &
              CTMP,DENCOR,CVTOT,CWORK,CVZERO,EXC,XCENC,XCSIF,.TRUE.)

      ENDIF

! Take the charge density back to reciprocal space
      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO
!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3D(CVTOT(1,ISP),GRIDC,-1)
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
      ENDDO

! Flip CVTOT from (up,down) to (total,magnetization) storage mode
! This is 1._q because the routine FORTAU needs the "total" component
      IF (WDES%LNONCOLLINEAR.AND.XC%LDOMETAGGA) CALL RC_FLIP_POTENTIAL(CVTOT,GRIDC,WDES%NCDIJ,.FALSE.)

!$ACC EXIT DATA DELETE(CWORK) 
      DEALLOCATE(CWORK)

      IF (XC%LMU) THEN
!-----------------------------------------------------------------------
! The same procedure with MUTOT
!-----------------------------------------------------------------------
         DO  ISP=1,WDES%NCDIJ
            CALL RL_ADD(MUTOT_(1,ISP),RINPL,MUTOT_(1,ISP),0._q,MUTOT_(1,ISP),GRIDC)
            CALL FFT3D(MUTOT_(1,ISP),GRIDC,-1)
            CALL SETUNB_COMPAT(MUTOT_(1,ISP),GRIDC)
         ENDDO
! Flip MUTOT from (up,down) to (total,magnetization) storage mode
! This is 1._q because the routine FORTAU needs the "total" component
         CALL RC_FLIP_POTENTIAL(MUTOT_,GRIDC,WDES%NCDIJ,.FALSE.)
!$ACC EXIT DATA COPYOUT(MUTOT_) 
      ELSE
!$ACC EXIT DATA DELETE(MUTOT_) 
         DEALLOCATE(MUTOT_)
      ENDIF
      NULLIFY(MUTOT_)

!$ACC EXIT DATA COPYOUT(CVTOT) 

!$ACC EXIT DATA DELETE(CHTOT,CTMP) 
      IF (XC%LNO_AUG_XC) THEN
!$ACC EXIT DATA DELETE(CHDEN) 
         DEALLOCATE(CTMP)
      ENDIF

      

      RETURN
      END SUBROUTINE POTXC


!> compute individual contributions to the potential and
!> write potential information to OUTCAR, LOCPOT, and HDF5 file
     SUBROUTINE WRITE_POTENTIAL(IO, KINEDEN, GRID, GRIDC, GRID_SOFT, WDES, &
            INFO, XC_DATA, P, T_INFO, COULOMB_POT, E, LATT_CUR, CHDEN, CHTOT, CSTRF, &
            CVTOT, DENCOR, SV, HAMILTONIAN, SOFT_TO_C, XCSIF, DYN, &
            LMDIM, N_MIX_PAW, CDIJ)

         USE base, ONLY: in_struct, info_struct, energy
         USE fileio, ONLY: OUTPOT
         USE hamil_struct_def, ONLY: ham_handle
         USE mdipol, ONLY: DIP, WRITE_VACUUM_LEVEL
         USE mgrid, ONLY: grid_3d, transit
         USE poscar, ONLY: latt, dynamics, type_info, OUTPOS
         USE pseudo_struct_def, ONLY: potcar
         USE setexm, ONLY: xc_info, XC_DATA_NOXC, PUSH_XC_TYPE
         USE tau_mu, ONLY: tau_handle
         USE wave_struct_def, ONLY: wavedes
         USE tutor, ONLY: vtutor
         USE pot_struct_def

         USE vhdf5, ONLY: GRP_RESULTS, GRP_POTENTIAL, IH5OUTFILEID, VH5_WRITE_POTENTIAL


         IMPLICIT NONE

         TYPE(in_struct), INTENT(IN)  :: IO        !< input output unit handle
! begin required for POTLOK
         TYPE(tau_handle), INTENT(IN) :: KINEDEN   !< kinetic energy density handle
         TYPE(grid_3d), INTENT(IN)    :: GRID      !< grid for wavefunctions
         TYPE(grid_3d), INTENT(IN)    :: GRIDC     !< grid for potentials/charge
         TYPE(grid_3d), INTENT(IN)    :: GRID_SOFT !< soft grid for potentials/charge
         TYPE(wavedes), INTENT(IN)    :: WDES      !< wavefunction descriptor
         TYPE(info_struct), INTENT(IN):: INFO      !< info about algo etc.
         TYPE(xc_info), INTENT(INOUT) :: XC_DATA   !< xc settings
         TYPE(potcar), INTENT(IN)     :: P(:)      !< structure to support radial grid
         TYPE(type_info), INTENT(IN)  :: T_INFO    !< ionic positions in superlattice
         TYPE(coulomb_potential), INTENT(IN) :: COULOMB_POT
         TYPE(energy), INTENT(INOUT)  :: E         !< energy handle for double counting, band energies, etc.
         TYPE(latt), INTENT(IN)       :: LATT_CUR  !< lattice parameter
         COMPLEX(q), INTENT(IN)       :: CHDEN(:,:)!< pseudo charge density on soft grid
         COMPLEX(q), INTENT(IN)       :: CHTOT(:,:)!< total charge density on fine grid
         COMPLEX(q), INTENT(IN)       :: CSTRF(:,:)!< structure factor
         COMPLEX(q), INTENT(INOUT)    :: CVTOT(:,:)!< local potential
         COMPLEX(q), INTENT(INOUT)         :: DENCOR(:) !< partial core charge
         COMPLEX(q), INTENT(INOUT)         :: SV(:,:)   !< total local potential on the soft grid
         TYPE(ham_handle), INTENT(IN) :: HAMILTONIAN !< handle for MGGA potential MU and MUTOT
         TYPE(transit), INTENT(IN)    :: SOFT_TO_C !< index table between GRID_SOFT and GRIDC
         REAL(q), INTENT(INOUT)       :: XCSIF(3,3)!< stress tensor from XC
! end required for POTLOK
         TYPE(dynamics), INTENT(IN)   :: DYN       !< contains position of ions
         INTEGER, INTENT(IN)          :: LMDIM     !< no. of augmentation channels
         INTEGER, INTENT(IN)          :: N_MIX_PAW !< number of elements of the augmentation occupancies mixed on the local node
         COMPLEX(q), INTENT(INOUT)       :: CDIJ(:,:,:,:) !< strength parameters
! local
         REAL(q) DUMMY(1)

         IF (WDES%COMM_KINTER%NODE_ME /= 1) RETURN

         IF (SHOULD_WRITE_POT(IO%LH5, IO%WRT_POTENTIAL)) THEN
            CALL OUTPOT_TO_POT("POT", WRITE_PAW=.TRUE.)
         ELSE IF (SHOULD_VASPWAVE_HAVE_POT(IO%LH5, IO%WRT_POTENTIAL)) THEN
            CALL OUTPOT_TO_VASPWAVE
         END IF

         CALL OUTPOT_TO_VASPOUT

         IF (SHOULD_WRITE_LOCPOT(IO%WRT_POTENTIAL)) &
            CALL OUTPOT_TO_LOCPOT

         IF (DIP%IDIPCO > 0 .AND. DIP%LCOR_DIP) &
            CALL OUTPUT_VACUUM_LEVEL

      CONTAINS


!> compute the local potential write it to LOCPOT
      SUBROUTINE OUTPOT_TO_LOCPOT
         IF (IO%WRT_POTENTIAL%LVHAR) CALL PUSH_XC_TYPE(XC_DATA,XC_DATA_NOXC)
         CALL POTLOK(KINEDEN,GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  INFO,XC_DATA,P,T_INFO,E,LATT_CUR, &
                  CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF,COULOMB_POT)
         CALL OUTPOT_TO_POT("LOCPOT", WRITE_PAW=.FALSE.)
      END SUBROUTINE OUTPOT_TO_LOCPOT


      SUBROUTINE OUTPUT_VACUUM_LEVEL
         IF (DIP%VACPOTAV .and. IO%WRT_POTENTIAL%LVTOT) THEN
            CALL vtutor%alert("You have requested that vacuum potentials are determined &
            &using the total potential (LVTOT = .True.), as opposed to just the &
            &hartree and ionic contributions. This selection requires large vacuum in your &
            &cell and a tightly converged calculations to work. Consider using &
            &LVHAR = .True. for a more rapid determination of the vacuum potentials.")
         ENDIF
! call the dipol routine without changing the potential
! It only makes sense to write out the vacuum potential
! if corrections to the potential are requested
         DIP%LCOR_DIP = .FALSE.
         CALL CDIPOL_CHTOT_REC(GRIDC, LATT_CUR,P,T_INFO, &
            CHTOT,CSTRF,CVTOT, WDES%NCDIJ, INFO%NELECT, E%PSCENC )
         IF (DIP%IDIPCO < 4) THEN
            CALL WRITE_VACUUM_LEVEL(IO%IU6)
         END IF
         DIP%LCOR_DIP = .TRUE.
      END SUBROUTINE OUTPUT_VACUUM_LEVEL


!> write the information to restart from a potential to vaspwave
      SUBROUTINE OUTPOT_TO_VASPWAVE

         USE pawm, ONLY: SET_RHO_PAW
         USE vhdf5, ONLY: IH5WAVEFILEID, VH5_WRITE_PAW_OCCUPANCIES, GRP_POTENTIAL
         REAL(q) DLM_EXX(N_MIX_PAW,WDES%NCDIJ)
         INTEGER SPIN
!
         CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)
!
         CALL VH5_WRITE_POTENTIAL(IH5WAVEFILEID, "total", GRIDC, CVTOT, GROUP="potential")
         CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
              CDIJ, DLM_EXX)
         DO SPIN = 1, WDES%NCDIJ
            CALL VH5_WRITE_PAW_OCCUPANCIES(IH5WAVEFILEID, GRP_POTENTIAL, P, T_INFO, &
               INFO%LOVERL, DLM_EXX(:,SPIN), GRIDC%COMM, SPIN)
         END DO
!
         CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)

      END SUBROUTINE OUTPOT_TO_VASPWAVE


!> compute the local potentials and write them to the HDF5 file
      SUBROUTINE OUTPOT_TO_VASPOUT

         COMPLEX(q), ALLOCATABLE :: CHARGE(:), WORK(:)
         INTEGER SPIN
!
         IF (SHOULD_VASPOUT_HAVE_TOTAL_POT(IO%WRT_POTENTIAL)) THEN
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "total", GRIDC, CVTOT)
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)
         END IF
!
         IF (SHOULD_VASPOUT_HAVE_XC_POT(IO%WRT_POTENTIAL)) THEN
            CALL POTXC(KINEDEN, GRIDC, GRID_SOFT, SOFT_TO_C, INFO, XC_DATA, WDES, &
               LATT_CUR, CHDEN, CHTOT, DENCOR, CVTOT, HAMILTONIAN%MUTOT, 0, DUMMY, 1)
            DO SPIN = 1, SIZE(CVTOT, 2)
               CALL SETUNB_COMPAT(CVTOT(:,SPIN), GRIDC)
               CALL FFT3D(CVTOT(:,SPIN), GRIDC, 1)
            END DO
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "xc", GRIDC, CVTOT)
         END IF
!
         IF (SHOULD_VASPOUT_HAVE_HARTREE_POT(IO%WRT_POTENTIAL)) THEN
            ALLOCATE(CHARGE, SOURCE=CHTOT(:,1))
            CALL FFT3D(CHARGE, GRIDC, 1)
            CALL FFT_RC_SCALE(CHARGE, CHARGE, GRIDC)
            CALL SETUNB_COMPAT(CHARGE, GRIDC)
            CALL POTHAR(GRIDC, LATT_CUR, COULOMB_POT, CHARGE, CVTOT, DUMMY(1))
            CALL SETUNB_COMPAT(CVTOT(:,1), GRIDC)
            CALL FFT3D(CVTOT(:,1), GRIDC, 1)
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "hartree", GRIDC, CVTOT(:,1:1))
         END IF
!
         IF (SHOULD_VASPOUT_HAVE_IONIC_POT(IO%WRT_POTENTIAL)) THEN
            ALLOCATE(WORK(MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV)))
            CALL POTION(GRIDC, P, LATT_CUR, T_INFO, COULOMB_POT, CVTOT, CSTRF, DUMMY(1))
            CALL SETUNB_COMPAT(CVTOT(:,1), GRIDC)
            CALL FFT3D(CVTOT(:,1), GRIDC, 1)
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "ionic", GRIDC, CVTOT(:,1:1))
         END IF

      END SUBROUTINE OUTPOT_TO_VASPOUT


!> small routine to write the potential in POT format
      SUBROUTINE OUTPOT_TO_POT(FILENAME, WRITE_PAW)
         USE main_mpi, ONLY: DIR_APP, DIR_LEN
         USE pawm, ONLY: SET_RHO_PAW, WRT_RHO_PAW
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: FILENAME
         LOGICAL, INTENT(IN) :: WRITE_PAW
         REAL(q) DLM_EXX(N_MIX_PAW,WDES%NCDIJ)
         INTEGER :: ISP, I, NODE_ME, IONODE

         NODE_ME = WDES%COMM%NODE_ME
         IONODE = WDES%COMM%IONODE

! the potential is complex for noncollinear calculations unless represented
! in terms of the Pauli matrices
         IF (WDES%NCDIJ == 4) &
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)

         IF (NODE_ME==IONODE) THEN

         IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE=DIR_APP(1:DIR_LEN) // FILENAME,STATUS='UNKNOWN')
         REWIND IO%IUVTOT
! write lattice parameters and positions to specified unit
         CALL OUTPOS(IO%IUVTOT,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSION)
         ENDIF

! at the moment the spin up and down potential is written to the file
         CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
              CDIJ, DLM_EXX)
         CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,CVTOT(:,1))
         IF (WRITE_PAW) &
            CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM_EXX(:,1), GRIDC%COMM, IO%IUVTOT)

         DO ISP = 2, WDES%NCDIJ
            IF (NODE_ME==IONODE) WRITE( IO%IUVTOT,'(5E20.12)') (T_INFO%ATOMOM(I),I=1,T_INFO%NIONS)
            CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,CVTOT(:,ISP))
            IF (WRITE_PAW) &
               CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM_EXX(:,ISP), GRIDC%COMM, IO%IUVTOT )
         END DO
         IF (NODE_ME==IONODE) THEN
         CLOSE(IO%IUVTOT)
         ENDIF

         IF (WDES%NCDIJ == 4) &
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)

      END SUBROUTINE OUTPOT_TO_POT

      END SUBROUTINE WRITE_POTENTIAL

      FUNCTION WRT_POTENTIAL_TO_STRING(WRT_POTENTIAL) RESULT (RES)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         CHARACTER(LEN=:), ALLOCATABLE :: RES
!
         IF (SHOULD_VASPOUT_HAVE_TOTAL_POT(WRT_POTENTIAL)) THEN
            RES = "total"
         ELSE
            RES = ""
         END IF
         IF (SHOULD_VASPOUT_HAVE_HARTREE_POT(WRT_POTENTIAL)) &
            RES = TRIM(RES // " hartree")
         IF (SHOULD_VASPOUT_HAVE_IONIC_POT(WRT_POTENTIAL)) &
            RES = TRIM(RES // " ionic")
         IF (SHOULD_VASPOUT_HAVE_XC_POT(WRT_POTENTIAL)) &
            RES = TRIM(RES // " xc")
         IF (WRT_POTENTIAL%PAW) &
            RES = TRIM(RES // " paw")
         IF (RES == "") &
            RES = "false"
      END FUNCTION WRT_POTENTIAL_TO_STRING

      LOGICAL FUNCTION SHOULD_WRITE_POT(LH5, WRT_POTENTIAL)
         USE chi_glb, ONLY: LOEP
         USE fock_glb, ONLY: EXXOEP
         LOGICAL, INTENT(IN) :: LH5
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_WRITE_POT = .NOT.LH5 .AND. WRT_POTENTIAL%LVTOT .AND. (LOEP .OR. EXXOEP > 0)
      END FUNCTION SHOULD_WRITE_POT

      LOGICAL FUNCTION SHOULD_VASPWAVE_HAVE_POT(LH5, WRT_POTENTIAL)
         USE chi_glb, ONLY: LOEP
         USE fock_glb, ONLY: EXXOEP
         LOGICAL, INTENT(IN) :: LH5
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPWAVE_HAVE_POT = LH5 .AND. WRT_POTENTIAL%LVTOT .AND. (LOEP .OR. EXXOEP > 0)
      END FUNCTION SHOULD_VASPWAVE_HAVE_POT

      LOGICAL FUNCTION SHOULD_WRITE_LOCPOT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_WRITE_LOCPOT = WRT_POTENTIAL%TOTAL .OR. WRT_POTENTIAL%LVTOT .OR. WRT_POTENTIAL%LVHAR
      END FUNCTION SHOULD_WRITE_LOCPOT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_TOTAL_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_TOTAL_POT = WRT_POTENTIAL%TOTAL .OR. WRT_POTENTIAL%LVTOT
      END FUNCTION SHOULD_VASPOUT_HAVE_TOTAL_POT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_HARTREE_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_HARTREE_POT = WRT_POTENTIAL%HARTREE .OR. WRT_POTENTIAL%LVHAR
      END FUNCTION SHOULD_VASPOUT_HAVE_HARTREE_POT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_IONIC_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_IONIC_POT = WRT_POTENTIAL%IONIC .OR. WRT_POTENTIAL%LVHAR
      END FUNCTION SHOULD_VASPOUT_HAVE_IONIC_POT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_XC_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_XC_POT = WRT_POTENTIAL%XC
      END FUNCTION SHOULD_VASPOUT_HAVE_XC_POT

    END MODULE


!************************ SUBROUTINE  MAG_DIRECTION  *******************
!
!> On entry CVTOT must contain the v_xc(up) and v_xc(down).
!>
!> On return CVTOT contains \n
!> Collinear case: \n
!>  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2 \n
!>  CVTOT(:,2) =  (v_xc(up) - v_xc(down))/2 \n
!> Non collinear case: \n
!>  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2 \n
!>  CVTOT(:,2) =  hat m_x (v_xc(up) - v_xc(down))/2 \n
!>  CVTOT(:,3) =  hat m_y (v_xc(up) - v_xc(down))/2 \n
!>  CVTOT(:,4) =  hat m_z (v_xc(up) - v_xc(down))/2 \n
!> where hat m is the unit vector of the local magnetization density.
!
!***********************************************************************



      SUBROUTINE MAG_DIRECTION(CHTOT, CVTOT, GRID, NCDIJ)

!! USE moffload_struct_def

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ

      COMPLEX(q) CHTOT(GRID%MPLWV, NCDIJ), &
            CVTOT(GRID%MPLWV, NCDIJ)
! local variables
      INTEGER K
      REAL(q) :: NORM2,DELTAV,V0

      IF (NCDIJ==2) THEN

!$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(V0,DELTAV) 
# 1347

         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
            CVTOT(K,1) = V0

            CVTOT(K,2) = DELTAV
# 1356

         ENDDO
      ELSE IF (NCDIJ==4) THEN
!$ACC PARALLEL LOOP PRESENT(CVTOT,CHTOT) PRIVATE(V0,DELTAV,NORM2) 
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
# 1365

            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CONJG(CHTOT(K,2))+ &
           &           CHTOT(K,3)*CONJG(CHTOT(K,3)) + CHTOT(K,4)*CONJG(CHTOT(K,4)))),1.E-20_q)


            CVTOT(K,1) = V0
            CVTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            CVTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            CVTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2
         ENDDO

      ELSE
         CALL vtutor%bug("internal error: MAG_DIRECTION called with NCDIJ= " // str(NCDIJ), "pot.F", 1377)
      ENDIF

      END SUBROUTINE MAG_DIRECTION


      SUBROUTINE MAG_DIRECTION_KINDENS(CHTOT, TAU, MUTOT, CVTOT, GRID, NCDIJ, LATT_CUR)

!! USE moffload_struct_def

      USE prec
      USE constant
      USE lattice
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      INTEGER NCDIJ

      COMPLEX(q) CHTOT(GRID%MPLWV, NCDIJ), TAU(GRID%MPLWV, NCDIJ), &
            MUTOT(GRID%MPLWV, NCDIJ), CVTOT(GRID%MPLWV, NCDIJ)
! local variables
      INTEGER K
      REAL(q) :: NORM2,DELTAV,V0,TAUPROJ,FACT
!#define debug3
# 1406


      IF (NCDIJ==2) THEN
! local potential

!$ACC PARALLEL PRESENT(GRID,CVTOT,MUTOT) 
# 1414

!$ACC LOOP GANG VECTOR PRIVATE(V0,DELTAV)
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
            CVTOT(K,1) = V0

            CVTOT(K,2) = DELTAV
# 1424

         ENDDO
! \mu = dE_xc / d\tau
!$ACC LOOP GANG VECTOR PRIVATE(V0,DELTAV)
         DO K=1,GRID%RL%NP
            V0    =(MUTOT(K,1)+MUTOT(K,2))/2
            DELTAV=(MUTOT(K,1)-MUTOT(K,2))/2
            MUTOT(K,1) = V0

            MUTOT(K,2) = DELTAV
# 1436

         ENDDO
!$ACC END PARALLEL
      ELSE IF (NCDIJ==4) THEN
! local potential and \mu = dE_xc / d\tau
!$ACC PARALLEL PRESENT(GRID,CVTOT,CHTOT,MUTOT,TAU) 
!$ACC LOOP GANG VECTOR PRIVATE(V0,DELTAV,NORM2)
         DO K=1,GRID%RL%NP
            V0    =REAL((CVTOT(K,1)+CVTOT(K,2))/2,KIND=q)
            DELTAV=REAL((CVTOT(K,1)-CVTOT(K,2))/2,KIND=q)
# 1448

!           NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CONJG(CHTOT(K,2))+ &
!          &           CHTOT(K,3)*CONJG(CHTOT(K,3)) + CHTOT(K,4)*CONJG(CHTOT(K,4)))),1.E-20_q)
            NORM2 = MAX(SQRT(ABS(REAL(CHTOT(K,2))*REAL(CHTOT(K,2))+ &
          &           REAL(CHTOT(K,3))*REAL(CHTOT(K,3))+REAL(CHTOT(K,4))*REAL(CHTOT(K,4)))),1.E-20_q)


            CVTOT(K,1) = V0
            CVTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            CVTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            CVTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2
         ENDDO
! \mu = dE_xc / d\tau
!$ACC LOOP GANG VECTOR PRIVATE(V0,DELTAV,NORM2,TAUPROJ,FACT)
         DO K=1,GRID%RL%NP
            V0    =REAL((MUTOT(K,1)+MUTOT(K,2))/2,KIND=q)
            DELTAV=REAL((MUTOT(K,1)-MUTOT(K,2))/2,KIND=q)
# 1467

!           NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CONJG(CHTOT(K,2))+ &
!          &           CHTOT(K,3)*CONJG(CHTOT(K,3)) + CHTOT(K,4)*CONJG(CHTOT(K,4)))),1.E-20_q)
            NORM2 = MAX(SQRT(ABS(REAL(CHTOT(K,2))*REAL(CHTOT(K,2))+ &
           &           REAL(CHTOT(K,3))*REAL(CHTOT(K,3))+REAL(CHTOT(K,4))*REAL(CHTOT(K,4)))),1.E-20_q)


            MUTOT(K,1) = V0
            MUTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            MUTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            MUTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2

# 1481

            TAUPROJ = (REAL(TAU(K,2),KIND=q)*REAL(CHTOT(K,2),KIND=q)+ &
           &            REAL(TAU(K,3),KIND=q)*REAL(CHTOT(K,3),KIND=q)+ &
           &             REAL(TAU(K,4),KIND=q)*REAL(CHTOT(K,4),KIND=q)) /NORM2/NORM2


            FACT=1._q/HSQDTM!*2._q
!           FACT=0._q
            IF (NORM2>1E-3_q) THEN
               CVTOT(K,2)=CVTOT(K,2)+DELTAV*REAL(TAU(K,2)-TAUPROJ*CHTOT(K,2),KIND=q)/NORM2*FACT
               CVTOT(K,3)=CVTOT(K,3)+DELTAV*REAL(TAU(K,3)-TAUPROJ*CHTOT(K,3),KIND=q)/NORM2*FACT
               CVTOT(K,4)=CVTOT(K,4)+DELTAV*REAL(TAU(K,4)-TAUPROJ*CHTOT(K,4),KIND=q)/NORM2*FACT
            ENDIF
# 1520

         ENDDO
!$ACC END PARALLEL
# 1529

      ELSE
         CALL vtutor%bug("internal error: MAG_DIRECTION called with NCDIJ= " // str(NCDIJ), "pot.F", 1531)
      ENDIF

      END SUBROUTINE MAG_DIRECTION_KINDENS


!************************ SUBROUTINE MAG_DENSITY ***********************
!
!> This subroutine calculates the total charge density and the
!> absolute magnitude of the magnetization density.
!>
!> On entry: \n
!>  CHTOT  rho, m_x, m_y, m_z \n
!> On exit: \n
!>  CWORK  rho, sqrt(m_x^2 + m_y^2 + m_z^2) \n
!> in the collinear case, it this means a simple copy CHTOT to CWORK.
!
!***********************************************************************

      SUBROUTINE MAG_DENSITY(CHTOT, CWORK, GRID, NCDIJ)

!! USE moffload_struct_def

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ

      COMPLEX(q) CHTOT(GRID%MPLWV, NCDIJ), &
            CWORK(GRID%MPLWV, NCDIJ)
! local
      INTEGER K


      IF (NCDIJ==2) THEN
!$ACC PARALLEL LOOP PRESENT(CWORK,CHTOT) 
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)

            CWORK(K,2)=CHTOT(K,2)
# 1575

         ENDDO
      ELSE IF (NCDIJ==4) THEN
!$ACC PARALLEL LOOP PRESENT(CWORK,CHTOT) 
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)
            CWORK(K,2)=SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3) + CHTOT(K,4)*CHTOT(K,4)))
         ENDDO
      ELSE
         CALL vtutor%bug("internal error: MAG_DENSITY called with NCDIJ= " // str(NCDIJ), "pot.F", 1584)
      ENDIF
      END SUBROUTINE MAG_DENSITY

!******************************** REMOVE_SOURCES_FROM_BXC  **********************************
!> @brief Remove sources from the XC B-field as proposed in DOI:10.1021/acs.jctc.7b01049
! CVTOT must be given in (charge,magnetization) convention in real space
! NB: Do this before the constraining potential and the external potential is applied.
!********************************************************************************************
      SUBROUTINE REMOVE_SOURCES_FROM_BXC(E, CVTOT, CHTOT, GRIDC, WDES, LATT_CUR)
        USE wave_struct_def, only: wavedes
        USE mgrid, only: grid_3d
        USE base, only: energy
        USE lattice
        IMPLICIT NONE
        TYPE (latt),   intent(in) :: LATT_CUR !< lattice parameter
        TYPE(wavedes), intent(in) :: WDES     !< wave-function descriptor in supercell
        TYPE(grid_3d), intent(in) :: GRIDC    !< real space grid
        COMPLEX(q),    intent(in) :: CHTOT(GRIDC%MPLWV,WDES%NCDIJ) !< charge density in (rho, mag)
        COMPLEX(q),    intent(inout) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ) !< xc potential in real space
        TYPE (energy), intent(inout) :: E        !< energy
! local
        REAL(q) :: XCENC_SPT, XCENC_BEFORE
 
! XCENC_BEFORE = -<rho*v0>-<m*Bxc>
        CALL SPT_DCC(CVTOT, CHTOT, GRIDC, WDES, LATT_CUR, XCENC_BEFORE)
! update CVTOT
        CALL REMOVE_SOURCES(CVTOT, GRIDC, WDES, LATT_CUR)
! XCENC_SPT = -<rho*v0>-<m*Bsfxc> - XCENC_BEFORE
        CALL SPT_DCC(CVTOT, CHTOT, GRIDC, WDES, LATT_CUR, XCENC_SPT)
        XCENC_SPT = XCENC_SPT-XCENC_BEFORE
! add double counting corrections due to spin torque to E%XCENC
        E%XCENC = E%XCENC + XCENC_SPT

      END SUBROUTINE REMOVE_SOURCES_FROM_BXC

!***************************    SPT_DCC    *****************************
! Compute double-counting corrections for a general noncollinear run.
!
!                     1  / rho+m_z  m_x-im_y \   / v_0+B_z  B_x-iBy \
! - 0.5 * Tr[n*v] = - – |                    | * |                  |
!                     2  \ m_x+im_y rho-m_z  /   \ B_x+iBy  v_0-B_z /
!
!                 = - rho*v_0 - m_x*B_x - m_y*B_y - m_z*B_z
!
! Specifically Bxc || m is not assumed, such that spin torque is accounted for.
!
!> @brief: Return double-counting corrections -v0*rho + Bxc*m for ncl case in XCENC_SPT
!***********************************************************************
      SUBROUTINE SPT_DCC(CVTOT, CHTOT, GRIDC, WDES, LATT_CUR, XCENC_SPT)
!! USE moffload_struct_def
      USE prec
      USE lattice
      USE mpimy
      USE mgrid, only : grid_3d
      USE wave_struct_def, only: wavedes
      USE constant

      IMPLICIT NONE
      TYPE(grid_3d), intent(in) :: GRIDC    !< real space grid
      TYPE(wavedes), intent(in) :: WDES     !< wave-function descriptor in supercell
      TYPE (latt),   intent(in) :: LATT_CUR !< lattice parameter
      COMPLEX(q),    intent(in) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ) !< xc potential in real space
      COMPLEX(q),    intent(in) :: CHTOT(GRIDC%MPLWV,WDES%NCDIJ) !< charge density in (rho, mag)
      REAL(q),    intent(inout) :: XCENC_SPT  !< -v0*rho + Bxc*m (double counting correction)
!-----------------------------------------------------------------
      INTEGER :: NODE_ME, IONODE, IDUMP, I
      REAL(q) :: CHGMIN=1E-10
      REAL(q) :: RINPL
      REAL(q) :: RHO, MX, MY, MZ
      REAL(q) :: V0, VX, VY, VZ
! set to (1._q,0._q) for error-dumps

      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
      IDUMP=0
# 1662

# 1668


      XCENC_SPT = 0._q
! Loop over local number of grid points in real space
!$ACC PARALLEL LOOP PRESENT(CHTOT,CVTOT) PRIVATE(RHO,MX,MY,MZ,V0,VX,VY,VZ) 
      DO I=1,GRIDC%RL%NP
        RHO = CHTOT(I,1)
        MX  = CHTOT(I,2)
        MY  = CHTOT(I,3)
        MZ  = CHTOT(I,4)

        V0  = CVTOT(I,1)
        VX  = CVTOT(I,2)
        VY  = CVTOT(I,3)
        VZ  = CVTOT(I,4)

        XCENC_SPT = XCENC_SPT - RHO*V0 - VX*MX - VY*MY - VZ*MZ
      ENDDO

      XCENC_SPT = XCENC_SPT/GRIDC%NPLWV
      CALL M_sum_d(GRIDC%COMM, XCENC_SPT, 1)

# 1695


      END SUBROUTINE SPT_DCC

!********************* SUBROUTINE REMOVE_SOURCES ***********************
!> @brief: Remove the sources from the xc B field
!   Assuming the Helmholtz decomposition reads
!
!   B(r) = - (1/4pi) grad PHI(r) + rot A(r).
!
!   Then, the magnetic potential is obtained solving the scalar
!   Poisson equation:
!
!   nabla^2 PHI(r) = 4pi div B(r)
!           PHI(G) = 4pi*imag (G dot B)/G^2
!
!   On the plane-wave grid, we find the source free B field using FFTs
!
!   Bsf(r) = B(r) + (1/4pi) grad PHI(r)
!   Bsf(G) = B(G) -          G * (G*B(G))/G^2
!
!   Finally, CVTOT is updated.
!***********************************************************************
      SUBROUTINE REMOVE_SOURCES(CVTOT, GRIDC, WDES, LATT_CUR)
!! USE moffload_struct_def
      USE prec
      USE lattice
      USE mpimy
      USE mgrid, only : grid_3d
      USE tutor, only: vtutor
      USE wave_struct_def, only: wavedes
      USE constant
      USE setexm

      IMPLICIT NONE
      TYPE(grid_3d), intent(in) :: GRIDC    !< real space grid
      TYPE(wavedes), intent(in) :: WDES     !< wave-function descriptor in supercell
      TYPE (latt),   intent(in) :: LATT_CUR !< lattice parameter
      COMPLEX(q), intent(inout) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ) !< xc potential in real space
!-----------------------------------------------------------------
      COMPLEX(q)   :: CWORK_B(GRIDC%MPLWV,WDES%NCDIJ-1) !< components of B
      COMPLEX(q)   :: CWORK1(GRIDC%MPLWV)   !< [G*B(G)]/G^2 (mag pot= 4pi*imag*CWORK1)
      INTEGER :: ISP, N1, N2, N3, NC, NI, I
      INTEGER :: NODE_ME, IONODE, IDUMP
      REAL(q) :: GX, GY, GZ, GSQU, GXYZ, SCALE

      
!$ACC ENTER DATA CREATE(CWORK_B,CWORK1) 
!
! compute div(B); inspired by xcspin.F FEXCGS_NONCOL_ddsc_
!
!$ACC KERNELS PRESENT(CWORK1) 
      CWORK1 = REAL(0.0, KIND=q)
!$ACC END KERNELS
! loop over Bx, By, Bz
      DO ISP=2, WDES%NCDIJ
! set CWORK_B to total real potential in real space
!$ACC PARALLEL LOOP PRESENT(CWORK_B,CVTOT) 
        DO I=1,GRIDC%MPLWV
          CWORK_B(I,ISP-1) = CVTOT(I,ISP)
        END DO
! set CWORK_B to total real potential in reciprocal space
        CALL FFT_RC_SCALE(CWORK_B(1,ISP-1),CWORK_B(1,ISP-1),GRIDC)
! set contribution of unbalanced lattic-vectors in potential to (0._q,0._q)
        CALL SETUNB_COMPAT(CWORK_B,GRIDC)
        CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK_B(1,ISP-1))

! loop over local number of grid points of reciprocal grid
!$ACC PARALLEL LOOP PRESENT(CWORK_B,CWORK1,LATT_CUR,GRIDC) PRIVATE(N1,NC,N2,N3,GXYZ) 
        DO I=1,GRIDC%RC%NP
! Obtain indices for reciprocal vector
          N1= MOD((I-1),GRIDC%RC%NROW) +1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
! Compute reciprocal vector GX, GY, GZ in Cartesian coordinates
! Note: LATT_CUR%B contains base of reciprocal lattice vectors without factor 2pi
          GXYZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(ISP-1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(ISP-1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(ISP-1,3))
! Compute div(B) = GX*Bx + GY*By + GZ*Bz
! Note: here we account for the factor TPI=2pi
          CWORK1(I) = CWORK1(I) + CWORK_B(I,ISP-1)*GXYZ*TPI
        END DO

      END DO

! inspired by POTHAR
      SCALE=1.0/TPI**2
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,NI,GX,GY,GZ,GSQU) &
!$ACC PRESENT(GRIDC,CWORK1,LATT_CUR) 
      DO NC=1,GRIDC%RC%NCOL
   N2= GRIDC%RC%I2(NC)
   N3= GRIDC%RC%I3(NC)
        DO N1=1,GRIDC%RC%NROW
!!     N2= GRIDC%RC%I2(NC)
!!     N3= GRIDC%RC%I3(NC)
! compute G^2
          NI=(NC-1)*GRIDC%RC%NROW+N1

          GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
          GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
          GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
          GSQU=GX**2+GY**2+GZ**2

! compute CWORK1 = G*B(G)/G^2
! if G^2 not 0 divide by G^2
          IF ((GRIDC%LPCTX(N1)==0).AND.(GRIDC%LPCTY(N2)==0).AND.(GRIDC%LPCTZ(N3)==0)) THEN
            CWORK1(NI)=(0.0_q,0.0_q)
          ELSE
            CWORK1(NI)=CWORK1(NI)/GSQU*SCALE
          ENDIF
        END DO
      END DO

! loop over Bx, By, Bz
      DO ISP=2, WDES%NCDIJ

! loop over reciprocal grid
!$ACC PARALLEL LOOP PRESENT(CWORK_B) PRIVATE(N1,NC,N2,N3,GXYZ) 
        DO I=1,GRIDC%RC%NP
! Compute reciprocal vector GX, GY, GZ
          N1= MOD((I-1),GRIDC%RC%NROW) +1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GXYZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(ISP-1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(ISP-1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(ISP-1,3))

! B_sf(G) = B(G) - G * [G*B(G)]/G^2
          CWORK_B(I,ISP-1) = CWORK_B(I,ISP-1)  - TPI*GXYZ*CWORK1(I)
        END DO
! Set CWORK to source free potential in real space
        CALL FFT3D(CWORK_B(1,ISP-1),GRIDC,1)
! Update CVTOT with source free potential
!$ACC PARALLEL LOOP PRESENT(CWORK_B,CVTOT) 
        DO I=1,GRIDC%MPLWV
           CVTOT(I,ISP) = CWORK_B(I,ISP-1)
        END DO
      END DO

!$ACC EXIT DATA DELETE(CWORK1,CWORK_B) 
      

      END SUBROUTINE REMOVE_SOURCES

!************************ SUBROUTINE POT_FLIP **************************
!
!> Rearranges the storage mode for spin components of potentials.
!>
!> for the collinear case calculate: \n
!>  v0 1 + v_z \n
!>  v0 1 - v_z \n
!> for the non collinear case calculate \n
!>  v  = v0 1 + sigma_x v_x + simga_y v_y + sigma_z v_z
!
!***********************************************************************

      SUBROUTINE POT_FLIP(CVTOT, GRID, NCDIJ)

!! USE moffload_struct_def

      USE prec
      USE mgrid
      IMPLICIT NONE
      INTEGER NCDIJ
      TYPE (grid_3d)     GRID
      COMPLEX(q) :: CVTOT(GRID%MPLWV, NCDIJ)

! local
      COMPLEX(q) :: C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER K

      IF (NCDIJ==2) THEN
         FAC=1.0_q
!$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CZ) 
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CZ =CVTOT(K,2)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         FAC=1.0_q
!$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CX,CY,CZ) 
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CX =CVTOT(K,2)
            CY =CVTOT(K,3)
            CZ =CVTOT(K,4)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CVTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CVTOT(K,4)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==1) THEN
      ENDIF

    END SUBROUTINE POT_FLIP

    SUBROUTINE POT_FLIP_RL(CVTOT, GRID, NCDIJ)

!! USE moffload_struct_def

      USE prec
      USE mgrid
      IMPLICIT NONE
      INTEGER NCDIJ
      TYPE (grid_3d)     GRID

      COMPLEX(q) CVTOT(GRID%MPLWV, NCDIJ)
! local
      COMPLEX(q) :: C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER K

      IF (NCDIJ==2) THEN
         FAC=1.0_q
!$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CZ) 
         DO K=1,GRID%RL%NP
            C00=CVTOT(K,1)
            CZ =CVTOT(K,2)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         FAC=1.0_q
!$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CX,CY,CZ) 
         DO K=1,GRID%RL%NP
            C00=CVTOT(K,1)
            CX =CVTOT(K,2)
            CY =CVTOT(K,3)
            CZ =CVTOT(K,4)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CVTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CVTOT(K,4)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==1) THEN
      ENDIF

    END SUBROUTINE POT_FLIP_RL


!************************ SUBROUTINE EXTERNAL_POT **********************
!
!> This subroutine can be used to add an external potential.
!>
!> The units of the potential are eV.
!
!***********************************************************************

      SUBROUTINE EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT)

!! USE moffload_struct_def

      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q)      CVTOT(GRIDC%MPLWV)

      RETURN

      IF (GRIDC%RL%NFAST==3) THEN
! mpi version: x-> N2, y-> N3, z-> N1
         N2MAX=GRIDC%NGX
         N3MAX=GRIDC%NGY
         N1MAX=GRIDC%NGZ

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,NG) &
!$ACC& PRESENT(GRIDC,CVTOT) 
         DO NC=1,GRIDC%RL%NCOL
       N2= GRIDC%RL%I2(NC)
       N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
!!          N2= GRIDC%RL%I2(NC)
!!          N3= GRIDC%RL%I3(NC)
               NG=(NC-1)*GRIDC%RL%NROW+N1
               CVTOT(NG)=CVTOT(NG)
               IF (N1<=GRIDC%NGZ/2) THEN
                  CVTOT(NG)=CVTOT(NG)+(REAL(N1,q)-1)/GRIDC%NGZ*2*0.1
               ELSE
                  CVTOT(NG)=CVTOT(NG)+(GRIDC%NGZ-REAL(N1,q))/GRIDC%NGZ*2*0.1
               ENDIF
            ENDDO
         ENDDO
      ELSE
! conventional version: x-> N1, y-> N2, z-> N3
         N1MAX=GRIDC%NGX
         N2MAX=GRIDC%NGY
         N3MAX=GRIDC%NGZ

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,NG) &
!$ACC& PRESENT(GRIDC,CVTOT) 
         DO NC=1,GRIDC%RL%NCOL
       N2= GRIDC%RL%I2(NC)
       N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
!!          N2= GRIDC%RL%I2(NC)
!!          N3= GRIDC%RL%I3(NC)
               NG=(NC-1)*GRIDC%RL%NROW+N1
               CVTOT(NG)=CVTOT(NG)
               IF (N3<=GRIDC%NGZ/2) THEN
                  CVTOT(NG)=CVTOT(NG)+(REAL(N3,q)-1)/GRIDC%NGZ*2*0.1
               ELSE
                  CVTOT(NG)=CVTOT(NG)+(GRIDC%NGZ-REAL(N3,q))/GRIDC%NGZ*2*0.1
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE
