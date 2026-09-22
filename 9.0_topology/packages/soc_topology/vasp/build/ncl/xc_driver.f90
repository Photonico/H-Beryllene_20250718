# 1 "xc_driver.F"
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


# 2 "xc_driver.F" 2 
# 4

    MODULE xc_driver

    CONTAINS
!************************ SUBROUTINE GGA_DRIVER **********************************
!
!> This subroutine calls the subroutines that calculate xc quantities
!> (non-spin-polarized GGA) for the energy, potential, forces, stress tensor, etc.
!>
!> LLDA allows to include the LDA contribution directly in this
!> routine (not implemented for all subroutines, see below).
!> For Libxc, the LDA contributions are directly included.
!
!*********************************************************************************

      SUBROUTINE GGA_DRIVER(XC,D,DD,EXC,EXCD,EXCDD,LLDA)
!$ACC ROUTINE SEQ
      USE prec
      USE constant
      USE setexm_struct_def
      USE xc_name
      USE xc_family_name
      USE ldalib
      USE ggalib
# 30


      IMPLICIT NONE

      TYPE (xc_info) XC

      INTEGER :: ider, IXC

      REAL(q) THRD, DTHRD
      REAL(q) D, DD, D1, D2, DD1, DD2, DELTA, DDELTA
      REAL(q) FK, RS, S, SIGMA, SK, T, ukfactor
      REAL(q) RHO1, RHO2, DRHO1, DRHO2
      REAL(q) EC, ECD, ECDD, ECLDA, ECDLDA, EXLDA, EXDLDA
      REAL(q) EXDLDA_LR, EXDLDA_SR, EXLDA_LR, EXLDA_SR, EXWPBE_SR, EXWPBED_LR, EXWPBE_LR, EXWPBED_SR, EXWPBEDD_LR, EXWPBEDD_SR
      REAL(q) EXC, EXC1, EXC2, EXCD, EXCD1, EXCD2, EXCDD, EXCDD1, EXCDD2, EXCDDA
      REAL(q) FXC, FXCLDA, FXCD1, FXCD2, FXCDD1, FXCDD2
      REAL(q) VXC, VXCLDA1, VXCLDA2
      REAL(q) EXC_TMP, EXCD_TMP, EXCDD_TMP
      REAL(q) ALDAX,ALDAC,AGGAX,AGGAC,LDASCREEN
# 52

# 55

      PARAMETER (DDELTA=1E-4_q)
      PARAMETER (THRD=1._q/3._q)
      LOGICAL LLDA,LUSE_LONGRANGE_HF,LUSE_THOMAS_FERMI,LUSE_MODEL_HF
      LOGICAL:: FORCE_PBE=.FALSE.
!$ACC ROUTINE (GGA91_WB) SEQ
!$ACC ROUTINE (EXCHPBE) SEQ
!$ACC ROUTINE (CALC_EXCHWPBE_SP) SEQ
!$ACC ROUTINE (CORunspPBE) SEQ
!$ACC ROUTINE (B3LYPXCS) SEQ
!$ACC ROUTINE (AM05) SEQ
!$ACC ROUTINE (EXCHPBESOL) SEQ
!$ACC ROUTINE (CALC_EXCHWPBEsol_SP) SEQ
!$ACC ROUTINE (CORunspPBESOL) SEQ
!$ACC ROUTINE (GGAEALL) SEQ
# 75


      ALDAX=XC%ALDAX
      ALDAC=XC%ALDAC
      AGGAX=XC%AGGAX
      AGGAC=XC%AGGAC
      LDASCREEN=XC%LDASCREEN
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LUSE_THOMAS_FERMI=XC%LUSE_THOMAS_FERMI
      LUSE_MODEL_HF=XC%LUSE_MODEL_HF

      EXC_TMP   = 0._q
      EXCD_TMP  = 0._q
      EXCDD_TMP = 0._q

      NXCLOOP1: DO IXC=1,XC%NXC

      IF ((XC%FAMILY(IXC)==XC_FAM_GGA).OR.((XC%ID(IXC)==ID_XC_LIBXC).AND.(XC%FAMILY(IXC)==XC_FAM_LDA))) THEN

      EXC   = 0._q
      EXCD  = 0._q
      EXCDD = 0._q

      IF (XC%ID(IXC)==ID_XC_PW91) THEN

! PW91 using the routines of Bird and White

        CALL GGA91_WB(XC,D,DD,EXC,EXCD,EXCDD)
        EXC = 2*EXC / D
        EXCD =2*EXCD
        EXCDD=2*EXCDD
!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
!vdw jk
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE).OR. &
             & (XC%ID(IXC)==ID_XC_RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_REVPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB88PBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB86BPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PW86RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_CX)) THEN
!vdw jk

! Perdew Burke Ernzerhof and revised functional

        IF (XC%ID(IXC)==ID_XC_PBE) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
           T=DD/(D*SK*2._q)
        ELSE
           S=0.0_q
           T=0.0_q
        ENDIF
! Iann Gerber: range separating in GGA exchange
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD, &
          &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
        ELSE
           CALL CALC_EXCHWPBE_SP(XC,D,S, &
          &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
          &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
        ENDIF
! Iann Gerber: end modification

        CALL CORunspPBE(RS,ECLDA,ECDLDA,SK, &
             T,EC,ECD,ECDD,.TRUE.)

!        WRITE(*,'(10F14.7)') D,S,(EXC-EXLDA)/D,EXCD-EXDLDA,EXCDD
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD,ECDD

        IF (LLDA) THEN
! Add LDA contributions as well
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX+EXLDA/D*ALDAX+ECLDA  ! eventuell noch mit ALDAC multi
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX+EXDLDA*ALDAX +ECDLDA ! eventuell noch mit ALDAC multi
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
! gk replaced ECLD by ECLDA, 18.09.2018
                 EXC  =EC  *AGGAC+ EXWPBE_LR/D*AGGAX +EXWPBE_SR/D +ECLDA! ist EXWPBE nur diff
                 EXCD =ECD *AGGAC+ EXWPBED_LR*AGGAX+EXWPBED_SR+ECDLDA   !
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR !
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC +EXWPBE_LR/D*AGGAX +ECLDA
                 EXCD =ECD *AGGAC +EXWPBED_LR*AGGAX+ECDLDA
                 EXCDD=ECDD*AGGAC +EXWPBEDD_LR*AGGAX
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+ EXWPBE_SR/D*AGGAX +EXWPBE_LR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_SR*AGGAX+EXWPBED_LR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ELSE
! Do not add LDA contributions
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX  +(EXWPBE_SR-EXLDA_SR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX+(EXWPBED_SR-EXDLDA_SR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR-EXLDA_SR)/D*AGGAX  +(EXWPBE_LR-EXLDA_LR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_SR-EXDLDA_SR)*AGGAX+(EXWPBED_LR-EXDLDA_LR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ENDIF

! Hartree -> Rydberg conversion
        EXC = 2._q*EXC
        EXCD =2._q*EXCD
        EXCDD=2._q*EXCDD

!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_X) THEN

! exchange-only PBE-type functionals

        IF (XC%ID(IXC)==ID_XC_PBE_X) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
        ELSE
           S=0.0_q
        ENDIF
        CALL EXCHPBE(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD,ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

        IF (LLDA) THEN
           EXC  =(EXC-EXLDA)/D*AGGAX+EXLDA/D*ALDAX
           EXCD =(EXCD-EXDLDA)*AGGAX+EXDLDA*ALDAX
           EXCDD=EXCDD*AGGAX
        ELSE
           EXC  =(EXC-EXLDA)/D*AGGAX
           EXCD =(EXCD-EXDLDA)*AGGAX
           EXCDD=EXCDD*AGGAX
        ENDIF

! Hartree -> Rydberg conversion
        EXC = 2._q*EXC
        EXCD =2._q*EXCD
        EXCDD=2._q*EXCDD

      ELSE IF (XC%ID(IXC)==ID_XC_PBE_C) THEN

! correlation-only PBE-type functionals

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           T=DD/(D*SK*2._q)
        ELSE
           T=0.0_q
        ENDIF

        CALL CORunspPBE(RS,ECLDA,ECDLDA,SK,T,EC,ECD,ECDD,.TRUE.)

        IF (LLDA) THEN
           EXC  =EC  *AGGAC + ECLDA
           EXCD =ECD *AGGAC + ECDLDA
           EXCDD=ECDD*AGGAC
        ELSE
           EXC  =EC  *AGGAC
           EXCD =ECD *AGGAC
           EXCDD=ECDD*AGGAC
        ENDIF

! Hartree -> Rydberg conversion
        EXC = 2._q*EXC
        EXCD =2._q*EXCD
        EXCDD=2._q*EXCDD

      ELSE IF ((XC%ID(IXC)==ID_XC_B3).OR.(XC%ID(IXC)==ID_XC_B5)) THEN

        RHO1=D/2.0_q
        RHO2=D/2.0_q
        DRHO1=DD/2.0_q
        DRHO2=DD/2.0_q

        CALL  B3LYPXCS(XC,RHO1,RHO2,DRHO1,DRHO2,DD,EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA)


        EXC     = 2._q* EXC/D
        EXCD1   = 2._q* EXCD1
        EXCD2   = 2._q* EXCD2
        EXCDD1  = 2._q* EXCDD1
        EXCDD2  = 2._q* EXCDD2
        EXCDDA  = 2._q* EXCDDA

        EXCD    = EXCD1
        EXCDD   = EXCDD1+EXCDDA

!       WRITE(91,'(7F14.7)') RHO1, RHO2, DRHO1, DRHO2, EXC, EXCD, EXCDD

! jP
!aem AM05 added

      ELSE IF (XC%ID(IXC)==ID_XC_AM05) THEN

!aem DD sometimes comes in negative.
!aem Since AM05 assumes that DD actually IS |grad rho|,
!aem and thus DD always should be positive,
!aem we need to use ABS(DD) as input to AM05.
!aem However, my routine needs to give the same symmetries out as
!aem PBE and PW91, therefor the sign-compensation of EXCDD.

        D1=D/2._q
        D2=D/2._q
        DD1=DD/2._q
        DD2=DD/2._q
        ider = 1
        CALL AM05(XC,D1,D2,ABS(DD1),ABS(DD2),ider, &
             FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
             FXCDD1,FXCDD2)

        IF (LLDA) THEN
          EXC = FXC/D
          EXCD = FXCD1
        ELSE
          EXC = (FXC - FXCLDA)/D
          EXCD = FXCD1 - VXCLDA1
        ENDIF
        EXCDD = SIGN(1.0_q,DD)*FXCDD1

! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD

      ELSE IF (XC%ID(IXC)==ID_XC_PBESOL) THEN

! Perdew Burke Ernzerhof and revised functional

           ukfactor=1.0_q

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
           T=DD/(D*SK*2._q)
        ELSE
           S=0.0_q
           T=0.0_q
        ENDIF
! Iann Gerber: range separating in GGA exchange
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(XC,D,S, &
          &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
          &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
        ENDIF
! Iann Gerber: end modification

        CALL CORunspPBESOL(RS,ECLDA,ECDLDA,SK, &
             T,EC,ECD,ECDD,.TRUE.)

!        WRITE(*,'(10F14.7)') D,S,(EXC-EXLDA)/D,EXCD-EXDLDA,EXCDD
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD,ECDD

        IF (LLDA) THEN
! Add LDA contributions as well
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX+EXLDA/D*ALDAX+ECLDA
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX+EXDLDA*ALDAX +ECDLDA
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+ EXWPBE_LR/D*AGGAX +EXWPBE_SR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_LR*AGGAX+EXWPBED_SR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC +EXWPBE_LR/D*AGGAX +ECLDA
                 EXCD =ECD *AGGAC +EXWPBED_LR*AGGAX+ECDLDA
                 EXCDD=ECDD*AGGAC +EXWPBEDD_LR*AGGAX
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+ EXWPBE_SR/D*AGGAX +EXWPBE_LR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_SR*AGGAX+EXWPBED_LR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ELSE
! Do not add LDA contributions
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX  +(EXWPBE_SR-EXLDA_SR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX+(EXWPBED_SR-EXDLDA_SR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX
                 EXCDD=ECDD*AGGAC+EXWPBEDD_LR*AGGAX
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR-EXLDA_SR)/D*AGGAX  +(EXWPBE_LR-EXLDA_LR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_SR-EXDLDA_SR)*AGGAX+(EXWPBED_LR-EXDLDA_LR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ENDIF

! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD

!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
! jP: end of adding PBEsol

! BEEF
# 464


      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_03).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_05).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_10).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_20)) THEN
! functionals for range-separated ACFDT (LDA - short range RPA)
! dummy stub since these functionals are not gradient corrected
         EXC=0._q ; EXCD=0._q ; EXCDD=0._q

      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_RA).OR.(XC%ID(IXC)==ID_XC_ACFDT_PL)) THEN
         EXC=0._q ; EXCD=0._q ; EXCDD=0._q

# 534

      ELSE

! for all other functionals
! we use finite differences to calculate the required
! quantities
! presently no other functional are supported
! but in case (1._q,0._q) needs these routines

        DELTA=MIN(DDELTA,ABS(D)/100)
        D1=D-DELTA
        D2=D+DELTA
        CALL GGAEALL(D1,DD,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(D2,DD,1._q,1._q,VXC,EXC2)
        EXCD=(EXC2*D2-EXC1*D1)/MAX((D2-D1),1E-10_q)

        DELTA=MIN(DDELTA,ABS(DD)/100)
        DD1=DD-DELTA
        DD2=DD+DELTA
        CALL GGAEALL(D,DD1,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(D,DD2,1._q,1._q,VXC,EXC2)
        EXCDD=(EXC2*D-EXC1*D)/MAX((DD2-DD1),1E-10_q)
        CALL GGAEALL(D,DD,1._q,1._q,VXC,EXC)

!aem added in case somebody (like me) might want to use this routine
! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD
!aem end addition

      ENDIF

      EXC_TMP   = EXC_TMP   + XC%COEFF(IXC)*EXC
      EXCD_TMP  = EXCD_TMP  + XC%COEFF(IXC)*EXCD
      EXCDD_TMP = EXCDD_TMP + XC%COEFF(IXC)*EXCDD

      ENDIF

      ENDDO NXCLOOP1

      EXC   = EXC_TMP
      EXCD  = EXCD_TMP
      EXCDD = EXCDD_TMP

      RETURN
      END SUBROUTINE GGA_DRIVER

# 926



!************************ SUBROUTINE GGASPIN_DRIVER ****************************
!
!> This subroutine calls the subroutines that calculate xc quantities
!> (spin-polarized GGA) for the energy, potential, forces, stress tensor, etc.
!>
!> LLDA allows to include the LDA contribution directly in this
!> routine (not implemented for all subroutines, see below).
!> For Libxc, the LDA contributions are directly included.
!
!*******************************************************************************

      SUBROUTINE GGASPIN_DRIVER(XC,D1,D2,DD1,DD2,DDA,EXC, &
         excd1,excd2,excq1,excq2,ecq,LLDA)
!$ACC ROUTINE SEQ

!     D1   density up
!     D2   density down
!     DD1  |gradient of density up|
!     DD2  |gradient of density down|
!     DDA  |gradient of the total density|
!     LLDA add lda contributions

      USE prec
      USE constant
      USE setexm_struct_def
      USE xc_name
      USE xc_family_name
      USE tutor, ONLY : vtutor
      USE ggalib
# 960

      IMPLICIT NONE

      TYPE (xc_info) XC

      INTEGER :: ider, IXC

      REAL(q) THRD, DTHRD, THRD4
      REAL(q) D, D1, D2, DD1, DD2, DDA
      REAL(q) G, FK, RS, S, SK, T, ukfactor, ZETA
      REAL(q) ECQ, EC, ECD1, ECD1LDA, ECD2, ECD2LDA, ECLDA, EXDLDA1, EXDLDA2, EXLDA1, EXLDA2
      REAL(q) EXDLDA_LR1, EXDLDA_LR2, EXDLDA_SR1, EXDLDA_SR2, EXLDA_LR1, EXLDA_LR2, EXLDA_SR1, EXLDA_SR2
      REAL(q) EXWPBE_LR1, EXWPBE_LR2, EXWPBE_SR1, EXWPBE_SR2, EXWPBED_LR1, EXWPBED_LR2
      REAL(q) EXWPBED_SR1, EXWPBED_SR2, EXWPBEDD_LR1, EXWPBEDD_LR2, EXWPBEDD_SR1, EXWPBEDD_SR2
      REAL(q) EXC, EXCD1, EXCD2, EXCQ1, EXCQ2, EXC1, EXC2, EXCDD1, EXCDD2, EXCDDA
      REAL(q) FXC, FXCD1, FXCD2, FXCDD1, FXCDD2, FXCLDA
      REAL(q) VXCLDA1, VXCLDA2
      REAL(q) EXC_TMP, EXCD1_TMP, EXCD2_TMP, EXCQ1_TMP, EXCQ2_TMP, ECQ_TMP
      REAL(q) ALDAX,ALDAC,AGGAX,AGGAC,LDASCREEN
# 982

# 985

      LOGICAL LLDA
      LOGICAL LUSE_LONGRANGE_HF,LUSE_THOMAS_FERMI,LUSE_MODEL_HF
      PARAMETER (THRD=1._q/3._q,THRD4=4._q/3._q)
      LOGICAL:: FORCE_PBE=.FALSE.
!$ACC ROUTINE (GGAXCS) SEQ
!$ACC ROUTINE (EXCHPBE) SEQ
!$ACC ROUTINE (CALC_EXCHWPBE_SP) SEQ
!$ACC ROUTINE (CORPBE) SEQ
!$ACC ROUTINE (B3LYPXCS) SEQ
!$ACC ROUTINE (AM05) SEQ
!$ACC ROUTINE (EXCHPBESOL) SEQ
!$ACC ROUTINE (EXCHPBESOL) SEQ
!$ACC ROUTINE (CALC_EXCHWPBEsol_SP) SEQ
!$ACC ROUTINE (corpbesol) SEQ
# 1005


      ALDAX=XC%ALDAX
      ALDAC=XC%ALDAC
      AGGAX=XC%AGGAX
      AGGAC=XC%AGGAC
      LDASCREEN=XC%LDASCREEN
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LUSE_THOMAS_FERMI=XC%LUSE_THOMAS_FERMI
      LUSE_MODEL_HF=XC%LUSE_MODEL_HF

      EXC_TMP   = 0._q
      EXCD1_TMP = 0._q
      EXCD2_TMP = 0._q
      EXCQ1_TMP = 0._q
      EXCQ2_TMP = 0._q
      ECQ_TMP   = 0._q

      NXCLOOP4: DO IXC=1,XC%NXC

      IF ((XC%FAMILY(IXC)==XC_FAM_GGA).OR.((XC%ID(IXC)==ID_XC_LIBXC).AND.(XC%FAMILY(IXC)==XC_FAM_LDA))) THEN

      EXC   = 0._q
      EXCD1 = 0._q
      EXCD2 = 0._q
      EXCQ1 = 0._q
      EXCQ2 = 0._q
      ECQ   = 0._q

      IF (XC%ID(IXC)==ID_XC_PW91) THEN
         CALL GGAXCS(D1,D2,DD1,DD2,EXC,EXCD1,EXCD2,EXCQ1,EXCQ2)

         D = D1 + D2
         EXC    = 2._q* EXC / D
         EXCD1  = 2._q* EXCD1
         EXCD2  = 2._q* EXCD2
         EXCQ1  = 2._q* EXCQ1
         EXCQ2  = 2._q* EXCQ2
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE).OR. &
             & (XC%ID(IXC)==ID_XC_RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_REVPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB88PBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB86BPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PW86RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_CX).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_RA).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_PL).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_03).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_05).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_10).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_20)) THEN

        IF (XC%ID(IXC)==ID_XC_PBE) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        D=2*D1



# 1070

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
# 1078


        ELSE
           S=0.0_q
        ENDIF
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
          &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
        ELSE
           CALL CALC_EXCHWPBE_SP(XC,D,S, &
          &  EXLDA1,EXDLDA1,EXLDA_SR1,EXLDA_LR1,EXDLDA_SR1,EXDLDA_LR1, &
          &  EXWPBE_SR1,EXWPBE_LR1,EXWPBED_SR1,EXWPBED_LR1,EXWPBEDD_SR1,EXWPBEDD_LR1)
        ENDIF
!        WRITE(*,'(10F14.7)') D,S,(EXC1-EXLDA1)/D,EXCD1-EXDLDA1,EXCQ1
!        EXLDA1=0; EXC1=0; EXDLDA1=0; EXCD1=0; EXCQ1=0

        D=2*D2
# 1097


        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
# 1105

        ELSE
           S=0.0_q
        ENDIF

        IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
          &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
        ELSE
           CALL CALC_EXCHWPBE_SP(XC,D,S, &
          &  EXLDA2,EXDLDA2,EXLDA_SR2,EXLDA_LR2,EXDLDA_SR2,EXDLDA_LR2, &
          &  EXWPBE_SR2,EXWPBE_LR2,EXWPBED_SR2,EXWPBED_LR2,EXWPBEDD_SR2,EXWPBEDD_LR2)
        ENDIF
!        EXLDA2=0; EXC2=0; EXDLDA2=0; EXCD2=0; EXCQ2=0
!        WRITE(*,'(10F14.7)') D,S,(EXC2-EXLDA2)/D,EXCD2-EXDLDA2,EXCQ2

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)
!        ECLDA=0 ; ECD1LDA=0 ; ECD2LDA=0; EC=0 ; ECD1=0 ; ECD2=0 ; ECQ=0
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD1,ECD2,ECQ

        IF (LLDA) THEN
! Add LDA contributions as well
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*ALDAX+ECLDA
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*ALDAX+ECD1LDA
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*ALDAX+ECD2LDA
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX +(EXWPBE_SR1+EXWPBE_SR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX +EXWPBED_SR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX +EXWPBED_SR2+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ*AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ*AGGAC
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1+EXWPBE_SR2)/(2*D)*AGGAX +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_SR1*AGGAX +EXWPBED_LR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_SR2*AGGAX +EXWPBED_LR2+ECD2LDA
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ*AGGAC
              ENDIF
           ENDIF
        ELSE
! Do not add LDA contribution
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC*  AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX+(EXWPBED_SR1-EXDLDA_SR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX+(EXWPBED_SR2-EXDLDA_SR2)
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ *AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ *AGGAC
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_SR1-EXDLDA_SR1)*AGGAX+(EXWPBED_LR1-EXDLDA_LR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_SR2-EXDLDA_SR2)*AGGAX+(EXWPBED_LR2-EXDLDA_LR2)
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ *AGGAC
              ENDIF
           ENDIF
        ENDIF
! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 2*ECQ
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_X) THEN

        IF (XC%ID(IXC)==ID_XC_PBE_X) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        D=2*D1



# 1233

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
# 1241


        ELSE
           S=0.0_q
        ENDIF
        CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
       &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

        D=2*D2
# 1252


        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
# 1260

        ELSE
           S=0.0_q
        ENDIF

        CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
       &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

        D=D1+D2

        IF (LLDA) THEN
! Add LDA contributions as well
           EXC  =(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*ALDAX
           EXCD1=(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*ALDAX
           EXCD2=(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*ALDAX
           EXCQ1=EXCQ1 *AGGAX
           EXCQ2=EXCQ2 *AGGAX
        ELSE
! Do not add LDA contribution
           EXC  =(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
           EXCD1=(EXCD1-EXDLDA1)*AGGAX
           EXCD2=(EXCD2-EXDLDA2)*AGGAX
           EXCQ1=EXCQ1 *AGGAX
           EXCQ2=EXCQ2 *AGGAX
        ENDIF
! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 0._q
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_C) THEN

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)

        IF (LLDA) THEN
! Add LDA contributions as well
           EXC  =EC  *AGGAC + ECLDA
           EXCD1=ECD1*AGGAC + ECD1LDA
           EXCD2=ECD2*AGGAC + ECD2LDA
           ECQ  =ECQ *AGGAC
        ELSE
! Do not add LDA contribution
           EXC  =EC*  AGGAC
           EXCD1=ECD1*AGGAC
           EXCD2=ECD2*AGGAC
           ECQ  =ECQ *AGGAC
        ENDIF
! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 0._q
        EXCQ2  = 0._q
        ECQ    = 2*ECQ
! jP
      ELSE IF ((XC%ID(IXC)==ID_XC_B3).OR.(XC%ID(IXC)==ID_XC_B5)) THEN

        CALL  B3LYPXCS(XC,D1,D2,DD1,DD2,DDA,EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA)

        D = D1 + D2
        EXC    = 2.0_q*EXC/D
        EXCD1  = 2.0_q*EXCD1
        EXCD2  = 2.0_q*EXCD2
        EXCQ1  = 2.0_q*EXCDD1
        EXCQ2  = 2.0_q*EXCDD2
        ECQ    = 2.0_q*EXCDDA

!       WRITE(91,'(7F12.7)') D1, D2, EXC, EXCD1, EXCD2, EXCDD1, EXCDD2
! jP
!aem AM05 spin formulation added
      ELSE IF (XC%ID(IXC)==ID_XC_AM05) THEN

!aem DD sometimes comes in negative.
!aem Since AM05 assumes that DD actually IS |grad rho|,
!aem and thus DD always should be positive,
!aem we need to use ABS(DD) as input to AM05.
!aem However, my routine needs to give the same symmetries out as
!aem PBE and PW91, therefor the sign-compensation of EXCDD.

        ider = 1
        CALL AM05(XC,D1,D2,ABS(DD1),ABS(DD2),ider, &
                  FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
                  FXCDD1,FXCDD2)

        D = D1 + D2
        IF (LLDA) THEN
          EXC = FXC/D
          EXCD1 = FXCD1
          EXCD2 = FXCD2
        ELSE
          EXC = (FXC - FXCLDA)/D
          EXCD1 = FXCD1 - VXCLDA1
          EXCD2 = FXCD2 - VXCLDA2
        ENDIF
        EXCQ1 = SIGN(1.0_q,DD1)*FXCDD1
        EXCQ2 = SIGN(1.0_q,DD2)*FXCDD2
        ECQ = 0._q

! Hartree -> Rydberg conversion
        EXC    = 2._q* EXC
        EXCD1  = 2._q* EXCD1
        EXCD2  = 2._q* EXCD2
        EXCQ1  = 2._q* EXCQ1
        EXCQ2  = 2._q* EXCQ2
        ECQ    = 2._q* ECQ
!aem AM05 spin formulation added
! jP: adding PBEsol
!     Note that this is a plain copy-paste - but it offers
!     the possibility to add a screened PBEsol version later,
!     for a possible PBEsol hybrid-functional  - closely related to HSE
!     Further note: necessary condition: fitting the PBEsol-exchange hole!!!
      ELSE IF (XC%ID(IXC)==ID_XC_PBESOL) THEN
           ukfactor=1.0_q


        D=2*D1



# 1397

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
# 1405


        ELSE
           S=0.0_q
        ENDIF
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(XC,D,S, &
          &  EXLDA1,EXDLDA1,EXLDA_SR1,EXLDA_LR1,EXDLDA_SR1,EXDLDA_LR1, &
          &  EXWPBE_SR1,EXWPBE_LR1,EXWPBED_SR1,EXWPBED_LR1,EXWPBEDD_SR1,EXWPBEDD_LR1)
        ENDIF
!        WRITE(*,'(10F14.7)') D,S,(EXC1-EXLDA1)/D,EXCD1-EXDLDA1,EXCQ1
!        EXLDA1=0; EXC1=0; EXDLDA1=0; EXCD1=0; EXCQ1=0

        D=2*D2
# 1424


        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
# 1432

        ELSE
           S=0.0_q
        ENDIF

        IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(XC,D,S, &
          &  EXLDA2,EXDLDA2,EXLDA_SR2,EXLDA_LR2,EXDLDA_SR2,EXDLDA_LR2, &
          &  EXWPBE_SR2,EXWPBE_LR2,EXWPBED_SR2,EXWPBED_LR2,EXWPBEDD_SR2,EXWPBEDD_LR2)
        ENDIF
!        EXLDA2=0; EXC2=0; EXDLDA2=0; EXCD2=0; EXCQ2=0
!        WRITE(*,'(10F14.7)') D,S,(EXC2-EXLDA2)/D,EXCD2-EXDLDA2,EXCQ2

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbesol(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)
!        ECLDA=0 ; ECD1LDA=0 ; ECD2LDA=0; EC=0 ; ECD1=0 ; ECD2=0 ; ECQ=0
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD1,ECD2,ECQ

        IF (LLDA) THEN
! Add LDA contributions as well
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*ALDAX+ECLDA
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*ALDAX+ECD1LDA
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*ALDAX+ECD2LDA
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX +(EXWPBE_SR1+EXWPBE_SR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX +EXWPBED_SR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX +EXWPBED_SR2+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ*AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC*AGGAC  +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ*AGGAC
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1+EXWPBE_SR2)/(2*D)*AGGAX +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_SR1*AGGAX +EXWPBED_LR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_SR2*AGGAX +EXWPBED_LR2+ECD2LDA
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ*AGGAC
              ENDIF
           ENDIF
        ELSE
! Do not add LDA contribution
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC*  AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX+(EXWPBED_SR1-EXDLDA_SR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX+(EXWPBED_SR2-EXDLDA_SR2)
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ *AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC=EC*AGGAC+ (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ *AGGAC
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_SR1-EXDLDA_SR1)*AGGAX+(EXWPBED_LR1-EXDLDA_LR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_SR2-EXDLDA_SR2)*AGGAX+(EXWPBED_LR2-EXDLDA_LR2)
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ *AGGAC
              ENDIF
           ENDIF
        ENDIF
! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 2*ECQ
! jH-new otherwise spin-polarized not possible
      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_RA).OR.(XC%ID(IXC)==ID_XC_ACFDT_PL)) THEN
        EXC    = 0._q
        EXCD1  = 0._q
        EXCD2  = 0._q
        EXCQ1  = 0._q
        EXCQ2  = 0._q
        ECQ    = 0._q
! jP: adding PBEsol

!BEEF
# 1603


      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_03).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_05).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_10).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_20)) THEN
! functionals for range-separated ACFDT (LDA - short range RPA)
! N.B. even if these functionals would support spin polarization,
! this would be a dummy stub since these functionals are not
! gradient corrected.
    CALL vtutor%error('GGA = 03 | 05 | 10 | 20  does not support spin polarization!')
# 1690

      ELSE
    CALL vtutor%bug('internal ERROR GGASPIN_DRIVER: Wrong functional, scheme not implemented!', "xc_driver.F", 1692)
      ENDIF

      EXC_TMP   = EXC_TMP   + XC%COEFF(IXC)*EXC
      EXCD1_TMP = EXCD1_TMP + XC%COEFF(IXC)*EXCD1
      EXCD2_TMP = EXCD2_TMP + XC%COEFF(IXC)*EXCD2
      EXCQ1_TMP = EXCQ1_TMP + XC%COEFF(IXC)*EXCQ1
      EXCQ2_TMP = EXCQ2_TMP + XC%COEFF(IXC)*EXCQ2
      ECQ_TMP   = ECQ_TMP   + XC%COEFF(IXC)*ECQ

      ENDIF

      ENDDO NXCLOOP4

      EXC   = EXC_TMP
      EXCD1 = EXCD1_TMP
      EXCD2 = EXCD2_TMP
      EXCQ1 = EXCQ1_TMP
      EXCQ2 = EXCQ2_TMP
      ECQ   = ECQ_TMP

      RETURN
      END SUBROUTINE GGASPIN_DRIVER


!************************ SUBROUTINE METAGGA_DRIVER ***********************************
!
!> This subroutine calls the subroutines that calculate xc quantities
!> (spin-polarized meta-GGA) for the energy, potential, forces, stress tensor, etc.
!>
!> Unlike in GGA_DRIVER and GGASPIN_DRIVER, the LDA contribution is included directly.
!
!**************************************************************************************

      SUBROUTINE METAGGA_DRIVER(&
     &   XC,RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,LAPLUP,LAPLDW,TAUUP,TAUDW, &
     &   EXC,dEXCdRHOup,dEXCdRHOdw,dEXCdABSNABup,dEXCdABSNABdw,dEXCdABSNAB, &
     &   dEXCdTAUup,dEXCdTAUdw,dEXCdLAPup,dEXCdLAPdw,CVMBJ)
!$ACC ROUTINE SEQ
      USE prec
      USE constant
      USE setexm_struct_def, only : xc_info
      USE ldalib, only : ex
      USE ggalib, only : corpbe,exchpbe
      USE mggalib
      USE xc_name
      USE xc_family_name
# 1742

      IMPLICIT NONE

      TYPE (xc_info) XC

      REAL(q) RHOUP,RHODW
      REAL(q) ABSNABUP,ABSNABDW,ABSNAB
      REAL(q) LAPLUP,LAPLDW
      REAL(q) TAUUP,TAUDW

      REAL(q) EXC
      REAL(q) dEXCdRHOup,dEXCdRHOdw
      REAL(q) dEXCdABSNABup,dEXCdABSNABdw,dEXCdABSNAB
      REAL(q) dEXCdTAUup,dEXCdTAUdw
      REAL(q) dEXCdLAPup,dEXCdLAPdw

      REAL(q) EXC_TMP
      REAL(q) dEXCdRHOup_TMP,dEXCdRHOdw_TMP
      REAL(q) dEXCdABSNABup_TMP,dEXCdABSNABdw_TMP,dEXCdABSNAB_TMP
      REAL(q) dEXCdTAUup_TMP,dEXCdTAUdw_TMP
      REAL(q) dEXCdLAPup_TMP,dEXCdLAPdw_TMP

      REAL(q) CVMBJ

! local variables
      REAL(q) Ex_TPSS,Ec_TPSS
      REAL(q) Ex_revTPSS,Ec_revTPSS
      REAL(q) Ex_SCAN,Ec_SCAN
      REAL(q) Ex_M06,Ec_M06
      REAL(q) Ex_MSX,Ec_MSX
      REAL(q) Ex_SRTM,Ec_RTM
      REAL(q) Ex_TASK,Ec_CC,Ex_LAK,Ec_LAK
      REAL(q) VXD1,VXD2,VXDD1,VXDD2,AMUXD1,AMUXD2
      REAL(q) VCD1,VCD2,VCDD1,VCDD2,AMUCD1,AMUCD2
      REAL(q) VCDD12

! VMBJ related variables
      REAL(q), PARAMETER :: THRD=1._q/3._q
      REAL(q) D,DTHRD,RS,ZETA
      REAL(q) ECLDA,ECD1LDA,ECD2LDA,EC,ECD1,ECD2,ECQ

! for SREGTM2L
      REAL(q) QQ,TAU1,TAU2,TAUTOT,TD1,TD2,TDD1,TDD2,TDL1,TDL2
      REAL(q) ECDT,ECL1,ECL2,ECQ1,ECQ2,EXCL1,EXCL2

! for PBETEST
      REAL(q) FK,G,S,SK,T
      REAL(q) EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2

      REAL(q) AMGGAX,AMGGAC

      INTEGER :: IXC

# 1798


!$ACC ROUTINE(VrevTPSSx) SEQ
!$ACC ROUTINE(VrevTPSSc) SEQ
!$ACC ROUTINE(VTPSSx) SEQ
!$ACC ROUTINE(VTPSSc) SEQ
!$ACC ROUTINE(VMBJ_PROYNOV) SEQ
!$ACC ROUTINE(CORPBE) SEQ
!$ACC ROUTINE(VM06x) SEQ
!$ACC ROUTINE(VM06c) SEQ
!$ACC ROUTINE(VMSXx) SEQ
!$ACC ROUTINE(VMSXc) SEQ
!$ACC ROUTINE(VSCANx) SEQ
!$ACC ROUTINE(VSCANc) SEQ
!$ACC ROUTINE(VR2SCANc) SEQ
!$ACC ROUTINE(VR2SCANx) SEQ
!$ACC ROUTINE(V1SRTMx) SEQ
!$ACC ROUTINE(V2SRTMx) SEQ
!$ACC ROUTINE(V3SRTMx) SEQ
!$ACC ROUTINE(VRTMc) SEQ
!$ACC ROUTINE(VTASKx) SEQ
!$ACC ROUTINE(VLAKx) SEQ
!$ACC ROUTINE(VLAKc) SEQ
!$ACC ROUTINE(scan_orb_free_xc) SEQ
!$ACC ROUTINE(EXCHPBE) SEQ
!$ACC ROUTINE(XV2SRTML) SEQ
!$ACC ROUTINE(PCREPMGGA) SEQ
!$ACC ROUTINE(CRTML) SEQ
# 1833


      AMGGAX=XC%AMGGAX
      AMGGAC=XC%AMGGAC

      EXC=0._q
      dEXCdRHOup=0._q
      dEXCdRHOdw=0._q
      dEXCdABSNABup=0._q
      dEXCdABSNABdw=0._q
      dEXCdABSNAB=0._q
      dEXCdTAUup=0._q
      dEXCdTAUdw=0._q
      dEXCdLAPup=0._q
      dEXCdLAPdw=0._q

      NXCLOOP: DO IXC=1,XC%NXC

      IF (XC%FAMILY(IXC)/=XC_FAM_MGGA) CYCLE

      EXC_TMP           = 0._q
      dEXCdRHOup_TMP    = 0._q
      dEXCdRHOdw_TMP    = 0._q
      dEXCdABSNABup_TMP = 0._q
      dEXCdABSNABdw_TMP = 0._q
      dEXCdABSNAB_TMP   = 0._q
      dEXCdTAUup_TMP    = 0._q
      dEXCdTAUdw_TMP    = 0._q
      dEXCdLAPup_TMP    = 0._q
      dEXCdLAPdw_TMP    = 0._q

      IF (XC%ID(IXC)==ID_XC_RTPSS) THEN
! revTPSS
! Exchange
         CALL VrevTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_revTPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

! Correlation
         CALL VrevTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_revTPSS+AMGGAC*Ec_revTPSS)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1+AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2+AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_RTPSS_X) THEN
! revTPSS
! Exchange
         CALL VrevTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_revTPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

         EXC_TMP=AMGGAX*Ex_revTPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_RTPSS_C) THEN
! revTPSS
! Correlation
         CALL VrevTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_revTPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_TPSS) THEN
! TPSS
! Exchange
         CALL VTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_TPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

! Correlation
         CALL VTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_TPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_TPSS+AMGGAC*Ec_TPSS)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1+AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2+AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_TPSS_X) THEN
! TPSS
! Exchange
         CALL VTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_TPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

         EXC_TMP=AMGGAX*Ex_TPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_TPSS_C) THEN
! TPSS
! Correlation
         CALL VTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_TPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_TPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MBJ).OR. &
              (XC%ID(IXC)==ID_XC_LMBJ)) THEN
! MBJ

         CALL VMBJ_PROYNOV(RHOUP,ABSNABUP,LAPLUP,2._q*TAUUP,CVMBJ,dEXCdRHOup_TMP)
         CALL VMBJ_PROYNOV(RHODW,ABSNABDW,LAPLDW,2._q*TAUDW,CVMBJ,dEXCdRHOdw_TMP)

! get the LDA correlation potential (using a call to CORPBE)
         D=RHOUP+RHODW
         DTHRD=exp(log(D)*THRD)
         RS=(0.75_q/PI)**THRD/DTHRD
         ZETA=(RHOUP-RHODW)/D
         ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

         CALL CORPBE(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,0._q,0._q,0._q,EC,ECD1,ECD2,ECQ,.FALSE.)

! add the LDA correlation potential
         dEXCdRHOup_TMP=dEXCdRHOup_TMP+ECD1LDA
         dEXCdRHOdw_TMP=dEXCdRHOdw_TMP+ECD2LDA

! LDA exchange-correlation energy
         EXC_TMP=ECLDA+0.25_q*EX((0.75_q/PI)**THRD/RHOUP**THRD,2,.TRUE.)+0.25_q*EX((0.75_q/PI)**THRD/RHODW**THRD,2,.TRUE.)
      ELSEIF (XC%ID(IXC)==ID_XC_M06L) THEN
! M06-L
! Exchange
         CALL VM06x(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_M06,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2,1)

! Correlation
         CALL VM06c(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_M06,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2,1)

         EXC_TMP=(AMGGAX*Ex_M06+AMGGAC*Ec_M06)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1+AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2+AMGGAC*VCDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_M06L_X) THEN
! M06-L
! Exchange
         CALL VM06x(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_M06,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2,1)

         EXC_TMP=AMGGAX*Ex_M06/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_M06L_C) THEN
! M06-L
! Correlation
         CALL VM06c(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_M06,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2,1)

         EXC_TMP=AMGGAC*Ec_M06/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAC*VCDD2
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MS0).OR. &
            & (XC%ID(IXC)==ID_XC_MS1).OR. &
            & (XC%ID(IXC)==ID_XC_MS2).OR. &
            & (XC%ID(IXC)==ID_XC_MSPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_MSRPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_MSB86BL).OR. &
            & (XC%ID(IXC)==ID_XC_RMSRPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_RMSPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_RMSB86BL)) THEN
! Exchange
         CALL VMSXx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_MSX,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%PARAM(IXC,3),XC%PARAM(IXC,4),XC%ID(IXC))

! Correlation
         CALL VMSXc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_MSX,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_MSX+AMGGAC*Ec_MSX)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MS0_X).OR. &
            & (XC%ID(IXC)==ID_XC_MS1_X).OR. &
            & (XC%ID(IXC)==ID_XC_MS2_X)) THEN
! Exchange
         CALL VMSXx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_MSX,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%PARAM(IXC,3),XC%PARAM(IXC,4),XC%ID(IXC))

         EXC_TMP=AMGGAX*Ex_MSX/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MS0_C).OR. &
            & (XC%ID(IXC)==ID_XC_MS1_C).OR. &
            & (XC%ID(IXC)==ID_XC_MS2_C)) THEN
! Correlation
         CALL VMSXc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_MSX,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_MSX/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCAN).OR. &
              (XC%ID(IXC)==ID_XC_RSCAN)) THEN
! SCAN
! Exchange
         CALL VSCANx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_SCAN,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,3),XC%ID(IXC))

! Correlation
         CALL VSCANc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_SCAN,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2, &
       &   XC%PARAM(IXC,2),XC%PARAM(IXC,3),XC%ID(IXC))

! Sum everything
         EXC_TMP=(AMGGAX*Ex_SCAN+AMGGAC*Ec_SCAN)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCAN_X).OR. &
              (XC%ID(IXC)==ID_XC_RSCAN_X)) THEN
! Exchange
         CALL VSCANx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_SCAN,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%ID(IXC))

         EXC_TMP=AMGGAX*Ex_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCAN_C).OR. &
              (XC%ID(IXC)==ID_XC_RSCAN_C)) THEN
! Correlation
         CALL VSCANc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_SCAN,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2,XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%ID(IXC))

         EXC_TMP=AMGGAC*Ec_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_R2SCAN) THEN
! r^2SCAN
! Exchange
         Ex_SCAN = 0.0_q
         CALL VR2SCANx(RHOUP,ABSNABUP,TAUUP,EX_SCAN,VXD1,VXDD1,AMUXD1)
         CALL VR2SCANx(RHODW,ABSNABDW,TAUDW,EX_SCAN,VXD2,VXDD2,AMUXD2)

! Correlation
         CALL VR2SCANc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_SCAN,VCD1,VCD2,&
      &                VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_SCAN+AMGGAC*Ec_SCAN)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_R2SCAN_X) THEN
! r^2SCAN
! Exchange
         Ex_SCAN = 0.0_q
         CALL VR2SCANx(RHOUP,ABSNABUP,TAUUP,EX_SCAN,VXD1,VXDD1,AMUXD1)
         CALL VR2SCANx(RHODW,ABSNABDW,TAUDW,EX_SCAN,VXD2,VXDD2,AMUXD2)

         EXC_TMP=AMGGAX*Ex_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_R2SCAN_C) THEN
! r^2SCAN
! Correlation
         CALL VR2SCANc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_SCAN,VCD1,VCD2,&
      &                VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM1) THEN
! Exchange
         Ex_SRTM = 0.0_q
         CALL V1SRTMx(RHOUP,ABSNABUP,TAUUP,Ex_SRTM,VXD1,VXDD1,AMUXD1)
         CALL V1SRTMx(RHODW,ABSNABDW,TAUDW,Ex_SRTM,VXD2,VXDD2,AMUXD2)
!Correlation rrTM
         CALL VRTMc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_RTM,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)
!
!        ! Sum everything
         EXC_TMP=(AMGGAX*Ex_SRTM+AMGGAC*Ec_RTM)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM2) THEN
! Exchange
         Ex_SRTM = 0.0_q
         CALL V2SRTMx(RHOUP,ABSNABUP,TAUUP,Ex_SRTM,VXD1,VXDD1,AMUXD1)
         CALL V2SRTMx(RHODW,ABSNABDW,TAUDW,Ex_SRTM,VXD2,VXDD2,AMUXD2)
!Correlation rrTM
         CALL VRTMc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_RTM,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)
!
!        ! Sum everything
         EXC_TMP=(AMGGAX*Ex_SRTM+AMGGAC*Ec_RTM)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM3) THEN
! Exchange
         Ex_SRTM = 0.0_q
         CALL V3SRTMx(RHOUP,ABSNABUP,TAUUP,Ex_SRTM,VXD1,VXDD1,AMUXD1)
         CALL V3SRTMx(RHODW,ABSNABDW,TAUDW,Ex_SRTM,VXD2,VXDD2,AMUXD2)
!Correlation rrTM
         CALL VRTMc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_RTM,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)
!
!        ! Sum everything
         EXC_TMP=(AMGGAX*Ex_SRTM+AMGGAC*Ec_RTM)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_TASK_X) THEN
! Exchange
         Ex_TASK = 0.0_q
         CALL VTASKx(RHOUP,ABSNABUP,TAUUP,EX_TASK,VXD1,VXDD1,AMUXD1)
         CALL VTASKx(RHODW,ABSNABDW,TAUDW,EX_TASK,VXD2,VXDD2,AMUXD2)

         EXC_TMP=AMGGAX*Ex_TASK/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_CC_C) THEN
!Correlation
         CALL VCCc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_CC,VCD1,VCD2,&
       &   VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_CC/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_LAK) THEN
! Exchange
         Ex_LAK = 0.0_q
         CALL VLAKx(RHOUP,ABSNABUP,TAUUP,EX_LAK,VXD1,VXDD1,AMUXD1)
         CALL VLAKx(RHODW,ABSNABDW,TAUDW,EX_LAK,VXD2,VXDD2,AMUXD2)
!Correlation
         CALL VLAKc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_LAK,VCD1,VCD2,&
       &   VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_LAK+AMGGAC*Ec_LAK)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_LAK_X) THEN
! Exchange
         Ex_LAK = 0.0_q
         CALL VLAKx(RHOUP,ABSNABUP,TAUUP,EX_LAK,VXD1,VXDD1,AMUXD1)
         CALL VLAKx(RHODW,ABSNABDW,TAUDW,EX_LAK,VXD2,VXDD2,AMUXD2)

         EXC_TMP=AMGGAX*Ex_LAK/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_LAK_C) THEN
!Correlation
         CALL VLAKc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_LAK,VCD1,VCD2,&
       &   VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_LAK/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCANL).OR. &
              (XC%ID(IXC)==ID_XC_RSCANL).OR. &
              (XC%ID(IXC)==ID_XC_R2SCANL).OR. &
              (XC%ID(IXC)==ID_XC_OFR2)) THEN
        call scan_orb_free_xc(rhoup,rhodw,absnabup,absnabdw,absnab, &
        &  LAPLUP,LAPLDW,EXC_TMP,dEXCdRHOup_TMP,dEXCdRHOdw_TMP,dEXCdABSNABup_TMP,dEXCdABSNABdw_TMP, &
        &  dEXCdABSNAB_TMP,dEXCdLAPup_TMP,dEXCdLAPdw_TMP,XC%ID(IXC),XC%AMGGAX,XC%AMGGAC)
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM2L) THEN

         D=2*RHOUP
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABUP/(D*FK)
            QQ=LAPLUP/(2.0_q*FK*FK*D)
         ELSE
            S=0.0_q
            QQ=0.0_q
         ENDIF

         CALL XV2SRTML(D,DTHRD,S,QQ,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1,EXCL1)
         CALL PCREPMGGA(D,DTHRD,S,QQ,TAU1,TD1,TDD1,TDL1)
         TAU1=0.5_q*TAU1

         D=2*RHODW
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABDW/(D*FK)
            QQ=LAPLDW/(2.0_q*FK*FK*D)
         ELSE
            S=0.0_q
            QQ=0.0_q
         ENDIF

         CALL XV2SRTML(D,DTHRD,S,QQ,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2,EXCL2)
         CALL PCREPMGGA(D,DTHRD,S,QQ,TAU2,TD2,TDD2,TDL2)
         TAU2=0.5_q*TAU2

         D=RHOUP+RHODW
         DTHRD=exp(log(D)*THRD)
         RS=(0.75_q/PI)**THRD/DTHRD
         ZETA=(RHOUP-RHODW)/D
         ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
         FK=(3._q*PI*PI)**THRD*DTHRD
         G=(exp((2*THRD)*log(1._q+ZETA)) &
     &      +exp((2*THRD)*log(1._q-ZETA)))/2._q
         IF(D>1.E-10_q)THEN
            S=ABSNAB/(D*2._q*FK)
         ELSE
            S=0.0_q
         ENDIF
         TAUTOT=TAU1+TAU2

         CALL CRTML(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,FK,S,EC,ECD1,ECD2,ECQ,TAUTOT,ECDT)

         ECD1 = ECD1 + ECDT*TD1
         ECD2 = ECD2 + ECDT*TD2
         ECQ1 = ECDT*TDD1
         ECQ2 = ECDT*TDD2
         ECL1 = ECDT*TDL1
         ECL2 = ECDT*TDL2

         EXC_TMP = AMGGAX*(EXC1 + EXC2)/(2._q*D) + AMGGAC*(ECLDA+EC)
         dEXCdRHOup_TMP = AMGGAX*EXCD1 + AMGGAC*(ECD1LDA+ECD1)
         dEXCdRHOdw_TMP = AMGGAX*EXCD2 + AMGGAC*(ECD2LDA+ECD2)
         dEXCdABSNABup_TMP = AMGGAX*EXCQ1 + AMGGAC*ECQ1
         dEXCdABSNABdw_TMP = AMGGAX*EXCQ2 + AMGGAC*ECQ2
         dEXCdABSNAB_TMP = AMGGAC*ECQ
         dEXCdLAPup_TMP = AMGGAX*EXCL1 + AMGGAC*ECL1
         dEXCdLAPdw_TMP = AMGGAX*EXCL2 + AMGGAC*ECL2

# 2421

      ELSEIF (XC%ID(IXC)==ID_XC_PBETEST) THEN
! PBE, for testing mainly

         D=2*RHOUP
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABUP/(D*FK)
         ELSE
            S=0.0_q
         ENDIF
         CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1,1.0_q,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
         D=2*RHODW
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABDW/(D*FK)
         ELSE
            S=0.0_q
         ENDIF
         CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2,1.0_q,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

         D=RHOUP+RHODW
         DTHRD=exp(log(D)*THRD)
         RS=(0.75_q/PI)**THRD/DTHRD
         ZETA=(RHOUP-RHODW)/D
         ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
         FK=(3._q*PI*PI)**THRD*DTHRD
         SK = SQRT(4.0_q*FK/PI)
         G = (exp((2*THRD)*log(1._q+ZETA)) &
                 +exp((2*THRD)*log(1._q-ZETA)))/2._q
         T = ABSNAB/(D*2._q*SK*G)
         CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK,T,EC,ECD1,ECD2,ECQ,.TRUE.)

         EXC_TMP           = AMGGAX*(EXC1+EXC2)/(2._q*D)+AMGGAC*(ECLDA+EC)
         dEXCdRHOup_TMP    = AMGGAX*EXCD1+AMGGAC*(ECD1LDA+ECD1)
         dEXCdRHOdw_TMP    = AMGGAX*EXCD2+AMGGAC*(ECD2LDA+ECD2)
         dEXCdABSNABup_TMP = AMGGAX*EXCQ1
         dEXCdABSNABdw_TMP = AMGGAX*EXCQ2
         dEXCdABSNAB_TMP   = AMGGAC*ECQ
      ENDIF

! factor 2 for Hartree to Rydberg
      EXC           = EXC           + 2._q*XC%COEFF(IXC)*EXC_TMP
      dEXCdRHOup    = dEXCdRHOup    + 2._q*XC%COEFF(IXC)*dEXCdRHOup_TMP
      dEXCdRHOdw    = dEXCdRHOdw    + 2._q*XC%COEFF(IXC)*dEXCdRHOdw_TMP
      dEXCdABSNABup = dEXCdABSNABup + 2._q*XC%COEFF(IXC)*dEXCdABSNABup_TMP
      dEXCdABSNABdw = dEXCdABSNABdw + 2._q*XC%COEFF(IXC)*dEXCdABSNABdw_TMP
      dEXCdABSNAB   = dEXCdABSNAB   + 2._q*XC%COEFF(IXC)*dEXCdABSNAB_TMP
      dEXCdTAUup    = dEXCdTAUup    + 2._q*XC%COEFF(IXC)*dEXCdTAUup_TMP
      dEXCdTAUdw    = dEXCdTAUdw    + 2._q*XC%COEFF(IXC)*dEXCdTAUdw_TMP
      dEXCdLAPup    = dEXCdLAPup    + 2._q*XC%COEFF(IXC)*dEXCdLAPup_TMP
      dEXCdLAPdw    = dEXCdLAPdw    + 2._q*XC%COEFF(IXC)*dEXCdLAPdw_TMP

      ENDDO NXCLOOP

      RETURN
      END SUBROUTINE METAGGA_DRIVER

    END MODULE xc_driver
