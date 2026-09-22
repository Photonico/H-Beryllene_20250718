# 1 "fock_frc.F"
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


# 2 "fock_frc.F" 2 
!***********************************************************************
!
!***********************************************************************

      MODULE fock_frc

!! USE moffload

      USE prec
      USE fock

      IMPLICIT NONE

      PUBLIC :: FOCK_FORCE

      PRIVATE

# 21

      CONTAINS

!************************ SUBROUTINE FOCK_FORCE ************************
!
!> Calculate the contribution to the Hellmann-Feynman forces and
!> the stress tensors (if LSIF=.TRUE.) due to Fock exchange
!>
!> @details @ref openmp :
!> the loop over the occupied orbitals owned by the particular
!> 1-rank (label: mband) is distributed over all available OpenMP
!> threads.
!
!***********************************************************************

      SUBROUTINE FOCK_FORCE( W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,FORHF,SIFHF,LSIF,IU0,IU6 )

!! USE moffload_fft

      USE wave_high
      USE nonl_high
      USE lattice
      USE pseudo
      USE full_kpoints
      USE sym_prec

!$ACC ROUTINE(NI_GLOBAL) SEQ

      TYPE (wavespin)     :: W
      TYPE (latt)         :: LATT_CUR
      TYPE (nonlr_struct) :: NONLR_S
      TYPE (nonl_struct)  :: NONL_S
      TYPE (potcar)       :: P(NONLR_S%NTYP)

      INTEGER :: LMDIM

      REAL(q) :: FORHF(3,NONLR_S%NIONS)
      REAL(q) :: SIFHF(3,3)

      LOGICAL :: LSIF

      INTEGER :: IU0,IU6

! local variables
      TYPE (wavespin) :: WHF
      TYPE (wavefun1) :: WQ
      TYPE (wavefun1), ALLOCATABLE :: W1(:)
      TYPE (wavefun1), ALLOCATABLE :: WIN(:)
      TYPE (wavedes1), TARGET :: WDESK,WDESQ,WDESQ_IRZ

      TYPE (nonlr_struct), ALLOCATABLE :: FAST_AUG(:)

      COMPLEX(q), ALLOCATABLE :: CPROJKXYZ(:,:,:),CPROJXYZ(:,:,:)
      COMPLEX(q), ALLOCATABLE :: CRHOLM(:,:)
      COMPLEX(q), ALLOCATABLE :: GWORK(:,:)

      COMPLEX(q), ALLOCATABLE :: CDIJ(:,:,:,:,:)

      COMPLEX(q), ALLOCATABLE, TARGET:: CDLM0(:,:),CDLM(:,:)

      REAL(q), ALLOCATABLE :: POTFAK(:,:)
      REAL(q), ALLOCATABLE :: FSG(:)

      REAL(q) :: WEIGHT_Q,WEIGHT
      REAL(q) :: RTMP

      REAL(q), ALLOCATABLE :: ENL(:,:)
      REAL(q) :: SIF(0:6)

      INTEGER :: ISP,NK,NB,NQ,MQ,ISP_IRZ,N,ISPINOR,NI,NT,LMMAXC,NPRO,NIP
      INTEGER :: NB_TOT,NB_TOTK,NBLK,NBLOCK,NBLOCK_ACT,NPOS,NGLB
      INTEGER :: NDIR,IDIR,I,J
      LOGICAL :: LSKIP,LSHIFT
      REAL(q) :: FD

      TYPE (rotation_handle), POINTER :: ROT_HANDLE

# 101

      

! early exit if possible
      IF ((.NOT.W%WDES%LOVERL.AND..NOT.LSIF).OR.MODEL_GW>0.OR.AEXX==0) THEN
         FORHF=0; SIFHF=0
         
         RETURN
      ENDIF

# 124


      CALL CHECK_FULL_KPOINTS
      NULLIFY(ROT_HANDLE)

      WHF=W
      WHF%WDES=>WDES_FOCK
# 133


      NDIR=3; IF (LSIF) NDIR=9

! determine the number of (partially) occupied states
      NB_TOT=0
      DO ISP=1,WHF%WDES%ISPIN
         DO NK=1,WHF%WDES%NKPTS
            DO NB=1,WHF%WDES%NB_TOT
               IF (ABS(W%FERTOT(NB,NK,ISP))>1E-8_q) NB_TOT=MAX(NB_TOT,NB)
            ENDDO
         ENDDO
      ENDDO

      NBLK=NBLOCK_FOCK
      NBLK=MIN(NBLK,NB_TOT)

      CALL WRK_ALLOCATE
!$ACC WAIT IF(OFFLOAD_ON)

!$ACC KERNELS PRESENT(FORHF,SIF) 
      FORHF=0; SIF=0
!$ACC END KERNELS

      spn: DO ISP=1,WHF%WDES%ISPIN
      kpt: DO NK=1,WHF%WDES%NKPTS

! set all NK dependent stuff
         CALL PREAMBLE_K ; IF (LSKIP) CYCLE kpt

! run over all bands at NK in chunks of size NBLOCK_ACT
         band: DO NPOS=1,NB_TOTK,NBLOCK
            NBLOCK_ACT=MIN(NB_TOTK-NPOS+1,NBLOCK)

# 169

            CALL FFT_AND_GATHER
            qpt: DO NQ=1,KPOINTS_FULL%NKPTS

! set all NQ dependent stuff
               CALL PREAMBLE_Q ; IF (LSKIP) CYCLE qpt


!$            CALL SET_ROT_HANDLE(P,LATT_CUR,WDESQ,ROT_HANDLE,KPOINTS_FULL%ISYMOP(:,:,WDESQ%NK))

!!!$OMP PARALLEL PRIVATE(WQ) REDUCTION(+:SIF) REDUCTION(-:FORHF)
!$            CALL NEWWAV(WQ,WDESQ,.TRUE.)
!!!$OMP DO SCHEDULE(STATIC) FIRSTPRIVATE(GWORK,CRHOLM,AUG_DES,POTFAK,W1,CDLM0,CDLM,CDIJ) &
!!!$OMP PRIVATE(MQ,LSHIFT,N,NGLB,WEIGHT,IDIR,ENL,ISPINOR,NI,NT,LMMAXC,NPRO,NIP,RTMP)

               mband: DO MQ=1,WHF%WDES%NBANDS
                  IF (ABS(WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))<=1E-10_q .OR. &
                      (MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW<NBANDSGWLOW_FOCK) CYCLE mband

                  IF (NQ<=WHF%WDES%NKPTS) THEN

                     CALL W1_COPY(ELEMENT(WHF, WDESQ, MQ, ISP), WQ)
                     CALL FFTWAV_W1(WQ)
                  ELSE

!
! symmetry must be considered if the wavefunctions for this
! k-point NQ (containing all k-points in the entire BZ)
! are not stored in W
!
                     LSHIFT=.FALSE.
                     IF ((ABS(KPOINTS_FULL%TRANS(1,NQ))>TINY) .OR. &
                         (ABS(KPOINTS_FULL%TRANS(2,NQ))>TINY) .OR. &
                         (ABS(KPOINTS_FULL%TRANS(3,NQ))>TINY)) LSHIFT=.TRUE.

                     CALL W1_ROTATE_AND_FFT(WQ, ELEMENT(WHF, WDESQ_IRZ, MQ, ISP_IRZ), ROT_HANDLE, P, LATT_CUR, LSHIFT)

                  ENDIF

! calculate charge phi_q nq(r) phi_k nk(r)
                  CALL FOCK_CHARGE_MU(WIN(1:NBLOCK_ACT),WQ,GWORK,CRHOLM)

# 248

                  nband: DO N=1,NBLOCK_ACT

                     NGLB=NPOS+N-1
! fft to reciprocal space
                     CALL FFT3D(GWORK(1,N),GRIDHF,-1)

                     WEIGHT=WHF%WDES%RSPIN*WHF%WDES%WTKPT(NK)*WHF%FERTOT(NGLB,NK,ISP)*  &
                          WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)*WEIGHT_Q

                     IF (LSIF) CALL APPLY_GFAC_DER(GRIDHF,GWORK(1,N),POTFAK(1,0),SIF(0),WEIGHT)

                     IF (WHF%WDES%LOVERL) THEN
! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
! and multiply in the k-point and Fermi-weight (gK 19.06.2020)
! this is cleaner and more transparent
                        CALL APPLY_GFAC_WEIGHT(GRIDHF,GWORK(1,N),POTFAK(1,0),WEIGHT)
! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
                        CALL FFT3D(GWORK(1,N),GRIDHF,1)
                     ENDIF
                  ENDDO nband

                  IF (WHF%WDES%LOVERL) THEN
! multiplicative factor used in RPROMU_HF_*
                     AUG_DES%RINPL=1.0_q/GRIDHF%NPLWV
!$ACC UPDATE DEVICE(AUG_DES%RINPL) 
# 278

! workspace used by RPROMU_HF
                        DO N=1,NBLOCK_ACT
                           W1(N)%CPROJ => CDLM0(:,N)
                        ENDDO
                        CALL RPROMU_HF(FAST_AUG_FOCK,AUG_DES,W1,NBLOCK_ACT,GWORK(1,1),SIZE(GWORK,1))
# 286

                     IF (WHF%WDES%NRSPINORS==2) THEN
                        DO N=1,NBLOCK_ACT
                           CALL ZCOPY(AUG_DES%NPRO,CDLM0(1,N),1,CDLM0(AUG_DES%NPRO+1,N),1)
                        ENDDO
                     ENDIF

                     DO IDIR=1,NDIR
# 298

! workspace used by RPROMU_HF
                           DO N=1,NBLOCK_ACT
                              W1(N)%CPROJ => CDLM(:,N)
                           ENDDO
                           CALL RPROMU_HF(FAST_AUG(IDIR),AUG_DES,W1,NBLOCK_ACT,GWORK(1,1),SIZE(GWORK,1))
# 306

                        IF (WHF%WDES%NRSPINORS==2) THEN
                           DO N=1,NBLOCK_ACT
                              CALL ZCOPY(AUG_DES%NPRO,CDLM(1,N),1,CDLM(AUG_DES%NPRO+1,N),1)
                           ENDDO
                        ENDIF

                        CALL DLLMM_TRANS_ECCP_NL_FOCK

                        

                        IF (IDIR<=3) THEN
!$ACC PARALLEL LOOP GANG PRIVATE(RTMP,NIP) PRESENT(WHF,WHF%WDES,FORHF,ENL) 
                           DO NI=1,WHF%WDES%NIONS
                              NIP=NI_GLOBAL(NI,WHF%WDES%COMM_INB)
                              RTMP=0
!$ACC LOOP VECTOR REDUCTION(+:RTMP)
                              DO N=1,NBLOCK_ACT
                                 RTMP=RTMP-ENL(N,NI)
                              ENDDO
!$ACC ATOMIC UPDATE
                              FORHF(IDIR,NIP)=FORHF(IDIR,NIP)+RTMP
                           ENDDO
                        ELSE
!$ACC PARALLEL LOOP GANG PRIVATE(RTMP) PRESENT(WHF,WHF%WDES,SIF,ENL) 
                           DO NI=1,WHF%WDES%NIONS
                              RTMP=0
!$ACC LOOP VECTOR REDUCTION(+:RTMP)
                              DO N=1,NBLOCK_ACT
                                 RTMP=RTMP+ENL(N,NI)
                              ENDDO
!$ACC ATOMIC UPDATE
                              SIF(IDIR-3)=SIF(IDIR-3)+RTMP
                           ENDDO
                        ENDIF

                        

                     ENDDO
                  ENDIF

               ENDDO mband

!!!$OMP END DO
!$             CALL DELWAV(WQ,.TRUE.)
!!!$OMP END PARALLEL

            ENDDO qpt
         ENDDO band

      ENDDO kpt
      ENDDO spn

      CALL M_sum_d(WHF%WDES%COMM_KINTER,FORHF(1,1),NONLR_S%NIONS*3)
      CALL M_sum_d(WHF%WDES%COMM_KINTER,SIF,7)

      CALL M_sum_d(WDESK%COMM_KIN,FORHF(1,1),NONLR_S%NIONS*3)
      CALL M_sum_d(WDESK%COMM_KIN,SIF,7)

!$ACC UPDATE SELF(FORHF,SIF) 
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

      FORHF=-FORHF

      IDIR=0
      DO I=1,3
         DO J=1,I
            IDIR=IDIR+1
            SIFHF(I,J)=SIF(IDIR)
            SIFHF(J,I)=SIF(IDIR)
         ENDDO
      ENDDO

      CALL WRK_DEALLOCATE
# 382

      CALL WNULLIFY(WHF)

# 398


      

!***********************************************************************
!***********************************************************************
!
! Internal subroutines: begin
!
!***********************************************************************
!***********************************************************************
      CONTAINS

!************************ SUBROUTINE WRK_ALLOCATE **********************
!
!> Allocate workspace and setup the derivatives of the augmentation
!> charges w.r.t. the ionic positions and the lattice vectors
!> (if LSIF=.TRUE.).
!
!***********************************************************************

      SUBROUTINE WRK_ALLOCATE
      USE tutor, ONLY: vtutor, isAlert, FockForce
! local variables
      TYPE (latt) :: LATT_FIN1,LATT_FIN2

      REAL(q) :: DIS,DISPL1(3,NONLR_S%NIONS),DISPL2(3,NONLR_S%NIONS)
      INTEGER :: N,IDIR,I,J

      INTEGER :: ISTATUS,ISTT

      

      ALLOCATE(GWORK( GRIDHF%MPLWV,NBLK),STAT=ISTATUS)
!$ACC ENTER DATA CREATE(GWORK) 

      ALLOCATE(POTFAK(GRIDHF%MPLWV,0:NDIR-3),FSG(0:NDIR-3),STAT=ISTT)
      ISTATUS=ISTATUS+ISTT

!$ACC ENTER DATA CREATE(POTFAK,FORHF,SIF) 

! average electrostatic potential prefactor for k=k' and n=n'
      IF (LSIF) THEN
         CALL SET_FSG_DER(GRIDHF,LATT_CUR,FSG)
      ELSE
         FSG(0)=SET_FSG(GRIDHF,LATT_CUR)
      ENDIF

!$ACC ENTER DATA CREATE(WDESQ) 
      CALL SETWDES(WHF%WDES,WDESQ,0)

!$ACC ENTER DATA CREATE(WQ)  
      CALL NEWWAV(WQ,WDESQ,.TRUE.,ISTT)
      ISTATUS=ISTATUS+ISTT


!$ACC ENTER DATA CREATE(WDESK) 
      CALL SETWDES(WHF%WDES,WDESK,0)

      ALLOCATE(WIN(NBLK),W1(NBLK))
!$ACC ENTER DATA CREATE(WIN(:)) 
      DO N=1,NBLK
         CALL NEWWAV(WIN(N),WDESK,.TRUE.,ISTT)
         IF (ISTT/=0) EXIT
      ENDDO
      ISTATUS=ISTATUS+ISTT

      IF (WHF%WDES%LOVERL) THEN
         ALLOCATE(CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
                  CDLM0(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
                  CDLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
                  CPROJXYZ(WHF%WDES%NPROD,WHF%WDES%NBANDS,NDIR), &
                  CPROJKXYZ(WHF%WDES%NPROD,NBLK,NDIR), &
                  ENL(NBLK,NONLR_S%NIONS),STAT=ISTT)
         ISTATUS=ISTATUS+ISTT

         ALLOCATE(CDIJ(LMDIM,LMDIM,WHF%WDES%NIONS,WHF%WDES%NRSPINORS,2),STAT=ISTT)
         ISTATUS=ISTATUS+ISTT

         ALLOCATE(FAST_AUG(NDIR),STAT=ISTT)
!!$ACC ENTER DATA CREATE(FAST_AUG(:)) 
         ISTATUS=ISTATUS+ISTT

! setup first derivative of augmentation charges with respect to ionic positions
!!    

         DIS=fd_displacement

         DO IDIR=1,3
            CALL COPY_FASTAUG(FAST_AUG_FOCK,FAST_AUG(IDIR))

            DISPL1=0
            DISPL1(IDIR,:)=-DIS

            DISPL2=0
            DISPL2(IDIR,:)= DIS

            CALL RSPHER_ALL(GRIDHF,FAST_AUG(IDIR),LATT_CUR,LATT_CUR,LATT_CUR,DISPL1,DISPL2,1)

! phase factor identical to FAST_AUG_FOCK
            IF (ASSOCIATED(FAST_AUG(IDIR)%CRREXP)) DEALLOCATE(FAST_AUG(IDIR)%CRREXP)
            FAST_AUG(IDIR)%CRREXP=>FAST_AUG_FOCK%CRREXP

# 504

            FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*SQRT(LATT_CUR%OMEGA)*(1._q/(2._q*DIS))

         ENDDO

! setup first derivative of augmentation charges with respect to lattice vectors
         IF (LSIF) THEN
            IDIR=3
            DO I=1,3
               DO J=1,I
                  IDIR=IDIR+1
                  LATT_FIN1%A=LATT_CUR%A
                  LATT_FIN2%A=LATT_CUR%A
                  LATT_FIN1%A(I,:)=LATT_CUR%A(I,:)+DIS*LATT_CUR%A(J,:)
                  LATT_FIN2%A(I,:)=LATT_CUR%A(I,:)-DIS*LATT_CUR%A(J,:)

                  CALL LATTIC(LATT_FIN1)
                  CALL LATTIC(LATT_FIN2)

                  CALL COPY_FASTAUG(FAST_AUG_FOCK,FAST_AUG(IDIR))

                  DISPL1=0
                  CALL RSPHER_ALL(GRIDHF,FAST_AUG(IDIR),LATT_FIN2,LATT_FIN1,LATT_CUR,DISPL1,DISPL1,1,LOMEGA=.TRUE.)

! phase factor identical to FAST_AUG_FOCK
                  IF (ASSOCIATED(FAST_AUG(IDIR)%CRREXP)) DEALLOCATE(FAST_AUG(IDIR)%CRREXP)
                  FAST_AUG(IDIR)%CRREXP=>FAST_AUG_FOCK%CRREXP

# 535

                 FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*(1._q/(2._q*DIS))

               ENDDO
            ENDDO
         ENDIF

# 548

      ENDIF

      CALL M_sum_i(WHF%WDES%COMM,ISTATUS,1)
      IF (ISTATUS/=0) THEN
! return if there was not enough memory: the forces will be bogus
! but we might still end up with a usable WAVECAR file.
         CALL vtutor%write(isAlert, FockForce)
         
         RETURN
      ENDIF

# 564

      

      RETURN
      END SUBROUTINE WRK_ALLOCATE


!************************ SUBROUTINE WRK_DEALLOCATE ********************
!
!> Deallocate workspace
!
!***********************************************************************

      SUBROUTINE WRK_DEALLOCATE
! local variables
      INTEGER :: N,IDIR

      

# 587

      IF (WHF%WDES%LOVERL) THEN
         DO N=1,NBLK
            NULLIFY(W1(N)%CPROJ)
         ENDDO

!$ACC EXIT DATA DELETE(CRHOLM,CDLM0,CDLM,CPROJXYZ,CPROJKXYZ,ENL) 
         DEALLOCATE(CRHOLM,CDLM0,CDLM,CPROJXYZ,CPROJKXYZ,ENL)

         DEALLOCATE(CDIJ)

# 602

!!    

         DO IDIR=1,NDIR
            NULLIFY(FAST_AUG(IDIR)%CRREXP)
            CALL NONLR_DEALLOC(FAST_AUG(IDIR))
         ENDDO

!!    

!!!$ACC EXIT DATA DELETE(FAST_AUG) 
         DEALLOCATE(FAST_AUG)
      ENDIF
# 618


      CALL DELWAV(WQ,.TRUE.)
!$ACC EXIT DATA DELETE(WQ) 

      DO N=1,NBLK
         CALL DELWAV(WIN(N),.TRUE.)
      ENDDO
!$ACC EXIT DATA DELETE(WIN(:)) 
      DEALLOCATE(WIN,W1)

!$ACC EXIT DATA DELETE(GWORK) 
      DEALLOCATE(GWORK)

!$ACC EXIT DATA DELETE(POTFAK,FORHF,SIF) 
      DEALLOCATE(POTFAK)

      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)

      

      RETURN
      END SUBROUTINE WRK_DEALLOCATE


!************************ SUBROUTINE PREAMBLE_K ************************
!
!> Initialize all k-point dependent quantities
!
!***********************************************************************

      SUBROUTINE PREAMBLE_K


      IF (MOD(NK-1,WHF%WDES%COMM_KINTER%NCPU).NE.WHF%WDES%COMM_KINTER%NODE_ME-1) THEN
         LSKIP=.TRUE. ; RETURN
      ENDIF

      LSKIP=.FALSE.

      NB_TOTK=MIN(NB_TOT, W%WDES%NB_TOTK(NK,ISP))
      NBLOCK =MIN(NB_TOTK,NBLK)

      CALL SETWDES(WHF%WDES,WDESK,NK)

! first derivative of wavefunction character with respect to all ionic positions
      IF (WHF%WDES%LOVERL) THEN
         IF (NONLR_S%LREAL) THEN
            CALL PHASER(W%WDES%GRID,LATT_CUR,NONLR_S,NK,W%WDES)
            CALL RPROXYZ(W%WDES%GRID,NONLR_S,P,LATT_CUR,W,W%WDES,ISP,NK,CPROJXYZ)
            IF (LSIF) CALL RPROLAT_DER(W%WDES%GRID,NONLR_S,P,LATT_CUR,W,W%WDES,ISP,NK,CPROJXYZ(1,1,4))
         ELSE
            CALL PHASE(W%WDES,NONL_S,NK)
            CALL PROJXYZ(NONL_S,W%WDES,W,LATT_CUR,ISP,NK,CPROJXYZ(:,:,1:3))
            IF (LSIF) CALL PROJLAT_DER(P,NONL_S,W%WDES,W,LATT_CUR,ISP,NK,CPROJXYZ(1,1,4))
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE PREAMBLE_K

# 756

!!#define interleave_communication
# 864

      SUBROUTINE FFT_AND_GATHER
! local variables
      INTEGER :: NI,N,NB_LOCAL,IDIR,NP

      INTEGER, ALLOCATABLE :: requests(:)
      INTEGER :: nrequests

# 874


      

# 880


      DO N=NPOS,NPOS+NBLOCK_ACT-1
         IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
            NI=N-NPOS+1 ; NB_LOCAL=1+(N-1)/WHF%WDES%NB_PAR
            CALL W1_COPY(ELEMENT(WHF,WIN(NI)%WDES1,NB_LOCAL,ISP),WIN(NI))
            CALL FFTWAV_W1(WIN(NI))
         ENDIF
      ENDDO

! copy the derivatives of the wave function characters into CPROJK array
      IF (WHF%WDES%LOVERL) THEN
! copy the derivatives of the wave function characters into CPROJKXYZ
!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WHF,CPROJKXYZ,CPROJXYZ) PRIVATE(NI,NB_LOCAL) 
         DO IDIR=1,NDIR
            DO N=NPOS,NPOS+NBLOCK_ACT-1
               IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
                  NI=N-NPOS+1 ; NB_LOCAL=1+(N-1)/WHF%WDES%NB_PAR
!$ACC LOOP VECTOR
                  DO NP=1,WHF%WDES%NPROD
                     CPROJKXYZ(NP,NI,IDIR)=CPROJXYZ(NP,NB_LOCAL,IDIR)
                  ENDDO
               ENDIF
           ENDDO
        ENDDO
      ENDIF

! distribute WIN and CPROJKXYZ to all nodes
      redis: IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
         ALLOCATE(requests((NDIR+2)*NBLOCK_ACT))
         nrequests=0
# 914

         DO N=NPOS,NPOS+NBLOCK_ACT-1
            NI=N-NPOS+1

            nrequests=nrequests+1
            CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CR(1), &
           &     SIZE(WIN(NI)%CR),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))

            IF (WHF%WDES%LOVERL) THEN
               nrequests=nrequests+1

               CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CPROJ(1), &
              &     SIZE(WIN(NI)%CPROJ),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))

               DO IDIR=1,NDIR
                  nrequests=nrequests+1
                  CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,CPROJKXYZ(:,NI,IDIR), &
                 &     SIZE(CPROJKXYZ,1),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
               ENDDO
# 942

            ENDIF
         ENDDO
# 950

            CALL M_waitall(nrequests,requests(1))
# 954

         DEALLOCATE(requests)
      ENDIF redis

      

      RETURN
      END SUBROUTINE FFT_AND_GATHER


!************************ SUBROUTINE PREAMBLE_Q ************************
!
!> Initialize all q-point dependent quantities
!
!***********************************************************************

      SUBROUTINE PREAMBLE_Q

      IF (KPOINTS_FULL%WTKPT(NQ)==0.OR.(HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,NQ))) .OR. &
          (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-WHF%WDES%VKPT(:,NK)))) THEN
          LSKIP=.TRUE. ; RETURN
      ENDIF

      IF (ALLOCATED(WEIGHT_K_POINT_PAIR_SMALL_GROUP).AND.LSYMGRAD) THEN
         IF (WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)==0) THEN
            LSKIP=.TRUE. ; RETURN
         ENDIF
         WEIGHT_Q=WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)
      ELSE
         WEIGHT_Q=1
      ENDIF

      LSKIP=.FALSE.

      CALL SETWDES(WHF%WDES,WDESQ,NQ)

      CALL SETWDES(WHF%WDES,WDESQ_IRZ,KPOINTS_FULL%NEQUIV(NQ))

      ISP_IRZ=ISP; IF (KPOINTS_FULL%SPINFLIP(NQ)==1) ISP_IRZ=3-ISP

! set POTFAK for this q and k point
      IF (LSIF) THEN
         CALL SET_GFAC_DER(GRIDHF,LATT_CUR,NK,NQ,FSG,POTFAK)
      ELSE
         CALL SET_GFAC(GRIDHF,LATT_CUR,NK,NQ,FSG(0),POTFAK(1,0))
      ENDIF

      RETURN
      END SUBROUTINE PREAMBLE_Q


!************************ SUBROUTINE CALC_GFAC_DER_MU ******************
!
!***********************************************************************
# 1025


!************************ SUBROUTINE DLLMM_TRANS_ECCP_NL_FOCK **********
!
!***********************************************************************

      SUBROUTINE DLLMM_TRANS_ECCP_NL_FOCK
      USE mopenmp_struct_def, ONLY : omp_nthreads
! local variables
      COMPLEX(q) :: CDIJ0,CDIJ1,GTMP
      INTEGER :: ISPINOR,N,NI,NT,LMMAXC,NPRO,NAUG,L,LP,LM,NGLB

      

# 1092


      ENL=0
      DO N=1,NBLOCK_ACT

         CALL CALC_DLLMM_TRANS(WHF%WDES,AUG_DES,TRANS_MATRIX_FOCK,CDIJ(:,:,:,:,1),CDLM0(:,N))

         CALL CALC_DLLMM_TRANS(WHF%WDES,AUG_DES,TRANS_MATRIX_FOCK,CDIJ(:,:,:,:,2),CDLM(:,N))

         NGLB=N+NPOS-1

         DO ISPINOR=0,WHF%WDES%NRSPINORS-1
!$OMP PARALLEL DO PRIVATE(NI,NT,LMMAXC,NPRO)
            DO NI=1,WHF%WDES%NIONS
               NT=WHF%WDES%ITYP(NI)

               LMMAXC=WHF%WDES%LMMAX(NT)
               IF (LMMAXC==0) CYCLE

               NPRO=WHF%WDES%LMBASE(NI)+ISPINOR*WHF%WDES%NPRO/2

               CALL ECCP_NL_FOCK(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR,1),CDIJ(1,1,NI,1+ISPINOR,2), &
                    WQ%CPROJ(NPRO+1),CPROJKXYZ(NPRO+1,N,IDIR),WIN(N)%CPROJ(NPRO+1),ENL(N,NI),1.0_q)
            ENDDO
!$OMP END PARALLEL DO
         ENDDO

      ENDDO

      

      RETURN
      END SUBROUTINE DLLMM_TRANS_ECCP_NL_FOCK

!***********************************************************************
!***********************************************************************
!
! Internal subroutines: end
!
!***********************************************************************
!***********************************************************************

      END SUBROUTINE FOCK_FORCE

      END MODULE fock_frc
