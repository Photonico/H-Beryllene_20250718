# 1 "fock_dbl.F"
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


# 2 "fock_dbl.F" 2 
!***********************************************************************
!
!> This module implements the evaluation of the action of the Fock
!> operator on the orbitals.
!>
!> Communication is 1._q by means of a double buffer technique:
!> While non-blocking asynchronous communication is filling (1._q,0._q)
!> buffer, computation can proceed on data in the other buffer.
!> (to hide the communication behind computation).
!
!***********************************************************************

      MODULE fock_dbl

!! USE moffload

      USE prec
      USE fock
      USE wave

      IMPLICIT NONE

      PUBLIC :: FOCK_ALL_DBLBUF,FOCK_ACC_COPY

      PRIVATE

      TYPE buffer
         TYPE (wavefun1), ALLOCATABLE :: WIN(:)
         TYPE (wavefun1), ALLOCATABLE :: WXI(:)

         INTEGER :: NBSTRT_win,NBSTOP_win
         INTEGER :: NBSTRT_wxi,NBSTOP_wxi

         INTEGER, ALLOCATABLE :: requests_gth(:)
         INTEGER              :: nrequests_gth

         INTEGER, ALLOCATABLE :: requests_red(:)
         INTEGER              :: nrequests_red
# 42


         LOGICAL :: LDO_FWD,LDO_BCK
      END TYPE buffer

      TYPE (wavespin) :: WXI

      LOGICAL :: LALLOCATED=.FALSE.

      COMPLEX(q) :: CDCHF
      REAL(q)    :: EXHF_ACFDT

# 56


      CONTAINS

!************************ SUBROUTINE FOCK_ACC_COPY *********************
!
!***********************************************************************

      SUBROUTINE FOCK_ACC_COPY(NB1,NSTRIP,NK,ISP,CH,EX,EX_ACFDT)
      COMPLEX(q) :: CH(:,:) 
      INTEGER :: NB1,NSTRIP,NK,ISP
      COMPLEX(q) :: EX
      REAL(q), OPTIONAL :: EX_ACFDT
! local variables
      INTEGER :: NB,NS
      INTEGER :: NODE_ME

! sanity check
      IF (.NOT.LALLOCATED) THEN
         CALL vtutor%error("FOCK_ACC_COPY: ERROR: accelaration have not been calculated yet.")
      ENDIF

      NODE_ME=WXI%WDES%COMM%NODE_ME
# 81

      DO NB=NB1,NB1+NSTRIP-1
         IF (MOD(NB-1,WXI%WDES%NB_PAR)+1/=WXI%WDES%NB_LOW) CYCLE
         NS=(NB-1)/WXI%WDES%NB_PAR+1
         CH(:,NS)=WXI%CPTWFP(:,NS,NK,ISP)
      ENDDO

      EX=0; IF (NODE_ME==1) EX=CDCHF
      
      IF (PRESENT(EX_ACFDT)) THEN
         EX_ACFDT=0; IF (NODE_ME==1) EX_ACFDT=EXHF_ACFDT
      ENDIF
         

      RETURN
      END SUBROUTINE FOCK_ACC_COPY


!************************ SUBROUTINE FOCK_ALL_DBLBUF *******************
!
!> Computes the action of the Fock potential, generated by the orbitals
!> in W, on all orbitals in W (or alternatively WP).
!>
!> The result is stored in WXI\%CPTWFP (that points to XI if the latter is
!> provided).
!>
!> This routine uses a double buffer techniques and non-blocking
!> communication calls to maximize the overlap of communication and
!> computation.
!
!***********************************************************************

      SUBROUTINE FOCK_ALL_DBLBUF( &
         WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
         EX,EX_ACFDT,NBMAX,NKSTART,NKSTOP,LSYMGRAD,XI,WP)

!! USE moffload_fft

      USE pseudo
      USE lattice
      USE poscar
      USE wave_high
      USE nonl_high
      USE full_kpoints
      USE mopenmp_struct_def, ONLY : omp_nthreads,omp_nthreads_acc

      TYPE (wavedes)      :: WDES
      TYPE (wavespin)     :: W
      TYPE (latt)         :: LATT_CUR
      TYPE (nonlr_struct) :: NONLR_S
      TYPE (nonl_struct)  :: NONL_S
      TYPE (potcar)       :: P(:)
      INTEGER             :: LMDIM
      COMPLEX(q)             :: CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      REAL(q), OPTIONAL   :: EX
      REAL(q), OPTIONAL   :: EX_ACFDT
      INTEGER, OPTIONAL   :: NBMAX
      INTEGER, OPTIONAL   :: NKSTART,NKSTOP
      LOGICAL, OPTIONAL   :: LSYMGRAD
!> space to store the action
      TYPE (wavespin), OPTIONAL :: XI
!> compute the action onto an alternative set of orbitals
      TYPE (wavespin), OPTIONAL :: WP

! local variables
      TYPE (wavespin) :: WHF,WHP
      TYPE (wavedes1), TARGET      :: WDESK,WDESQ
      TYPE (wavefun1), TARGET      :: WQ
      TYPE (wavefun1), ALLOCATABLE :: W1(:)

      COMPLEX(q),      ALLOCATABLE :: GWORK(:,:)               ! fock pot in real sp
      COMPLEX(q),      ALLOCATABLE :: CRHOLM(:,:)              ! augmentation occupancy matrix

      COMPLEX(q),      ALLOCATABLE :: CDIJ(:,:,:,:)            ! D_lml'm'

      REAL(q),   ALLOCATABLE :: POTFAK(: )   ! 1/(G+dk)**2 (G)

      COMPLEX(q),ALLOCATABLE, TARGET :: CXI(:,:)         ! acc. in real space
      COMPLEX(q),      ALLOCATABLE, TARGET :: CDLM(:,:)        ! D_LM
      COMPLEX(q),      ALLOCATABLE, TARGET :: CKAPPA(:,:)      ! stores NL accelerations

      COMPLEX(q), ALLOCATABLE :: CWORK(:)

      TYPE (buffer) :: BUF1,BUF2

      INTEGER :: ISP,NK,NB
      INTEGER :: NB_TOT,NB_TOTK,NBLK
      INTEGER :: NBLK1,NBLK2,NBLK1P,NBLK2P
      INTEGER :: NSTR1,NSTP1,NSTR2,NSTP2
      INTEGER :: NKSTR,NKSTP

      REAL(q) :: FSG,FSG_AEXX ! singularity correction setting V(G=0) to a finite value

      TYPE( rotation_handle), POINTER :: ROT_HANDLE

      COMPLEX(q), ALLOCATABLE :: CTMP1(:),CTMP2(:)
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      

! early exit if possible
      IF (AEXX==0) THEN
         IF (PRESENT(EX)) EX=0
         IF (PRESENT(EX_ACFDT)) EX_ACFDT=0
         
         RETURN
      ENDIF

! test_
# 206

! test_

      CALL CHECK_FULL_KPOINTS
      NULLIFY(ROT_HANDLE)

      WHF=W
      WHF%WDES=>WDES_FOCK
# 216


      WHP=W; IF (PRESENT(WP)) WHP=WP
      WHP%WDES=>WDES_FOCK
# 222


! start at the first k-point
      NKSTR=1
! unless otherwise specified
      IF (PRESENT(NKSTART)) NKSTR=NKSTART

! number of k-points in the IBZ
      NKSTP=WDES%NKPTS
! unless otherwise specified
      IF (PRESENT(NKSTOP)) NKSTP=NKSTOP

      NB_TOT=WDES%NB_TOT
      IF (PRESENT(NBMAX)) NB_TOT=MIN(NB_TOT,NBMAX)

      NBLK=NBLOCK_FOCK
      NBLK=MIN(NBLK,NB_TOT)

      NBLK=MIN(NBLK, MINVAL(WDES%NB_TOTK(NKSTR:NKSTP,1:WDES%ISPIN)))

      CDCHF=0
      EXHF_ACFDT=0

      IF (PRESENT(XI)) THEN
! just point to the storage space passed down in the call
         WXI%CPTWFP=>XI%CPTWFP; WXI%WDES=>XI%WDES
         LALLOCATED=.TRUE.
# 253

      ELSE
! allocate space to store the action (if not already 1._q before)
         CALL ALLOCWXI(WDES)
      ENDIF

! allocate workspace
      CALL WRK_ALLOCATE

      spn: DO ISP=1,WDES%ISPIN
      kpt: DO NK=NKSTR,NKSTP

         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE kpt

         NB_TOTK=MIN(NB_TOT, WDES%NB_TOTK(NK,ISP))

! set all NK dependent stuff
         CALL PREAMBLE

         BUF1%LDO_FWD=.FALSE. ; BUF1%LDO_BCK=.FALSE.
         BUF2%LDO_FWD=.FALSE. ; BUF2%LDO_BCK=.FALSE.
# 276


         NB=1

! loop over bands on which the Fock potential acts
! (effectively runs over all NB_TOTK bands)
!$ACC WAIT IF(OFFLOAD_ON)
         bnd: DO

            NBLK1=MAX(MIN(NB_TOTK-(NB-1)*NBLK,NBLK),0)
            NSTR1=(NB-1)*NBLK+1
            NSTP1=NSTR1+NBLK1-1

            NBLK2=MAX(MIN(NB_TOTK- NB   *NBLK,NBLK),0)
            NSTR2= NB   *NBLK+1
            NSTP2=NSTR2+NBLK2-1

            IF (NBLK1>0) &
               CALL FFT_AND_GATHER(WHP,NSTR1,NSTP1,ISP,BUF1)

            IF (BUF1%LDO_BCK) &
               CALL REDUCE(BUF1)

            IF (BUF2%LDO_BCK) THEN
               CALL WAITRED(BUF2)
               CALL FOCK_BCK(BUF2)
            ENDIF

            IF (BUF2%LDO_FWD) THEN
               CALL WAITGTH(BUF2)
               CALL FOCK_FWD(BUF2)
            ENDIF

            IF (NBLK2>0) &
               CALL FFT_AND_GATHER(WHP,NSTR2,NSTP2,ISP,BUF2)

            IF (BUF2%LDO_BCK) &
               CALL REDUCE(BUF2)

            IF (BUF1%LDO_BCK) THEN
               CALL WAITRED(BUF1)
               CALL FOCK_BCK(BUF1)
            ENDIF

            IF (BUF1%LDO_FWD) THEN
               CALL WAITGTH(BUF1)
               CALL FOCK_FWD(BUF1)
            ENDIF

            IF ((.NOT.BUF1%LDO_FWD).AND.(.NOT.BUF1%LDO_BCK).AND. &
           &    (.NOT.BUF2%LDO_FWD).AND.(.NOT.BUF2%LDO_BCK)) EXIT bnd

            NB=NB+2

         ENDDO bnd
!$ACC WAIT IF(OFFLOAD_ON)
# 334


!         bnd: DO
!
!            NBLK1=MAX(MIN(NB_TOTK-(NB-1)*NBLK,NBLK),0)
!            NSTR1=(NB-1)*NBLK+1
!            NSTP1=NSTR1+NBLK1-1
!
!            NBLK2=MAX(MIN(NB_TOTK- NB   *NBLK,NBLK),0)
!            NSTR2= NB   *NBLK+1
!            NSTP2=NSTR2+NBLK2-1
!
!            IF (BUF2%LDO_BCK) &
!               CALL REDUCE(BUF2)
!
!            IF (BUF1%LDO_BCK) THEN
!               CALL WAITRED(BUF1)
!               CALL FOCK_BCK(BUF1)
!            ENDIF
!
!            IF (NBLK1>0) &
!               CALL FFT_AND_GATHER(WHF,NSTR1,NSTP1,ISP,BUF1)
!
!            IF (BUF2%LDO_BCK) THEN
!               CALL WAITRED(BUF2)
!               CALL FOCK_BCK(BUF2)
!            ENDIF
!
!            IF (NBLK2>0) &
!               CALL FFT_AND_GATHER(WHF,NSTR2,NSTP2,ISP,BUF2)
!
!            IF (BUF1%LDO_FWD) THEN
!               CALL WAITGTH(BUF1)
!               CALL FOCK_FWD(BUF1)
!               CALL REDUCE(BUF1)
!            ENDIF
!
!            IF (BUF2%LDO_FWD) THEN
!               CALL WAITGTH(BUF2)
!               CALL FOCK_FWD(BUF2)
!            ENDIF
!
!            IF ((.NOT.BUF1%LDO_FWD).AND.(.NOT.BUF1%LDO_BCK).AND. &
!           &    (.NOT.BUF2%LDO_FWD).AND.(.NOT.BUF2%LDO_BCK)) EXIT bnd
!
!            NB=NB+2
!         ENDDO bnd


! wait until all ranks working on the same k-point
! are finished before we move to the next
         CALL M_barrier(WDES%COMM_KIN)

      ENDDO kpt
      ENDDO spn

! deallocate workspace
      CALL WRK_DEALLOCATE

# 403

      CALL WNULLIFY(WHP)
      CALL WNULLIFY(WHF)

! test_
# 424

! test_

      CALL M_sum_z(WDES%COMM_KIN,    CDCHF, 1)
      CALL M_sum_z(WDES%COMM_KINTER, CDCHF, 1)
      IF (PRESENT(EX)) EX=CDCHF

      CALL M_sum_d(WDES%COMM_KIN,    EXHF_ACFDT, 1)
      CALL M_sum_d(WDES%COMM_KINTER, EXHF_ACFDT, 1)
      IF (PRESENT(EX_ACFDT)) EX_ACFDT=EXHF_ACFDT
 
      

!***********************************************************************
!***********************************************************************
!
! Internal subroutines: begin
!
!***********************************************************************
!***********************************************************************
      CONTAINS

!************************ SUBROUTINE WASSIGN ***************************
!
!***********************************************************************

      SUBROUTINE WASSIGN(W_LHS,W_RHS)
      TYPE (wavespin) :: W_LHS,W_RHS
! local variables
      INTEGER :: NB_LOW,NB_TOT,NB_PAR

      NB_LOW=W_RHS%WDES%NB_LOW
      NB_TOT=W_RHS%WDES%NB_TOT
      NB_PAR=W_RHS%WDES%NB_PAR

      W_LHS%WDES   => W_RHS%WDES
      W_LHS%CPTWFP => W_RHS%CPTWFP
      W_LHS%CPROJ  => W_RHS%CPROJ
      W_LHS%FERTOT => W_RHS%FERTOT
      W_LHS%CELTOT => W_RHS%CELTOT
      W_LHS%AUXTOT => W_RHS%AUXTOT

      W_LHS%FERWE  => W_RHS%FERTOT(NB_LOW:NB_TOT:NB_PAR,:,:)
      W_LHS%CELEN  => W_RHS%CELTOT(NB_LOW:NB_TOT:NB_PAR,:,:)
      W_LHS%AUX    => W_RHS%AUXTOT(NB_LOW:NB_TOT:NB_PAR,:,:)

      W_LHS%OVER_BAND=W_RHS%OVER_BAND

      RETURN
      END SUBROUTINE WASSIGN


!************************ SUBROUTINE WRK_ALLOCATE **********************
!
!***********************************************************************

      SUBROUTINE WRK_ALLOCATE
      USE ini
! local variables
      INTEGER :: ISTT,ISTATUS,N
      LOGICAL, SAVE :: LFIRST=.TRUE.
# 488

! Register (and write out) the memory demands of FOCK_ACC before the actual allocation.
! We include only the part that scales with NBLK. This is 1._q only the first time FOCK_ACC is called.
      IF (LFIRST) THEN
         CALL REGISTER_ALLOCATE( &
            16._q*      W%WDES%GRID%MPLWV*W%WDES%NRSPINORS+ &     ! CWORK
           ( 8._q*2* GRIDHF%MPLWV+ &                        ! GWORK
             8._q*2*2*AUG_DES%NPROD+ &                        ! CRHOLM+CDLM
            32._q*      GRID_FOCK%MPLWV+ &                        ! CXI
            16._q*2*WHF%WDES%NPROD+ &                         ! CKAPPA
           (32._q*      WHF%WDES%NRPLWV+ &                        ! WIN%CPTWFP
            32._q*      WHF%WDES%GRID%MPLWV*WHF%WDES%NRSPINORS+ & ! WIN%CR
            16._q*2*WHF%WDES%NPROD) &                         ! WIN%CPROJ
# 503

             )*NBLK,'fock_wrk')

         LFIRST=.FALSE.
      ENDIF

      ALLOCATE( &
           GWORK( GRIDHF%MPLWV,NBLK), &
           CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
           W1(NBLK), &

           CDIJ(LMDIM,LMDIM,WHF%WDES%NIONS,WHF%WDES%NRSPINORS), &

           CDLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
           CWORK(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS), &
           STAT=ISTATUS)

!$    IF (MODEL_GW/=2) THEN
!$       ALLOCATE(POTFAK(GRIDHF%MPLWV,omp_nthreads),STAT=ISTT)
!$    ELSE
         ALLOCATE(POTFAK(GRIDHF%MPLWV ),STAT=ISTT)
!$    ENDIF
      ISTATUS=ISTATUS+ISTT

      ALLOCATE(CTMP1(WHF%WDES%NPROD),CTMP2(WHF%WDES%NPROD),STAT=ISTT)
! CTMP1=0; CTMP2=0 ! Beaks bulk_InP_SOC_G0W0_nosym (gnu), may be uninitialized when passed to OVERL1.
      ISTATUS=ISTATUS+ISTT
!$ACC ENTER DATA CREATE(CTMP1,CTMP2,CWORK) IF(OFFLOAD_ON)

!$ACC ENTER DATA CREATE(WDESQ,WQ) IF(OFFLOAD_ON)
      CALL SETWDES(WHF%WDES,WDESQ,0)
      CALL NEWWAV(WQ, WDESQ, .TRUE., ISTT)
      ISTATUS=ISTATUS+ISTT
# 581

!$ACC ENTER DATA CREATE(WDESK) IF(OFFLOAD_ON)
      CALL SETWDES(WHF%WDES,WDESK,0)

      ALLOCATE(BUF1%WIN(NBLK))
!$ACC ENTER DATA CREATE(BUF1%WIN(:)) IF(OFFLOAD_ON)
      DO N=1,NBLK
         CALL NEWWAV(BUF1%WIN(N),WDESK,.TRUE.,ISTT)
         IF (ISTT/=0) EXIT
      ENDDO
      ISTATUS=ISTATUS+ISTT

      ALLOCATE(BUF2%WIN(NBLK))
!$ACC ENTER DATA CREATE(BUF2%WIN(:)) IF(OFFLOAD_ON)
      DO N=1,NBLK
         CALL NEWWAV(BUF2%WIN(N),WDESK,.TRUE.,ISTT)
         IF (ISTT/=0) EXIT
      ENDDO
      ISTATUS=ISTATUS+ISTT

      ALLOCATE( &
           CXI(GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS,2*NBLK), &
           CKAPPA(WHF%WDES%NPROD,2*NBLK), &
           STAT=ISTT)
      ISTATUS=ISTATUS+ISTT
!$ACC ENTER DATA CREATE(CXI,CKAPPA) IF(OFFLOAD_ON)

      ALLOCATE(BUF1%WXI(NBLK))
!$ACC ENTER DATA CREATE(BUF1) IF(OFFLOAD_ON)
!$ACC ENTER DATA CREATE(BUF1%WXI(:)) IF(OFFLOAD_ON)
      DO N=1,NBLK
         BUF1%WXI(N)%CR   =>   CXI(:,N)
         BUF1%WXI(N)%CPROJ=>CKAPPA(:,N)
!$ACC ENTER DATA CREATE(BUF1%WXI(N)%CR,BUF1%WXI(N)%CPROJ) IF(OFFLOAD_ON)
      ENDDO
      ALLOCATE(BUF2%WXI(NBLK))
!$ACC ENTER DATA CREATE(BUF2) IF(OFFLOAD_ON)
!$ACC ENTER DATA CREATE(BUF2%WXI(:)) IF(OFFLOAD_ON)
      DO N=1,NBLK
         BUF2%WXI(N)%CR   =>   CXI(:,N+NBLK)
         BUF2%WXI(N)%CPROJ=>CKAPPA(:,N+NBLK)
!$ACC ENTER DATA CREATE(BUF2%WXI(N)%CR,BUF2%WXI(N)%CPROJ) IF(OFFLOAD_ON)
      ENDDO
 
      ALLOCATE(BUF1%requests_gth(3*NBLK),BUF2%requests_gth(3*NBLK))
      ALLOCATE(BUF1%requests_red(2*NBLK),BUF2%requests_red(2*NBLK))
      BUF1%nrequests_gth=0; BUF2%nrequests_gth=0
      BUF1%nrequests_red=0; BUF2%nrequests_red=0

      CALL M_sum_i(WHF%WDES%COMM, ISTATUS, 1)

      IF (ISTATUS/=0) THEN
         CALL VTUTOR%ERROR("WRK_ALLOCATE: ERROR: could not allocate enough workspace. Try reducing NBLOCK_FOCK.")
      ENDIF

!$ACC ENTER DATA CREATE(CRHOLM,CDLM,POTFAK) 
# 640

      RETURN
      END SUBROUTINE WRK_ALLOCATE


!************************ SUBROUTINE WRK_DEALLOCATE ********************
!
!***********************************************************************

      SUBROUTINE WRK_DEALLOCATE
      USE ini
! local variables
      INTEGER :: N
      LOGICAL, SAVE :: LFIRST=.TRUE.

! Register deallocation of FOCK_ACC workspace (only on first call).
! When a shared memory segment is used for WIN this will not be deallocated.
      IF (LFIRST) THEN
         CALL DEREGISTER_ALLOCATE( &
            16._q*      W%WDES%GRID%MPLWV*W%WDES%NRSPINORS+ &     ! CWORK
           ( 8._q*2* GRIDHF%MPLWV+ &                        ! GWORK
             8._q*2*2*AUG_DES%NPROD+ &                        ! CRHOLM+CDLM
            32._q*      GRID_FOCK%MPLWV+ &                        ! CXI
            16._q*2*WHF%WDES%NPROD  &                         ! CKAPPA
             )*NBLK,'fock_wrk')

         LFIRST=.FALSE.
      ENDIF

# 685

      DEALLOCATE(CDIJ)

      DEALLOCATE(CXI,CKAPPA,CRHOLM,CDLM,W1,POTFAK,CWORK,CTMP1,CTMP2)

      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
      CALL DELWAV(WQ,.TRUE.)
!$ACC EXIT DATA DELETE(WQ) IF(OFFLOAD_ON)
# 704

      DO N=1,NBLK
         CALL DELWAV(BUF1%WIN(N) ,.TRUE.)
         CALL DELWAV(BUF2%WIN(N) ,.TRUE.)
      ENDDO
!$ACC EXIT DATA DELETE(BUF1%WIN(:)) IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(BUF2%WIN(:)) IF(OFFLOAD_ON)

      DEALLOCATE(BUF1%WIN,BUF2%WIN)

      DO N=1,NBLK
         NULLIFY(BUF1%WXI(N)%CR,BUF1%WXI(N)%CPROJ,BUF1%WXI(N)%WDES1)
         NULLIFY(BUF2%WXI(N)%CR,BUF2%WXI(N)%CPROJ,BUF2%WXI(N)%WDES1)
      ENDDO
      DEALLOCATE(BUF1%WXI,BUF2%WXI)

      DEALLOCATE(BUF1%requests_gth,BUF2%requests_gth)
      DEALLOCATE(BUF1%requests_red,BUF2%requests_red)

      RETURN
      END SUBROUTINE WRK_DEALLOCATE


!************************ SUBROUTINE ALLOCWXI **************************
!
!***********************************************************************

      SUBROUTINE ALLOCWXI(WDES)
      TYPE (wavedes) :: WDES
! local variables
      TYPE (wavefun) :: WTMP

      IF (LALLOCATED) RETURN

      CALL ALLOCW(WDES,WXI,WTMP,WTMP)
      LALLOCATED=.TRUE.

# 745


      RETURN
      END SUBROUTINE ALLOCWXI


!************************ SUBROUTINE PREAMBLE **************************
!
!***********************************************************************

      SUBROUTINE PREAMBLE
! local variables
      INTEGER :: N

! average electrostatic potential for k=k' and n=n'
      FSG=FSG_STORE(NK)
      CALL FOCKCORR_SELECT(MCALPHA,FSG,FSG_AEXX)
      CALL SETWDES(WHF%WDES,WDESK,NK)

      IF (NONLR_S%LREAL) THEN
         CALL PHASER(WDES%GRID,LATT_CUR,NONLR_S,NK,WDES)
      ELSE
         CALL PHASE(WDES,NONL_S,NK)
      ENDIF

      DO N=1,NBLK
         BUF1%WIN(N)%WDES1=>WDESK
         BUF1%WXI(N)%WDES1=>WDESK

         BUF2%WIN(N)%WDES1=>WDESK
         BUF2%WXI(N)%WDES1=>WDESK
      ENDDO

      RETURN
      END SUBROUTINE PREAMBLE


# 866

!************************ SUBROUTINE FFT_AND_GATHER ********************
!
!> Take the wavefunction NB1,...,NB2 to real space, store them in BUF
!> and broadcast them to all nodes in COMM_INTER.
!>
!> NB1 and NB2 are global band indices.
!>
!> We use non-blocking broadcasts and do NOT wait for them to end.
!
!***********************************************************************
      SUBROUTINE FFT_AND_GATHER(W,NB1,NB2,ISP,BUF)
# 880

      TYPE (wavespin) :: W
      INTEGER :: NB1
      INTEGER :: NB2
      INTEGER :: ISP
      TYPE (buffer) :: BUF

! local variables
      INTEGER :: NI, N, NB_LOCAL
# 891


      
# 896


      BUF%NBSTRT_win=NB1
      BUF%NBSTOP_win=NB2

      DO N=NB1,NB2
         NI=N-NB1+1
         IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
            NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
            CALL W1_COPY(ELEMENT(W,BUF%WIN(NI)%WDES1,NB_LOCAL,ISP),BUF%WIN(NI))
            CALL FFTWAV_W1(BUF%WIN(NI))
         ENDIF
      ENDDO

! distribute WIN to all nodes
      IF (W%WDES%COMM_INTER%NCPU>1) THEN

# 916

# 919

         BUF%nrequests_gth=0
         DO N=NB1,NB2
            NI=N-NB1+1
!!            BUF%nrequests_gth=BUF%nrequests_gth+1
!!            CALL M_ibcast_z_from(W%WDES%COMM_INTER,BUF%WIN(NI)%CPTWFP(1), &
!!           &     SIZE(BUF%WIN(NI)%CPTWFP),MOD(N-1,W%WDES%NB_PAR)+1,BUF%requests_gth(BUF%nrequests_gth))
            BUF%nrequests_gth=BUF%nrequests_gth+1
            CALL M_ibcast_z_from(W%WDES%COMM_INTER,BUF%WIN(NI)%CR(1), &
           &     SIZE(BUF%WIN(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,BUF%requests_gth(BUF%nrequests_gth))
            IF (W%WDES%LOVERL) THEN
               BUF%nrequests_gth=BUF%nrequests_gth+1

               CALL M_ibcast_z_from(W%WDES%COMM_INTER,BUF%WIN(NI)%CPROJ(1), &
              &     SIZE(BUF%WIN(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,BUF%requests_gth(BUF%nrequests_gth))
# 937

            ENDIF
         ENDDO
# 942

 
!        CALL M_waitall(BUF%nrequests_gth,BUF%requests_gth(1))

      ENDIF

      BUF%LDO_FWD=.TRUE.

      

      RETURN
      END SUBROUTINE FFT_AND_GATHER


!************************ SUBROUTINE REDUCE ****************************
!
!> Reduce WXI(i)\%CR i=1,...,NB2-NB1+1 over all nodes in COMM_inter_node,
!> and store the result only on the node that owns band i+NB1-1.
!>
!> NB1 and NB2 are global band indices.
!>
!> We use non-blocking reduction_to and do NOT wait for them to end.
!
!***********************************************************************

      SUBROUTINE REDUCE(BUF)
# 970

      TYPE (buffer) :: BUF

! local variables
      INTEGER :: NB1,NB2,N,NI
# 977


      
# 983


! reduce WXI

      IF (BUF%WXI(1)%WDES1%DO_REDIS) THEN
 
         NB1=BUF%NBSTRT_wxi
         NB2=BUF%NBSTOP_wxi

# 994

         BUF%nrequests_red=0
         DO N=NB1,NB2
            NI=N-NB1+1

            BUF%nrequests_red=BUF%nrequests_red+1
            CALL M_ireduce_z_to(BUF%WXI(NI)%WDES1%COMM_INTER,BUF%WXI(NI)%CR(1), &
           &     SIZE(BUF%WXI(NI)%CR),MOD(N-1,BUF%WXI(NI)%WDES1%NB_PAR)+1,BUF%requests_red(BUF%nrequests_red))

            IF (BUF%WXI(1)%WDES1%LOVERL) THEN
               BUF%nrequests_red=BUF%nrequests_red+1

               CALL M_ireduce_z_to(BUF%WXI(NI)%WDES1%COMM_INTER,BUF%WXI(NI)%CPROJ(1), &
              &     SIZE(BUF%WXI(NI)%CPROJ),MOD(N-1,BUF%WXI(NI)%WDES1%NB_PAR)+1,BUF%requests_red(BUF%nrequests_red))
# 1011

            ENDIF
         ENDDO
# 1016

 
!        CALL M_waitall(BUF%nrequests_red,BUF%requests_red(1))
 
      ENDIF

      

      RETURN
      END SUBROUTINE REDUCE


!************************ SUBROUTINE WAITGTH ***************************
!
!> Wait for the broadcast operations in BUF, started by a call
!> to #FFT_AND_GATHER, to end.
!
!***********************************************************************

      SUBROUTINE WAITGTH(BUF)
      TYPE (buffer) :: BUF

      
# 1048


      IF (BUF%nrequests_gth>0) THEN
         CALL M_waitall(BUF%nrequests_gth,BUF%requests_gth(1))
         BUF%nrequests_gth=0
      ENDIF
# 1056


# 1060

      

      RETURN
      END SUBROUTINE WAITGTH
 

!************************ SUBROUTINE WAITRED ***************************
!
!> Wait for the reduction operations in BUF, started by a call
!> to #REDUCE, to end.
!
!***********************************************************************

      SUBROUTINE WAITRED(BUF)
      TYPE (buffer) :: BUF
      
# 1086


      IF (BUF%nrequests_red>0) THEN
         CALL M_waitall(BUF%nrequests_red,BUF%requests_red(1))
         BUF%nrequests_red=0
      ENDIF

# 1095

!     CALL M_barrier(W%WDES%COMM_INTER)

      

      RETURN
      END SUBROUTINE WAITRED
 

!************************ SUBROUTINE FOCK_FWD **************************
!
!> Compute the action of the Fock potential on \phi_i:
!> ~~~
!>   X_i(r) = \sum_j \phi_j(r)\int \phi^*_j(r')\phi_i(r')/|r-r'| d3r',
!> ~~~
!> where i labels the orbitals gathered into BUF\%WIN, and j runs
!> over the orbitals owned locally by the 1-ranks.
!>
!> The result is stored in BUF\%WXI.
!
!***********************************************************************

      SUBROUTINE FOCK_FWD(BUF)
      USE sym_prec

      TYPE (buffer) :: BUF

! local variables
      TYPE (wavedes1), TARGET :: WDESQ_IRZ

      INTEGER :: NBSTART,NBSTOP,NDO
      INTEGER :: NQ_USED,NQ,MQ
      INTEGER :: N,M,NS
      INTEGER :: ISP_IRZ

      REAL(q) :: WEIGHT_Q,WEIGHT
      REAL(q) :: EXHF

      REAL(q) :: EXX
# 1136

      REAL(q) :: FD

      LOGICAL :: LSHIFT
# 1142


!$    INTEGER :: i,NSTRIP,NSTRIP_ACT
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      

# 1151


      NBSTART=BUF%NBSTRT_win
      NBSTOP =BUF%NBSTOP_win
      NDO    =NBSTOP-NBSTART+1

      EXHF=0


      CDIJ=0
# 1167


      IF (WHF%WDES%LOVERL) THEN
!$ACC KERNELS DEFAULT(PRESENT) 
         CDLM=0
         DO N=1,NDO
            BUF%WXI(N)%CR=0; BUF%WXI(N)%CPROJ=0
         ENDDO
!$ACC END KERNELS
      ELSE
!$ACC KERNELS DEFAULT(PRESENT) 
         DO N=1,NDO
            BUF%WXI(N)%CR=0
         ENDDO
!$ACC END KERNELS
      ENDIF

      NQ_USED=0
!==========================================================================
!  loop over all q-points (index NQ)
!  sum_nq phi_nq mq (r') \int phi_nq mq(r) phi_nk mk(r) / (r-r') d3r
!==========================================================================
      qpoints: DO NQ=1,KPOINTS_FULL%NKPTS
         IF( KPOINTS_FULL%WTKPT(NQ)==0 .OR. &
            (HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,NQ))) .OR. &
            (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-WHF%WDES%VKPT(:,NK)))) CYCLE qpoints

         NQ_USED=NQ_USED+1
         WEIGHT_Q=1

         IF (ALLOCATED(WEIGHT_K_POINT_PAIR_SMALL_GROUP) .AND. PRESENT(LSYMGRAD) ) THEN
            IF (LSYMGRAD) THEN
               IF (WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)==0) CYCLE qpoints
               WEIGHT_Q=WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)
            ENDIF
         ENDIF

         CALL SETWDES(WHF%WDES,WDESQ,NQ)
         CALL SETWDES(WHF%WDES,WDESQ_IRZ,KPOINTS_FULL%NEQUIV(NQ))

         ISP_IRZ=ISP
         IF (KPOINTS_FULL%SPINFLIP(NQ)==1) THEN
            ISP_IRZ=3-ISP
         ENDIF

! set POTFAK for this q and k point
         CALL SET_GFAC(GRIDHF,LATT_CUR,NK,NQ,FSG,POTFAK)

! loop over bands mq (occupied bands for present q-point on the local CPU)
         mband: DO MQ=1,WHF%WDES%NBANDS
            IF (ABS(WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))<=1E-10) CYCLE mband
            IF ((MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW<NBANDSGWLOW_FOCK) CYCLE mband

            IF (NQ<=WHF%WDES%NKPTS) THEN
               CALL W1_COPY(ELEMENT(WHF, WDESQ, MQ, ISP), WQ)
               CALL FFTWAV_W1(WQ)
            ELSE

! symmetry must be considered if the wavefunctions for this
! k-point NQ (containing all k-points in the entire BZ)
! are not stored in W
               LSHIFT=.FALSE.
               IF ((ABS(KPOINTS_FULL%TRANS(1,NQ))>TINY) .OR. &
                   (ABS(KPOINTS_FULL%TRANS(2,NQ))>TINY) .OR. &
                   (ABS(KPOINTS_FULL%TRANS(3,NQ))>TINY)) LSHIFT=.TRUE.

               CALL W1_ROTATE_AND_FFT(WQ, ELEMENT(WHF, WDESQ_IRZ, MQ, ISP_IRZ), ROT_HANDLE, P, LATT_CUR, LSHIFT)

            ENDIF
!-----------------------------------------------------------------------------
! calculate fock potential and add to accelerations
!-----------------------------------------------------------------------------
! calculate charge phi_q nq(r) phi_k nk(r)
            CALL FOCK_CHARGE_MU( BUF%WIN(1:NDO), WQ, GWORK, CRHOLM)

            WEIGHT=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)/GRIDHF%NPLWV*WEIGHT_Q


!!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N,EXX,NS,FD) REDUCTION(+:EXHF) IF (MOD(NBLK,omp_nthreads)==0)
            nband: DO N=1,NDO
!$             i=1
               NS=NBSTART+N-1        ! global storage index of present band
! fft to reciprocal space
               CALL FFT3D(GWORK(1,N),GRIDHF,-1)
! model-GW set GFAC state dependent
               IF (MODEL_GW==2) THEN
!$                i=OMP_GET_THREAD_NUM()+1
                  CALL MODEL_GW_SET_GFAC(GRIDHF, LATT_CUR, NK, NQ, KPOINTS_FULL%NEQUIV(NQ), &
                     NS, (MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW, ISP, ISP_IRZ, FSG, POTFAK(1 ))
               ENDIF

               IF (MCALPHA/=0) THEN
! with finite size corrections:
                  CALL APPLY_GFAC_MULTIPOLE(GRIDHF, GWORK(1,N), POTFAK(1 ))
               ELSE
                  CALL APPLY_GFAC_EXCHANGE(GRIDHF, GWORK(1,N), POTFAK(1 ), EXX)
                  EXX=EXX*(0.5_q/GRIDHF%NPLWV)  ! divide by grid points
! correct for self-interaction
                  IF (NS==(MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW .AND. NQ==NK) THEN

                     IF (W%WDES%COMM_INB%NODE_ME==1) THEN
                        EXX=EXX+FSG_AEXX*0.5_q  ! (1._q,0._q) node adds corrections
                     ENDIF
# 1272

                  ENDIF
! use smaller occupancy
! FD=MIN(WHF%FERTOT(NS,NK,ISP),WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))
! use (1._q,0._q)-electron occupancy of state at greater energy
                  IF (REAL(WHF%CELTOT(NS,NK,ISP),q) > REAL(WHF%CELEN(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ),q)) THEN
                     FD=WHF%FERTOT(NS,NK,ISP)
                  ELSE
                     FD=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)
                  ENDIF

                  EXHF=EXHF-EXX &
                        *WHF%WDES%RSPIN*WHF%WDES%WTKPT(NK)*WEIGHT_Q* &
                        (FD- & 
                         WHF%FERTOT(NS,NK,ISP)*WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))
               ENDIF

! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
               CALL FFT3D(GWORK(1,N),GRIDHF,1)

! add to acceleration xi in real space
               CALL VHAMIL_TRACE(WDESK, GRID_FOCK, GWORK(1,N), WQ%CR(1), BUF%WXI(N)%CR(1), WEIGHT)
            ENDDO nband
!!!$OMP END PARALLEL DO
# 1388


            IF (WHF%WDES%LOVERL) THEN
! calculate D_LM
! build the descriptor for RPRO1
               DO N=1,NDO
                  W1(N)%CPROJ => CDLM(:,N)
               ENDDO
               AUG_DES%RINPL=WEIGHT ! multiplicator for RPRO1
!$ACC UPDATE DEVICE(AUG_DES%RINPL) 
# 1402

               CALL RPROMU_HF(FAST_AUG_FOCK, AUG_DES, W1, NDO, GWORK(1,1), SIZE(GWORK,1))
# 1413

# 1422

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N,CDIJ) !! IF(MOD(NBLK,omp_nthreads)==0)
               DO N=1,NDO
                  IF (WHF%WDES%NRSPINORS==2) CALL ZCOPY(AUG_DES%NPRO,CDLM(1,N),1,CDLM(AUG_DES%NPRO+1,N),1)
! transform D_LM -> D_lml'm'
                  CALL CALC_DLLMM_TRANS(WHF%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ, CDLM(:,N))
! add D_lml'm' to kappa_lm_N (sum over l'm')
                  CALL OVERL_FOCK(WHF%WDES, LMDIM, CDIJ(1,1,1,1), WQ%CPROJ(1), BUF%WXI(N)%CPROJ(1), .TRUE.)
               ENDDO
!$OMP END PARALLEL DO


            ENDIF
         ENDDO mband
      ENDDO qpoints

# 1441


! sanity check
      IF (((ODDONLY.OR. EVENONLY) .AND. NQ_USED*2 /=KPOINTS_FULL%NKPTS_NON_ZERO) .OR. &
          (.NOT. (ODDONLY.OR. EVENONLY) .AND. NQ_USED*NKREDX*NKREDY*NKREDZ /=KPOINTS_FULL%NKPTS_NON_ZERO)) THEN
         CALL vtutor%bug("internal error in FOCK_FCC: number of k-points incorrect " // str(NQ_USED) &
            // " " // str(KPOINTS_FULL%NKPTS_NON_ZERO) // " " // str(NKREDX) // " " // str(NKREDY) // &
            " " // str(NKREDZ), "fock_dbl.F", 1448)
      ENDIF

      EXHF_ACFDT=EXHF_ACFDT+EXHF

      BUF%LDO_FWD=.FALSE.
      BUF%LDO_BCK=.TRUE.

      BUF%NBSTRT_wxi=BUF%NBSTRT_win
      BUF%NBSTOP_wxi=BUF%NBSTOP_win

      

      RETURN
      END SUBROUTINE FOCK_FWD


!************************ SUBROUTINE FOCK_BCK **************************
!
!> Take the action of the Fock potential on \phi_i
!> ~~~
!>   X_i(r) = \sum_j \phi_j(r)\int \phi^*_j(r')\phi_i(r')/|r-r'| d3r',
!> ~~~
!> from real to reciprocal space (FFT).
!>
!> Each 1-rank works on those X_i that correspond to orbitals it
!> owns locally
!>
!> The results are stored in WXI\%CPTWFP(:,i,NK,ISP).
!
!***********************************************************************
      SUBROUTINE FOCK_BCK(BUF)
      TYPE (buffer) :: BUF

! local variables
      INTEGER :: NBSTART,NBSTOP,NDO
      INTEGER :: NS,N,NP,MM,ISPINOR

      REAL(q) :: WEIGHT

      
# 1491


      NBSTART=BUF%NBSTRT_wxi
      NBSTOP =BUF%NBSTOP_wxi
      NDO    =NBSTOP-NBSTART+1

! take the descriptor from the W%WDES:
! only difference compared to WHF%WDES is the FFT mesh (GRID and not GRIDHF)
      CALL SETWDES(W%WDES,WDESQ,NK)

! fourier transform local accelerations xi (only own bands)
!!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(N,NS,CTMP1,CTMP2,CWORK,WEIGHT,ISPINOR,NP,MM) REDUCTION(+:CDCHF)
      fft_back:DO N=1,NDO
         IF (MOD(NBSTART+N-1-1,W%WDES%NB_PAR)+1/=W%WDES%NB_LOW) CYCLE fft_back

! local storage index
         NS=(NBSTART+N-1-1)/W%WDES%NB_PAR+1

!!!$ IF (MODULO(NS,OMP_GET_NUM_THREADS())/=OMP_GET_THREAD_NUM()) CYCLE fft_back

!$ACC KERNELS PRESENT(WXI,WXI%CPTWFP) 
         WXI%CPTWFP(:,NS,NK,ISP)=0
!$ACC END KERNELS

! add CKAPPA to CXI (full acceleration on band N now in CXI)
         IF (WHF%WDES%LOVERL) THEN
! convergence correction non local part
!$ACC KERNELS PRESENT(CTMP1,WHF,WHF%CPROJ,WHF%FERWE) 
            CTMP1(:)=FSG_AEXX*WHF%CPROJ(:,NS,NK,ISP)*WHF%FERWE(NS,NK,ISP)
!$ACC END KERNELS
            CALL OVERL1(WDESK, LMDIM, CQIJ(1,1,1,1), CQIJ(1,1,1,1), 0.0_q, CTMP1(1), CTMP2(1))
!$ACC PARALLEL LOOP PRESENT(BUF,BUF%WXI,CTMP2) 
            DO NP=1,SIZE(BUF%WXI(N)%CPROJ)
               BUF%WXI(N)%CPROJ(NP)=BUF%WXI(N)%CPROJ(NP)+CTMP2(NP)
            ENDDO
            IF (NONLR_S%LREAL) THEN
!$ACC KERNELS PRESENT(CWORK) 
               CWORK=0
!$ACC END KERNELS
# 1534

                  CALL RACC0(NONLR_S, WDESQ, BUF%WXI(N)%CPROJ(1), CWORK(1))
# 1538

               DO ISPINOR=0,WDESQ%NRSPINORS-1
                  CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                       CWORK(1+ISPINOR*WDESQ%GRID%MPLWV), &
                       WXI%CPTWFP(1+ISPINOR*WDESQ%NGVECTOR,NS,NK,ISP),WDESQ%GRID,.TRUE.)
               ENDDO
            ELSE
               CALL VNLAC0(NONL_S, WDESK, BUF%WXI(N)%CPROJ(1), WXI%CPTWFP(1,NS,NK,ISP))
            ENDIF
         ENDIF

! convergence correction non local part
         CALL ZAXPY(WXI%WDES%NRPLWV,CMPLX(FSG_AEXX*WHF%FERWE(NS,NK,ISP),0._q,KIND=q),WHF%CPTWFP(1,NS,NK,ISP),1,WXI%CPTWFP(1,NS,NK,ISP),1)

! double counting hence subtract half the self energy
! and change sign since we have use e^2 to calculate the potential
         WEIGHT=WHF%FERWE(NS,NK,ISP)*WHF%WDES%WTKPT(NK)*0.5_q*WHF%WDES%RSPIN

         DO ISPINOR=0,WDESK%NRSPINORS-1
            CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW(1), &
                 BUF%WXI(N)%CR(1+ISPINOR*GRID_FOCK%MPLWV), &
                 WXI%CPTWFP(1+ISPINOR*WDESK%NGVECTOR,NS,NK,ISP),GRID_FOCK,.TRUE.)

!$ACC PARALLEL LOOP PRESENT(WXI,WXI%CPTWFP,WDESK,WHF,WHF%CPTWFP) PRIVATE(MM) REDUCTION(-:CDCHF) &
!$ACC 
            DO NP=1,WDESK%NGVECTOR
               MM=NP+ISPINOR*WDESK%NGVECTOR
               WXI%CPTWFP(MM,NS,NK,ISP)=-WXI%CPTWFP(MM,NS,NK,ISP)
               CDCHF=CDCHF-CONJG(WXI%CPTWFP(MM,NS,NK,ISP))*WHF%CPTWFP(MM,NS,NK,ISP)*WEIGHT
            ENDDO
         ENDDO
      ENDDO fft_back
!!!$OMP END PARALLEL

      BUF%LDO_BCK=.FALSE.

      

      RETURN
      END SUBROUTINE FOCK_BCK

!***********************************************************************
!***********************************************************************
!
! Internal subroutines: end
!
!***********************************************************************
!***********************************************************************

      END SUBROUTINE FOCK_ALL_DBLBUF

      END MODULE
