# 1 "subrot_lr.F"
!#define dotiming
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


# 4 "subrot_lr.F" 2 

MODULE subrot_lr
  USE prec
  USE dfast
CONTAINS
!************************ SUBROUTINE EDDIAG_LR *************************
!> RCS:  $Id: subrot.F,v 1.10 2003/06/27 13:22:23 kresse Exp kresse $
!> this subroutine performs Loewdin perturbation theory  in the
!> sub space spanned by the calculated orbitals
!> more specifically it recalculates the wavefunctions W
!>
!> ~~~
!>                            |phi(0)_j> <phi(0)_j| xi_i> +  K_ji
!> phi_i(1) ->phi_i(1) +sum_j -------------------------------------------
!>                                      e(0)_i - e(0)_j
!> ~~~
!>
!> where
!> ~~~
!>  K_ji =  <phi(0)_j| (H(0) - e_i(0)S(0))  |phi_i(1)>
!>          K_ji =  (e_j - e_i(0))  <phi(0)_j| S(0) |phi_i(1)>
!> ~~~
!>
!> this leads to the update equation
!>
!> ~~~
!>                             |phi(0)_j> <phi(0)_j| xi_i>
!> phi_i(1) ->phi_i(1) +sum_j ---------------------------  - O_ji
!>                                     e(0)_i - e(0)_j
!> ~~~
!>
!> with
!>  ~~~
!>  O_ji =  <phi(0)_j| S(0) |phi(1)_i>
!>  ~~~
!>
!> upon exit of the k-points and spin loop CHAM contains
!> ~~~
!>                    <phi(0)_j| xi_i>
!>  CHAM(j,i) =   ---------------------------
!>                    e(0)_i - e(0)_j
!> ~~~
!
!***********************************************************************

  SUBROUTINE EDDIAG_LR(W,W0,WXI,LMDIM,CQIJ,LOVERL,LHERM,CSHIFT,IU0,DEG_CLUSTER,RESOLVE_DEG)

!! USE moffload

      USE prec
      USE wave_mpi
      USE wave_high
      USE mpimy
      USE nonl_high
      USE hamil
      USE constant
      USE main_mpi
      USE subrot_cluster

      IMPLICIT NONE

      TYPE (wavespin)    W
      TYPE (wavespin)    W0
      TYPE (wavespin)    WXI

      INTEGER LMDIM
      COMPLEX(q) CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
!> debug output
      INTEGER IU0
!> overlap matrix used (i.e. US-PP, PAW)
      LOGICAL LOVERL
!> <phi(0)_j| xi_i> hermitian or not
      LOGICAL LHERM
      TYPE (eigenf_cluster_pointer), OPTIONAL:: DEG_CLUSTER(:,:)
      LOGICAL , OPTIONAL:: RESOLVE_DEG
      REAL(q) CSHIFT
      
! local variables
      COMPLEX(q)    CSUM
      INTEGER ISP, NK, N_TO, NSTRIP
      INTEGER NPOS, NSTRIP_ACT, N
      INTEGER NB_TOT, NB_TOT_W0, NDONE, NP, ISPINOR, M, MM, NPOS_RED, NSTRIP_RED
      INTEGER N1, N2
      REAL(q) DIFCEL, FAKT
      COMPLEX(q) CROT
      LOGICAL :: LSETCHAM=.TRUE.        ! reset the matrix

      TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
      TYPE (wavedes1)    WDES1_W0       ! descriptor for (1._q,0._q) k-point for W0
      TYPE (wavefuna)    WA             ! array to store wavefunction
      TYPE (wavefuna)    WAXI           ! array to store wavefunction
      TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character

      COMPLEX(q),ALLOCATABLE,TARGET::  CHAM(:,:),COVL(:,:)

      

# 111


!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      NB_TOT   =W %WDES%NB_TOT
      NB_TOT_W0=W0%WDES%NB_TOT
      IF (PRESENT(DEG_CLUSTER) .AND. PRESENT(RESOLVE_DEG) .AND. NB_TOT/=NB_TOT_W0) THEN
         CALL vtutor%bug("internal error in EDDIAG_LR: resolving degeneracies is presently not " &
            // "tested for non square matrices", "subrot_lr.F", 120)
      ENDIF

! set NSTRIP between [1 and 32]
      NSTRIP=NSTRIP_STANDARD_GLOBAL

!$ACC ENTER DATA CREATE(WNONL,WDES1,WDES1_W0) 
      CALL SETWDES(W%WDES,WDES1,0)
      CALL NEWWAVA_PROJ(WNONL, WDES1)

      ALLOCATE(CHAM(NB_TOT_W0,NB_TOT),COVL(NB_TOT_W0,NB_TOT))
!$ACC ENTER DATA CREATE(CHAM,COVL) 
!=======================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoint: DO NK=1,W0%WDES%NKPTS

      IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
      CALL SETWDES(W%WDES ,WDES1,NK)
      CALL SETWDES(W0%WDES,WDES1_W0,NK)
      WA  =ELEMENTS(W, WDES1, ISP)
      WAXI=ELEMENTS(WXI, WDES1, ISP)
!=======================================================================
!  calculate Hamiltonian CHAM
!                <phi(0)_j| xi_i>
!  CHAM(j,i) =   ----------------
!                e(0)_i - e(0)_j
!=======================================================================
   IF (LSETCHAM) THEN

      CALL OVERL(WDES1, .TRUE.,LMDIM,CQIJ(1,1,1,ISP),W%CPROJ(1,1,NK,ISP),WNONL%CPROJ(1,1))
! redistribute the projected wavefunctions
! wavefunctions are still required at this point
      IF (W0%WDES%DO_REDIS) THEN
         CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, WNONL%CPROJ(1,1))
         CALL REDIS_PW  (WDES1_W0, W0%WDES%NBANDS, W0%CPTWFP   (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1_W0, W0%WDES%NBANDS, W0%CPROJ(1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, W%WDES%NBANDS, W%CPTWFP    (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, W%CPROJ (1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, W%WDES%NBANDS, WXI%CPTWFP  (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, WXI%CPROJ(1,1,NK,ISP))
      ENDIF

!$ACC PARALLEL PRESENT(CHAM) 
      CHAM=0
!$ACC END PARALLEL
      NDONE=0

! first calculate CHAM(j,i) =  <phi(0)_j| xi_i>
      DO NPOS=1,NB_TOT,NSTRIP
         NSTRIP_ACT=MIN(NB_TOT+1-NPOS, NSTRIP)

         IF (LOVERL .OR. .NOT. LHERM) THEN
! this matrix is Hermitian only if no overlap matrix exists
            CALL ORTH2( &
                 W0%CPTWFP(1,1,NK,ISP),WAXI%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
                 WNONL%CPROJ_RED(1,NPOS),NB_TOT_W0, &
                 NPOS,NSTRIP_ACT,WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1))
         ELSE
            CALL ORTH1("U", &
                 W0%CPTWFP(1,1,NK,ISP),WAXI%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
                 WNONL%CPROJ_RED(1,NPOS),NB_TOT_W0, &
                 NPOS,NSTRIP_ACT,WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1))
         ENDIF

      ENDDO
      CALL M_sum_z(W0%WDES%COMM_KIN,CHAM(1,1),NB_TOT_W0*NB_TOT)
# 190

1     FORMAT(1I2,3X,20F9.5)
2     FORMAT(1I2,3X,20E9.1)


! add Hermitian elements
      IF (.NOT. (LOVERL .OR. .NOT. LHERM)) THEN
!$ACC PARALLEL LOOP PRESENT(CHAM) 
         DO N2=1,NB_TOT_W0
         DO N1=1,N2-1
            CHAM(N2,N1)=CONJG(CHAM(N1,N2))
         ENDDO
         ENDDO
      ENDIF
!
! resolve degeneragy problem
! this should be 1._q via the DEG_CLUSTER structure which
! stores a list of degenerated eigenvalue/eigenfunction pairs
!
      IF (PRESENT(DEG_CLUSTER).AND.PRESENT(RESOLVE_DEG)) THEN
         CALL SETUP_DEG_CLUSTERS(DEG_CLUSTER(NK,ISP)%DEG_CLUSTER, CHAM,  &
              WXI%CELTOT(:,NK,ISP), W0%CELTOT(:,NK,ISP))
      ELSE IF(PRESENT(DEG_CLUSTER).AND.CSHIFT==0) THEN
         CALL ZERO_HAM_DEG_CLUSTERS(DEG_CLUSTER(NK,ISP)%DEG_CLUSTER, CHAM)
      ENDIF
! now set CHAM(j,i) =  <phi(0)_j| xi_i> / ((e0(i)-e0(j))
!  (Loewdin perturbation theory)
!$ACC PARALLEL LOOP PRESENT(W0,CHAM) PRIVATE(DIFCEL) 
      DO N2=1,NB_TOT
      DO N1=1,NB_TOT_W0
         DIFCEL= REAL( W0%CELTOT(N2,NK,ISP)-W0%CELTOT(N1,NK,ISP) ,KIND=q)
         IF (ABS(DIFCEL)<1E-10_q .AND. CSHIFT==0) THEN
            CHAM(N1,N2)=0
         ELSE
! it is inclear whether 2*CSHIFT or CSHIFT is correct here
! at least for ionic contributions of elastic tensor in SiO2 only first version works
            CHAM(N1,N2)=CHAM(N1,N2)/(DIFCEL+CMPLX(0.0_q,2*CSHIFT,q))
!            CHAM(N1,N2)=CHAM(N1,N2)/(DIFCEL+CMPLX(0.0_q,CSHIFT,q))
         ENDIF
      ENDDO
      ENDDO
   END IF
!-----------------------------------------------------------------------
! calculate the overlap matrix
! OVERL(j,i) =  <phi(0)_j| S(0) |phi(1)_i>
!-----------------------------------------------------------------------
!$ACC PARALLEL PRESENT(COVL) 
      COVL=(0._q,0._q)
!$ACC END PARALLEL

      DO NPOS=1,NB_TOT,NSTRIP
         NSTRIP_ACT=MIN(NB_TOT+1-NPOS, NSTRIP)

         CALL ORTH2( &
              W0%CPTWFP(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
              WNONL%CPROJ_RED(1,NPOS),NB_TOT_W0, &
              NPOS,NSTRIP_ACT,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
      ENDDO

      CALL M_sum_z(W0%WDES%COMM_KIN,COVL(1,1),NB_TOT_W0*NB_TOT)
# 252

!=======================================================================
! now add COVL to CHAM
!=======================================================================
!$ACC PARALLEL LOOP PRESENT(CHAM,COVL) 
      DO N2=1,NB_TOT
      DO N1=1,NB_TOT_W0
         CHAM(N1,N2)=CHAM(N1,N2)-COVL(N1,N2)
      ENDDO
      ENDDO
# 264

!=======================================================================
! rotate wavefunctions
!=======================================================================
      IF (WDES1%NPL_RED/=0) &
      CALL ZGEMM('N', 'N',  WDES1%NPL_RED, NB_TOT, NB_TOT_W0, (1._q,0._q), &
     &               W0%CPTWFP(1,1,NK,ISP),  WDES1%NRPLWV_RED, CHAM(1,1), &
     &               NB_TOT_W0, (1._q,0._q), W%CPTWFP(1,1,NK,ISP) ,  WDES1%NRPLWV_RED)
      IF (WDES1%NPRO_RED/=0) &
      CALL ZGEMM('N', 'N', WDES1%NPRO_RED, NB_TOT, NB_TOT_W0, (1._q,0._q), &
     &               W0%CPROJ(1,1,NK,ISP), WDES1%NPROD_RED, CHAM(1,1), &
     &               NB_TOT_W0, (1._q,0._q), W%CPROJ(1,1,NK,ISP) , WDES1%NPROD_RED)

      ! "lincom ok"

      IF (PRESENT(DEG_CLUSTER) .AND. PRESENT(RESOLVE_DEG)) THEN
         IF (ASSOCIATED(DEG_CLUSTER(NK,ISP)%DEG_CLUSTER)) THEN
            CALL SUBROT_DEG_CLUSTERS(W0%WDES, WDES1%NPL_RED, WDES1%NPRO_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, &
                W0%CPTWFP(:,:,NK,ISP), W0%CPROJ(:,:,NK,ISP), DEG_CLUSTER(NK,ISP)%DEG_CLUSTER, .FALSE., .FALSE.)
            CALL SUBROT_DEG_CLUSTERS(W0%WDES, WDES1%NPL_RED, WDES1%NPRO_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, &
                W%CPTWFP(:,:,NK,ISP),  W%CPROJ(:,:,NK,ISP), DEG_CLUSTER(NK,ISP)%DEG_CLUSTER, .FALSE., .FALSE.)
            CALL SUBROT_DEG_CLUSTERS(W0%WDES, WDES1%NPL_RED, WDES1%NPRO_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, &
                WXI%CPTWFP(:,:,NK,ISP), WXI%CPROJ(:,:,NK,ISP), DEG_CLUSTER(NK,ISP)%DEG_CLUSTER, .FALSE., .FALSE.)
         ENDIF
      ENDIF
      
      IF (W0%WDES%DO_REDIS) THEN
         CALL REDIS_PW  (WDES1_W0, W0%WDES%NBANDS, W0%CPTWFP   (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1_W0, W0%WDES%NBANDS, W0%CPROJ(1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, W%WDES%NBANDS, W%CPTWFP    (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, W%CPROJ (1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, W%WDES%NBANDS, WXI%CPTWFP  (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, WXI%CPROJ(1,1,NK,ISP))
         ! "redis ok"
      ENDIF

!=======================================================================
! restore CHAM
!=======================================================================
!$ACC PARALLEL LOOP PRESENT(CHAM, COVL) 
      DO N2=1,NB_TOT
      DO N1=1,NB_TOT_W0
         CHAM(N1,N2)=CHAM(N1,N2)+COVL(N1,N2)
      ENDDO
      ENDDO

!=======================================================================
      ENDDO kpoint
      ENDDO spin
!=======================================================================
# 317

      CALL DELWAVA_PROJ(WNONL)
!$ACC EXIT DATA DELETE(CHAM,COVL,WNONL) 
      DEALLOCATE(CHAM,COVL)

# 336


      

      RETURN
      END SUBROUTINE EDDIAG_LR


!************************ SUBROUTINE HRPSI *****************************
!
!***********************************************************************

      SUBROUTINE HRPSI(W,W0,WXI,LMDIM,CQIJ,LOVERL,LHERM,IU0)

      USE prec
      USE wave_mpi
      USE wave_high
      USE mpimy
      USE nonl_high
      USE hamil
      USE constant
      USE main_mpi
      USE subrot_cluster

      IMPLICIT NONE

      TYPE (wavespin)    W
      TYPE (wavespin)    W0
      TYPE (wavespin)    WXI

      INTEGER LMDIM
      COMPLEX(q) CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
!> debug output
      INTEGER IU0
!> overlap matrix used (i.e. US-PP, PAW)
      LOGICAL LOVERL
!> <phi(0)_j| xi_i> hermitian or not
      LOGICAL LHERM

! local variables
      INTEGER ISP, NK
      INTEGER NPOS, NSTRIP, NSTRIP_ACT
      INTEGER NB_TOT, NB_TOT_W0, NP, NPOS_RED, NSTRIP_RED
      INTEGER N1, N2
      REAL(q) DIFCEL

      TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
      TYPE (wavedes1)    WDES1_W0       ! descriptor for (1._q,0._q) k-point for W0
      TYPE (wavefuna)    WA             ! array to store wavefunction
      TYPE (wavefuna)    WAXI           ! array to store wavefunction

      COMPLEX(q), ALLOCATABLE, TARGET::  CH1(:,:),CH2(:,:)
! redistributed plane wave coefficients
      INTEGER NCPU


      NCPU   =W0%WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 395

!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      NB_TOT   =W %WDES%NB_TOT
      NB_TOT_W0=W0%WDES%NB_TOT

! set NSTRIP between [1 and 32]
      NSTRIP=NSTRIP_STANDARD_GLOBAL

      CALL SETWDES(W%WDES,WDES1,0)

      ALLOCATE(CH1(NB_TOT_W0,NB_TOT),CH2(NB_TOT_W0,NB_TOT))
!=======================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoint: DO NK=1,W0%WDES%NKPTS

      IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
      CALL SETWDES(W%WDES ,WDES1,NK)
      CALL SETWDES(W0%WDES,WDES1_W0,NK)
      WA  =ELEMENTS(W, WDES1, ISP)
      WAXI=ELEMENTS(WXI, WDES1, ISP)

! redistribute the projected wavefunctions
! wavefunctions are still required at this point
      IF (W0%WDES%DO_REDIS) THEN
        CALL REDIS_PW  (WDES1_W0, W0%WDES%NBANDS, W0%CPTWFP   (1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1_W0, W0%WDES%NBANDS, W0%CPROJ(1,1,NK,ISP))
        CALL REDIS_PW  (WDES1, W%WDES%NBANDS, W%CPTWFP    (1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, W%CPROJ (1,1,NK,ISP))
        CALL REDIS_PW  (WDES1, W%WDES%NBANDS, WXI%CPTWFP  (1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, WXI%CPROJ(1,1,NK,ISP))
      ENDIF

      CH1=0

! first calculate CHAM(j,i)= \sum_k < phi(0)_j|p_k> (\sum_l Q_kl <d/dk p_l| - tau_kl <p_l|) |phi(0)_i>
      DO NPOS=1,NB_TOT,NSTRIP
         NSTRIP_ACT=MIN(NB_TOT+1-NPOS, NSTRIP)

         IF (LOVERL .OR. .NOT. LHERM) THEN
! this matrix is Hermitian only if no overlap matrix exists
            CALL ORTH2( &
                 W0%CPTWFP(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
                 WA%CPROJ_RED(1,NPOS),NB_TOT_W0, &
                 NPOS,NSTRIP_ACT,0,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CH1(1,1))
         ELSE
            CALL ORTH1("U", &
                 W0%CPTWFP(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
                 WA%CPROJ_RED(1,NPOS),NB_TOT_W0, &
                 NPOS,NSTRIP_ACT,0,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CH1(1,1))
         ENDIF
      ENDDO
      CALL M_sum_z(W0%WDES%COMM_KIN,CH1(1,1),NB_TOT_W0*NB_TOT)
# 453

      DO N2=1,NB_TOT
      DO N1=1,NB_TOT_W0
         DIFCEL= REAL( W0%CELTOT(N2,NK,ISP)-W0%CELTOT(N1,NK,ISP) ,KIND=q)
         CH1(N1,N2)=CH1(N1,N2)*DIFCEL
      ENDDO
      ENDDO
# 462


      WA%CPROJ_RED=0

      CH2 = 0

! CH2(j,i) = <phi(0)_j| xi_i>
      DO NPOS=1,NB_TOT,NSTRIP
         NSTRIP_ACT=MIN(NB_TOT+1-NPOS, NSTRIP)

         IF (LOVERL .OR. .NOT. LHERM) THEN
! this matrix is Hermitian only if no overlap matrix exists
            CALL ORTH2( &
                 W0%CPTWFP(1,1,NK,ISP),WAXI%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
                 WAXI%CPROJ_RED(1,NPOS),NB_TOT_W0, &
                 NPOS,NSTRIP_ACT,WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CH2(1,1))
         ELSE
            CALL ORTH1("U", &
                 W0%CPTWFP(1,1,NK,ISP),WAXI%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
                 WAXI%CPROJ_RED(1,NPOS),NB_TOT_W0, &
                 NPOS,NSTRIP_ACT,WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CH2(1,1))
         ENDIF

      ENDDO
      CALL M_sum_z(W0%WDES%COMM_KIN,CH2(1,1),NB_TOT_W0*NB_TOT)
# 489

1     FORMAT(1I2,3X,20F9.5)
2     FORMAT(1I2,3X,20E9.1)

      CH1=CH1+CH2 

! test_
! add trace elements WXI%CELTOT
      DO N1=1,MIN(NB_TOT_W0,NB_TOT)
         CH1(N1,N1)=CH1(N1,N1)+WXI%CELTOT(N1,NK,ISP)
      ENDDO
! test_

! add Hermitian elements
      IF (.NOT. (LOVERL .OR. .NOT. LHERM)) THEN
         DO N2=1,NB_TOT_W0
         DO N1=1,N2-1
            CH1(N2,N1)=CONJG(CH1(N1,N2))
         ENDDO
         ENDDO
      ENDIF

# 513

!=======================================================================
! rotate wavefunctions
!=======================================================================
      IF (WDES1%NPL_RED/=0) &
      CALL ZGEMM('N', 'N',  WDES1%NPL_RED, NB_TOT, NB_TOT_W0, (1._q,0._q), &
     &               W0%CPTWFP(1,1,NK,ISP),  WDES1%NRPLWV_RED, CH1(1,1), &
     &               NB_TOT_W0, (1._q,0._q), W%CPTWFP(1,1,NK,ISP) ,  WDES1%NRPLWV_RED)
      IF (WDES1%NPRO_RED/=0) &
      CALL ZGEMM('N', 'N', WDES1%NPRO_RED, NB_TOT, NB_TOT_W0, (1._q,0._q), &
     &               W0%CPROJ(1,1,NK,ISP), WDES1%NPROD_RED, CH1(1,1), &
     &               NB_TOT_W0, (1._q,0._q), W%CPROJ(1,1,NK,ISP) , WDES1%NPROD_RED)

      ! "lincom ok"

      IF (W0%WDES%DO_REDIS) THEN
         CALL REDIS_PW  (WDES1_W0, W0%WDES%NBANDS, W0%CPTWFP   (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1_W0, W0%WDES%NBANDS, W0%CPROJ(1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, W%WDES%NBANDS, W%CPTWFP    (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, W%CPROJ (1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, W%WDES%NBANDS, WXI%CPTWFP  (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, W%WDES%NBANDS, WXI%CPROJ(1,1,NK,ISP))
         ! "redis ok"
      ENDIF

!=======================================================================
      ENDDO kpoint
      ENDDO spin
!=======================================================================
      DEALLOCATE(CH1,CH2)

      RETURN
      END SUBROUTINE HRPSI


!************************ SUBROUTINE ORTHO_LR **************************
!
!> this subroutine removes any hermitian or anti-hermitian (unitary)
!> part from the first order change of wave functions phi(1)
!> upon calling
!>
!> W%CPTWFP must store     <G | phi(1)>
!> W%CPROJ must store  <p(1) |  phi(0)> + <p(0) | phi(1)>
!>   (i.e. total change of wave function character)
!>
!> after the call the first order change of the orbitals is
!> essentially S orthogonal to the occupied bands
!> i.e. after calling the routine the following equation holds:
!> anti-hermitian case:
!> ~~~
!>     <phi(0)_j | S(0) | phi(1)_i> + <phi(0)_i | S(0) | phi(1)_j>*
!>   + <phi(0)_j | S(1) | phi(0)_i>   = 0
!> ~~~
!> where S(1) is the first order change of the overlap operator
!> hermitian case:
!> ~~~
!>     <phi(0)_j | S(0) | phi(1)_i> - <phi(0)_i | S(0) | phi(1)_j>*
!>   + <phi(0)_j | S(1) | phi(0)_i>   = 0
!> ~~~
!>
!> furthermore the routine also removes the diagonal part
!> ~~~
!>     <phi(0)_i | S(0) | phi(1)_i> + <phi(1)_i | S(0) | phi(0)_i>
!>   + <phi(0)_i | S(1) | phi(0)_i>   = 0
!> ~~~
!> This is required since the orbital must remain S orthogonal
!> ~~~
!>   <phi(0)_i+delta phi(1)_i| S(0)+ delta S(1) |phi(0)_i+delta phi(1)_i>=1
!> ~~~
!>
!> note that all linear response routine yields a response compatible with
!> ~~~
!>   <phi(0)_i | S(0) | phi(1)_i> = 0
!> ~~~
!> which is obviously not quite right and should rather be
!> ~~~
!>   <phi(0)_i | S(0) | phi(1)_i> + <phi(0)_i | S(1) | phi(0)_i> = 0
!> ~~~
!
!***********************************************************************

  SUBROUTINE ORTHO_LR(W,W0,WDES,LMDIM,CQIJ,LOVERL,IU0,LHERM)
      USE prec
      USE wave_mpi
      USE wave
      USE wave_high
      USE mpimy

      IMPLICIT NONE

      TYPE (wavespin)    W
      TYPE (wavespin)    W0
      TYPE (wavedes)     WDES

      INTEGER LMDIM
      COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!> debug output
      INTEGER IU0
!> overlap matrix used (i.e. US-PP, PAW)
      LOGICAL LOVERL
      LOGICAL LHERM
! local variables
      INTEGER ISP, NK, NBANDS, NSTRIP
      INTEGER NPOS,  N
      INTEGER NB_TOT
      INTEGER N1, N2
! work arrays (do max of 16 strips simultaneously)
      TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point

      COMPLEX(q),ALLOCATABLE,TARGET::  COVL(:,:)
! redistributed plane wave coefficients
      TYPE (wavefuna)    WA             ! array to store wavefunction
      TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
      INTEGER NCPU
      COMPLEX(q) :: C

      


      NCPU   =WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 635

!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS

! set NSTRIP between [1 and 32]
      NSTRIP=NSTRIP_STANDARD_GLOBAL

      CALL SETWDES(WDES,WDES1,0)
      CALL NEWWAVA_PROJ(WNONL, WDES1)

! allocate work space
      ALLOCATE(COVL(NB_TOT,NB_TOT))

!=======================================================================
      spin:  DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS

      IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
      CALL SETWDES(WDES,WDES1,NK)
      WA=ELEMENTS(W, WDES1, ISP)

!-----------------------------------------------------------------------
! calculate the overlap matrix
! OVERL(j,i) =  <phi(0)_j| S(0) |phi(1)_i>
!-----------------------------------------------------------------------
      CALL OVERL(WDES1, .TRUE.,LMDIM,CQIJ(1,1,1,ISP),W%CPROJ(1,1,NK,ISP),WNONL%CPROJ(1,1))
! redistribute the projected wavefunctions
! wavefunctions are still required at this point
      IF (WDES%DO_REDIS) THEN
        CALL REDISTRIBUTE_PROJ(WNONL)
        CALL REDIS_PW  (WDES1, NBANDS, W0%CPTWFP   (1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1, NBANDS, W0%CPROJ(1,1,NK,ISP))
        CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP    (1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ (1,1,NK,ISP))
      ENDIF

      COVL=(0._q,0._q)

      DO NPOS=1,NB_TOT-NSTRIP,NSTRIP
         CALL ORTH2( &
              W0%CPTWFP(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
              WNONL%CPROJ_RED(1,NPOS),NB_TOT, &
              NPOS,NSTRIP,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
      ENDDO

      CALL ORTH2( &
           W0%CPTWFP(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
           WNONL%CPROJ_RED(1,NPOS),NB_TOT, &
           NPOS,NB_TOT-NPOS+1,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))

      CALL M_sum_z(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT)
# 693

1     FORMAT(1I2,3X,20F9.5)
2     FORMAT(1I2,3X,20E9.1)
!-----------------------------------------------------------------------
! now remove unitary part or Hermitian part of the rotation matrix
! (unitary part for Hermitian perturbation, Hermitian part for
!              anti-Hermitian perturbation)
! after reconsideration of the equations
! Anything in the occupied manyfold must be removed since
! non unitary contributions would violate the orthogonality constraint
! and can be only the result of noise indicative of numerical errors
! the errors are a result of huge rotations among states close in energy
!-----------------------------------------------------------------------


      DO N2=1,NB_TOT
      DO N1=1,NB_TOT
         IF (FILLED(W0%FERTOT(N1,NK,ISP)) .AND. FILLED(W0%FERTOT(N2,NK,ISP)) ) THEN
            COVL(N1,N2)=-COVL(N1,N2)
         ELSE
            COVL(N1,N2)=0
         ENDIF
      ENDDO
      ENDDO
# 754

!=======================================================================
! rotate wavefunctions
!=======================================================================
      IF (WDES1%NPL_RED/=0) &
      CALL ZGEMM('N', 'N',  WDES1%NPL_RED, NB_TOT, NB_TOT, (1._q,0._q), &
     &               W0%CPTWFP(1,1,NK,ISP),  WDES1%NRPLWV_RED, COVL(1,1), &
     &               NB_TOT, (1._q,0._q), W%CPTWFP(1,1,NK,ISP) ,  WDES1%NRPLWV_RED)
      IF (WDES1%NPRO_RED/=0) &
      CALL ZGEMM('N', 'N', WDES1%NPRO_RED, NB_TOT, NB_TOT, (1._q,0._q), &
     &               W0%CPROJ(1,1,NK,ISP), WDES1%NPROD_RED, COVL(1,1), &
     &               NB_TOT, (1._q,0._q), W%CPROJ(1,1,NK,ISP) , WDES1%NPROD_RED)

      ! "lincom ok"

      IF (WDES%DO_REDIS) THEN
         CALL REDIS_PW  (WDES1, NBANDS, W0%CPTWFP   (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, NBANDS, W0%CPROJ(1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP    (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ (1,1,NK,ISP))
         ! "redis ok"
      ENDIF
!=======================================================================
      ENDDO kpoint
      ENDDO spin
!=======================================================================
      CALL DELWAVA_PROJ(WNONL)
      DEALLOCATE(COVL)

      

      RETURN
    END SUBROUTINE ORTHO_LR

!************************ SUBROUTINE ORTHO_LR_TEST *********************
!
!> this subroutine is similar to #ORTHO_LR
!>
!> but only dumps the "overlap" matrix
!> it also assumes the W%CPROJ has already been multiplied by
!> CQIJ
!
!***********************************************************************

  SUBROUTINE ORTHO_LR_TEST(W,W0,WDES,LMDIM,CQIJ,LOVERL,IU0,LHERM)
      USE prec
      USE wave_mpi
      USE wave
      USE wave_high
      USE mpimy

      IMPLICIT NONE

      TYPE (wavespin)    W
      TYPE (wavespin)    W0
      TYPE (wavedes)     WDES

      INTEGER LMDIM
      COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!> debug output
      INTEGER IU0
!> overlap matrix used (i.e. US-PP, PAW)
      LOGICAL LOVERL
      LOGICAL LHERM
! local variables
      INTEGER ISP, NK, NBANDS, NSTRIP
      INTEGER NPOS,  N
      INTEGER NB_TOT
      INTEGER N1, N2
! work arrays (do max of 16 strips simultaneously)
      TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point

      COMPLEX(q),ALLOCATABLE,TARGET::  COVL(:,:)
! redistributed plane wave coefficients
      TYPE (wavefuna)    WA             ! array to store wavefunction
      INTEGER NCPU
      COMPLEX(q) :: C
      COMPLEX(q) CRESUL(WDES%NPRO)


      NCPU   =WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 837

!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS

! set NSTRIP between [1 and 32]
      NSTRIP=NSTRIP_STANDARD_GLOBAL

      CALL SETWDES(WDES,WDES1,0)

! allocate work space
      ALLOCATE(COVL(NB_TOT,NB_TOT))

!=======================================================================
      spin:  DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS

      IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
      CALL SETWDES(WDES,WDES1,NK)
      WA=ELEMENTS(W, WDES1, ISP)

!-----------------------------------------------------------------------
! calculate the overlap matrix
! OVERL(j,i) =  <phi(0)_j| S(0) |phi(1)_i>
!----------------------------------------------------------------------
! redistribute the projected wavefunctions
! wavefunctions are still required at this point
      IF (WDES%DO_REDIS) THEN
        CALL REDIS_PW  (WDES1, NBANDS, W0%CPTWFP   (1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1, NBANDS, W0%CPROJ(1,1,NK,ISP))
        CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP    (1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ (1,1,NK,ISP))
      ENDIF

      COVL=(0._q,0._q)

      DO NPOS=1,NB_TOT-NSTRIP,NSTRIP
         CALL ORTH2( &
              W0%CPTWFP(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
              WA%CPROJ_RED(1,NPOS),NB_TOT, &
              NPOS,NSTRIP,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
      ENDDO

      CALL ORTH2( &
           W0%CPTWFP(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
           WA%CPROJ_RED(1,NPOS),NB_TOT, &
           NPOS,NB_TOT-NPOS+1,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))

      CALL M_sum_z(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT)
!#ifdef debug
      IF (IU0>=0) CALL DUMP_HAM( "ORTHO_LR: <phi(0)_j| S(0) |phi_i(1)>",WDES, COVL)
!#endif
1     FORMAT(1I2,3X,20F9.5)
2     FORMAT(1I2,3X,20E9.1)
!-----------------------------------------------------------------------
! now remove the hermitian or unitary part of the rotation matrix
!-----------------------------------------------------------------------
      IF (LHERM) THEN
      DO N2=1,NB_TOT
      DO N1=1,N2-1
         IF (FILLED(W0%FERTOT(N1,NK,ISP)) .AND. FILLED(W0%FERTOT(N2,NK,ISP)) ) THEN
            C=(COVL(N1,N2)+CONJG(COVL(N2,N1)))/2
            COVL(N1,N2)=-C
            COVL(N2,N1)=-CONJG(C)
         ELSE
            COVL(N1,N2)=0
            COVL(N2,N1)=0
         ENDIF
      ENDDO
      ENDDO
      ELSE
      DO N2=1,NB_TOT
      DO N1=1,N2-1
         IF (FILLED(W0%FERTOT(N1,NK,ISP)) .AND. FILLED(W0%FERTOT(N2,NK,ISP)) ) THEN
            C=(COVL(N1,N2)-CONJG(COVL(N2,N1)))/2
            COVL(N1,N2)=-C
            COVL(N2,N1)= CONJG(C)
         ELSE
            COVL(N1,N2)=0
            COVL(N2,N1)=0
         ENDIF
      ENDDO
      ENDDO
      ENDIF
!=======================================================================
! rotate wavefunctions
!=======================================================================
      ! "lincom ok"

      IF (WDES%DO_REDIS) THEN
         CALL REDIS_PW  (WDES1, NBANDS, W0%CPTWFP   (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, NBANDS, W0%CPROJ(1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP    (1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ (1,1,NK,ISP))
         ! "redis ok"
      ENDIF

!=======================================================================
      ENDDO kpoint
      ENDDO spin
!=======================================================================
      DEALLOCATE(COVL)

      RETURN
    END SUBROUTINE ORTHO_LR_TEST

!***********************************************************************
!
!> function that yields .TRUE. if a state is fully occupied
!>
!> this is a little bit ambiguous the threshold is presently set
!> to something 0.99 electrons
!
!***********************************************************************

  FUNCTION FILLED(F)
    USE prec
    IMPLICIT NONE
    LOGICAL FILLED
    REAL(q) :: F
    IF (ABS(F-1)<0.1) THEN
       FILLED=.TRUE.
    ELSE
       FILLED=.FALSE.
    ENDIF
  END FUNCTION FILLED

!***********************************************************************
!
!> add the first order change of the orbitals parallel to phi(0)
!>
!> this involves adding the term
!> ~~~
!>  phi(1) +=  -<phi(0)| S(1) |phi(0)> phi(0) /2
!> ~~~
!> it is essential for stability reasons to do this at the very end,
!> since outherwise the norm of the residual vector (rmm-diis_lr)
!> converges to a small non (0._q,0._q) value, which decreases the stability
!> the routine EDDIAG_LR therefore removes this element
!> (i.e.
!> ~~~
!> <phi(0)_i| S(0) |phi(1)_i> = 0
!> ~~~
!>)
!
!***********************************************************************

    SUBROUTINE NORM_LR( W1, W0, WXI,WDES)

!! USE moffload

      USE wave_mpi
      USE wave
      
      IMPLICIT NONE
      
      TYPE (wavespin)    W1
      TYPE (wavespin)    W0
      TYPE (wavespin)    WXI
      TYPE (wavedes)     WDES

      INTEGER ISP, NK, N, NPL, NPR, I

      
# 1012

      NPL=SIZE(W1%CPTWFP   ,1)
      NPR=SIZE(W1%CPROJ,1)

      DO ISP=1,WDES%ISPIN
         DO NK=1,WDES%NKPTS

            IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W1,W0,WXI) 
            DO N=1,WDES%NBANDS
               DO I=1,NPL
                  W1%CPTWFP   (I,N,NK,ISP) = W1%CPTWFP   (I,N,NK,ISP) - W0%CPTWFP   (I,N,NK,ISP) *WXI%FERWE(N,NK,ISP)/2
               ENDDO
            ENDDO
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W1,W0,WXI) 
            DO N=1,WDES%NBANDS
               DO I=1,NPR
                  W1%CPROJ(I,N,NK,ISP) = W1%CPROJ(I,N,NK,ISP) - W0%CPROJ(I,N,NK,ISP) *WXI%FERWE(N,NK,ISP)/2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
# 1043

      

    END SUBROUTINE NORM_LR

END MODULE subrot_lr
