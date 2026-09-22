# 1 "subrot.F"
!#define timing
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


# 3 "subrot.F" 2 




MODULE subrot
  USE prec
  USE dfast
  USE hamil
  USE pead, ONLY : LUSEPEAD,LPEAD_NO_SCF,W_STORE,PEAD_EIGENVALUES,PEAD_ACC_ADD_CPROJ,PEAD_ACC_ADD_PW

CONTAINS
!************************ SUBROUTINE EDDIAG ****************************
!> RCS:  $Id: subrot.F,v 1.10 2003/06/27 13:22:23 kresse Exp kresse $
!> this subroutine calculates the electronic eigenvalues and
!> optionally performes a  sub-space diagonalisation
!> i.e. unitary transforms the wavefunctions so that the Hamiltonian
!>  becomes diagonal in the subspace spanned by the wavefunctions
!>
!> IFLAG:
!>  0 only eigenvalues (without  diagonalisation no sub-space matrix)
!>  1 only eigenvalues and sub-space matrix (no diagonalisation)
!>  2 eigenvalues using diagonalisation of sub-space matrix
!>    do not rotate wavefunctions
!> 12 eigenvalues using diagonalisation of sub-space matrix
!>    no kinetic energy
!>  3 eigenvalues and sub-space diagonalisation rotate wavefunctions
!> 13 (for 13 no Jacobi algorithm is allowed)
!> 23 eigenvalues and sub-space diagonaliation in the occupied
!>    manifold (occupancies > 1-1E-3)
!>  4 eigenvalues and sub-space diagonalisation rotate wavefunctions
!>    using Loewdin pertubation theory (conserves ordering)
!>  5 eigenvalues and sub-space diagonalisation + orthogonalization
!>    unfortunately this option turns out to slow down the
!>    convergence of IALGO=48
!>
!> notes on the matrix CHAM used internally
!> it is defined as
!>  ~~~
!>  chaM(n2,n1)   = < phi_n2 | H | phi_n1 >
!>  ~~~
!> however the returned matrix CHAMHF is defined as the complex conjugated
!>  ~~~
!>  CHAMHF(n2,n1) = < phi_n2 | H | phi_n1 >* = < phi_n1 | H | phi_n2 >
!>  ~~~
!> since this is what rot.F uses
!>
!>  written by gK
!>  last update Sep 24 2004 massive cleanup
!>  optional arguments:
!>  NBANDS_MAX   maximum number of bands for which eigenvalues are
!>               calculated (only supported for IFLAG=0)
!>  \todo: Martijn new arguments
!> @details @ref openacc :
!
!***********************************************************************

  SUBROUTINE EDDIAG(HAMILTONIAN, &
       GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
       LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IU0,EXHF, & 
       CHAMHF,LFIRST,LLAST,NKSTART,NKSTOP,NBANDS_MAX,EXHF_ACFDT)

!! USE moffload

    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_struct_def
    USE constant
    USE jacobi
    USE scala
    USE main_mpi
    USE fock
    USE pseudo
    USE ini
    USE sym_grad
    USE fileio
    USE mopenmp_struct_def, ONLY : omp_nthreads

    USE fock_dbl


    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry) :: SYMM
    INTEGER LMDIM
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!> determines mode of diagonalisation
    INTEGER            IFLAG            ! determines mode of diagonalisation
!> local potential
    COMPLEX(q)   SV(GRID%MPLWV,WDES%NCDIJ) ! local potential
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    INTEGER IU0
    REAL(q) EXHF
!> if CHAMHF is present, CHAMHF is added to the subspace matrix at each iteration
!> in addition calculation of Fock part is bypassed
    COMPLEX(q), OPTIONAL :: CHAMHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
!> if LFIRST is supplied CHAMHF is set to CONJG(CHAMHF)-CHAM
    LOGICAL, OPTIONAL :: LFIRST
!> if LLAST  is supplied CHAMHF is set to CONJG(CHAM)
    LOGICAL, OPTIONAL :: LLAST
!> start k-point
    INTEGER, OPTIONAL :: NKSTART
!> stop k-point
    INTEGER, OPTIONAL :: NKSTOP 
    INTEGER, OPTIONAL :: NBANDS_MAX  ! maximum band index
    LOGICAL :: LSCAAWARE_LOCAL
    REAL(q), OPTIONAL :: EXHF_ACFDT
    
! local
! work arrays for ZHEEV (blocksize times number of bands)
    INTEGER, PARAMETER :: LWORK=32
    COMPLEX(q)       CWRK(LWORK*WDES%NB_TOT)
    REAL(q)    R(WDES%NB_TOT)
# 128

    REAL(q)    RWORK(7*WDES%NB_TOT), ABSTOL, VL, VU
    INTEGER IWORK(5*WDES%NB_TOT), INFO(WDES%NB_TOT),IL, IU, NB_CALC

! work arrays (do max of 16 strips simultaneously)

    TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
    TYPE (wavefun1)    W1             ! current wavefunction
    TYPE (wavefuna)    WA             ! array to store wavefunction
    TYPE (wavespin)    WFOCK          ! array to store the Fock (exchange contribution)
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WOVL           ! array to hold non local part Q * wave function character
    TYPE (wavefuna)    WHAM           ! array to store accelerations for a selected block
    COMPLEX(q)         CDCHF          ! HF double counting energy
    INTEGER :: ICALL=0                ! number of calls
    COMPLEX(q),ALLOCATABLE,TARGET::  CHAM(:,:),COVL(:,:) ! Hamiltonian and overlap matrix
    TYPE (wavefun1)    WTMP(NSTRIP_STANDARD)

! for asynchronous data redistribution
    TYPE (REDIS_PW_CTR),POINTER :: H_PW1, H_PW2
    LOGICAL LASYNC_LOCAL

    LOGICAL, EXTERNAL :: USEFOCK_CONTRIBUTION

    INTEGER :: NB_TOT, NBANDS, NSTRIP, ISP, NK, N, NP, NPOS, NSTRIP_ACT, NSTRIP_NEXT, &
               NPOS_RED, NSTRIP_RED, IFAIL, MY_NKSTART, MY_NKSTOP, NSIM_LOCAL
!$  INTEGER NSIM_FOCK

    

! input checks
    IF (IFLAG<0) THEN
      CALL vtutor%bug("ERROR in EDDIAG: IFLAG not properly set.", "subrot.F", 160)
    ENDIF

# 180


    LscaAWARE_LOCAL=LscaAWARE.AND.(.NOT.(IFLAG==4))

# 188


    IF (PRESENT(CHAMHF) .OR. IFLAG==1 .OR. IFLAG==23) LscaAWARE_LOCAL=.FALSE.

    LASYNC_LOCAL = LASYNC

    CDCHF=0         ! double counting HF
    IF (PRESENT(EXHF_ACFDT)) EXHF_ACFDT=0
    ICALL=ICALL+1

    NB_TOT=WDES%NB_TOT
    NBANDS=WDES%NBANDS
    IF (PRESENT(NBANDS_MAX) .AND. IFLAG == 0) NBANDS=NBANDS_MAX

    NSTRIP=NSTRIP_STANDARD
!$  NSTRIP=NSTRIP*omp_nthreads
 
    NSTRIP=MIN(NSTRIP,NBANDS)

! allocate work space
    IF (.NOT. LscaAWARE_LOCAL) THEN
       ALLOCATE(CHAM(NB_TOT,NB_TOT))
    ELSE
       CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
       ALLOCATE(CHAM(SCALA_NP(),SCALA_NQ()))
    ENDIF
!$ACC ENTER DATA CREATE(CHAM) IF(OFFLOAD_ON)

!$ACC ENTER DATA CREATE(WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(WDES,WDES1,0)

!$ACC ENTER DATA CREATE(WHAM,WNONL) IF(OFFLOAD_ON)
    CALL NEWWAVA(WHAM, WDES1, NSTRIP)
    CALL NEWWAVA_PROJ(WNONL, WDES1)

!$ACC ENTER DATA CREATE(W1) IF(OFFLOAD_ON)
    CALL NEWWAV_R(W1, WDES1)

    IF (IFLAG==5) THEN
       ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL,WOVL) IF(OFFLOAD_ON)
       CALL NEWWAVA_PROJ(WOVL, WDES1)
    ENDIF

    IF (PRESENT(NKSTART)) THEN
       MY_NKSTART=NKSTART
    ELSE
       MY_NKSTART=1
    ENDIF
    IF (PRESENT(NKSTOP)) THEN
       MY_NKSTOP=NKSTOP
    ELSE
       MY_NKSTOP=WDES%NKPTS
    ENDIF

!=======================================================================
! start with HF part and store results WFOCK
!=======================================================================
    IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN

       CALL ALLOCW(WDES,WFOCK)
# 251

# 351

       IF (LPEAD_NO_SCF()) THEN
          IF (PRESENT(EXHF_ACFDT)) THEN
             CALL FOCK_ALL_DBLBUF(WDES,W_STORE,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,EX_ACFDT=EXHF_ACFDT,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK,WP=W)
          ELSE
             CALL FOCK_ALL_DBLBUF(WDES,W_STORE,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK,WP=W)
          ENDIF
       ELSE
          IF (PRESENT(EXHF_ACFDT)) THEN
             CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,EX_ACFDT=EXHF_ACFDT,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK)
          ELSE
             CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK)
          ENDIF
       ENDIF

       IF (LSYMGRAD) &
            CALL APPLY_SMALL_SPACE_GROUP_OP( W, WFOCK, NONLR_S, NONL_S,P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1, MY_NKSTART)
    ENDIF
!=======================================================================
    spin:  DO ISP=1,WDES%ISPIN
    kpoint: DO NK=MY_NKSTART,MY_NKSTOP

       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
       IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))

       CALL SETWDES(WDES,WDES1,NK)
!=======================================================================
!  IFLAG=0 calculate eigenvalues  only
!=======================================================================
       IF (IFLAG==0) THEN
          W%CELEN(:,NK,ISP)=0
          DO N=1,NBANDS
! transform wavefunction to real space
! and calculate eigenvalues calling ECCP, no redistribution !
             CALL SETWAV(W, W1, WDES1, N, ISP) ! allocation for W1%CR 1._q above
             CALL FFTWAV_W1(W1)
!            IF (ASSOCIATED(HAMILTONIAN%AVEC)) THEN
!               CALL ECCP_VEC(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),HAMILTONIAN%AVEC, W1%CELEN)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL ECCP_TAU(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),W1%CELEN)    
             ELSE
                CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W1%CELEN)
             ENDIF
             IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
                W1%CELEN=W1%CELEN+W1_DOT( ELEMENT( W, WDES1, N, ISP), ELEMENT (WFOCK, WDES1, N, ISP))
             ENDIF

             W%CELEN(N,NK,ISP)=W1%CELEN
          ENDDO
          CALL PEAD_EIGENVALUES(W,NK,ISP)
          CYCLE kpoint
       ENDIF

       WA=ELEMENTS(W, WDES1, ISP)
!=======================================================================
!  IFLAG /= 0 calculate Hamiltonian CHAM
!=======================================================================
!  caclulate D |cfin_n> (D = non local strength of PP)
       IF (WDES%DO_REDIS .AND. LASYNC_LOCAL) THEN
          CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW1)
          CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW2)
# 427

          CALL REDIS_PW_STRIP_START(WDES, WA%CPTWFP(1,1), NSTRIP, 1, H_PW1)
# 433

       ENDIF

       CALL OVERL(WDES1,.TRUE.,LMDIM,CDIJ(1,1,1,ISP),WA%CPROJ(1,1),WNONL%CPROJ(1,1))
       DO N=1,NBANDS
          CALL PEAD_ACC_ADD_CPROJ(WNONL%CPROJ(:,N),N,NK,ISP)
       ENDDO

       IF (IFLAG==5) THEN
          CALL OVERL(WDES1,WDES1%LOVERL,LMDIM,CQIJ(1,1,1,ISP),WA%CPROJ(1,1),WOVL%CPROJ(1,1))
       ENDIF

! redistribute the wavefunction characters
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PROJ(WNONL)
       IF (IFLAG==5) CALL REDISTRIBUTE_PROJ(WOVL)
!$ACC KERNELS PRESENT(CHAM) 
       CHAM=0
!$ACC END KERNELS
       strip: DO NPOS=1,NBANDS,NSTRIP
          NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSTRIP)

!  calculate V_{local} |phi> + T | phi >
!  for a block containing NSTRIP wavefunctions

! set Fock contribution
          IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
!$ACC KERNELS PRESENT(WHAM,WFOCK) 
             WHAM%CPTWFP(:,1:NSTRIP_ACT)=WFOCK%CPTWFP(:,NPOS:NPOS+NSTRIP_ACT-1,NK,ISP)
!$ACC END KERNELS
          ENDIF

          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1
             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, HAMILTONIAN%MU, ISP, WHAM%CPTWFP(:,NP), &
               &   USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12) 
             ELSE
                CALL HAMILT_LOCAL(W1, SV, ISP, WHAM%CPTWFP(:,NP), USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12)
             ENDIF
             CALL PEAD_ACC_ADD_PW(WHAM%CPTWFP(:,NP),N,NK,ISP)
             IF (WDES%DO_REDIS.AND. LASYNC_LOCAL) THEN
# 482

                CALL REDIS_PW_START(WDES, WHAM%CPTWFP(1,NP), N, H_PW2)
# 488

             ENDIF
          ENDDO
! redistribute wavefunctions
! after this redistributed up to and including 1...NPOS+NSTRIP_ACT
          IF (WDES%DO_REDIS) THEN
             IF (LASYNC_LOCAL) THEN
                CALL REDIS_PW_STRIP_STOP(WDES, WA%CPTWFP(1,NPOS), NSTRIP_ACT, NPOS, H_PW1)
                CALL REDIS_PW_STRIP_STOP(WDES, WHAM%CPTWFP(1,1), NSTRIP_ACT, NPOS, H_PW2)
# 502

                NSTRIP_NEXT = MIN(NSTRIP, NBANDS-NPOS-NSTRIP_ACT+1)
                IF (NSTRIP_NEXT>0) CALL REDIS_PW_STRIP_START(WDES, WA%CPTWFP(1,NPOS+NSTRIP_ACT), NSTRIP_NEXT, NPOS+NSTRIP_ACT, H_PW1)
# 509

             ELSE
                CALL REDISTRIBUTE_PW( ELEMENTS( WA, NPOS, NPOS-1+NSTRIP_ACT))
                CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
             ENDIF
          ENDIF

          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

          IF (.NOT. LscaAWARE_LOCAL) THEN
             CALL ORTH1('U', &
               WA%CW_RED(1,1),WHAM%CPTWFP(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1))
          ELSE
             CALL ORTH1_DISTRI('U', &
               WA%CW_RED(1,1),WHAM%CPTWFP(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1), & 
               WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
          ENDIF
       ENDDO strip

       IF (WDES%DO_REDIS .AND. LASYNC_LOCAL) THEN
          CALL REDIS_PW_DEALLOC(H_PW1)
          CALL REDIS_PW_DEALLOC(H_PW2)
       ENDIF

       IF (.NOT. LscaAWARE_LOCAL) THEN
          CALL M_sum_z8(WDES%COMM_KIN,CHAM(1,1),SIZE(CHAM,KIND=qi8))
! add lower triangle
!$ACC PARALLEL LOOP GANG PRESENT(CHAM) 
          DO N=1,NB_TOT
!$ACC LOOP VECTOR
             DO NP=N+1,NB_TOT
                CHAM(NP,N)=CONJG(CHAM(N,NP))
             ENDDO
          ENDDO
          IF (IFLAG==23) CALL RESTRICT_TO_OCCUPIED_ONLY(WDES%NB_TOTK(NK,ISP), CHAM, W%FERTOT(:,NK, ISP))
       ENDIF

       IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
!$ACC KERNELS PRESENT(CHAMHF,CHAM) 
          IF (LFIRST) THEN
             CHAMHF(:,:,NK,ISP)=CONJG(CHAMHF(:,:,NK,ISP))-CHAM(:,:)
          ENDIF
          CHAM(:,:)=CHAM(:,:)+CHAMHF(:,:,NK,ISP)
!$ACC END KERNELS
       ENDIF
# 563

!-----------------------------------------------------------------------
! calculate the overlap matrix
!-----------------------------------------------------------------------
       IF (IFLAG==5) THEN
!$ACC KERNELS PRESENT(COVL) 
          COVL=(0._q,0._q)
!$ACC END KERNELS
          DO NPOS=1,NB_TOT-NSTRIP_STANDARD_GLOBAL,NSTRIP_STANDARD_GLOBAL
             CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
                  WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
                  NPOS,NSTRIP_STANDARD_GLOBAL,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
          ENDDO

          CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
               WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
               NPOS,NB_TOT-NPOS+1,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
          CALL M_sum_z8(WDES%COMM_KIN,COVL(1,1),SIZE(COVL,KIND=qi8))
# 585

       ENDIF
!=======================================================================
! IFLAG =2
!  simply copy eigenvalues
!=======================================================================
       IF (IFLAG==12) THEN
!$ACC UPDATE SELF(CHAM) 
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
          IF (IU0>=0) CALL DUMP_HAM( "Hamilton matrix",WDES, CHAM)
          CALL vtutor%stopCode()
       ENDIF

       IF (.NOT. LscaAWARE_LOCAL) THEN
!$ACC PARALLEL LOOP PRESENT(WDES,W,CHAM) 
          DO N=1,WDES%NB_TOTK(NK,ISP)
             W%CELTOT(N,NK,ISP)=CHAM(N,N)
          ENDDO
!$ACC UPDATE SELF(W%CELTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)) 
       ENDIF
!=======================================================================
! IFLAG =4 use Loewdin perturbation to get rotation matrix
! this preserves the ordering of the eigenvalues
! MIND: does not work for real matrices
!=======================================================================
       IF (IFLAG==4) THEN
          CALL LOEWDIN_DIAG(WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM)

          CALL ORSP(WDES%NB_TOTK(NK,ISP), NB_TOT, NB_TOT, CHAM)
!test writegamma
!          IF (IU0>=0) CALL WRITEGAMMA(NK, WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM, .TRUE.)
# 618

       ELSE
!=======================================================================
! IFLAG > 1 and IFLAG <4
! diagonalization of CHAM
! we have lots of choices for the parallel version
! this  makes things rather complicated
! to allow for reasonable simple programming, once the diagonalisation
! has been 1._q I jump to line 100
!=======================================================================
          IF (IFLAG==1) GOTO 1000


          IF (.NOT. LscaAWARE_LOCAL) THEN
 !$OMP PARALLEL DO DEFAULT(SHARED) &
 !$OMP PRIVATE(N)
!$ACC KERNELS PRESENT(WDES,CHAM) 
             DO N=1,WDES%NB_TOTK(NK,ISP)
                IF (ABS(AIMAG(CHAM(N,N)))>1E-2_q .AND. IU0>=0) THEN
 !$OMP CRITICAL (omp_wrt_stdout)
              WRITE(IU0,'(A,I5,E14.3)')'EDDIAG: WARNING: subspace matrix is not hermitian',N,AIMAG(CHAM(N,N))
 !$OMP END CRITICAL (omp_wrt_stdout)
                ENDIF
                CHAM(N,N)= REAL( CHAM(N,N) ,KIND=q)
             ENDDO
!$ACC END KERNELS
 !$OMP END PARALLEL DO
          ELSE
             CALL BG_CHANGE_DIAGONALE(WDES%NB_TOTK(NK,ISP),CHAM(1,1),IU0)
          ENDIF


!
! parallel versions
! if fast Jacobi method exists use it (T3D, T3E only)
! use the first line in that case
          CALL SHIFT_BANDS_BETWEEN_LOW_HIGH( WDES, WDES%NB_TOTK(NK,ISP), LscaAWARE_LOCAL, CHAM )

          IFAIL=0


          IF ( LJACOBI .AND. IFLAG /=13 ) THEN
             IF (IU0>=0) WRITE(IU0,*)'jacobi called'
             CALL jacDSSYEV(WDES%COMM_KIN, CHAM(1,1), R, NB_TOT)
             CALL M_sum_z8(WDES%COMM_KIN, CHAM(1,1), SIZE(CHAM,KIND=qi8))
             CALL M_sum_z(WDES%COMM_KIN, R , NB_TOT)

             GOTO 100
          ENDIF


! when IFLAG/=5 we may use ScaLAPACK
          IF (IFLAG/=5) THEN
! Is the Hamiltonian distributed block-cyclically?
! If so we will definitely use ScaLAPACK
             IF (LscaAWARE_LOCAL) THEN

              CALL BG_pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  WDES%NB_TOTK(NK,ISP))
!  alternatively use divide and conquer, somewhat slower though
!CALL BG_pDSYEV_ZHEEVD(WDES%COMM_KIN, CHAM(1,1), R,  WDES%NB_TOTK(NK,ISP))

                
                GOTO 100

! If the Hamiltonian is not block-cyclically distributed
! we may still use ScaLAPACK ... (but not in the OpenACC version)
             ELSEIF(LscaLAPACK) THEN

                CALL pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  NB_TOT, WDES%NB_TOTK(NK,ISP))
!  alternatively use divide and conquer, somewhat slower though
!                CALL pDSYEV_ZHEEVD(WDES%COMM_KIN, CHAM(1,1), R,  NB_TOT, WDES%NB_TOTK(NK,ISP))
                CALL M_sum_z8(WDES%COMM_KIN, CHAM(1,1), SIZE(CHAM,KIND=qi8))

                
                GOTO 100

             ENDIF
          ENDIF

!
!  serial codes
!
# 725

          IF (IFLAG == 5) THEN
             CALL ZHEGV &
                  (1,'V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,COVL(1,1),NB_TOT, &
                  R,CWRK,LWORK*NB_TOT, RWORK, IFAIL)
          ELSE
# 735

             ABSTOL=1E-10_q
             VL=0 ; VU=0 ; IL=0 ; IU=0
             ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL) 
             CALL ZHEEVX &
                  ( 'V', 'A', 'U', WDES%NB_TOTK(NK,ISP), CHAM(1,1) , NB_TOT, VL, VU, IL, IU, &
                  ABSTOL , NB_CALC , R, COVL(1,1), NB_TOT, CWRK, &
                  LWORK*NB_TOT, RWORK, IWORK, INFO, IFAIL)
             CALL ZCOPY(NB_TOT*NB_TOT,COVL,1,CHAM,1)
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(COVL) IF(OFFLOAD_ON)
             DEALLOCATE(COVL)

          ENDIF


100       CONTINUE
          

          IF (IFAIL/=0) THEN
             CALL vtutor%error("ERROR in EDDIAG: call to ZHEEV/ZHEEVX/DSYEV/DSYEVX failed! error code &
                &was " // str(IFAIL))
          ENDIF

! shift eigenvalues back
          IF (IFLAG==23) CALL SHIFT_UNOCCUPIED_BACK(WDES%NB_TOTK(NK,ISP), R, W%FERTOT(:,NK, ISP))
          CALL SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH(WDES, WDES%NB_TOTK(NK,ISP), R )

          DO N=1,WDES%NB_TOTK(NK,ISP)
             W%CELTOT(N,NK,ISP)=R(N)
          ENDDO
       ENDIF
!=======================================================================
! IFLAG > 2
! rotate wavefunctions
!=======================================================================
       IF (IFLAG==2) GOTO 1000

       IF (.NOT. LscaAWARE_LOCAL) THEN
          CALL LINCOM('F',WA%CW_RED(:,:),WA%CPROJ_RED(:,:),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
            WA%CW_RED(:,:),WA%CPROJ_RED(:,:))
       ELSE
          CALL LINCOM_DISTRI('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
            WDES%COMM_KIN, NBLK )
       ENDIF

1000   CONTINUE
       
!  back redistribution over bands
       IF (WDES%DO_REDIS) THEN
          CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
          IF (LASYNC_LOCAL) THEN
             W%OVER_BAND=.TRUE.
          ELSE
             CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))
          ENDIF
          ! "redis ok"
       ENDIF

! return Hamiltonian (IFLAG==1) or rotation matrix (IFLAG=3,4)
       IF (PRESENT(LLAST)) THEN       
          IF (LLAST) THEN
!$ACC KERNELS PRESENT(CHAMHF,CHAM) 
             IF (IFLAG==1) CHAM(:,:)=CONJG(CHAM(:,:))
             CHAMHF(:,:,NK,ISP)=CHAM(:,:)
!$ACC END KERNELS
          ENDIF
       ENDIF

    ENDDO kpoint
    ENDDO spin

!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

!    CALL M_sum_z(WDES%COMM_KIN,CDCHF,1)
!    CALL M_sum_z(WDES%COMM_KINTER,CDCHF,1)
!    EXHF=CDCHF
!
!    IF (PRESENT(EXHF_ACFDT)) THEN
!       CALL M_sum_d(WDES%COMM_KIN,EXHF_ACFDT,1)
!       CALL M_sum_d(WDES%COMM_KINTER,EXHF_ACFDT,1)
!    ENDIF

! need to correct CELTOT at this point so that it is correct on all nodes
    IF (IFLAG==0) THEN
        CALL MRG_CEL(WDES,W)
    ENDIF

    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       CALL KPAR_SYNC_CELTOT(WDES,W)
    ENDIF

    IF (IFLAG>=3.AND.(LHFCALC.OR.LUSEPEAD()).AND.WDES%COMM_KINTER%NCPU.GT.1) THEN
       CALL KPAR_SYNC_FERTOT(WDES,W)
       CALL KPAR_SYNC_CELTOT(WDES,W)
       CALL KPAR_SYNC_AUXTOT(WDES,W)
       CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
    ENDIF

!$ACC UPDATE DEVICE(W%CELTOT) 

! deallocation ...
!$ACC WAIT IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CHAM) 
    DEALLOCATE(CHAM)

    CALL DELWAV_R(W1)
!$ACC EXIT DATA DELETE(W1) 

    CALL DELWAVA(WHAM)
    CALL DELWAVA_PROJ(WNONL)
!$ACC EXIT DATA DELETE(WHAM,WNONL) 

    IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN 
# 856

       CALL DEALLOCW(WFOCK)
# 864

    ENDIF

    IF (IFLAG==5) THEN
!$ACC EXIT DATA DELETE(COVL) 
       DEALLOCATE(COVL)
       CALL DELWAVA_PROJ(WOVL)
!$ACC EXIT DATA DELETE(WOVL) 
    ENDIF

    IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
       LFIRST=.FALSE.
    ENDIF

# 880


# 899


    

    RETURN
    
    CONTAINS

!***********************************************************************
!
!> The following subroutines are used to constrain the
!> diagonalization to the occupied manifold
!> here orbitals with an occupancy close to 1
!
!***********************************************************************

    SUBROUTINE RESTRICT_TO_OCCUPIED_ONLY(NB_TOT, CHAM, FERTOT)

!! USE moffload_struct_def

      INTEGER :: NB_TOT
      COMPLEX(q)    :: CHAM(:,:)
      REAL(q) :: FERTOT(:)
      INTEGER :: NB_OCC, I, J

! seek first occupancy that differs from 1.00
      DO NB_OCC=1, NB_TOT
         IF (ABS((FERTOT(NB_OCC))-1.0_q)>1E-5_q) EXIT
      ENDDO

! now set the lower (upper) triangle starting with row NB_OCC to (0._q,0._q)
! loop over row index
!$ACC PARALLEL LOOP GANG PRESENT(CHAM) 
      DO I=NB_OCC, NB_TOT
! loop over column index
!$ACC LOOP VECTOR
         DO J=1,I-1
            CHAM(I,J)=0
            CHAM(J,I)=0
         ENDDO
! shift diagonal elements by 10 eV
! so that the diagonalization does not mix occupied and unoccupied manyfold
! this shift needs to be removed by  SHIFT_UNOCCUPIED_BACK
         CHAM(I,I)=CHAM(I,I)+10.0_q
      ENDDO

    END SUBROUTINE RESTRICT_TO_OCCUPIED_ONLY

    SUBROUTINE SHIFT_UNOCCUPIED_BACK(NB_TOT, R, FERTOT)
      INTEGER :: NB_TOT
      REAL(q) :: R(:)
      REAL(q) :: FERTOT(:)
      INTEGER :: NB_OCC, I

! seek first occupancy that differs from 1.00
      DO NB_OCC=1, NB_TOT
         IF (ABS((FERTOT(NB_OCC)-1.0_q))>1E-5_q) EXIT
      ENDDO

! now set the lower (upper) triangle starting with row NB_OCC to (0._q,0._q)
! loop over row index
      DO I=NB_OCC, NB_TOT
         R(I)=R(I)-10.0_q
      ENDDO
    END SUBROUTINE SHIFT_UNOCCUPIED_BACK

  END SUBROUTINE EDDIAG

!***********************************************************************
!
!>  for a selected k-points and spin component
!>  calculate onsite ((1._q,0._q) centre) contribution to the Hamilton matrix
!>
!>  if LOVERL is .TRUE.
!>    ~~~
!>    CCORR(m,n) += <psi_k,n| beta_i> D_ij <beta_j | psi_k,m>
!>    ~~~
!>
!>  if the local potential SV is passed  down the Hamilton matrix
!>
!>   ~~~
!>    CCORR(m,n) += <psi_k,n|  SV | psi_k,m>
!>   ~~~
!>  is added as well
!>  the calling routine must initialize CCORR to (0._q,0._q)
!
!***********************************************************************

  SUBROUTINE ONE_CENTER_BETWEEN_STATES(HAMILTONIAN, LATT_CUR, LOVERL, WDES, W, NK, ISP, LMDIM, &
       CDIJC, CCORR, SV)
    USE prec
    USE wave_high
    USE dfast
    USE lattice
    USE hamil_struct_def

    TYPE (ham_handle)  HAMILTONIAN
    TYPE (latt)        LATT_CUR
    LOGICAL LOVERL
    INTEGER NK, ISP, LMDIM
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
!> (1._q,0._q) centre correction
    COMPLEX(q) CDIJC(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)  
    COMPLEX(q) ::  CCORR(WDES%NB_TOT,WDES%NB_TOT)
!> local potential
    COMPLEX(q), OPTIONAL  :: SV(WDES%GRID%MPLWV,WDES%NCDIJ)
! local
    COMPLEX(q)   , POINTER :: CPROJ_RED(:,:)
    INTEGER NPOS
    INTEGER NSTRIP, NSTRIP_ACT, NPOS_RED, NSTRIP_RED
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WOVL           ! array to hold non local part
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1             ! current wavefunction
    TYPE (wavefuna)    WHAM           ! store Hamiltonian times wavefunction
    TYPE (wavefuna)    WA             ! array to store wavefunction

    NSTRIP=NSTRIP_STANDARD
!-----------------------------------------------------------------------
! non local part
!-----------------------------------------------------------------------
    IF (LOVERL) THEN
       CALL SETWDES(WDES,WDES1,NK)

       CALL NEWWAVA_PROJ(WNONL, WDES1)
       IF (WDES%DO_REDIS) THEN
          CALL NEWWAVA_PROJ(WOVL, WDES1)
          CALL WA_COPY_CPROJ(ELEMENTS(W, WDES1, ISP), WOVL)
       ELSE
          WOVL=ELEMENTS(W, WDES1, ISP)
       ENDIF
       
       CALL OVERL(WDES1, .TRUE., LMDIM, CDIJC(1,1,1,ISP), WOVL%CPROJ(1,1), WNONL%CPROJ(1,1))
       
       IF (WDES%DO_REDIS) THEN
          CALL REDIS_PROJ(WDES1, WDES%NBANDS, WNONL%CPROJ(1,1))
          CALL REDIS_PROJ(WDES1, WDES%NBANDS, WOVL%CPROJ (1,1))
       ENDIF
       
       DO NPOS=1,WDES%NBANDS,NSTRIP
          NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)
          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR
          
          CALL ORTH1('U', &
               W%CPTWFP(1,1,NK,ISP),W%CPTWFP(1,1,NK,ISP),WOVL%CPROJ(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, 0,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CCORR(1,1))
!             attention              -

          
       ENDDO
       CALL DELWAVA_PROJ(WNONL)
       IF (WDES%DO_REDIS) CALL DELWAVA_PROJ(WOVL)
    ENDIF
!-----------------------------------------------------------------------
! local part
!-----------------------------------------------------------------------
    IF (PRESENT(SV)) THEN
! allocate work space
       ALLOCATE(W1%CR(WDES%GRID%MPLWV*WDES%NRSPINORS))

       CALL SETWDES(WDES,WDES1,NK)
       CALL NEWWAVA(WHAM, WDES1, NSTRIP)

       WA=ELEMENTS(W, WDES1, ISP)

       strip: DO NPOS=1,WDES%NBANDS,NSTRIP
          NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)

!  calculate V_{local} |phi> + T | phi >
!  for a block containing NSTRIP wavefunctions
          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1

             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, HAMILTONIAN%MU, ISP,  WHAM%CPTWFP(:,NP), .FALSE., .FALSE.)
             ELSE
                CALL HAMILT_LOCAL(W1, SV, ISP,  WHAM%CPTWFP(:,NP), .FALSE., .FALSE.)
             ENDIF
          ENDDO
! redistribute wavefunctions
! after this redistributed up to and including 1...NPOS+NSTRIP_ACT
          IF (WDES%DO_REDIS) THEN
             CALL REDISTRIBUTE_PW( ELEMENTS( WA, NPOS, NPOS-1+NSTRIP_ACT))
             CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
          ENDIF

          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

          CALL ORTH1('U', &
               WA%CW_RED(1,1),WHAM%CPTWFP(1,1),WA%CPROJ_RED(1,1), &
               WA%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CCORR(1,1))
!           attention                   ---

       ENDDO strip

       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

! deallocation ...
       DEALLOCATE(W1%CR)
       CALL DELWAVA(WHAM)

    ENDIF

    CALL M_sum_z8(WDES%COMM_KIN,CCORR(1,1),SIZE(CCORR,KIND=qi8))

  END SUBROUTINE ONE_CENTER_BETWEEN_STATES


!***********************************************************************
!
!> read in the density matrix or Hamiltonian from the GAMMA or vaspgamma.h5
!> file to the diagonal density matrix F ((1._q,0._q)-electron occupancies) or
!> Hamiltonian ((1._q,0._q)-electron eigenvalues)
!> then diagonalize the resulting matrix and rotate the
!> (1._q,0._q) electron orbitals accordingly
!
!> @details @ref openacc :
!> On entering this subroutine W%CPTWFP and W%CPROJ are updated on the host
!> (if present on the device) and then OpenACC execution is switched off.
!> On exit the OpenACC execution mode reverts to its original status,
!> and when moffload_struct_def::offload_on = .TRUE. the arrays
!> W%CPTWFP and W%CPROJ are updated on the device (if present).
!
!***********************************************************************
  SUBROUTINE ADD_GAMMA_FROM_FILE( WDES, W, KPOINTS, NELECT, NUP_DOWN, LWRITE, IO )

!! USE moffload

    USE prec
    USE wave_high
    USE mkpoints
    USE dfast
    USE base
    USE fileio

    USE vhdf5_base
    USE vhdf5
    USE string

    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
!> number of electrons
    REAL (q)  NELECT
!> number total spin
    REAL (q)  NUP_DOWN
!> k-points structure
    TYPE (kpoints_struct) KPOINTS
!> write final eigenvalues to OUTCAR
    LOGICAL :: LWRITE
    TYPE (in_struct)   IO
! local
    INTEGER NK, ISP, NB, LMDIM
    LOGICAL FILE_PRESENT
    TYPE (wavedes1)    WDES1
    TYPE (wavefuna)    WA             ! array to store wavefunction
    COMPLEX(q), ALLOCATABLE :: CHAM(:,:,:,:)
    REAL(q), ALLOCATABLE :: CELTOT(:,:,:)
    INTEGER :: IMODE                  ! add to Hamiltonian (IMODE=2) or density matrix (IMODE=1)
! LAPACK
    REAL(q)    R(WDES%NB_TOT)
    INTEGER :: IFAIL
    INTEGER, PARAMETER :: LWORK=32
    COMPLEX(q)       CWRK(LWORK*WDES%NB_TOT)
    REAL(q)    RWORK(3*WDES%NB_TOT)

    INTEGER N1, N2, dim_c, bnd_c_max, bnd_c_min
    COMPLEX(q), ALLOCATABLE :: deltaN_k(:,:)
    REAL(q), ALLOCATABLE :: triqs_band_window(:,:)
    CHARACTER(len=:), allocatable :: gamma_h5
    INTEGER(HID_T):: th5filouteid
    INTEGER(HID_T):: th5gammagid


    

!$ACC UPDATE SELF(W%CPTWFP,W%CPROJ) IF_PRESENT IF(OFFLOAD_ON)
!! 
   
! read CHAM from file
    ALLOCATE(CHAM(WDES%NB_TOT, WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN))
    CHAM=0
    
! if the text file GAMMA exists it takes priority
    CALL FILE_EXISTS(WDES%COMM, 'GAMMA', FILE_PRESENT)
    IF (FILE_PRESENT) THEN 
      CALL OPENGAMMA
      CALL READGAMMA_HEAD( IMODE, WDES%NKPTS , WDES%NB_TOT, IO)

! for IMODE = 2 we save the CELTOT locally
      IF (IMODE==2) THEN
         ALLOCATE(CELTOT(SIZE(W%CELTOT,1),SIZE(W%CELTOT,2),SIZE(W%CELTOT,3)))
         CELTOT=W%CELTOT
      ENDIF
      DO ISP=1,WDES%ISPIN
        DO NK=1,WDES%NKPTS
          CALL READGAMMA(NK, WDES%NB_TOTK(NK,ISP), SIZE(CHAM,1), CHAM(:,:,NK,ISP),  IO)

          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

! add (1._q,0._q)-particle occupancies
          IF (IMODE==1) THEN
            DO NB=1, WDES%NB_TOTK(NK,ISP)
                CHAM(NB,NB,NK,ISP)=CHAM(NB,NB,NK,ISP)+W%FERTOT(NB,NK,ISP)
            ENDDO
! add (1._q,0._q)-particle eigenvalues
          ELSE
            DO NB=1, WDES%NB_TOTK(NK,ISP)
                CHAM(NB,NB,NK,ISP)=CHAM(NB,NB,NK,ISP)+W%CELTOT(NB,NK,ISP)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      
! change sign of density matrix to sort occupied states as lowest states
      IF (IMODE==1) CHAM=-CHAM
      CALL CLOSEGAMMA
    ELSE 
      CALL FILE_EXISTS(WDES%COMM, 'vaspgamma.h5', FILE_PRESENT)
      IF (.NOT.(FILE_PRESENT)) THEN 
          CALL vtutor%error("ERROR in ADD_GAMMA_FROM_FILE: no GAMMA or vaspgamma.h5 file found")
      ENDIF

      gamma_h5= 'vaspgamma.h5'
      ALLOCATE(triqs_band_window(2, WDES%NKPTS))
! open h5 file
      call vh5_error(vh5_file_open_read(gamma_h5, th5filouteid),"subrot.F",1232)
      IMODE = 1 ! only density matrix is stored in the hdf5 file. No option for Hamiltonian

      DO ISP=1,WDES%ISPIN
! load deltaN / GAMMA matrix and band_window
        triqs_band_window = 0
        IF (ISP == 1) THEN
          call vh5_error(vh5_read(th5filouteid, "band_window/0", triqs_band_window(:,:)),"subrot.F",1239)
          call vh5_error(vh5_group_open(th5filouteid, "deltaN/up", th5gammagid),"subrot.F",1240)
        ELSE
          call vh5_error(vh5_read(th5filouteid, "band_window/1", triqs_band_window(:,:)),"subrot.F",1242)
          call vh5_error(vh5_group_open(th5filouteid, "deltaN/down", th5gammagid),"subrot.F",1243)
        ENDIF
        DO NK=1,WDES%NKPTS
! read CHAM from h5
          bnd_c_min = INT(triqs_band_window(1, NK))
          bnd_c_max = INT(triqs_band_window(2, NK))
          dim_c = bnd_c_max - bnd_c_min + 1
! deltaN can have different dim each k point
          ALLOCATE(deltaN_k(dim_c, dim_c))
          deltaN_k = 0
! load deltaN at NK from h5
          call vh5_error(vh5_read(th5gammagid, trim(str(NK-1)), deltaN_k),"subrot.F",1254)
! write into CHAM
          DO N1= bnd_c_min, bnd_c_max
            DO N2= bnd_c_min, bnd_c_max
! switch indices due to different indexing in C and Fortran
                CHAM(N2,N1,NK,ISP) = deltaN_k(N1-bnd_c_min+1,N2-bnd_c_min+1)
            ENDDO
          ENDDO
          DEALLOCATE(deltaN_k)

          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

! add (1._q,0._q)-particle occupancies
          DO NB=1, WDES%NB_TOTK(NK,ISP)
            CHAM(NB,NB,NK,ISP)=CHAM(NB,NB,NK,ISP)+W%FERTOT(NB,NK,ISP)
          ENDDO
        ENDDO
      ENDDO
! change sign of density matrix to sort occupied states as lowest states
      CHAM=-CHAM
! close groups and file
    call vh5_error(vh5_group_close(th5gammagid),"subrot.F",1275)
    call vh5_error(vh5_file_close(th5filouteid),"subrot.F",1276)
# 1279

    ENDIF

    DO ISP=1,WDES%ISPIN
      DO NK=1,WDES%NKPTS
# 1288

       CALL ZHEEV &
            ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1,NK,ISP),WDES%NB_TOT, &
            R,CWRK,LWORK*WDES%NB_TOT, RWORK,  IFAIL)

# 1295

! change sign of eigenvalues
       IF (IMODE==1) R=-R

       IF (IFAIL/=0) THEN
          CALL vtutor%error("ERROR in ADD_GAMMA_FROM_FILE: call to ZHEEV/ DSYEV failed! error code was &
             &" // str(IFAIL))
       ENDIF

       CALL SETWDES(WDES,WDES1,NK)

       WA=ELEMENTS(W, WDES1, ISP)

!  distribution over plane wave coefficients
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

       CALL LINCOM('F',WA%CW_RED(:,:),WA%CPROJ_RED(:,:),CHAM(1,1,NK,ISP), &
            WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,WDES%NB_TOT, &
            WA%CW_RED(:,:),WA%CPROJ_RED(:,:))

!  back redistribution over bands
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

! updated (1._q,0._q)-electron occupancies
       IF (IMODE==1) THEN
          W%FERTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)=R(1:WDES%NB_TOTK(NK,ISP))
       ELSE
! copy R to CELTOT
          W%CELTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)=R(1:WDES%NB_TOTK(NK,ISP))
       ENDIF
      ENDDO
    ENDDO

! just in case sync everything back to all nodes if KPAR is used
    CALL KPAR_SYNC_ALL(WDES,W)

! IMODE==2, update (1._q,0._q)-electron occupancies now
    IF (IMODE==2) THEN
! update fermi-weights
# 1339

       CALL DENSTA_SIMPLE(W, KPOINTS, NELECT, NUP_DOWN )
       IF (LWRITE) THEN
          IF (IO%IU6>=0) THEN
             WRITE(IO%IU6,*) 'eigenvalues after inclusion of GAMMA file'
          ENDIF
          CALL WRITE_EIGENVAL( W%WDES, W, IO%IU6)
       ENDIF
! restore old (1._q,0._q)-electron eigenvalues
       W%CELTOT=CELTOT
    ENDIF

    DEALLOCATE(CHAM)

!! 
!$ACC UPDATE DEVICE(W%CPTWFP,W%CPROJ) IF_PRESENT IF(OFFLOAD_ON)

  

  END SUBROUTINE ADD_GAMMA_FROM_FILE


!************************ SUBROUTINE EDDIAG_EXACT **********************
!
!> this subroutine performs a full diagonalization of the Hamiltonian
!
!> right now it is stupidly implemented since
!> the number of bands is increased for all k-points
!> then EDDIAG is called
!> and finally the number of bands is set back to the original
!> value
!> this requires a lot of storage but is still convenient for
!> GW and RPA calculations
!
!> @details @ref openacc :
!> On entering this subroutine W%CPTWFP and W%CPROJ are updated on the host
!> (if present on the device) and then OpenACC execution is switched off.
!> On exit the OpenACC execution mode reverts to its original status,
!> and when moffload_struct_def::offload_on = .TRUE. the arrays
!> W%CPTWFP and W%CPROJ are updated on the device (if present).
!
!***********************************************************************

  SUBROUTINE EDDIAG_EXACT(HAMILTONIAN, &
       GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
       LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IU0,IU6,EXHF,EXHF_ACFDT,NKSTART,NKSTOP)

!! USE moffload

    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_struct_def
    USE main_mpi
    USE pseudo
    USE poscar
    USE ini
    USE choleski
    USE fock
    USE scala
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry) ::   SYMM      
    INTEGER LMDIM
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!> determines mode of diagonalisation
    INTEGER            IFLAG 
!> local potential
    COMPLEX(q)   SV(GRID%MPLWV,WDES%NCDIJ)
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    INTEGER IU0, IU6
    REAL(q) EXHF
    REAL(q) EXHF_ACFDT
    INTEGER, OPTIONAL :: NKSTART, NKSTOP     ! start k-point
! local
    INTEGER NB_TOT    ! maximum number of plane wave coefficients = number of bands
    TYPE (wavedes)     WDES_TMP
    TYPE (wavespin)    W_TMP
    INTEGER NK, DEGREES_OF_FREEDOM
# 1430


    

! just make sure that data distribution is over bands
    CALL REDIS_PW_OVER_BANDS(WDES, W)
! are all bands calculated anyway
    DEGREES_OF_FREEDOM=MAXVAL(WDES%NPLWKP_TOT)
    IF (WDES%LGAMMA) THEN
       DEGREES_OF_FREEDOM=DEGREES_OF_FREEDOM*2-1
    ENDIF

    IF (DEGREES_OF_FREEDOM<=WDES%NB_TOT) THEN
       call vtutor%warning('The number of requested bands '//str(WDES%NB_TOT)//&
                           ' is larger than the basis size '//str(DEGREES_OF_FREEDOM)//&
                          '. You might want to increase ENCUT.')
       IFLAG=3
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IU0,EXHF,EXHF_ACFDT=EXHF_ACFDT,&
            NKSTART=NKSTART,NKSTOP=NKSTOP)
    ELSE

       NB_TOT=((DEGREES_OF_FREEDOM+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    
       WDES_TMP=WDES
       WDES_TMP%NB_TOT=NB_TOT
       WDES_TMP%NBANDS=NB_TOT/WDES%NB_PAR
       CALL INIT_SCALAAWARE( WDES_TMP%NB_TOT, WDES_TMP%NRPLWV, WDES_TMP%COMM_KIN )
       
       NULLIFY(WDES_TMP%NB_TOTK)
       ALLOCATE(WDES_TMP%NB_TOTK(WDES%NKDIM,2))
! set the maximum number of bands k-point dependent
       DO NK=1,WDES_TMP%NKPTS
          IF (WDES_TMP%LGAMMA) THEN
             WDES_TMP%NB_TOTK(NK,:)=MIN(WDES_TMP%NB_TOT,WDES_TMP%NPLWKP_TOT(NK)*2-1)
          ELSE
             WDES_TMP%NB_TOTK(NK,:)=MIN(WDES_TMP%NB_TOT,WDES_TMP%NPLWKP_TOT(NK))
          ENDIF
       ENDDO
       CALL RESETUP_FOCK_WDES(WDES_TMP, LATT_CUR, LATT_CUR, -1)

       CALL ALLOCW(WDES_TMP,W_TMP)

# 1487


       CALL DUMP_ALLOCATE(IU6)

       W_TMP%FERTOT=0
       W_TMP%CELTOT=0
! copy data back to work array
!$ACC KERNELS PRESENT(W,W_TMP) 
       W_TMP%CPTWFP(:,1:WDES%NBANDS,:,:)    =W%CPTWFP(:,1:WDES%NBANDS,:,:)
       W_TMP%CPROJ(:,1:WDES%NBANDS,:,:) =W%CPROJ(:,1:WDES%NBANDS,:,:)
!$ACC END KERNELS
       W_TMP%CELTOT(1:WDES%NB_TOT,:,:)=W%CELTOT(1:WDES%NB_TOT,:,:)
       W_TMP%FERTOT(1:WDES%NB_TOT,:,:)=W%FERTOT(1:WDES%NB_TOT,:,:)

! random initialization beyond WDES%NBANDS
       CALL WFINIT(WDES_TMP, W_TMP, 1E10_q, WDES%NB_TOT+1) ! ENINI=1E10 not cutoff restriction
!$ACC UPDATE DEVICE(W_TMP%CPTWFP(:,WDES%NBANDS+1:,:,:),W_TMP%CPROJ(:,WDES%NBANDS+1:,:,:)) 

! get characters
       CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP)
! orthogonalization
       CALL ORTHCH(WDES_TMP,W_TMP, WDES%LOVERL, LMDIM,CQIJ, NKSTART=NKSTART)
! and diagonalization
       IFLAG=3
!note(sm): Updating FERTOT because FOCK_ALL_DBLBUF needs FERWE in sync on GPU and CPYU
!$ACC UPDATE DEVICE(W_TMP%FERTOT) IF(OFFLOAD_ON .AND. USEFOCK_CONTRIBUTION()) ASYNC(ACC_ASYNC_Q)
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP,WDES_TMP,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IU0,EXHF,EXHF_ACFDT=EXHF_ACFDT, &
            NKSTART=NKSTART,NKSTOP=NKSTOP)
! copy data back to original array

!$ACC KERNELS PRESENT(W,W_TMP) 
       W%CPTWFP(:,1:WDES%NBANDS,:,:)    =W_TMP%CPTWFP(:,1:WDES%NBANDS,:,:)
       W%CPROJ(:,1:WDES%NBANDS,:,:) =W_TMP%CPROJ(:,1:WDES%NBANDS,:,:)
!$ACC END KERNELS
       W%CELTOT(1:WDES%NB_TOT,:,:)=W_TMP%CELTOT(1:WDES%NB_TOT,:,:)
       W%FERTOT(1:WDES%NB_TOT,:,:)=W_TMP%FERTOT(1:WDES%NB_TOT,:,:)

# 1528

       CALL DEALLOCW(W_TMP)
       DEALLOCATE(WDES_TMP%NB_TOTK)

       CALL RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_CUR, -1)

    ENDIF

    

  END SUBROUTINE EDDIAG_EXACT


!***********************************************************************
!
!> The following subroutine is used to constrain the diagonalization
!> to the set of orbitals between WDES\%NBANDSLOW and WDES%NBANDSHIGH
!>
!> see also #SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH
!
!***********************************************************************

    SUBROUTINE SHIFT_BANDS_BETWEEN_LOW_HIGH(WDES, NB_TOT, LscaAWARE_LOCAL, CHAM )
      USE prec
      USE wave_high
      USE scala

      TYPE (wavedes)     WDES
      INTEGER :: NB_TOT
!> use 1
      LOGICAL :: LscaAWARE_LOCAL  ! use 1
!> Hamiltonian (either distributed or global array)
      COMPLEX(q)    :: CHAM(:,:)
      
      INTEGER :: I
      REAL(q) :: RSHIFT
      

! electron volts by which states below WDES%NBANDSLOW-1 and above WDES%NBANDSHIGH+1 should be shifted
      RSHIFT = 2.0_q * WDES%ENMAX

! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSLOW > 0) THEN
         IF (LscaAWARE_LOCAL) THEN

            CALL SCA_SHIFT_BANDS_LOW_NOINT( NB_TOT, WDES%NBANDSLOW, CHAM(1,1),-RSHIFT, DESCSTD)

         ELSE
! set off diagonal elements below NBANDSLOW to 0
!$ACC KERNELS PRESENT(WDES,CHAM) 
            DO I=1, MIN(WDES%NBANDSLOW-1,NB_TOT)
               CHAM(I,I+1:)=0
               CHAM(I+1:,I)=0
! shift diagonal elements by - RSHIFT eV
! so that the diagonalization does not mix occupied and unoccupied manyfold
! this shift needs to be removed by  SHIFT_UNOCCUPIED_BACK
               CHAM(I,I)=CHAM(I,I) - RSHIFT
            ENDDO
!$ACC END KERNELS
         ENDIF
      ENDIF

! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSHIGH > 0) THEN
         IF (LscaAWARE_LOCAL) THEN
            WRITE(*,*) 'shift high bands lapack',NB_TOT

            CALL SCA_SHIFT_BANDS_HIGH_NOINT( NB_TOT, WDES%NBANDSHIGH, CHAM(1,1), RSHIFT, DESCSTD)

         ELSE
! set off diagonal elements below NBANDSLOW to 0
!$ACC KERNELS PRESENT(WDES,CHAM) 
            DO I=WDES%NBANDSHIGH+1, NB_TOT
               CHAM(I,:I-1)=0
               CHAM(:I-1,I)=0
! shift diagonal elements by +100 eV
! so that the diagonalization does not mix occupied and unoccupied manyfold
! this shift needs to be removed by  SHIFT_UNOCCUPIED_BACK
               CHAM(I,I)=CHAM(I,I)+RSHIFT
            ENDDO
!$ACC END KERNELS
         ENDIF
      ENDIF

      
    END SUBROUTINE SHIFT_BANDS_BETWEEN_LOW_HIGH


!***********************************************************************
!
!> The following subroutine is used to constrain the diagonalization
!> to the set of orbitals between WDES%NBANDSLOW and WDES%NBANDSHIGH
!> see also SHIFT_BANDS_BETWEEN_LOW_HIGH
!
!***********************************************************************

    SUBROUTINE SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH(WDES, NB_TOT, R )
      USE prec
      USE wave_high

      TYPE (wavedes)     WDES
      INTEGER :: NB_TOT
!> eigenvalues
      REAL(q) :: R(:)
      
      INTEGER :: I
      REAL(q) :: RSHIFT

! electron volts by which states below WDES%NBANDSLOW-1 and above WDES%NBANDSHIGH+1 should be shifted
      RSHIFT = 2.0_q * WDES%ENMAX

! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSLOW > 0) THEN
! set off diagonal elements below NBANDSLOW to 0
         DO I=1, MIN(WDES%NBANDSLOW-1, NB_TOT)
            R(I)=R(I)+RSHIFT
         ENDDO
      ENDIF

! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSHIGH > 0) THEN
! set off diagonal elements below NBANDSLOW to 0
         DO I=WDES%NBANDSHIGH+1, NB_TOT
            R(I)=R(I)-RSHIFT
         ENDDO
      ENDIF

    END SUBROUTINE SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH

END MODULE subrot




!***********************************************************************
!
!> dump a "Hamilton matrix" between the calculated states
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM( STRING, WDES, CHAM )
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    COMPLEX(q) ::  CHAM(WDES%NB_TOT,WDES%NB_TOT)
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING

    NPL2=MIN(12,NB_TOT)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)

    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)


1   FORMAT(1I2,3X,40F9.5)
!1   FORMAT(1I2,3X,40F14.9)
2   FORMAT(1I2,3X,40F9.5)

  END SUBROUTINE DUMP_HAM


!***********************************************************************
!
!> dump a "Hamilton matrix" between the calculated states
!> single precision version
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM_SINGLE( STRING, WDES, CHAM)
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    COMPLEX(qs) ::  CHAM(WDES%NB_TOT,WDES%NB_TOT)
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING
    NPL2=MIN(10,NB_TOT)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)

    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)

1   FORMAT(1I2,3X,40F9.5)
2   FORMAT(1I2,3X,40F9.5)
!2   FORMAT(1I2,3X,40E9.1)

  END SUBROUTINE DUMP_HAM_SINGLE


!***********************************************************************
!
!> dump a "Hamilton matrix" between the calculated states
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM_SELECTED( STRING, WDES, CHAM, NDIM, NBANDS)
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    INTEGER NDIM
    INTEGER NBANDS
    COMPLEX(q) ::  CHAM(NDIM,NDIM)
! local
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING
    NPL2=MIN(12,NBANDS)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)

    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)

1   FORMAT(1I2,3X,24F9.5)
2   FORMAT(1I2,3X,24F9.5)

  END SUBROUTINE DUMP_HAM_SELECTED


!=======================================================================
!
!> small routine to dump a distributed matrix
!>
!> the descriptor in DESCA must properly describe the matrix
!> since RECON_SLICE calls check
!
!=======================================================================


  SUBROUTINE DUMP_HAM_DISTRI( STRING, WDES, CHAM_DISTRI, NB_TOT, DESCA, IU)
    USE wave
    USE scala
    IMPLICIT NONE
!> string to dump
    CHARACTER (LEN=*) :: STRING
!> wave function descriptor
    TYPE (wavedes)    :: WDES
!> distributed matrix
    COMPLEX(q)              :: CHAM_DISTRI(*)
    INTEGER           :: NB_TOT
!> distributed matrix descriptor array
    INTEGER           :: DESCA(*)
!> unit to write to (not dump for IU<0)
    INTEGER           :: IU
! local
    INTEGER, PARAMETER :: NDUMP=16
    COMPLEX(q),ALLOCATABLE  :: CHAM(:,:)
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW

    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  DESCA, COLUMN_LOW, COLUMN_HIGH)

    CALL M_sum_z8(WDES%COMM_KIN, CHAM(1,1), SIZE(CHAM,KIND=qi8))

    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1,N2+COLUMN_LOW-1) ,KIND=q) ,N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)

    DO N1=1,NDUMP
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1,N2+COLUMN_LOW-1)),N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)

    ENDIF

    DEALLOCATE(CHAM)
!1   FORMAT(1I2,3X,40F9.5)
!2   FORMAT(1I2,3X,40F9.5)
1   FORMAT(1I2,3X,40F7.4)
2   FORMAT(1I2,3X,40F7.4)
!1   FORMAT(1I2,3X,40F14.9)

  END SUBROUTINE DUMP_HAM_DISTRI


!************************ SUBROUTINE ORSP   ****************************
!
!> this subroutine perfomes a gram-schmidt orthogonalistion
!>
!> of a set
!> of vectors (all elements on local node)
!> the subroutine uses BLAS 3 calls
!
!***********************************************************************

  SUBROUTINE ORSP(NBANDS, NPL, NRPLWV, CPTWFP)

!! USE moffload

    USE prec
    USE tutor, ONLY: vtutor
    IMPLICIT NONE

    INTEGER NBANDS
    INTEGER NPL
    INTEGER NRPLWV
    COMPLEX(q) CPTWFP(NRPLWV,NBANDS)
! local
    COMPLEX(q) CPRO(NBANDS)
    REAL(q) WFMAG
    INTEGER I,N

    COMPLEX(q), EXTERNAL :: ZDOTC

    IF (NBANDS> NRPLWV) THEN
       CALL vtutor%bug("internal error in ORSP: leading dimension of matrix too small", "subrot.F", 1870)
    ENDIF

!$ACC DATA CREATE(CPRO) 
!$ACC KERNELS PRESENT(CPRO) 
    CPRO=0
!$ACC END KERNELS
    DO N=1,NBANDS

! normalise the vector

       WFMAG= ZDOTC (NPL,CPTWFP(1,N),1,CPTWFP(1,N),1)
       CALL   ZDSCAL(NPL,1/SQRT(WFMAG),CPTWFP(1,N),1)

! now orthogonalise all higher vectors to the
! present vector

       IF (NBANDS/=N ) THEN
          CALL ZGEMV( 'C', NPL , NBANDS-N ,(1._q,0._q) , CPTWFP(1,N+1), &
               &             NRPLWV, CPTWFP(1,N), 1 , (0._q,0._q) ,  CPRO, 1)

!$ACC PARALLEL LOOP PRESENT(CPRO) 
          DO I=1,NBANDS
             CPRO(I)=CONJG(CPRO(I))
          ENDDO

          CALL ZGEMM( 'N', 'T' , NPL , NBANDS-N , 1 , -(1._q,0._q) , &
               &             CPTWFP(1,N), NRPLWV , CPRO , NBANDS , &
               &             (1._q,0._q) , CPTWFP(1,N+1) , NRPLWV )
       ENDIF
    ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!$ACC END DATA
    RETURN
  END SUBROUTINE ORSP


!***********************************************************************
!
!> use Loewdin perturbation to determine a rotation matrix
!>
!> this preserves the ordering of the eigenvalues
!> MIND: does not work for real matrices
!
!***********************************************************************

  SUBROUTINE LOEWDIN_DIAG(NB_TOT, NBDIM, CHAM)

!! USE moffload

    USE prec
    USE tutor, ONLY: vtutor
    IMPLICIT NONE
    INTEGER NB_TOT, NBDIM
    COMPLEX(q) :: CHAM(NBDIM, NB_TOT)
! local
    REAL(q), PARAMETER  :: DIFMAX=0.001_q
    REAL(q) DIFCEL
    INTEGER N1, N2
    COMPLEX(q) :: CROT
    REAL(q) :: FAKT

    IF (NB_TOT>NBDIM) THEN
       CALL vtutor%bug("internal error in LOEWDIN_DIAG: leading dimension of matrix too small", "subrot.F", 1933)
    ENDIF
!!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CHAM) PRIVATE(DIFCEL,CROT,FAKT) 
!$ACC PARALLEL LOOP PRESENT(CHAM) PRIVATE(N1,DIFCEL,CROT,FAKT) 
    DO N2=1,NB_TOT
       DO N1=1,N2-1
          DIFCEL= REAL( CHAM(N2,N2)-CHAM(N1,N1) ,KIND=q)
          IF (ABS(DIFCEL)<DIFMAX) THEN
             CROT  =0
          ELSE
             CROT  =CONJG(CHAM(N1,N2))/DIFCEL
             IF (ABS(CROT)>0.1_q) THEN
                FAKT= 0.1_q/ABS(CROT)
                CROT  = CROT*FAKT
             ENDIF
          ENDIF
          CHAM(N2,N1) =-CROT
          CHAM(N1,N2) =-CONJG(CROT)
       ENDDO
    ENDDO
!$ACC PARALLEL LOOP PRESENT(CHAM) 
    DO N1=1,NB_TOT
       CHAM(N1,N1)=1
    ENDDO
  END SUBROUTINE LOEWDIN_DIAG


!*******************************************************************
!>  calculate the matrix elements of a local potential
!>
!>  ~~~
!>  CHAM(i,j)= <psi_i,k| V |psi_j,k> = int psi_i,k*(r) V(r) psi_j,k(r)
!>  ~~~
!>
!> between states
!> the argument CVTOT must be in real space
!> the result is retured in CHAM
!
!*******************************************************************

  SUBROUTINE LOCAL_BETWEEN_STATES( HAMILTONIAN, W, LATT_CUR, P, T_INFO, IRDMAX, LMDIM, &
       GRID_SOFT, GRIDC, GRIDUS, SOFT_TO_C, C_TO_US, CVTOT, CHAM)
    USE prec
    USE wave_high
    USE lattice
    USE poscar
    USE pseudo
    USE pot
    USE pawm
    USE subrot
    USE hamil_struct_def
    USE us
    USE tutor, ONLY: vtutor
    
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (wavespin)    W
    
    INTEGER  IRDMAX      ! allocation required for augmentation
    TYPE (latt)        LATT_CUR
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
!> grid for potentials/charge
    TYPE (grid_3d)     GRIDC
!> grid for soft chargedensity
    TYPE (grid_3d)     GRID_SOFT
!> grid for augmentation
    TYPE (grid_3d)     GRIDUS
!> index table between GRID_SOFT and GRIDC
    TYPE (transit)     SOFT_TO_C
!> index table between GRID_SOFT and GRIDC
    TYPE (transit)     C_TO_US
!> local potential
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,W%WDES%NCDIJ)
    INTEGER LMDIM
    COMPLEX(q)       CHAM(W%WDES%NB_TOT,W%WDES%NB_TOT,W%WDES%NKPTS,W%WDES%ISPIN)
! local
    INTEGER ISP, NK
    COMPLEX(q) ::   SV(W%WDES%GRID%MPLWV,W%WDES%NCDIJ)   ! local potential
    COMPLEX(q) :: CDIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    COMPLEX(q) :: CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    INTEGER IRDMAA
    REAL(q)  DISPL(3,T_INFO%NIONS)


    IF (W%WDES%COMM_KINTER%NCPU.NE.1) THEN
!PK Trivial but callers must be adapted
       CALL vtutor%error("LOCAL_BETWEEN_STATES: KPAR>1 not tested (but seems ok), sorry.")
    END IF


    DISPL=0

! get  the non local strenght parameters
    CALL SETDIJ_(W%WDES, GRIDC, GRIDUS, C_TO_US, LATT_CUR, P, T_INFO, W%WDES%LOVERL, &
         LMDIM, CDIJ, CQIJ, CVTOT, .FALSE., IRDMAA, IRDMAX, DISPL)

! transform CVTOT to reciprocal space (required by SET_SV)
    DO ISP=1,W%WDES%NCDIJ
       CALL FFT_RC_SCALE(CVTOT(1,ISP),CVTOT(1,ISP),GRIDC)
    ENDDO

! now set SV from CVTOT
    CALL SET_SV( W%WDES%GRID, GRIDC, GRID_SOFT, W%WDES%COMM_INTER, SOFT_TO_C, W%WDES%NCDIJ, SV, CVTOT)
    
    CHAM=0

    DO ISP=1,W%WDES%NCDIJ
       DO NK=1,W%WDES%NKPTS
          CALL ONE_CENTER_BETWEEN_STATES( HAMILTONIAN, LATT_CUR, W%WDES%LOVERL, W%WDES, W, NK, ISP, LMDIM, &
               CDIJ, CHAM(1,1,NK,ISP), SV)
!         CALL DUMP_HAM( "Hamiltonian", W%WDES, CHAM(1,1,NK,ISP))

       ENDDO
    ENDDO
       

  END SUBROUTINE LOCAL_BETWEEN_STATES
