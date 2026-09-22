# 1 "choleski2.F"
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


# 4 "choleski2.F" 2 

MODULE choleski
  USE prec
  USE dfast
CONTAINS
!************************ SUBROUTINE ORTHCH ****************************
! RCS:  $Id: choleski2.F,v 1.5 2002/04/16 07:28:38 kresse Exp $
!
! this subroutine orthonormalizes a set of complex (wave-)functions
! using a Choleski-decomposition of the overlap matrix (O = L L^H)
! in conjunction with inversion of the result of this decomposition
! (U --> U^-1). If we have arbitrary functions {|cfin_i>} on input,
! we have first to build up the overlap matrix OVL_ij=<cfin_i|cfin_j>
! then we have to decompose it (OVL_ij=L_ik U_kj), have to invert
! U_ij and then to form the output set |cfout_i>=U^-1_ji|cfin_j>. As
! (1._q,0._q) can show easily it holds: <cfout_i|cfout_j>=delta_ij !! Voila!
!
!> @details @ref openmp :
!> NSTRIP is multiplied with mopenmp::omp_nthreads.
!> This improves the performance of the ::orth1 and ::orth1_distri
!> calls under OpenMP.
!
!***********************************************************************

      SUBROUTINE ORTHCH(WDES,W, LOVERL,LMDIM,CQIJ,NKSTART,NKSTOP)

# 33


      USE prec
      USE scala
      USE wave_high
      USE wave
      USE wave_mpi
      USE tutor, ONLY: ORTHONORMALIZATION_FAILED, newArgument, isError
      USE mopenmp_struct_def, ONLY : omp_nthreads
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin) W
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      INTEGER, OPTIONAL :: NKSTART,NKSTOP

      INTEGER :: MY_NKSTART,MY_NKSTOP
      LOGICAL   LOVERL
      COMPLEX(q)   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),ALLOCATABLE,TARGET :: CPROW(:,:),COVL(:,:)
! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW_RED(:,:)
      COMPLEX(q)   , POINTER :: CPROW_RED(:,:),CPROJ_RED(:,:)
      LOGICAL DO_REDIS
      TYPE (REDIS_PW_CTR),POINTER :: H_PW

      LOGICAL IS_ALREADY_REDIS   ! are the plane wave coefficients already redistributed
      LOGICAL LASYNC_LOCAL       ! ansyncronous exchange of pw coefficients

      LOGICAL LscaAWARE_LOCAL

      

# 76


      LscaAWARE_LOCAL = LscaAWARE
# 83



      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 93

!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      IF (NCPU /= 1) THEN

        DO_REDIS=.TRUE.
        NRPLWV_RED=WDES%NRPLWV/NCPU
        NPROD_RED =WDES%NPROD /NCPU
        LASYNC_LOCAL=LASYNC
! it is possible that the bands are distributed of plane wave coefficients upon entry
        IS_ALREADY_REDIS=W%OVER_BAND
! of course no need to do it in overlap with calculations
        IF (IS_ALREADY_REDIS) LASYNC_LOCAL=.FALSE.

      ELSE

        DO_REDIS=.FALSE.
        NRPLWV_RED=WDES%NRPLWV
        NPROD_RED =WDES%NPROD

      ENDIF
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS
      NSTRIP=NSTRIP_STANDARD
!$    NSTRIP=NSTRIP*omp_nthreads

      NSTRIP=MIN(NSTRIP,NBANDS)

! allocate work space
      ALLOCATE(CPROW(WDES%NPROD,NBANDS))

      IF (.NOT.LscaAWARE_LOCAL) THEN
         ALLOCATE(COVL(NB_TOT,NB_TOT))
      ELSE
         CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
         ALLOCATE(COVL(SCALA_NP(), SCALA_NQ()))
      ENDIF
!$ACC ENTER DATA CREATE(WDES1,CPROW,COVL) 

      MY_NKSTART=1
      MY_NKSTOP=WDES%NKPTS
      IF(PRESENT(NKSTART))MY_NKSTART=NKSTART
      IF(PRESENT(NKSTOP))MY_NKSTOP=NKSTOP
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=MY_NKSTART,MY_NKSTOP

      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
      IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES%COMM_KIN,WDES%NB_TOTK(NK,ISP))
      CALL SETWDES(WDES,WDES1,NK)
!   get pointers for redistributed wavefunctions
!   I can not guarantee that this works with all f90 compilers
!   please see comments in wave_mpi.F
      IF (NCPU /= 1) THEN
        CALL SET_WPOINTER(CW_RED,   NRPLWV_RED, NB_TOT, W%CPTWFP(1,1,NK,ISP))
        CALL SET_GPOINTER(CPROJ_RED, NPROD_RED, NB_TOT, W%CPROJ(1,1,NK,ISP))
        CALL SET_GPOINTER(CPROW_RED, NPROD_RED, NB_TOT, CPROW(1,1))
      ELSE
        CW_RED    => W%CPTWFP(:,:,NK,ISP)
        CPROJ_RED => W%CPROJ(:,:,NK,ISP)
        CPROW_RED => CPROW(:,:)
      ENDIF

!   set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO

      CALL SET_NPL_NPRO(WDES1, NPL, NPRO)

      NPRO_O=NPRO
      IF (.NOT. LOVERL) NPRO_O=0
!=======================================================================
!  calculate overlap matrix (only upper triangle is needed):
!=======================================================================
      IF (DO_REDIS .AND. LASYNC_LOCAL) THEN
         CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
# 177

         DO NPOS=1,NSTRIP
            CALL REDIS_PW_START(WDES, W%CPTWFP(1,NPOS,NK,ISP), NPOS, H_PW)
         ENDDO
# 185

      ENDIF

      CALL OVERL(WDES1,LOVERL,LMDIM,CQIJ(1,1,1,ISP),W%CPROJ(1,1,NK,ISP),CPROW(1,1))
! redistribute everything

      
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
        
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        
        IF (.NOT. LASYNC_LOCAL .AND. .NOT. IS_ALREADY_REDIS)  CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP   (1,1,NK,ISP))
        
      ENDIF

      IF (.NOT.LscaAWARE_LOCAL) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(COVL) 
         DO N=1,NB_TOT
            DO I=1,NB_TOT
               COVL(I,N)=(0._q,0._q)
            ENDDO
         ENDDO
       ENDIF
!
! there is a strange bug in the PII optimized blas DGEMM, which seems
! to access in certain instances data beyond the last array element
! if a matrix is multiplied with a vector (second call to ORTH1)
! to work around this I avoid calling ORTH1 with NB_TOT-NPOS+1=1
      DO NPOS=1,NBANDS-NSTRIP,NSTRIP
        IF (DO_REDIS .AND. LASYNC_LOCAL) THEN
           DO NPOS_=NPOS,NPOS+NSTRIP-1
              CALL REDIS_PW_STOP (WDES, W%CPTWFP(1,NPOS_,NK,ISP), NPOS_, H_PW)
           ENDDO
# 224

           DO NPOS_=NPOS+NSTRIP,MIN(NPOS+2*NSTRIP-1,NBANDS)
              CALL REDIS_PW_START(WDES, W%CPTWFP(1,NPOS_,NK,ISP), NPOS_, H_PW)
           ENDDO
# 232

        ENDIF

        NPOS_RED  =(NPOS-1)*NCPU+1
        NSTRIP_RED=NSTRIP*NCPU

        IF (.NOT.LscaAWARE_LOCAL) THEN
           CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
                CPROW_RED(1,NPOS_RED),NB_TOT, &
                NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
        ELSE
           CALL ORTH1_DISTRI('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
                CPROW_RED(1,NPOS_RED),NB_TOT, &
                NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1), & 
                WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
        ENDIF
      ENDDO

      IF (DO_REDIS .AND. LASYNC_LOCAL) THEN
         DO NPOS_=NPOS,NBANDS
            CALL REDIS_PW_STOP (WDES, W%CPTWFP(1,NPOS_,NK,ISP), NPOS_, H_PW)
         ENDDO
      ENDIF
      
      NPOS_RED  =(NPOS-1)*NCPU+1
      NSTRIP_RED=(NBANDS-NPOS+1)*NCPU

      IF (.NOT.LscaAWARE_LOCAL) THEN
         CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
              CPROW_RED(1,NPOS_RED),NB_TOT, &
              NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
      ELSE
         CALL ORTH1_DISTRI('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
              CPROW_RED(1,NPOS_RED),NB_TOT, &
              NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1), &  
              WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
      ENDIF

      IF (DO_REDIS .AND. LASYNC_LOCAL) CALL REDIS_PW_DEALLOC(H_PW)

      
      IF (.NOT.LscaAWARE_LOCAL) THEN
         CALL M_sum_z8(WDES%COMM_KIN,COVL(1,1),SIZE(COVL,KIND=qi8))
      ENDIF
      

# 282

!=======================================================================
! Choleski-decomposition of the overlap matrix + inversion of the result
! calling LAPACK-routines ZPOTRF (decomposition) and ZTRTRI (inversion):
!=======================================================================
      

! Is the overlap matrix distributed block-cyclically?
      IF (.NOT.LscaAWARE_LOCAL) THEN
! If the Hamiltonian is not block-cyclically distributed
! we may still use ScaLAPACK ... (but not in the OpenACC version)
         IF (LscaLAPACK.AND.LscaLU) THEN

            CALL pPOTRF_TRTRI(WDES%COMM_KIN, COVL(1,1),NB_TOT, WDES%NB_TOTK(NK,ISP))
            CALL M_sum_z8(WDES%COMM_KIN,COVL(1,1),SIZE(COVL,KIND=qi8))

! or call plain LAPACK routines
         ELSE

            INFO=0
# 304

            CALL ZPOTRF &

                 & ('U',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%write(isError, ORTHONORMALIZATION_FAILED, newArgument(ival=[1, NK, ISP]))
            ENDIF
# 313

            CALL ZTRTRI &

                 & ('U','N',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%write(isError, ORTHONORMALIZATION_FAILED, newArgument(ival=[2, NK, ISP]))
            ENDIF

         ENDIF

! If the Hamiltonian is block-cyclically distributed
! we will definitely use ScaLAPACK
      ELSE

         CALL BG_pPOTRF_TRTRI(COVL(1,1), WDES%NB_TOTK(NK,ISP), INFO)
         IF (INFO/=0) THEN
            CALL vtutor%write(isError, ORTHONORMALIZATION_FAILED, newArgument(ival=[3, NK, ISP]))
         ENDIF

      ENDIF

      

      
!=======================================================================
!  construct the orthogonal set:
!=======================================================================
      IF (.NOT.LscaAWARE_LOCAL) THEN
         CALL LINCOM('U',CW_RED,CPROJ_RED,COVL(1,1), &
              WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP),NPL,NPRO,NRPLWV_RED,NPROD_RED,NB_TOT, &
              CW_RED,CPROJ_RED)
      ELSE
         CALL LINCOM_DISTRI('U',CW_RED(1,1),CPROJ_RED(1,1),COVL(1,1), &
              WDES%NB_TOTK(NK,ISP),NPL,NPRO,NRPLWV_RED,NPROD_RED,NB_TOT, & 
              WDES%COMM_KIN, NBLK)
      ENDIF

      

!  back redistribution
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        IF (LASYNC_LOCAL .OR. IS_ALREADY_REDIS) THEN
           W%OVER_BAND=.TRUE.
        ELSE
! if the routine was entered with already redistributed wave-functions
! then no need to redistribute them back
           CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP   (1,1,NK,ISP))
        ENDIF
      ENDIF
      
!=======================================================================
      ENDDO kpoints
      ENDDO spin
!=======================================================================
!$ACC EXIT DATA DELETE(CPROW,COVL) 
      DEALLOCATE(CPROW,COVL)

# 382


      

      RETURN
    END SUBROUTINE ORTHCH


!************************ SUBROUTINE ORTHCH_DUAL ***********************
! RCS:  $Id: choleski2.F,v 1.5 2002/04/16 07:28:38 kresse Exp $
!
! this subroutine determines a dual set of orbitals with the
! property
!  < phi_i | phi_j> = delta_ij
! the dual set W_DUAL is passed down as a pointer to a wavespin variable
! if NC potential are used W_DUAL is simply linked to W
! otherwise if W_DUAL is not associated, the required arrays are
! allocated and determined (make sure to NULLIFY W_DUAL before
!   calling the routined)
! if it is already associated, the routine recalculates W_DUAL from W but
! assumes that the proper allocation was 1._q before
!
!***********************************************************************

    SUBROUTINE ORTHCH_DUAL(WDES,W,W_DUAL, LOVERL,LMDIM,CQIJ)
      USE prec
      USE scala
      USE wave_high
      USE wave
      USE wave_mpi
      IMPLICIT NONE

      TYPE (wavespin), TARGET  ::  W
      TYPE (wavespin), POINTER ::  W_DUAL
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1

      LOGICAL   LOVERL
      INTEGER   LMDIM
      COMPLEX(q)   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),ALLOCATABLE,TARGET :: COVL(:,:)
! redistributed plane wave coefficients
      TYPE (wavefuna)    WA             ! array pointers
      TYPE (wavefuna)    WA_DUAL        ! array pointers
      INTEGER   NODE_ME, IONODE, NCPU
      INTEGER   NB_TOT, NBANDS, NSTRIP, NSTRIP_ACT, ISP, NK, N, I, INFO, NPOS, NPOS_RED, NSTRIP_RED
! scaAWARE is often a little bit tricky to implement
! routine was tested and currenty LscaAWARE_LOCAL is set to LscaAWARE
! see below
      LOGICAL :: LscaAWARE_LOCAL=.FALSE.

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 440

      IF (.NOT. LOVERL) THEN
         W_DUAL => W
         RETURN
      ENDIF

! if there are any subspected problems with duality comment out
! the following line (or set LscaAWARE_LOCAL =.FALSE.)
      LscaAWARE_LOCAL=LscaAWARE

      IF (.NOT. ASSOCIATED(W_DUAL)) THEN
         ALLOCATE(W_DUAL)
         CALL ALLOCW( WDES, W_DUAL)
!#define test_w_dual
# 459

      ENDIF


      W_DUAL%CPROJ=W%CPROJ
      W_DUAL%CPTWFP   =W%CPTWFP
      W_DUAL%CELTOT=W%CELTOT
      W_DUAL%FERTOT=W%FERTOT

!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS
      NSTRIP=NSTRIP_STANDARD

      IF (.NOT. LscaAWARE_LOCAL) THEN
         ALLOCATE(COVL(NB_TOT,NB_TOT))
      ELSE
         CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
         ALLOCATE(COVL(SCALA_NP(), SCALA_NQ()))
      ENDIF
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS

      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

!=======================================================================
      IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES%COMM_KIN,WDES%NB_TOTK(NK,ISP))

      CALL SETWDES(WDES,WDES1,NK)
      WA     =ELEMENTS(W, WDES1, ISP)
      WA_DUAL=ELEMENTS(W_DUAL, WDES1, ISP)
!=======================================================================
!  calculate overlap matrix (only upper triangle is needed):
!=======================================================================

      
      IF (WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
         CALL REDISTRIBUTE_PROJ(WA_DUAL)
         CALL REDISTRIBUTE_PW(WA_DUAL)
      ENDIF

      IF (.NOT. LscaAWARE_LOCAL) THEN
         DO N=1,NB_TOT
            DO I=1,NB_TOT
               COVL(I,N)=(0._q,0._q)
            ENDDO
         ENDDO
      ENDIF
!
      DO NPOS=1,NBANDS,NSTRIP
        NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)
        NPOS_RED  =(NPOS-1)*NCPU+1
        NSTRIP_RED=NSTRIP_ACT*NCPU

        IF (.NOT. LscaAWARE_LOCAL) THEN
! W^+ x W_DUAL (W_DUAL=W usually)
           CALL ORTH1('U',WA%CW_RED(1,1),WA_DUAL%CW_RED(1,NPOS_RED),WA%CPROJ_RED(1,1), &
                WA_DUAL%CPROJ_RED(1,1),NB_TOT, &
                NPOS_RED,NSTRIP_RED, WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
        ELSE
           CALL ORTH1_DISTRI('U',WA%CW_RED(1,1),WA_DUAL%CW_RED(1,NPOS_RED),WA%CPROJ_RED(1,1), &
                WA_DUAL%CPROJ_RED(1,1),NB_TOT, &
                NPOS_RED,NSTRIP_RED, WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1), & 
                WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
        ENDIF
      ENDDO
      
      
      IF (.NOT. LscaAWARE_LOCAL) THEN
         CALL M_sum_z8(WDES%COMM_KIN,COVL(1,1),SIZE(COVL,KIND=qi8))
      ENDIF
      

# 541

!=======================================================================
! Choleski-decomposition of the overlap matrix + inversion of the result
! calling LAPACK-routines ZPOTRF (decomposition) and ZPOTRI (inversion):
!=======================================================================
      IF (.NOT. LscaAWARE_LOCAL) THEN
         IF (LscaLAPACK .AND. LscaLU ) THEN
            CALL pPOTRF_POTRI(WDES%COMM_KIN, COVL(1,1),NB_TOT, WDES%NB_TOTK(NK,ISP))
            CALL M_sum_z8(WDES%COMM_KIN,COVL(1,1),SIZE(COVL,KIND=qi8))
         ELSE
            INFO=0
# 554

            CALL ZPOTRF &

                 & ('U',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%error("LAPACK: Routine ZPOTRF failed! " // str(INFO) // " " // str(NK) // &
                  " " // str(ISP))
            ENDIF
# 564

            CALL ZPOTRI &

                 & ('U',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%error("LAPACK: Routine ZTRTRI failed! " // str(INFO) // " " // str(NK) // &
                  " " // str(ISP))
            ENDIF
         ENDIF
! add lower triangle not calculated by POTRI
         DO N=1,NB_TOT
            DO I=N+1,NB_TOT
               COVL(I,N)=CONJG(COVL(N,I))
            ENDDO
         ENDDO
      ELSE
         CALL BG_pPOTRF_POTRI(COVL(1,1), WDES%NB_TOTK(NK,ISP), INFO)
         IF (INFO/=0) THEN
            CALL vtutor%error("LAPACK: Routine ZTRTRI failed! " // str(INFO) // " " // str(NK) // " &
               &" // str(ISP))
         ENDIF
      ENDIF
# 590


      
!=======================================================================
!  construct the dual set by inplace transformation of W_DUAL
!=======================================================================
      IF (.NOT. LscaAWARE_LOCAL) THEN
         CALL LINCOM('F',WA_DUAL%CW_RED,WA_DUAL%CPROJ_RED,COVL(1,1), &
              WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
              WDES1%NPL_RED, WDES1%NPRO_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, NB_TOT, &
              WA%CW_RED, WA%CPROJ_RED)
      ELSE
! POTRF set only the upper triangle of the transformation matrix
! add the lower triangle using DISTRI_HERM
         CALL LINCOM_DISTRI_HERM('F',WA_DUAL%CW_RED(1,1),WA_DUAL%CPROJ_RED(1,1),COVL(1,1), &
              WDES%NB_TOTK(NK,ISP), & 
              WDES1%NPL_RED, WDES1%NPRO_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, NB_TOT, & 
              WDES%COMM_KIN, NBLK)
      ENDIF

      

!  back redistribution
      IF (WDES%DO_REDIS) THEN
        CALL REDISTRIBUTE_PROJ(WA)
        CALL REDISTRIBUTE_PW  (WA)
        CALL REDISTRIBUTE_PROJ(WA_DUAL)
        CALL REDISTRIBUTE_PW  (WA_DUAL)
      ENDIF
      
!=======================================================================
      ENDDO kpoints
      ENDDO spin
!=======================================================================
      DEALLOCATE(COVL)

      RETURN
      END SUBROUTINE ORTHCH_DUAL


!************************ SUBROUTINE DUMP_S ****************************
!
! dump the overlap matrix between the occupied states
! this is essentially identical to the previous routine
!
!***********************************************************************


      SUBROUTINE DUMP_S(WDES,W, LOVERL,LMDIM,CQIJ)
      USE prec
      USE scala
      USE wave_high
      USE wave
      USE wave_mpi
      USE dfast
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin) W
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1

      LOGICAL   LOVERL
      COMPLEX(q)   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),ALLOCATABLE,TARGET :: CPROW(:,:),COVL(:,:)

! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW_RED(:,:)
      COMPLEX(q)   , POINTER :: CPROW_RED(:,:),CPROJ_RED(:,:)
      LOGICAL DO_REDIS
      REAL(q) :: NELECT

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 669

      IF (NCPU /= 1) THEN

         DO_REDIS=.TRUE.
         NRPLWV_RED=WDES%NRPLWV/NCPU
         NPROD_RED =WDES%NPROD /NCPU
      ELSE
         DO_REDIS=.FALSE.
         NRPLWV_RED=WDES%NRPLWV
         NPROD_RED =WDES%NPROD
      ENDIF

      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS

      NSTRIP=NSTRIP_STANDARD

! allocate work space
      ALLOCATE(CPROW(WDES%NPROD,NBANDS),COVL(NB_TOT,NB_TOT))

      NELECT=0

      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS

      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

      CALL SETWDES(WDES,WDES1,NK)
! get pointers for redistributed wavefunctions
      IF (NCPU /= 1) THEN
         CALL SET_WPOINTER(CW_RED,   NRPLWV_RED, NB_TOT, W%CPTWFP(1,1,NK,ISP))
         CALL SET_GPOINTER(CPROJ_RED, NPROD_RED, NB_TOT, W%CPROJ(1,1,NK,ISP))
         CALL SET_GPOINTER(CPROW_RED, NPROD_RED, NB_TOT, CPROW(1,1))
      ELSE
         CW_RED    => W%CPTWFP(:,:,NK,ISP)
         CPROJ_RED => W%CPROJ(:,:,NK,ISP)
         CPROW_RED => CPROW(:,:)
      ENDIF

! set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO

      CALL SET_NPL_NPRO(WDES1, NPL, NPRO)

      NPRO_O=NPRO
      IF (.NOT. LOVERL) NPRO_O=0

      CALL OVERL(WDES1, LOVERL,LMDIM,CQIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP),CPROW(1,1))

! redistribute everything
      IF (DO_REDIS) THEN
         CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
         CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP   (1,1,NK,ISP))
      ENDIF

! calculate overlap
      DO N=1,NB_TOT
         DO I=1,NB_TOT
            COVL(I,N)=(0._q,0._q)
         ENDDO
      ENDDO

      DO NPOS=1,NBANDS-NSTRIP,NSTRIP
         NPOS_RED  =(NPOS-1)*NCPU+1
         NSTRIP_RED=NSTRIP*NCPU

         CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
              CPROW_RED(1,NPOS_RED),NB_TOT, &
              NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
      ENDDO

      
      NPOS_RED  =(NPOS-1)*NCPU+1
      NSTRIP_RED=(NBANDS-NPOS+1)*NCPU

      CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
           CPROW_RED(1,NPOS_RED),NB_TOT, &
           NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))

      DO N=1,NB_TOT
         NELECT=NELECT+ COVL(N, N)*W%FERWE(N,NK,ISP)*WDES%RSPIN*WDES%WTKPT(NK)
      ENDDO

      CALL M_sum_z8(WDES%COMM_KIN,COVL(1,1),SIZE(COVL,KIND=qi8))

! back redistribution
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP   (1,1,NK,ISP))
      ENDIF

      WRITE(*,*) 'NK=',NK,'ISP=',ISP
      NPL2=MIN(10,NB_TOT)
      DO  N1=1,NPL2
         WRITE(*,1)N1,(REAL( COVL(N1,N2) ,KIND=q) ,N2=1,NPL2)
      ENDDO
      WRITE(*,*)

      DO N1=1,NPL2
         WRITE(*,2)N1,(AIMAG(COVL(N1,N2)),N2=1,NPL2)
      ENDDO

      WRITE(*,*)
      WRITE(*,'(20E10.3)') 0,(REAL( COVL(N1,N1)-1 ,KIND=q) ,N1=1,NPL2)
      WRITE(*,*)
1     FORMAT(1I2,3X,20F9.5)
2     FORMAT(1I2,3X,20E9.1)

      ENDDO kpoints
      ENDDO spin

      WRITE(*,'(A,F14.8)') 'total number of electrons NELECT=',NELECT
      DEALLOCATE(CPROW,COVL)

      RETURN
      END SUBROUTINE

!                              __________________
!_____________________________/ CALC_PAW_OVERLAP \______________________________
!
!> @brief Calculates the S-overlap between two wavefunctions WA1 and WA2.
!>
!> The wavefunction descriptors associated with WA1 and WA2 have to be
!> compatible, preferably identical. WA1 and WA2 may refer to the same storage
!> location. This is detected and 1 redistribution is handled correctly for
!> this case.
!>
!> @param[in] WA1             Bra state
!> @param[in] WA2             Ket state
!> @param[in] CQIJ            Augmentation charges
!> @param[in,out] OVERLAP_MAT Calculated overlap matrix
!> @param[in] ADD             If true, result is added to OVERLAP_MAT instead
!>                            of overwriting it.
!> @param[in] ONLY_AUG        If true, only augmentation part is calculated.
!_______________________________________________________________________________
      SUBROUTINE CALC_PAW_OVERLAP(WA1, WA2, CQIJ, OVERLAP_MAT, ADD, ONLY_AUG)
         USE wave_struct_def, ONLY: wavefuna, &
                                    wavedes1
         USE wave_high,       ONLY: OVERL, &
                                    REDISTRIBUTE_PW, &
                                    REDISTRIBUTE_PROJ, &
                                    NEWWAVA_PROJ, &
                                    DELWAVA_PROJ
         USE wave_mpi,        ONLY: SET_WPOINTER
         USE tutor,           ONLY: vtutor
         USE string,          ONLY: str
         IMPLICIT NONE
   
         TYPE(wavefuna),      INTENT(IN)     :: WA1, WA2
         COMPLEX(q),             INTENT(IN)     :: CQIJ(:, :, :, :)
         COMPLEX(q),                INTENT(INOUT)  :: OVERLAP_MAT(:, :)
         LOGICAL, OPTIONAL,   INTENT(IN)     :: ADD
         LOGICAL, OPTIONAL,   INTENT(IN)     :: ONLY_AUG
   
         INTEGER        :: NB_TOT
         INTEGER        :: STRIP
         INTEGER        :: NPL
         LOGICAL        :: WA_SAME
         TYPE(wavefuna) :: WOVL
         TYPE(wavedes1) :: WDES1
         LOGICAL        :: MY_ADD, MY_ONLY_AUG
   
         MY_ONLY_AUG = .FALSE.; IF (PRESENT(ONLY_AUG)) MY_ONLY_AUG = ONLY_AUG
         MY_ADD = .FALSE.; IF (PRESENT(ADD)) MY_ADD = ADD
   
         WA_SAME = ASSOCIATED(WA1%CPTWFP, WA2%CPTWFP)
   
         WDES1 = WA1%WDES1
         NB_TOT = WDES1%NB_TOT
   
         IF (ANY(SHAPE(OVERLAP_MAT) /= NB_TOT)) &
            CALL vtutor%bug("CALC_PAW_OVERLAP: Shape mismatch: " // &
               str(NB_TOT) // ", " // str(SIZE(OVERLAP_MAT, 1)) // ", " // str(SIZE(OVERLAP_MAT, 2)), &
               "choleski2.F", 844)
   
! Allocate only projections
         CALL NEWWAVA_PROJ(WOVL, WDES1)
   
         IF (MY_ADD) THEN

            OVERLAP_MAT = OVERLAP_MAT / WDES1%COMM%NCPU

         ELSE
            OVERLAP_MAT = 0
         ENDIF
   
! Reuse plane-wave components of WA2
         WOVL%CPTWFP=>WA2%CPTWFP
         IF (WDES1%DO_REDIS) THEN
            CALL SET_WPOINTER(WOVL%CW_RED, WDES1%NRPLWV_RED, NB_TOT, WOVL%CPTWFP(1, 1))
         ELSE
            WOVL%CW_RED=>WOVL%CPTWFP
         ENDIF
         
! Compute \sum_j Q_ij < p_j | W_P >
         CALL OVERL(WDES1, WDES1%LOVERL, SIZE(CQIJ, 1), CQIJ(1, 1, 1, 1), WA2%CPROJ(1, 1), WOVL%CPROJ(1, 1))
   
! Redistribute plane-wave components
         IF (.NOT. MY_ONLY_AUG) THEN
            CALL REDISTRIBUTE_PW(WA1)
            IF (.NOT. WA_SAME) CALL REDISTRIBUTE_PW(WOVL)
         ENDIF
! Redistribute projections
         CALL REDISTRIBUTE_PROJ(WA1)
         CALL REDISTRIBUTE_PROJ(WOVL)
   
         NPL = WDES1%NPL_RED
! NPL == 0 skips calculation of plane-wave part
         IF (MY_ONLY_AUG) NPL = 0
   
! Compute < WA1_m | S | WA2_n >
         DO STRIP = 1, NB_TOT - NSTRIP_STANDARD_GLOBAL, NSTRIP_STANDARD_GLOBAL
            CALL ORTH2( &
               WA1%CW_RED(1, 1), WOVL%CW_RED(1, STRIP), WA1%CPROJ_RED(1, 1), WOVL%CPROJ_RED(1, STRIP), &
               NB_TOT, STRIP, NSTRIP_STANDARD_GLOBAL, NPL, &
               WDES1%NPRO_O_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, OVERLAP_MAT(1, 1))
         ENDDO
   
         CALL ORTH2( &
            WA1%CW_RED(1, 1), WOVL%CW_RED(1, STRIP), WA1%CPROJ_RED(1, 1), WOVL%CPROJ_RED(1, STRIP), &
            NB_TOT, STRIP, NB_TOT-STRIP+1, NPL, &
            WDES1%NPRO_O_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, OVERLAP_MAT(1, 1))
   
         CALL M_sum_z8(WDES1%COMM_KIN, OVERLAP_MAT(1, 1), INT(NB_TOT, kind=qi8) * INT(NB_TOT, kind=qi8))
   
! Return to original plane-wave distribution
         IF (.NOT. MY_ONLY_AUG) THEN
            CALL REDISTRIBUTE_PW(WA1)
            IF (.NOT. WA_SAME) CALL REDISTRIBUTE_PW(WA2)
         ENDIF
! Return to original projection distribution
         CALL REDISTRIBUTE_PROJ(WA1)
   
         NULLIFY(WOVL%CPTWFP)
         CALL DELWAVA_PROJ(WOVL)
      END SUBROUTINE

END MODULE
