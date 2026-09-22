# 1 "wave_mpi.F"
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


# 2 "wave_mpi.F" 2 
!************************************************************************
! RCS:  $Id: wave_mpi.F,v 1.3 2001/04/03 10:43:18 kresse Exp $
!
!>  this module contains the routines required to communicate
!>  wavefunctions and projected wavefunctions
!>
!>  there are also two quite tricky routines SET_WPOINTER
!>  which allow to generate pointer to F77-sequenced arrays
!>  There is no guarantee that this will work on all computers, but
!>  currently it seems to be ok
!
!***********************************************************************

  MODULE wave_mpi

!! USE moffload

      USE prec
      USE mpimy
!
!> I need an interface block here because I pass the
!> first element of a pointer to CPTWFP
!> most F90 comiler do not allow such a construct
!
      INTERFACE
      SUBROUTINE SET_WPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP
      END SUBROUTINE
      END INTERFACE

      INTERFACE
      SUBROUTINE SET_GPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP
      END SUBROUTINE
      END INTERFACE


!> each communcation package needs a unique identifier (ICOMM)
!> this is handled by the global variable ICOMM_BASE_HANDLE
!> which is incremented/decremented by ICOMM_INCREMENT
!> whenever a redistribution handle is allocated

      INTEGER, PARAMETER :: ICOMM_BASE=10000, ICOMM_INCREMENT=1000
      INTEGER, SAVE :: ICOMM_BASE_HANDLE=ICOMM_BASE

!> is assyncronous communication allowed or not
      LOGICAL:: LASYNC=.FALSE.

      TYPE REDIS_PW_CTR
         INTEGER :: NB                   !< number of bands to be redistributed
         INTEGER :: NBANDS               !< maximum band index
         INTEGER :: NEXT                 !< next vacant slot
         COMPLEX(q), POINTER :: CPTWFP(:,:)  !< storage for redistribution
         INTEGER,POINTER :: BAND(:)      !< bands that are currently redistributed
         INTEGER,POINTER :: SREQUEST(:,:)!< handles for outstanding send requests
         INTEGER,POINTER :: RREQUEST(:,:)!< handles for outstanding receive requests
         INTEGER :: ICOMM                !< identifier for communication
         TYPE(communic),POINTER :: COMM  !< communication handle
      END TYPE REDIS_PW_CTR

    CONTAINS
!************************ SUBROUTINE REDIS_PW_ALLOC ********************
!
!> This subroutine allocates the storage required for asynchronous
!> redistribution of wave function coeffients
!>
!> the redistribution can be initiated for up to NSIM bands
!> at the same time
!
!***********************************************************************

    SUBROUTINE REDIS_PW_ALLOC(WDES, NSIM, H)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)   :: WDES
      INTEGER :: NSIM                  !< number of bands 1._q simultaneausly
      TYPE (redis_pw_ctr),POINTER :: H !< handle

      INTEGER NRPLWV

      NRPLWV=WDES%NRPLWV

      IF (MOD(NRPLWV,WDES%NB_PAR) /= 0) THEN
        CALL vtutor%bug("REDIS_PW_ALLOC: internal error(1) " // str(NRPLWV) // " " // str(WDES%NB_PAR), "wave_mpi.F", 93)
      ENDIF

      ALLOCATE(H)
!$ACC ENTER DATA CREATE(H) 
      H%NB=NSIM
      H%NBANDS=WDES%NBANDS
      H%COMM=>WDES%COMM_INTER
      H%ICOMM=ICOMM_BASE_HANDLE

      ALLOCATE(H%SREQUEST(H%COMM%NCPU,NSIM), H%RREQUEST(H%COMM%NCPU,NSIM))
      ALLOCATE(H%CPTWFP(NRPLWV,NSIM),H%BAND(NSIM))
!$ACC ENTER DATA CREATE(H%CPTWFP) 
      H%BAND=0

      H%NEXT=1

      ICOMM_BASE_HANDLE=ICOMM_BASE_HANDLE+ICOMM_INCREMENT

    END SUBROUTINE REDIS_PW_ALLOC


    SUBROUTINE REDIS_PW_DEALLOC(H)
      USE prec
      USE wave
      IMPLICIT NONE

      INTEGER :: NSIM                  ! number of bands 1._q simultaneausly
      TYPE (redis_pw_ctr),POINTER :: H ! handle
      INTEGER :: IND

      IF (.NOT.ASSOCIATED(H)) THEN
        CALL vtutor%bug("REDIS_PW_DEALLOC: internal error(1)", "wave_mpi.F", 125)
      ENDIF

      DO IND=1, H%NB
         IF ( H%BAND(IND) /= 0) THEN
            CALL vtutor%bug("internal error: REDIS_PW_DEALLOC not all lines STOPED", "wave_mpi.F", 130)
         ENDIF
      ENDDO
!$ACC EXIT DATA DELETE(H%CPTWFP) 
      DEALLOCATE(H%CPTWFP,H%BAND)
      DEALLOCATE(H%SREQUEST, H%RREQUEST)
!$ACC EXIT DATA DELETE(H) 
      DEALLOCATE(H)

      ICOMM_BASE_HANDLE=ICOMM_BASE_HANDLE-ICOMM_INCREMENT
      IF (ICOMM_BASE_HANDLE < ICOMM_BASE) THEN
        CALL vtutor%bug("REDIS_PW_DEALLOC: internal error(2) " // str(ICOMM_BASE_HANDLE) // " " // &
           str(ICOMM_BASE), "wave_mpi.F", 142)
      ENDIF

    END SUBROUTINE REDIS_PW_DEALLOC

  END MODULE wave_mpi

!=======================================================================
!
!> this routine returns a pointer to an SEQUENCED F77 like array
!> with a given storage convention
!>
!> 1. dimension is N1 second (1._q,0._q) N2
!
!=======================================================================

      SUBROUTINE SET_WPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      IMPLICIT NONE
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP(N1,N2)
      CW_P => CPTWFP(:,:)

      END SUBROUTINE

      SUBROUTINE SET_GPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      IMPLICIT NONE
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP(N1,N2)
      CW_P => CPTWFP(:,:)

      END SUBROUTINE

!=======================================================================
!>  distribute projector part of wavefunctions over all nodes
!=======================================================================

      SUBROUTINE DIS_PROJ(WDES1,CPROJ,CPROJ_LOCAL)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPROJ(WDES1%NPRO_TOT)
      COMPLEX(q) CPROJ_LOCAL(WDES1%NPRO)

      INTEGER NPRO,NT,LMMAXC,NPRO_POS,L,NI,NPRO_TOT
!
! quick copy if possible
      IF (WDES1%COMM_INB%NCPU==1) THEN
!$ACC KERNELS PRESENT(CPROJ_LOCAL,WDES1,CPROJ) 
         CPROJ_LOCAL(1:WDES1%NPRO)=CPROJ(1:WDES1%NPRO)
!$ACC END KERNELS
         RETURN
      ENDIF
      CALL M_sum_z(WDES1%COMM_INB, CPROJ, WDES1%NPRO_TOT)

!$ACC PARALLEL LOOP PRESENT(WDES1,CPROJ,CPROJ_LOCAL) PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) 
      DO NI=1,WDES1%NIONS
         NT=WDES1%ITYP(NI)
         LMMAXC=WDES1%LMMAX(NT)
         IF (LMMAXC/=0) THEN
            NPRO_POS=WDES1%NPRO_POS(NI)
            NPRO=WDES1%LMBASE(NI)
!$ACC LOOP
            DO L=1,LMMAXC
               CPROJ_LOCAL(L+NPRO)=CPROJ(L+NPRO_POS)
            ENDDO
          ENDIF
      ENDDO

      IF (WDES1%LNONCOLLINEAR) THEN
         
         NPRO_TOT= WDES1%NPRO_TOT/2
!$ACC PARALLEL LOOP PRESENT(WDES1,CPROJ,CPROJ_LOCAL) PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) 
         DO NI=1,WDES1%NIONS
            NT=WDES1%ITYP(NI)
            LMMAXC=WDES1%LMMAX(NT)
            IF (LMMAXC/=0) THEN
               NPRO_POS=WDES1%NPRO_POS(NI)
               NPRO=WDES1%LMBASE(NI)+WDES1%LMBASE(WDES1%NIONS+1)
!$ACC LOOP
               DO L=1,LMMAXC
                  CPROJ_LOCAL(L+NPRO)=CPROJ(L+NPRO_POS+NPRO_TOT)
               ENDDO
             ENDIF
         ENDDO

      ENDIF
# 241

      END SUBROUTINE

!=======================================================================
!>  merge projector part wavefunctions from all nodes
!>  definitely not the most efficient implementation
!>  but it works
!=======================================================================

      SUBROUTINE MRG_PROJ(WDES1,CPROJ,CPROJ_LOCAL)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPROJ(WDES1%NPRO_TOT)
      COMPLEX(q) CPROJ_LOCAL(WDES1%NPRO)

      INTEGER NPRO,NT,LMMAXC,NPRO_POS,L,NI,NPRO_TOT
!
! quick copy if possible
      IF (WDES1%COMM_INB%NCPU==1) THEN
!$ACC KERNELS PRESENT(WDES1,CPROJ,CPROJ_LOCAL) 
         CPROJ(1:WDES1%NPRO)=CPROJ_LOCAL(1:WDES1%NPRO)
!$ACC END KERNELS
         RETURN
      ENDIF
!$ACC KERNELS PRESENT(WDES1,CPROJ) 
      CPROJ(1:WDES1%NPRO_TOT)=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP PRESENT(WDES1,WDES1%LMBASE,WDES1%NPRO_POS,WDES1%LMMAX,WDES1%ITYP,CPROJ,CPROJ_LOCAL) &
!$ACC& PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) 
      DO NI=1,WDES1%NIONS
         NT = WDES1%ITYP(NI)
         LMMAXC=WDES1%LMMAX(NT)
         IF (LMMAXC/=0) THEN
            NPRO = WDES1%LMBASE(NI)
            NPRO_POS=WDES1%NPRO_POS(NI)
!$ACC LOOP
            DO L=1,LMMAXC
               CPROJ(L+NPRO_POS)=CPROJ_LOCAL(L+NPRO)
            ENDDO
         ENDIF
      ENDDO

      IF (WDES1%LNONCOLLINEAR) THEN

         NPRO_TOT= WDES1%NPRO_TOT/2
!$ACC PARALLEL LOOP PRESENT(WDES1,CPROJ,CPROJ_LOCAL) PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) &
!$ACC 
         DO NI=1,WDES1%NIONS
            NT=WDES1%ITYP(NI)
            LMMAXC=WDES1%LMMAX(NT)
            IF (LMMAXC/=0) THEN
               NPRO = WDES1%LMBASE(NI)+WDES1%LMBASE(WDES1%NIONS+1)
               NPRO_POS=WDES1%NPRO_POS(NI)
!$ACC LOOP
               DO L=1,LMMAXC
                  CPROJ(L+NPRO_POS+NPRO_TOT)=CPROJ_LOCAL(L+NPRO)
               ENDDO
            ENDIF
         ENDDO

      ENDIF
      CALL M_sum_z(WDES1%COMM_INB,CPROJ, WDES1%NPRO_TOT)
# 314


      END SUBROUTINE

!=======================================================================
!>  distribute plane wave part of (1._q,0._q) wavefunction over all COMM_INB nodes
!>
!>  the supplied wavefunction must have the standard serial layout and
!>  must be defined on the node COMM_INB%IONODE
!=======================================================================

      SUBROUTINE DIS_PW(WDES1,CPTWFP,CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)

      INTEGER NC,IND,I
      INTEGER NPL_TOT,NPL

      CALL M_bcast_z(WDES1%COMM_INB, CPTWFP, WDES1%NPL_TOT)

# 341


      IND=0
      DO NC=1,WDES1%NCOL
        DO I=1,WDES1%PL_COL(NC)
          CW_LOCAL(IND+I)=CPTWFP(WDES1%PL_INDEX(NC)+I)
        ENDDO
        IND=IND+WDES1%PL_COL(NC)
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         IND=0
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO NC=1,WDES1%NCOL
            DO I=1,WDES1%PL_COL(NC)
               CW_LOCAL(IND+I+NPL)=CPTWFP(WDES1%PL_INDEX(NC)+I+NPL_TOT)
            ENDDO
            IND=IND+WDES1%PL_COL(NC)
         ENDDO
      ENDIF

      IF (IND*WDES1%NRSPINORS/=WDES1%NPL) THEN
        CALL vtutor%bug("internal error in DIS_PW: " // str(IND) // " " // str(WDES1%NPL), "wave_mpi.F", 364)
      ENDIF
# 368

      END SUBROUTINE


!=======================================================================
!>  merge plane wave part of wavefunctions from all COMM_INB nodes
!>
!>  using an index array
!>  only used when reading orbitals in Gamma only
!=======================================================================

      SUBROUTINE DIS_PW_GAMMA(WDES1,CPTWFP,CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)

      INTEGER NC,IND,I
      INTEGER NPL_TOT,NPL

      CALL M_bcast_z(WDES1%COMM_INB, CPTWFP, WDES1%NPL_TOT)

# 397


      DO IND=1,WDES1%NGVECTOR
         IF (LCONJG(IND)) THEN
            CW_LOCAL(IND)=CONJG(CPTWFP(INDEX(IND)))
         ELSE
            CW_LOCAL(IND)=CPTWFP(INDEX(IND))
         ENDIF
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO IND=1,WDES1%NGVECTOR
            IF (LCONJG(IND)) THEN
               CW_LOCAL(IND+NPL)=CONJG(CPTWFP(INDEX(IND)+NPL_TOT))
            ELSE
               CW_LOCAL(IND+NPL)=CPTWFP(INDEX(IND)+NPL_TOT)
            ENDIF
         ENDDO
      ENDIF

# 421

      END SUBROUTINE


!=======================================================================
!> distribute (1._q,0._q) specific band to nodes
!>
!> only required after i.e. reading wavefunctions
!> it is sufficient if wavefunctions are defined on master node
!> but mind, that all nodes must call this communication routine
!=======================================================================

      SUBROUTINE DIS_PW_BAND(WDES1, NB, CPTWFP, CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      INTEGER ierror


! first copy wavefunction to all nodes (it would be sufficient to
!     copy it to master nodes in each inband-communicator, but
!     I am not sure who the master node is)
      CALL M_bcast_z(WDES1%COMM_KIN, CPTWFP, WDES1%NPL_TOT)
# 450


      NB_L=NB_LOCAL(NB,WDES1)
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
        CALL DIS_PW(WDES1,CPTWFP,CW_LOCAL(1,NB_L))
      ENDIF
# 458

      END SUBROUTINE

!
! new version uses an index array and works even if parallel or seriell FFT is used
!

      SUBROUTINE DIS_PW_BAND_GAMMA(WDES1, NB, CPTWFP, CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)
      INTEGER ierror


! first copy wavefunction to all nodes (it would be sufficient to
!     copy it to master nodes in each inband-communicator, but
!     I am not sure who the master node is)
      CALL M_bcast_z(WDES1%COMM_KIN, CPTWFP, WDES1%NPL_TOT)
# 484


      NB_L=NB_LOCAL(NB,WDES1)
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
         CALL DIS_PW_GAMMA(WDES1,CPTWFP,CW_LOCAL(1,NB_L), INDEX, LCONJG)
      ENDIF
# 493

      END SUBROUTINE


!=======================================================================
!>  merge plane wave part of wavefunctions from all COMM_INB nodes
!=======================================================================

      SUBROUTINE MRG_PW(WDES1,CPTWFP,CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)

      INTEGER NC,IND,I
      INTEGER NPL_TOT,NPL

      CPTWFP=0
      IND=0
      DO NC=1,WDES1%NCOL
        DO I=1,WDES1%PL_COL(NC)
          CPTWFP(WDES1%PL_INDEX(NC)+I)=CW_LOCAL(IND+I)
        ENDDO
        IND=IND+WDES1%PL_COL(NC)
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         IND=0
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO NC=1,WDES1%NCOL
            DO I=1,WDES1%PL_COL(NC)
               CPTWFP(WDES1%PL_INDEX(NC)+I+NPL_TOT)=CW_LOCAL(IND+I+NPL)
            ENDDO
            IND=IND+WDES1%PL_COL(NC)
         ENDDO
      ENDIF

      IF (IND*WDES1%NRSPINORS/=WDES1%NPL) THEN
        CALL vtutor%bug("internal error MRG_PW: " // str(IND) // " " // str(WDES1%NPL), "wave_mpi.F", 535)
      ENDIF

      CALL M_sum_z(WDES1%COMM_INB, CPTWFP(1), WDES1%NPL_TOT)
# 541

      END SUBROUTINE


!=======================================================================
!>  merge plane wave part of wavefunctions from all COMM_INB nodes
!>
!>  using an index array
!>  only used when writing orbitals in Gamma only
!=======================================================================

      SUBROUTINE MRG_PW_GAMMA(WDES1, CPTWFP, CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)

      INTEGER IND
      INTEGER NPL_TOT,NPL

      CPTWFP=0
      DO IND=1,WDES1%NGVECTOR
         IF (LCONJG(IND)) THEN
            CPTWFP(INDEX(IND))=CONJG(CW_LOCAL(IND))
         ELSE
            CPTWFP(INDEX(IND))=CW_LOCAL(IND)
         ENDIF
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO IND=1,WDES1%NGVECTOR
            IF (LCONJG(IND)) THEN
               CPTWFP(INDEX(IND)+NPL_TOT)=CONJG(CW_LOCAL(IND+NPL))
            ELSE
               CPTWFP(INDEX(IND)+NPL_TOT)=CW_LOCAL(IND+NPL)
            ENDIF
         ENDDO
      ENDIF

      CALL M_sum_z(WDES1%COMM_INB, CPTWFP(1), WDES1%NPL_TOT)
# 590

      END SUBROUTINE

!=======================================================================
!> merge  band NB from all nodes into a local copy CPTWFP
!> only required for i.e. writing wavefunctions
!=======================================================================

      SUBROUTINE MRG_PW_BAND(WDES1, NB, CPTWFP, CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)

      CPTWFP=0
      NB_L=NB_LOCAL(NB,WDES1)
! if wavefunction is located on this node merge it to CPTWFP
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
        CALL MRG_PW(WDES1,CPTWFP,CW_LOCAL(1,NB_L))
      ENDIF
! then merge over all nodes using inter-band communication
      CALL M_sum_z(WDES1%COMM_INTER, CPTWFP(1), WDES1%NPL_TOT)
# 617

      END SUBROUTINE

!
!> new version uses an index array and works even if parallel or seriell FFT is used
!

      SUBROUTINE MRG_PW_BAND_GAMMA(WDES1, NB, CPTWFP, CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)

      CPTWFP=0
      NB_L=NB_LOCAL(NB,WDES1)
! if wavefunction is located on this node merge it to CPTWFP
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
         CALL MRG_PW_GAMMA(WDES1,CPTWFP,CW_LOCAL(1,NB_L), INDEX, LCONJG)
      ENDIF
! then merge over all nodes using inter-band communication
      CALL M_sum_z(WDES1%COMM_INTER, CPTWFP(1), WDES1%NPL_TOT)
# 647

      END SUBROUTINE

!=======================================================================
!> merge eigenvalues from all nodes
!=======================================================================

      SUBROUTINE MRG_CEL(WDES,W)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL

!! 

!
! first (0._q,0._q) out all components from other nodes
!
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%CELTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
         CALL SETWDES(WDES,WDES1,0)
         band: DO NB_GLOBAL=1,WDES%NB_TOT
            NB=NB_LOCAL(NB_GLOBAL,WDES1)
            IF (NB==0) THEN
               W%CELTOT(NB_GLOBAL,NK,I)=0
            ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin
!
! then merge CELTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALL M_sum_z( WDES%COMM_KINTER, W%CELTOT, NCEL)
      CALL M_sum_z( WDES%COMM_INTER,  W%CELTOT, NCEL)

! finally set eigenvalues of states beyond WDES%NB_TOTK to large values
      CALL WFSET_HIGH_CELEN(W)

!! 

      END SUBROUTINE


      SUBROUTINE MRG_CELTOT(WDES,W)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
!$ACC ROUTINE(NB_LOCAL) SEQ

!
! first (0._q,0._q) out all components from other nodes
!
!$ACC ENTER DATA CREATE(WDES1) 
      CALL SETWDES(WDES,WDES1,0)

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,WDES1,W) 
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!$ACC LOOP VECTOR
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%CELTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
!$ACC LOOP VECTOR
         band: DO NB_GLOBAL=1,WDES%NB_TOT
            NB=NB_LOCAL(NB_GLOBAL,WDES1)
            IF (NB==0) THEN
               W%CELTOT(NB_GLOBAL,NK,I)=0
            ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin

# 746

!
! then merge CELTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALL M_sum_z( WDES%COMM_KINTER, W%CELTOT(1,1,1), NCEL)
      CALL M_sum_z( WDES%COMM_INTER,  W%CELTOT(1,1,1), NCEL)


! finally set eigenvalues of states beyond WDES%NB_TOTK to large values
      CALL WFSET_HIGH_CELEN(W)

      END SUBROUTINE


!=======================================================================
!> merge partial occupancies from all nodes
!=======================================================================

      SUBROUTINE MRG_FERWE(WDES,W)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
!$ACC ROUTINE(NB_LOCAL) SEQ

!
! first (0._q,0._q) out all components from other nodes
!
!$ACC ENTER DATA CREATE(WDES1) 
      CALL SETWDES(WDES,WDES1,0)

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,WDES1,W) 
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!$ACC LOOP VECTOR
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%FERTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
!$ACC LOOP VECTOR
         band: DO NB_GLOBAL=1,WDES%NB_TOT
           NB=NB_LOCAL(NB_GLOBAL,WDES1)
           IF (NB==0) THEN
             W%FERTOT(NB_GLOBAL,NK,I)=0
           ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin

# 807

!
! then merge FERTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALL M_sum_d( WDES%COMM_KINTER, W%FERTOT(1,1,1), NCEL)
      CALL M_sum_d( WDES%COMM_INTER,  W%FERTOT(1,1,1), NCEL)

      END SUBROUTINE


!=======================================================================
!> merge auxilary array from all nodes
!=======================================================================

      SUBROUTINE MRG_AUX(WDES,W)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
!$ACC ROUTINE(NB_LOCAL) SEQ

!
! first (0._q,0._q) out all components from other nodes
!
!$ACC ENTER DATA CREATE(WDES1) 
      CALL SETWDES(WDES,WDES1,0)

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,WDES1,W) 
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!$ACC LOOP VECTOR
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%AUXTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
!$ACC LOOP VECTOR
         band: DO NB_GLOBAL=1,WDES%NB_TOT
           NB=NB_LOCAL(NB_GLOBAL,WDES1)
           IF (NB==0) THEN
             W%AUXTOT(NB_GLOBAL,NK,I)=0
           ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin

# 864

!
! then merge AUXTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALL M_sum_d( WDES%COMM_KINTER, W%AUXTOT(1,1,1), NCEL)
      CALL M_sum_d( WDES%COMM_INTER,  W%AUXTOT(1,1,1), NCEL)

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW **************************
!
!> redistribute plane wave coefficients from over band to over
!> plane wave coefficient or vice versa
!>
!> this operation is 1._q in place to reduce storage demands
!>
!> mind that if the routine is called twice the original distribution
!> is obtained
!
!***********************************************************************

      SUBROUTINE REDIS_PW(WDES1, NBANDS, CPTWFP)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NBANDS
      COMPLEX(q) :: CPTWFP(WDES1%NRPLWV*NBANDS)
# 902

      COMPLEX(q), ALLOCATABLE :: CWORK(:)

! local variables
      INTEGER :: NRPLWV,N,NB,INFO
# 910


      

! quick return if possible
      IF (WDES1%COMM_INTER%NCPU ==1) THEN
         
         RETURN
      ENDIF

      NRPLWV=WDES1%NRPLWV

      IF (MOD(NRPLWV,WDES1%NB_PAR) /= 0) THEN
        CALL vtutor%bug("REDIS_PW: internal error(1) " // str(NRPLWV) // " " // str(WDES1%NB_PAR), "wave_mpi.F", 923)
      ENDIF
# 929

      ALLOCATE( CWORK(NRPLWV))
! CWORK=(0.0_q, 0.0_q) ! No test failures reported.
!$ACC DATA CREATE(CWORK) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)

      DO NB=0,NBANDS-1
         CALL M_alltoall_z(WDES1%COMM_INTER, NRPLWV, CPTWFP(NB*NRPLWV+1), CWORK(1))
         CALL ZCOPY( NRPLWV, CWORK(1), 1, CPTWFP(NB*NRPLWV+1), 1)
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_ACTIVE)
!$ACC END DATA
      DEALLOCATE(CWORK)

      

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_OVER_BANDS ***************
!
!> redistribute all coefficients
!>
!> if and only if LOVER_BAND
!> is .TRUE.
!> after this call the bands are distributed over nodes
!> and LOVER_BAND is .FALSE.
!
!***********************************************************************

      SUBROUTINE REDIS_PW_OVER_BANDS(WDES, W)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W


      

      IF ( W%OVER_BAND ) THEN
         CALL REDIS_PW_ALL(WDES, W)
      ENDIF

      

      END SUBROUTINE REDIS_PW_OVER_BANDS

!************************ SUBROUTINE REDIS_PW_ALL **********************
!
!> redistribute all plane wave coeffcients,
!> change the order W%OVER_BAND
!> calling the routine twice restores the original order
!
!***********************************************************************

      SUBROUTINE REDIS_PW_ALL(WDES, W)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W


      TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
      INTEGER NCPU, ISP, NK

      
     
      NCPU   =WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
      IF (NCPU /= 1) THEN
         DO ISP=1,WDES%ISPIN
         DO NK=1,WDES%NKPTS
            IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
            CALL SETWDES(WDES,WDES1,NK)
            CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CPTWFP   (1,1,NK,ISP))
         ENDDO
         ENDDO

         W%OVER_BAND=.NOT. W%OVER_BAND
      ENDIF

      

      END SUBROUTINE REDIS_PW_ALL

!************************ SUBROUTINE REDIS_PW_START ********************
!
!> redistribute plane wave coefficients
!> from over band to over
!> plane wave coefficient or vice versa asynchronously
!> (1._q,0._q) band is initiated
!
!***********************************************************************

      SUBROUTINE REDIS_PW_START(WDES, CPTWFP, BANDINDEX, H)

!! USE moffload

      USE prec
      USE wave_mpi
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)  WDES
      INTEGER BANDINDEX
      COMPLEX(q) :: CPTWFP(WDES%NRPLWV)
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: NRPLWV
      INTEGER :: IND,N


      

      IF (BANDINDEX>H%NBANDS) &
         CALL vtutor%bug("REDIS_PW_START: non existing band " // str(BANDINDEX), "wave_mpi.F", 1046)

      NRPLWV=WDES%NRPLWV

      IND=H%NEXT
      H%NEXT=H%NEXT+1 ; IF (H%NEXT > H%NB) H%NEXT=1

      IF (H%BAND(IND)/=0) &
         CALL vtutor%bug("REDIS_PW_START: slot NEXT not empty " // str(H%NEXT), "wave_mpi.F", 1054)

! fill in band index
      H%BAND(IND)=BANDINDEX
# 1063

      CALL M_alltoall_d_async(H%COMM, 2*NRPLWV, CPTWFP(1), H%CPTWFP(1,IND), &
           H%ICOMM+IND, H%SREQUEST(1,IND), H%RREQUEST(1,IND) )

      

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_STRIP_START **************
!
!> start asynchrpnous redistribution plane wave coefficients from
!> over-band to over-plane-wave coefficient or vice versa,
!> for NSTRIP bands (starting at index NAMDINDEX)
!
!***********************************************************************

      SUBROUTINE REDIS_PW_STRIP_START(WDES, CPTWFP, NSTRIP, BANDINDEX, H)
# 1083

      USE prec
      USE wave_mpi
      USE wave
      IMPLICIT NONE

      TYPE (wavedes) :: WDES
      COMPLEX(q) :: CPTWFP(WDES%NRPLWV, NSTRIP)
      INTEGER :: NSTRIP, BANDINDEX
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: IND

# 1098



      

# 1110


      DO IND = 1, NSTRIP
         CALL REDIS_PW_START(WDES, CPTWFP(1,IND), BANDINDEX+IND-1, H)
      ENDDO

# 1123


      

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_STOP  ********************
!
!> finish redistribution of  plane wave coefficients
!> from over band to over plane wave coefficient or vice versa.
!
!***********************************************************************

      SUBROUTINE REDIS_PW_STOP(WDES, CPTWFP, BANDINDEX, H)

!! USE moffload_struct_def

      USE prec
      USE wave_mpi
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)  WDES
      INTEGER BANDINDEX
      COMPLEX(q) :: CPTWFP(WDES%NRPLWV)
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: NRPLWV
      INTEGER :: IND,N


      

      IF (BANDINDEX>H%NBANDS) &
         CALL vtutor%bug("REDIS_PW_STOP: non existing band " // str(BANDINDEX), "wave_mpi.F", 1157)

      NRPLWV=WDES%NRPLWV

! usually the required band is in the slot NEXT, however
! possibly we have to search for it

      IND=H%NEXT
      IF (H%BAND(IND) /= BANDINDEX) THEN
         DO IND=1, H%NB
            IF ( H%BAND(IND) == BANDINDEX) EXIT
         ENDDO
         IND=MIN(IND,H%NB)
         IF  (H%BAND(IND) /= BANDINDEX) &
            CALL vtutor%bug("REDIS_PW_STOP: band " // str(BANDINDEX) // " not found in " // str(H%BAND), "wave_mpi.F", 1171)
      ENDIF
! fill in (0._q,0._q) into the band index array
      H%BAND(IND)=0
# 1179

! wait for send and receive to finish
      CALL M_alltoall_wait(H%COMM, H%SREQUEST(1,IND), H%RREQUEST(1,IND) )

!$ACC PARALLEL LOOP PRESENT(CPTWFP,H%CPTWFP) 
      DO N=1,NRPLWV
         CPTWFP(N)=H%CPTWFP(N,IND)
      ENDDO

      

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_STRIP_STOP ***************
!
!> finish redistribution of plane wave coefficients
!> from over-band to over-plane-wave coefficient or vice versa,
!> for NSTRIP bands (starting at index BANDINDEX)
!
!***********************************************************************

      SUBROUTINE REDIS_PW_STRIP_STOP(WDES, CPTWFP, NSTRIP, BANDINDEX, H)

!! USE moffload_struct_def

      USE prec
      USE wave_mpi
      USE wave

      IMPLICIT NONE

      TYPE (wavedes) :: WDES
      COMPLEX(q) :: CPTWFP(WDES%NRPLWV, NSTRIP)
      INTEGER :: NSTRIP, BANDINDEX
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: IND, IB, NRPLWV, N, ISLOT


      

      ISLOT = H%NEXT - NSTRIP
      IF (ISLOT <= 0) ISLOT = ISLOT + H%NB
! usually the first required band will be in the slot ISLOT
! ... if not we have to search for it
      IF (H%BAND(ISLOT) /= BANDINDEX) THEN
         DO ISLOT = 1, H%NB
            IF (H%BAND(ISLOT) == BANDINDEX) EXIT
         ENDDO
! when BANDINDEX was not found in H%BAND, then ISLOT = H%NB+1
         IF (ISLOT > H%NB) &
            CALL vtutor%bug("REDIS_PW_STRIP_STOP: band " // str(BANDINDEX) // " not found in " // str(H%BAND), "wave_mpi.F", 1230)
      ENDIF

# 1239


      NRPLWV = WDES%NRPLWV

      DO IND = 1, NSTRIP

         IB = BANDINDEX+IND-1

! the required band should be in the slot ISLOT
         IF (H%BAND(ISLOT) /= IB) &
            CALL vtutor%bug("REDIS_PW_STRIP_STOP: band " // str(IB) // " not in slot " // str(ISLOT),  "wave_mpi.F", 1249)

# 1255

# 1260

! wait for send and receive to finish
         CALL M_alltoall_wait(H%COMM, H%SREQUEST(1,ISLOT), H%RREQUEST(1,ISLOT))

!$ACC PARALLEL LOOP PRESENT(CPTWFP,H%CPTWFP) 
         DO N = 1, NRPLWV
            CPTWFP(N,IND) = H%CPTWFP(N,ISLOT)
         ENDDO


! set entry in band index array to (0._q,0._q) to indicate
! this slot is now empty
         H%BAND(ISLOT) = 0

! the next band should be in the following slot
         ISLOT = ISLOT + 1 ; IF (ISLOT > H%NB) ISLOT = 1
      ENDDO

! Set H%NEXT to the first slot that was freed above
      ISLOT = ISLOT - NSTRIP ; IF (ISLOT <= 0) ISLOT = ISLOT + H%NB
      H%NEXT = ISLOT

      

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PROJ ************************
!
!> redistribute projector part of wave function from over band to over
!> plane wave coefficient or vice versa
!>
!> this operation is 1._q in place to reduce storage demands
!>
!> mind that if the routine is called twice the original distribution
!> is obtained
!
!***********************************************************************

      SUBROUTINE REDIS_PROJ(WDES1, NBANDS, CPROJ)

!! USE moffload

      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NBANDS
      COMPLEX(q) :: CPROJ(WDES1%NPROD*NBANDS)
# 1314

      COMPLEX(q), ALLOCATABLE :: CWORK(:)


# 1320

      INTEGER, PARAMETER :: MCOMP=2

! local variables

      INTEGER :: NPROD,N,NB,INFO
# 1329

      

      IF (WDES1%COMM_INTER%NCPU==1 .OR. WDES1%NPROD==0) THEN
         
         RETURN
      ENDIF

      NPROD=WDES1%NPROD

      IF (MOD(NPROD,WDES1%NB_PAR) /= 0) THEN
        CALL vtutor%bug("REDIS_PW: NPROD not dividable by NB_PAR:" // str(NPROD) // " " // str(WDES1%NB_PAR), "wave_mpi.F", 1340)
      ENDIF
# 1346

      ALLOCATE( CWORK(NPROD))
!$ACC ENTER DATA CREATE(CWORK) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
! CWORK=0 ! No test failures reported.

      DO NB=0,NBANDS-1
# 1354

         CALL M_alltoall_z(WDES1%COMM_INTER, NPROD, CPROJ(NB*NPROD+1), CWORK(1))

         CALL ZCOPY( NPROD, CWORK(1), 1, CPROJ(NB*NPROD+1), 1)
      ENDDO

!$ACC EXIT DATA DELETE(CWORK) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
      DEALLOCATE(CWORK)

      

      END SUBROUTINE


!************************ SUBROUTINE SET_NPL_NPRO **********************
!
!> set the local number of plane waves on a node after redistribution
!> of plane wave coefficients and projected wavefunctions
!> on entry the global number of plane waves must be given
!
!***********************************************************************

      SUBROUTINE SET_NPL_NPRO(WDES1, NPL, NPRO)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NPL,NPRO,NPL_REM,NPRO_REM,I
! local variables
      NPL_REM =NPL
      NPRO_REM=NPRO

      NPL =WDES1%NRPLWV/WDES1%COMM_INTER%NCPU
      NPRO=WDES1%NPROD /WDES1%COMM_INTER%NCPU
      DO I=1,WDES1%COMM_INTER%NODE_ME
        NPL =MIN(NPL,NPL_REM)
        NPRO=MIN(NPRO,NPRO_REM)
        NPL_REM  =NPL_REM -NPL
        NPRO_REM =NPRO_REM-NPRO
      ENDDO
# 1398


      END SUBROUTINE

!***********************************************************************
!
!> Copy wavefunction coefficients from each k-point communicator to
!> replicated copies in other communicators
!> Should use a gather operation
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      INTEGER :: ISP,K,inode


      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               inode=MOD(K-1,WDES%COMM_KINTER%NCPU)+1
! the product of the number of bands and number of orbitals is possible exceeding INTEGER (KIND=4)
               CALL M_bcast_z8_from(WDES%COMM_KINTER, W%CPTWFP(1,1,K,ISP),SIZE(W%CPTWFP,1,KIND=qi8)*SIZE(W%CPTWFP,2,KIND=qi8),inode)
# 1429

               CALL M_bcast_z8_from(WDES%COMM_KINTER, W%CPROJ(1,1,K,ISP),SIZE(W%CPROJ,1,KIND=qi8)*SIZE(W%CPROJ,2,KIND=qi8),inode)

            END DO
         END DO
      ENDIF

    END SUBROUTINE KPAR_SYNC_WAVEFUNCTIONS

!***********************************************************************
!
!> Copy eigenvalues from each k-point communicator to
!> replicated copies in other communicators
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_CELTOT(WDES,W)
! Not the same as a MRG_CEL since we do not sum over within bands
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
      INTEGER :: ISP,K


      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN

!!    

         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                  W%CELTOT(:,K,ISP)=0.
               END IF
            END DO
         END DO
! Reduce and broadcast for simplicity. Really an allgather.
         CALL M_sum_master_z( WDES%COMM_KINTER, W%CELTOT, SIZE(W%CELTOT) )
         CALL M_bcast_z( WDES%COMM_KINTER, W%CELTOT, SIZE(W%CELTOT) )

!!    

      ENDIF

    END SUBROUTINE KPAR_SYNC_CELTOT

!***********************************************************************
!
!> Copy occupancies from each k-point communicator to
!> replicated copies in other communicators
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_FERTOT(WDES,W)
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
      INTEGER :: ISP,K


      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN

!!    

         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                  W%FERTOT(:,K,ISP)=0.
               END IF
            END DO
         END DO
! Reduce and broadcast for simplicity. Really an allgather.
         CALL M_sum_master_d( WDES%COMM_KINTER, W%FERTOT, SIZE(W%FERTOT) )
         CALL M_bcast_d( WDES%COMM_KINTER, W%FERTOT, SIZE(W%FERTOT) )

!!    

      ENDIF

    END SUBROUTINE KPAR_SYNC_FERTOT

!***********************************************************************
!
!> Copy auxilary array from each k-point communicator to
!> replicated copies in other communicators
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_AUXTOT(WDES,W)
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
      INTEGER :: ISP,K


      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN

!!    

         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                  W%AUXTOT(:,K,ISP)=0.
               END IF
            END DO
         END DO
! Reduce and broadcast for simplicity. Really an allgather.
         CALL M_sum_master_d( WDES%COMM_KINTER, W%AUXTOT, SIZE(W%AUXTOT) )
         CALL M_bcast_d( WDES%COMM_KINTER, W%AUXTOT, SIZE(W%AUXTOT) )

!!    

      ENDIF

    END SUBROUTINE KPAR_SYNC_AUXTOT

!***********************************************************************
!
!> Sync orbitals, eigenvalues and occupancies
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_ALL(WDES,W)
      USE wave
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
      

      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         CALL KPAR_SYNC_FERTOT(WDES,W)
         CALL KPAR_SYNC_CELTOT(WDES,W)
         CALL KPAR_SYNC_AUXTOT(WDES,W)
         CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
      ENDIF

      
    END SUBROUTINE KPAR_SYNC_ALL
