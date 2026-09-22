# 1 "pawsym.F"
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


# 2 "pawsym.F" 2 
!***********************************************************************
!
!> PAW symmetry routines
! mostly written by gK (non linear case written by Martijn Marsmann)
!
!***********************************************************************

  MODULE pawsym
    CONTAINS
!***********************************************************************
!
!>   Routine AUGSYM symmetrizes arrays like RHOLM or DLM
!>   or other arrays with a similar structure.
!>
!>   The procedure is a quite
!>   simple: Rotate the input array and add it to the array at the
!>   rotated atomic position (for all space group operations ...).
!>   Finally divide the result by the number of symmetry operations!
!>
!>
!>   Input parameters:
!>   -----------------
!>   see routine FSYM in symlib.F
!>
!>      MAP(NAX,1,48,NP) contains a table connecting the rotated and
!>                   the unrotated positions (stored is the index of
!>                   the rotated position). \n
!>      S(3,3,48) contains the INTEGER rotation matrices. \n
!>      NR contains the number of given space group operations. \n
!>      NP contains the number of "primitive" translations in the cell. \n
!>
!>      NTYP is the number of atomic species. \n
!>      NITYP(NSP) contains the number of atoms per species.
!>
!>      A(:,I)  contains the direct (real space) lattice vectors. \n
!>      B(:,I)  contains the reciprocal lattice vectors.
!>
!>
!>   Output parameters:
!>   ------------------
!>
!>      Array MAT contains the symmetrized matrix on output (input is
!>      overwritten!).
!                                                                      *
!                                                                      *
!***********************************************************************


      SUBROUTINE AUGSYM(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
            NIOND,NR,NP,MAP,MAGROT,S,A,B,ISP)
      USE prec
      USE pseudo
      USE asa
      IMPLICIT NONE

! parameters required from the symmetry package
      INTEGER NIOND      ! dimension for number of ions
      INTEGER NIONS      ! number of ions
      INTEGER NP         ! number of primitive translations
      INTEGER NR         ! number of rotations
      INTEGER MAP(NIOND,NR,NP)
      INTEGER NTYP       ! number of species
      INTEGER NITYP(NTYP)! number of atoms per species
      INTEGER ISP        ! spin component
      INTEGER S(3,3,48)
      REAL(q) A(3,3),B(3,3),MAGROT(48,NP)
!
      INTEGER LMDIM      ! first dimension of MAT
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS)
      TYPE (potcar) P(NTYP)

! local varibales
      COMPLEX(q),ALLOCATABLE :: TMP(:,:,:),ROTMAT(:,:,:)
      REAL(q),ALLOCATABLE :: SL(:,:,:)
      INTEGER NROT,LMAX,LDIM,ITRANS,IA,IAP,MMAX,IROT,IS,ISTART
      INTEGER,EXTERNAL :: MAXL,MAXL1
      REAL(q) SCALE
      INTEGER L,LP,NI
!-------------------------------------------------------------------
! allocate work arrays
!-------------------------------------------------------------------
# 90


      LDIM=MAXL(NTYP,P  )         ! maximum l quantum number
      MMAX=(2*LDIM+1)             ! maximum m quantum number
      ALLOCATE (TMP(LMDIM,LMDIM,NIONS),SL(MMAX,MMAX,0:LDIM), &
                ROTMAT(LMDIM,LMDIM,NIONS))

      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
! loop over all species
      DO IS=1,NTYP
        LMAX=MAXL1(P(IS))      ! maximum l for this species
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
! loop over all rotations
        DO IROT=1,NR
! setup rotation matrices for L=0,...,LMAX
          CALL SETUP_SYM_LL(MMAX,LMAX,S(1,1,IROT),SL,A,B)
! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
! rotate the matrix and store result in ROTMAT
            CALL ROTATE_MATRIX(LMDIM,MAT(1,1,IA),ROTMAT(1,1,IA),MMAX,LMAX,SL,P(IS))
          ENDDO

! loop over all space group operations (translations+ rotations)
          DO IA=ISTART,NITYP(IS)+ISTART-1
            DO ITRANS=1,NP
! destination atom
               IAP=MAP(IA,IROT,ITRANS)
               SCALE=1._q
               IF (ISP==2) SCALE=MAGROT(IROT,ITRANS)
               TMP(:,:,IA)=TMP(:,:,IA)+ROTMAT(:,:,IAP)*SCALE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
! divide final result by the number of translations and rotations
      SCALE=1._q/(NP*NR)
      MAT=TMP*SCALE

      DEALLOCATE(TMP,ROTMAT,SL)

# 141

      END SUBROUTINE


!************** SUBROUTINE AUGSYM_NONCOL *******************************
!
!>   Routine AUGSYM_NONCOL symmetrizes arrays like RHOLM or DLM
!>   or other arrays with a similar structure.
!>
!>   Input parameters:
!>   -----------------
!>   see routine FSYM in symlib.F
!>
!>      MAP(NAX,1,48,NP) contains a table connecting the rotated and
!>                   the unrotated positions (stored is the index of
!>                   the rotated position). \n
!>      S(3,3,48) contains the INTEGER rotation matrices. \n
!>      NR contains the number of given space group operations. \n
!>      NP contains the number of "primitive" translations in the cell. \n
!>
!>      NTYP is the number of atomic species. \n
!>      NITYP(NSP) contains the number of atoms per species.
!>
!>      A(:,I)  contains the direct (real space) lattice vectors. \n
!>      B(:,I)  contains the reciprocal lattice vectors.
!>
!>
!>   Output parameters:
!>   ------------------
!>
!>      Array MAT contains the symmetrized matrix on output (input is
!>      overwritten!).
!
!
!***********************************************************************

      SUBROUTINE AUGSYM_NONCOL(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
     &      NIOND,NR,NP,MAP,MAGROT,SAXIS,S,INVMAP,A,B)
      USE prec
      USE pseudo
      USE relativistic
      USE paw
      USE asa
      IMPLICIT NONE

! parameters required from the symmetry package
      INTEGER NIOND      ! dimension for number of ions
      INTEGER NIONS      ! number of ions
      INTEGER NP         ! number of primitive translations
      INTEGER NR         ! number of rotations
      INTEGER MAP(NIOND,NR,NP)
      INTEGER NTYP       ! number of species
      INTEGER NITYP(NTYP)! number of atoms per species
      INTEGER ISP        ! spin component
      INTEGER S(3,3,48),INVMAP(48),I
      REAL(q) A(3,3),B(3,3),MAGROT(48,NP)
!
      INTEGER LMDIM      ! first dimension of MAT
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS,3)
      TYPE (potcar) P(NTYP)

! local varibales
      COMPLEX(q),ALLOCATABLE :: TMP(:,:,:,:),ROTMAT(:,:,:,:)
      COMPLEX(q),ALLOCATABLE :: ROTMAT_TEMP(:,:,:),TROTMAT(:,:,:)
      REAL(q),ALLOCATABLE :: SL(:,:,:)
      INTEGER NROT,LMAX,LDIM,ITRANS,IA,IAP,MMAX,IROT,IROTI,IS,ISTART,IDIR
      INTEGER,EXTERNAL :: MAXL,MAXL1
      REAL(q) SCALE,SAXIS(3),ALPHA,BETA
      INTEGER L,LP,NI

      INTEGER J,NROTK,DET,TMPM(3,3)

!-------------------------------------------------------------------
! allocate work arrays
!-------------------------------------------------------------------
# 225


      LDIM=MAXL(NTYP,P  )         ! maximum l quantum number
      MMAX=(2*LDIM+1)             ! maximum m quantum number
      ALLOCATE (TMP(LMDIM,LMDIM,NIONS,3),SL(MMAX,MMAX,0:LDIM), &
           ROTMAT(LMDIM,LMDIM,NIONS,3),ROTMAT_TEMP(LMDIM,LMDIM,3), &
           TROTMAT(LMDIM,LMDIM,3))

      CALL EULER(SAXIS,ALPHA,BETA)

      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
! loop over all species
      DO IS=1,NTYP
        LMAX=MAXL1(P(IS))      ! maximum l for this species
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
! loop over all rotations
        DO IROT=1,NR
! setup rotation matrices for L=0,...,LMAX
          CALL SETUP_SYM_LL(MMAX,LMAX,S(1,1,IROT),SL,A,B)
! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
! rotate the matrix and store result in ROTMAT
          DO IDIR=1,3
             CALL ROTATE_MATRIX(LMDIM,MAT(1,1,IA,IDIR),ROTMAT(1,1,IA,IDIR),MMAX,LMAX,SL,P(IS))
          ENDDO
          ENDDO
! loop over all space group operations (translations+ rotations)
          DO IA=ISTART,NITYP(IS)+ISTART-1
            DO ITRANS=1,NP
! destination atom
               IAP=MAP(IA,IROT,ITRANS)
               SCALE=1._q
! Transform from "SAXIS basis" to the system of cartesian axes
               ROTMAT_TEMP(:,:,1)=COS(BETA)*COS(ALPHA)*ROTMAT(:,:,IAP,1)- &
                    SIN(ALPHA)*ROTMAT(:,:,IAP,2)+ &
                    SIN(BETA)*COS(ALPHA)*ROTMAT(:,:,IAP,3)
               ROTMAT_TEMP(:,:,2)=COS(BETA)*SIN(ALPHA)*ROTMAT(:,:,IAP,1)+ &
                    COS(ALPHA)*ROTMAT(:,:,IAP,2)+ &
                    SIN(BETA)*SIN(ALPHA)*ROTMAT(:,:,IAP,3)
               ROTMAT_TEMP(:,:,3)=-SIN(BETA)*ROTMAT(:,:,IAP,1)+ &
                    COS(BETA)*ROTMAT(:,:,IAP,3)
! Bring to direct coordinates in which the integer rotation matrices are defined
               CALL MAT_KARDIR(LMDIM,ROTMAT_TEMP,B)

               TMPM=S(:,:,IROT)
               DET=TMPM(1,1)*TMPM(2,2)*TMPM(3,3) - &
                   TMPM(1,1)*TMPM(2,3)*TMPM(3,2) + &
                   TMPM(1,2)*TMPM(2,3)*TMPM(3,1) - &
                   TMPM(1,2)*TMPM(2,1)*TMPM(3,3) + &
                   TMPM(1,3)*TMPM(2,1)*TMPM(3,2) - &
                   TMPM(1,3)*TMPM(2,2)*TMPM(3,1)
               IF (DET<0) TMPM=-TMPM

               TROTMAT(:,:,1)= &
                    ROTMAT_TEMP(:,:,1)*TMPM(1,1)+ &
                    ROTMAT_TEMP(:,:,2)*TMPM(2,1)+ &
                    ROTMAT_TEMP(:,:,3)*TMPM(3,1)
               
               TROTMAT(:,:,2)= &
                    ROTMAT_TEMP(:,:,1)*TMPM(1,2)+ &
                    ROTMAT_TEMP(:,:,2)*TMPM(2,2)+ &
                    ROTMAT_TEMP(:,:,3)*TMPM(3,2)
               
               TROTMAT(:,:,3)= &
                    ROTMAT_TEMP(:,:,1)*TMPM(1,3)+ &
                    ROTMAT_TEMP(:,:,2)*TMPM(2,3)+ &
                    ROTMAT_TEMP(:,:,3)*TMPM(3,3)
# 309

! bring TROTMAT  to cartesian coordinates
               CALL MAT_DIRKAR(LMDIM,TROTMAT,A)
! And back to SAXIS representation
               ROTMAT_TEMP(:,:,1)=COS(BETA)*COS(ALPHA)*TROTMAT(:,:,1)+ &
                    COS(BETA)*SIN(ALPHA)*TROTMAT(:,:,2)- &
                    SIN(BETA)*TROTMAT(:,:,3)
               ROTMAT_TEMP(:,:,2)=-SIN(ALPHA)*TROTMAT(:,:,1)+ &
                    COS(ALPHA)*TROTMAT(:,:,2)
               ROTMAT_TEMP(:,:,3)=SIN(BETA)*COS(ALPHA)*TROTMAT(:,:,1)+ &
                    SIN(BETA)*SIN(ALPHA)*TROTMAT(:,:,2)+ &
                    COS(BETA)*TROTMAT(:,:,3)
! And then we sum
               TMP(:,:,IA,:)=TMP(:,:,IA,:)+ROTMAT_TEMP(:,:,:)*SCALE
# 331

            ENDDO
          ENDDO
        ENDDO
      ENDDO
! divide final result by the number of translations and rotations
      SCALE=1._q/(NP*NR)
      MAT=TMP*SCALE

      DEALLOCATE(TMP,ROTMAT,SL,ROTMAT_TEMP,TROTMAT)

# 351

      END SUBROUTINE


!**************** SUBROUTINE MAT_KARDIR ********************************
!
!> Transform a set of vectors from cartesian coordinates to
!> - direct lattice      (BASIS must be equal to B reciprocal lattice)
!> - reciprocal lattice  (BASIS must be equal to A direct lattice)
!
!***********************************************************************

      SUBROUTINE MAT_KARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(NMAX,NMAX,3),V1,V2,V3
      DIMENSION  BASIS(3,3)

      DO N=1,NMAX
      DO M=1,NMAX
        V1=V(N,M,1)*BASIS(1,1)+V(N,M,2)*BASIS(2,1)+V(N,M,3)*BASIS(3,1)
        V2=V(N,M,1)*BASIS(1,2)+V(N,M,2)*BASIS(2,2)+V(N,M,3)*BASIS(3,2)
        V3=V(N,M,1)*BASIS(1,3)+V(N,M,2)*BASIS(2,3)+V(N,M,3)*BASIS(3,3)
        V(N,M,1)=V1
        V(N,M,2)=V2
        V(N,M,3)=V3
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE


!**************** SUBROUTINE MAT_DIRKAR ********************************
!> Transform a set of vectors from
!> - direct lattice      (BASIS must be equal to A direct lattice)
!> - reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
!> to cartesian coordinates.
!***********************************************************************

      SUBROUTINE MAT_DIRKAR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(NMAX,NMAX,3),V1,V2,V3
      DIMENSION  BASIS(3,3)
      
      DO N=1,NMAX
      DO M=1,NMAX
        V1=V(N,M,1)*BASIS(1,1)+V(N,M,2)*BASIS(1,2)+V(N,M,3)*BASIS(1,3)
        V2=V(N,M,1)*BASIS(2,1)+V(N,M,2)*BASIS(2,2)+V(N,M,3)*BASIS(2,3)
        V3=V(N,M,1)*BASIS(3,1)+V(N,M,2)*BASIS(3,2)+V(N,M,3)*BASIS(3,3)
        V(N,M,1)=V1
        V(N,M,2)=V2
        V(N,M,3)=V3
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!-------------------------------------------------------------------
!> Subroutine to rotate (1._q,0._q) matrix MAT and store result in ROTMAT.
!-------------------------------------------------------------------

      SUBROUTINE ROTATE_MATRIX(LMDIM,MAT,ROTMAT,MMAX,LMAX,SL,P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER LMDIM,MMAX,LMAX

      COMPLEX(q) :: MAT(LMDIM,LMDIM)    ! initial matrix
      COMPLEX(q) :: ROTMAT(LMDIM,LMDIM) ! final matrix
      REAL(q) :: SL(MMAX,MMAX,0:LMAX)     ! rotation matrix (allways symmetric)
      TYPE (potcar) P
! local variables
      INTEGER CHANNEL,CHANNELS,IND,L,M,MP

      COMPLEX(q) :: TMP(LMDIM,LMDIM)
# 431


      CHANNELS=P%LMAX
! left hand transformation
      IND=0
      TMP=0

      DO CHANNEL=1,CHANNELS
! l-qantum number of this channel
        L=P%LPS(CHANNEL)
! rotate this l-block
        DO M=1,(2*L+1)
        DO MP=1,(2*L+1)
          TMP(IND+M,:)=TMP(IND+M,:)+SL(M,MP,L)*MAT(IND+MP,:)
        ENDDO
        ENDDO

        IND=IND+(2*L+1)
      ENDDO
! right hand transformation
      IND=0
      ROTMAT=0

      DO CHANNEL=1,CHANNELS
! l-qantum number of this channel
        L=P%LPS(CHANNEL)
! rotate this l-block
        DO M=1,(2*L+1)
        DO MP=1,(2*L+1)
          ROTMAT(:,IND+M)=ROTMAT(:,IND+M)+SL(M,MP,L)*TMP(:,IND+MP)
        ENDDO
        ENDDO

        IND=IND+(2*L+1)
      ENDDO
      END SUBROUTINE

   END MODULE pawsym


!************************************************************************
!
!> (Non explicit) interface for AUGSYM.
!
!************************************************************************
      SUBROUTINE AUGSYM_(P,LMDIM,NIONS,NIOND,NTYP,NITYP,MAT,ROTMAP,MAGROT,A,B,ISP)
      USE prec
      USE pawsym
      USE pseudo
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (potcar) P(NTYP)
      INTEGER LMDIM,NIONS,NIOND,NTYP,NITYP(NTYP)
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS)
      REAL(q) MAGROT(48,NPCELL)
      INTEGER ROTMAP(NIOND,NROTK,NPCELL)
      REAL(q) A(3,3),B(3,3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      CALL AUGSYM(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
            NIOND,NROTK,NPCELL,ROTMAP,MAGROT,ISYMOP,A,B,ISP)

      END SUBROUTINE
      
!************************************************************************
!
!> (Non explicit) interface for AUGSYM_NONCOL.
!
!************************************************************************
      SUBROUTINE AUGSYM_NONCOL_(P,LMDIM,NIONS,NIOND,NTYP,NITYP,MAT,ROTMAP,MAGROT,SAXIS,A,B)
      USE prec
      USE pawsym
      USE pseudo
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (potcar) P(NTYP)
      INTEGER LMDIM,NIONS,NIOND,NTYP,NITYP(NTYP)
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS,3)
      REAL(q) MAGROT(48,NPCELL),SAXIS(3)
      INTEGER ROTMAP(NIOND,NROTK,NPCELL)
      REAL(q) A(3,3),B(3,3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      CALL AUGSYM_NONCOL(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
     &      NIOND,NROTK,NPCELL,ROTMAP,MAGROT,SAXIS,ISYMOP,INVMAP,A,B)

      END SUBROUTINE
