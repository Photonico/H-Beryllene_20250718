#include "symbol.inc"
!=======================================================================
! RCS:  $Id: scala.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!=======================================================================
!> Module containing wrapper for scaLAPACK
!>
!> written by Gilles de Wijs (gD) and Georg Kresse (gK),
!> modified to run on any number of nodes by Dario Alfe,
!> BG_ routines were added by IBM fellow (bG)
!>
!> 10.09.2014   gK:
!> - removed the T3D support (broken anyway most likely)
!> - BG_ routines are now always compiled in
!> - added some comments
!> - added double check that BG_INIT_scala was called before calling
!>    compute routines (simply not sure this is done properly in VASP)
!> - removed a lot of unused variables
!> - some clean up, but a lot of off and creapy stuff is still around
!>   specifically the redundance between BG_INIT_SCALA and INIT_SCALA
!>   is odd
!>
!> @note This is a common note for routines which are called with descriptor
!> DESCA (in the source file starting with #DISTRI).
!> The routines can be called from any context or subgroup
!> of MPI processed as long as DESCA is properly set up.
!> They also do not refer to GSD. In that sense they are much safer to use, but
!> matter of fact the caller must "know" scalapack.
!> For the time being the check however that the passed matrix dimensions are
!> consistent with those set up in GSD by calling CHECK_scala. This should be
!> removed after more extensive testing.
!>
!=======================================================================

#ifndef scaLAPACK


 MODULE scala
   USE prec
   USE mpimy
   USE chi_glb
   USE minimax_struct

   IMPLICIT NONE

   LOGICAL, PUBLIC :: LscaLAPACK = .FALSE.      !< do not use scaLAPACK
   LOGICAL, PUBLIC :: LscaLU     = .FALSE.      !< do not perform parallel LU decomposition
   LOGICAL, PUBLIC :: LscaAWARE  = .FALSE.      !< do not store subspace matrix block cyclically
   LOGICAL, PUBLIC :: LscaAWARE_read  = .FALSE. !< LscaAWARE was read not from file

   LOGICAL, PUBLIC :: LELPA = .FALSE.           !< do not use ELPA

  CONTAINS
    SUBROUTINE  INIT_SCALAAWARE( NB_TOT, NRPLWV, COMM )
      INTEGER :: NB_TOT, NRPLWV
      TYPE (communic) COMM
    END SUBROUTINE INIT_SCALAAWARE

    SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN, NDIM, N)
      TYPE (communic) COMM
      INTEGER NDIM,N            !< NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N)    !< input/output matrix
    END SUBROUTINE pPOTRF_TRTRI

    SUBROUTINE pPOTRF_POTRI (COMM, AMATIN, NDIM, N)
      TYPE (communic) COMM
      INTEGER NDIM,N            !< NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N)    !< input/output matrix
    END SUBROUTINE pPOTRF_POTRI

    SUBROUTINE pDSSYEX_ZHEEVX(COMM,AMATIN,W, NDIM, N, ABSTOL, AMATIN_SCALA)
      TYPE (communic) COMM
      INTEGER NDIM,N                        !< NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N)                !< input/output matrix
      REAL(q) W(N)                          !< eigenvalues
      REAL(q), OPTIONAL ::  ABSTOL          !< tolerance
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(N*N) !< alternative input/output matrix (scaLAPACK format)
    END SUBROUTINE pDSSYEX_ZHEEVX

    SUBROUTINE BG_pDSSYEX_ZHEEVX(COMM,A, W, N)
      TYPE (communic) COMM
      INTEGER N                             !< NxN matrix to be distributed
      GDEF    A(1,1)                        !< input/output matrix
      REAL(q) W(N)                          !< eigenvalues
    END SUBROUTINE BG_pDSSYEX_ZHEEVX

    SUBROUTINE BG_pPOTRF_TRTRI(A, N, INFO)
      TYPE (communic) COMM
      INTEGER N          !< NxN matrix to be distributed
      GDEF    A(1,1)     !< input/output matrix
      INTEGER INFO
    END SUBROUTINE BG_pPOTRF_TRTRI

    SUBROUTINE BG_pPOTRF_POTRI(A, N, INFO)
      TYPE (communic) COMM
      INTEGER N          !< NxN matrix to be distributed
      GDEF    A(1,1)     !< input/output matrix
      INTEGER INFO
    END SUBROUTINE BG_pPOTRF_POTRI

    SUBROUTINE BG_CHANGE_DIAGONALE(N,A,IU0)
      INTEGER N,IU0
      GDEF    A(*)
    END SUBROUTINE BG_CHANGE_DIAGONALE

    SUBROUTINE INIT_scala(COMM, N)
      IMPLICIT NONE
      TYPE (communic) :: COMM
      INTEGER N
    END SUBROUTINE INIT_scala

    FUNCTION SCALA_NP()
      INTEGER SCALA_NP
      SCALA_NP=1
    END FUNCTION SCALA_NP

    FUNCTION SCALA_NQ()
      INTEGER SCALA_NQ
      SCALA_NQ=1
    END FUNCTION SCALA_NQ

    SUBROUTINE INIT_scala_t3d
    END SUBROUTINE INIT_scala_t3d

 END MODULE

#else

 MODULE scala
    USE scala_struct
    USE prec
    USE mpimy
#ifdef ELPA
    USE elpa, ONLY : elpa_t, ELPA_INIT, ELPA_ALLOCATE, ELPA_SOLVER_2STAGE, ELPA_OK, ELPA_DEALLOCATE, ELPA_UNINIT
#endif

    IMPLICIT NONE

#ifdef ELPA
    LOGICAL, PUBLIC :: LELPA      = .TRUE.       !< use ELPA
    INTEGER,SAVE,PRIVATE :: ELPA_MPI_COMM_ROWS,ELPA_MPI_COMM_COLS,ELPA_MPI_COMM_ALL
#else
    LOGICAL, PUBLIC :: LELPA      = .FALSE.      !< do not use ELPA
#endif

    LOGICAL, PUBLIC :: LscaLAPACK = .TRUE.       !< use scaLAPACK
    LOGICAL, PUBLIC :: LscaLU     = .FALSE.      !< perform parallel LU decomposition
    LOGICAL, PUBLIC :: LscaAWARE  = .FALSE.      !< distribute matrices block cyclically
    LOGICAL, PUBLIC :: LscaAWARE_read  = .FALSE. !< LscaAWARE was read from file

!
! customize if required
!
      !> blocking factor for distribution of matrices
      !>
      !> P4 optimal, larger values were even slightly better (160) but still
      !> slower on a Gigabit cluster than ZHEEVX
      INTEGER,PRIVATE :: NB=16

      ! customization of matrix diagonalization

      !> maximum cluster of eigenvector
      !>
      !>  see manpage for pSSYEVX pZHEEVX:
      !>    NCLUST is the number of eigenvalues in the largest cluster,
      !>    where a cluster is defined as a set of close eigenvalues:
      !>    {W(K),...,W(K+NCLUST-1)|W(J+1)<= W(J)+orfac*norm(A)}
      INTEGER,PRIVATE :: NCLUST=48
      REAL(q),PRIVATE :: ABSTOL_DEF=1e-10_q   !< specifies eigenvector orthogonality tolerance
      REAL(q),PRIVATE :: ORFAC=-1.e0_q        !< controls reorthogonalisation of eigenvectors
      REAL(qs),PRIVATE:: ORFAC_single=-1.e0_qs!< controls reorthogonalisation of eigenvectors

!
! end customization
!


!
! INIT_scala initializes the following variables
!
! usually VASP uses a "light-weighted" implementation of scalapack, in the
! sense that the matrices are stored in the conventional seriel manner
! and only when the scaLAPACK routines are called, the matrices are
! stored distributed and then treated by scaLAPACK
! this has some disadvantages, for instance the storage requirements
! are huge
! if the variable  LscaAWARE is set in the INCAR, however
! VASP uses distributed matrices from the outside
! however, even in this case, scaLAPACK is essentially completely
! "unvisible" to VASP, and all scaLAPACK related routines
! should go in here
!
! in several place (bse.F, chi_GG.F) that strategy did not work out,
! and the implementations is directly in the main VASP code
!
! most routines in this module use a global array descriptor DESCSTD
! which is set up by INIT_SCALA. INIT_SCALA can be called
! repeatedly if the matrix dimensions or the processor grid change
! however, all subroutines that have an argument DESCA
! are save and they can be used in different contexts
!
      TYPE scalapack_des
         INTEGER :: MAGIC=0        !< tells whether MODULE is initialized
         INTEGER :: N              !< dimension of matrix currently treated
                                   !< in all routines checked against the actual dimension
         INTEGER :: MPI_COMM       !< mpi communicator used for scaLAPACK
         INTEGER :: NPROCS         !< number of PEs
         INTEGER :: ICTXT          !< context handle of grid (usually MPI context)
         INTEGER :: NP             !< number of rows on the processor
         INTEGER :: NQ             !< number of cols on the processor
         INTEGER :: NPROW,NPCOL    !< processor grid dimensions
         INTEGER :: MYROW,MYCOL    !< processor coordinates in ps. grid
         INTEGER :: LWWORK,LIWORK,LRWORK,MALLOCPQ !< sizes of scalapack workarrays
      END TYPE scalapack_des

      !> global sca-lapack descriptor used as "default" in VASP
      TYPE (scalapack_des), SAVE, PRIVATE :: GSD
      INTEGER,SAVE   :: DESCSTD( DLEN_ ) !< distributed matrix descriptor array

! few notes: many of the variables are redundant with DESCSTS
!
! GSD%NP= NUMROC(GSD%N,DESCSTD(MB_),GSD%MYROW,0,GSD%NPROW)
! GSD%NQ= NUMROC(GSD%N,DESCSTD(NB_),GSD%MYCOL,0,GSD%NPCOL)
! CALL BLACS_GRIDINFO(DESCSTD(CTXT_),GSD%NPROW,GSD%NPCOL,GSD%MYROW,GSD%MYCOL) to determine GSD%NPROW,...GSD%MYCOL
!
! not sure how lightweighted the evaluation is
!

!> size of arrays defined with GDEF
#ifdef gammareal
      INTEGER, PARAMETER,PRIVATE :: MCOMP=1
#else
      INTEGER, PARAMETER,PRIVATE :: MCOMP=2
#endif

!> describes a distributed (fermionic) Green's function or self-energy matrix in orbital basis
   TYPE greens_mat_des
      TYPE (communic) :: COMM             !< communicator for communication between all nodes
      TYPE (communic) :: COMM_INTRA       !< communicator for nodes sharing one time/frequency
      TYPE (communic) :: COMM_INTER       !< communicator between time/frequency points
      INTEGER :: NROWS                    !< maximal matrix size (required for allocation)
      INTEGER :: MY_NROWS                 !< local row matrix size
      INTEGER :: MY_NCOLS                 !< local column matrix size 
      INTEGER :: ICTXT                    !< BLACS context defined on top of COMM_INTAU
      INTEGER :: DESC(DLEN_)              !< matrix descriptor

      INTEGER :: NPROW,NPCOL              !< processor grid dimensions
      INTEGER :: MYROW,MYCOL              !< processor coordinates in ps. grid
   END TYPE greens_mat_des

 CONTAINS
!=======================================================================
!
!> this subroutine sets LscaAWARE to a default value if not read from
!> INCAR
!
!=======================================================================
    SUBROUTINE  INIT_SCALAAWARE( NB_TOT, NRPLWV, COMM )
      USE reader_tags, ONLY: INCAR_F
      USE tutor, ONLY: vtutor
      USE incar_reader, ONLY: PROCESS_INCAR
      IMPLICIT NONE

      INTEGER :: NB_TOT, NRPLWV
      TYPE (communic) COMM
    ! local
      REAL(q) :: STORAGE_HAM, STORAGE_WAVE

      CHARACTER(LEN=256) :: ENVVAR_VALUE
      INTEGER :: ENVVAR_STAT, IERR

#ifdef gammareal
      STORAGE_HAM =8._q*NB_TOT*NB_TOT
      STORAGE_WAVE=16._q*NB_TOT*NRPLWV/ COMM%NCPU
#else
      STORAGE_HAM =16._q*NB_TOT*NB_TOT
      STORAGE_WAVE=16._q*NB_TOT*NRPLWV/ COMM%NCPU
#endif
      IF (.NOT. LscaAWARE_read ) THEN
         IF (STORAGE_HAM>STORAGE_WAVE) THEN
            LscaAWARE=LscaLAPACK
         ELSE
            LscaAWARE=.FALSE.
         ENDIF
      ENDIF

      ! defaults
#ifdef CUSOLVERMP
      NB=256
#else
      NB=16
#endif

      ! check and block size NB from environment
      CALL GET_ENVIRONMENT_VARIABLE("VASP_SCALA_NB", ENVVAR_VALUE, STATUS=ENVVAR_STAT)

      IF (ENVVAR_STAT == 0) THEN
          READ(ENVVAR_VALUE,*,iostat=ENVVAR_STAT) NB
          IF (ENVVAR_STAT /= 0) THEN
              CALL vtutor%error("ERROR: INIT_SCALAAWARE: could not parse VASP_SCALA_NB value from env.")
          ENDIF
      ENDIF

      ! read block size from INCAR
      CALL PROCESS_INCAR(INCAR_F, "NBSCALA", NB, IERR)

    END SUBROUTINE

!=======================================================================
!
!> function to determine the blocking factor used in scaLAPACK
!
!=======================================================================

    SUBROUTINE  QUERRY_SCALA_NB( NB_QUERRY)
      INTEGER, INTENT(INOUT) :: NB_QUERRY 
      NB_QUERRY=NB
    END SUBROUTINE QUERRY_SCALA_NB

!=======================================================================
!> helper function, factorizes NPROCS into NPROW*NPCOL
!>
!> with restriction that NPROW ~ NPCOL and NPROW >= NPCOL
!=======================================================================
    SUBROUTINE FERMAT_RAZOR( NPROCS, NPROW, NPCOL )
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: NPROCS
       INTEGER, INTENT(OUT) :: NPROW, NPCOL
       ! local
       INTEGER              :: NPRIMES
       INTEGER, POINTER     :: PRIMES(:) => NULL()
       INTEGER, POINTER     :: ALL_PRIMES(:) => NULL()
       INTEGER, POINTER     :: NEXP(:) => NULL()
       INTEGER              :: I,IP
       INTEGER              :: NODD
       INTEGER              :: A, B
       PROFILING_START( 'fermat_razor' )

       IF ( NPROCS == 1 ) THEN
          NPROW = 1
          NPCOL = 1
          PROFILING_STOP( 'fermat_razor' )
          RETURN
       ENDIF

       ! obtain all primes up to SQRT( NPROCS  )
       CALL ERATOSTHENES_SIEVE( NPROCS, PRIMES, ALL_PRIMES )
       NPRIMES = SIZE( PRIMES )
       ! obtain prime factorization of NPROCS
       CALL PRIME_FACTORIZATION( NPROCS, NEXP, PRIMES )

       ! first check if NPROCS is a perfect square
       IF ( IS_PERFECT_SQUARE( NPROCS ) ) THEN
          NPROW = PERFECT_SQUARE( NPROCS )
          NPCOL = PERFECT_SQUARE( NPROCS )
       ELSE
       ! the idea is to factorize NPROCS = A*B
       ! with restriction |A - B| -> 0
       ! fermat sieve works well if N is odd
       NODD = NPROCS/(PRIMES(1)**NEXP(1))

       ! in case NPROCS is even
       IF ( NODD == 1 ) THEN

          NPROW = PRIMES(1)**( NEXP(1)/2 )
          NPCOL = PRIMES(1)**( NEXP(1)/2+MOD(NEXP(1),2) )
       ! in case NPROCS is not even
       ! use fermats sieve on NODD
       ELSE

          ! perfect square nothing to do
          IF ( IS_PERFECT_SQUARE( NODD ) ) THEN
             NPROW = PERFECT_SQUARE( NODD )
             NPCOL = PERFECT_SQUARE( NODD )
          ELSE
             IF ( IS_PRIME( NODD, ALL_PRIMES ) ) THEN
                NPROW = 1
                NPCOL = NODD
             ELSE
                ! Fermat's way to factorize an odd integer:
                ! NODD = A**2- B**2 = ( A - B ) * ( A + B )
                ! and write A**2 - N = B**2
                ! and increase A by 1 until A**2 - N is a perfect square
                ! this gives B and finally the factorization
                A = CEILING( SQRT( REAL(NODD) ) )
                B = A*A - NODD
                DO
                   IF ( IS_PERFECT_SQUARE( B ) ) EXIT
                   A = A+1
                   B = A*A-NODD
                ENDDO
                NPROW = A+PERFECT_SQUARE( B )
                NPCOL = A-PERFECT_SQUARE( B )
             ENDIF
          ENDIF
          ! distribute remaining factors of 2
          ! such that NPROW is roughly NPCOL
          DO I = 1, NEXP(1)
             IF ( NPROW < NPCOL ) THEN
                NPROW = NPROW*2
             ELSE
                NPCOL = NPCOL*2
             ENDIF
          ENDDO
       ENDIF

       ENDIF ! perfect square

#ifdef CUSOLVERMP
       ! enforce NPCOL >= NPROW via check and swap (needed for e.g. NPROCS = 240 or 480)
       IF(NPROW > NPCOL) THEN
          A = NPCOL
          NPCOL = NPROW
          NPROW = A
       ENDIF
#endif

       IF ( NPROW*NPCOL /= NPROCS ) THEN
          CALL vtutor%error("ERROR, FERMAT_RAZOR failed " // str(NPROCS) // " " // str(NPROW) // " "&
              // str(NPCOL))
       ENDIF

       DEALLOCATE( PRIMES )
       DEALLOCATE( ALL_PRIMES )
       DEALLOCATE( NEXP )
       NULLIFY( PRIMES )
       NULLIFY( ALL_PRIMES )
       NULLIFY( NEXP )

       PROFILING_STOP( 'fermat_razor' )

       CONTAINS

       !> Eratosthenes_Sieve to obtain primes up to N
       SUBROUTINE ERATOSTHENES_SIEVE(N, PRIM_SQRT, PRIMES)
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: N
          INTEGER, POINTER     :: PRIM_SQRT(:)
          INTEGER              :: I, J, MAXPRIME, STAT
          LOGICAL              :: A(N)
          INTEGER, POINTER     :: PRIMES(:)

          MAXPRIME = FLOOR(SQRT(REAL(N)))
          A = .TRUE.

          DO I=2,MAXPRIME
             J = I*I
             DO
                A(J) = .FALSE.
                J = J + I ; IF ( J .GT. N ) EXIT
             ENDDO
          ENDDO !I

          ALLOCATE( PRIMES( COUNT(A)-1 ) )
          J = 1
          DO I=2,N ! SKIP 1
             IF ( .NOT. A(I) ) CYCLE
             PRIMES( J ) = I
             J = J + 1 ; IF ( J > SIZE(PRIMES) ) EXIT
          ENDDO


          J = 0
          DO I = 1, SIZE( PRIMES )
             IF ( PRIMES( I )*PRIMES( I ) > N ) THEN
                J = J+1
                EXIT
             ENDIF
             J = J+1
          ENDDO

          ALLOCATE( PRIM_SQRT( J ) )
          DO I = 1, J
             PRIM_SQRT( I ) = PRIMES( I )
          ENDDO
       END SUBROUTINE ERATOSTHENES_SIEVE

       !> obtains prime factoriztion of N
       SUBROUTINE PRIME_FACTORIZATION( N, NEXP, PRIMES )
          IMPLICIT NONE
          INTEGER          :: N
          INTEGER, POINTER :: NEXP(:)
          INTEGER          :: PRIMES(:)
          ! local
          INTEGER          :: NPRIMES
          INTEGER          :: IP
          INTEGER          :: M

          NPRIMES=SIZE( PRIMES )
          ALLOCATE( NEXP( NPRIMES ) )
          NEXP = 0
          ! 1 is not considered as prime here
          M = N
          !loop over all primes in PRIMES
          prime: DO IP = 1 , NPRIMES
             ! find exponent
             DO
                ! if this prime is a divisor
                IF( MOD( M, PRIMES( IP ) ) == 0 ) THEN
                   M = M/PRIMES( IP )
                   ! exponent of factor
                   NEXP( IP ) = NEXP( IP ) + 1
                ELSE
                   EXIT
                ENDIF
             ENDDO
          ENDDO prime
       END SUBROUTINE PRIME_FACTORIZATION

       !> returns true if B is prime
       FUNCTION IS_PRIME( B, PRIMES )
          LOGICAL :: IS_PRIME
          INTEGER :: PRIMES(:)
          INTEGER :: B
          ! local
          INTEGER :: I
          INTEGER :: ILIMIT

          ILIMIT = SIZE( PRIMES )

          DO I = 1, ILIMIT
             IF ( B == PRIMES( I ) ) EXIT
          ENDDO

          IF ( I > ILIMIT ) THEN
             IS_PRIME = .FALSE.
          ELSE
             IS_PRIME = .TRUE.
          ENDIF
       END FUNCTION IS_PRIME

       !> returns true if B is perfect square
       FUNCTION IS_PERFECT_SQUARE( B )
          LOGICAL :: IS_PERFECT_SQUARE
          INTEGER :: B
          ! local
          INTEGER :: I
          INTEGER :: ILIMIT

          ILIMIT =  CEILING( SQRT( REAL(B) ) )
          DO I = 1, ILIMIT
             IF ( I*I == B ) EXIT
          ENDDO
          IF ( I > ILIMIT ) THEN
             IS_PERFECT_SQUARE = .FALSE.
          ELSE
             IS_PERFECT_SQUARE = .TRUE.
          ENDIF
       END FUNCTION IS_PERFECT_SQUARE

       !> returns perfect square of B
       FUNCTION PERFECT_SQUARE( B )
          INTEGER :: PERFECT_SQUARE
          INTEGER :: B
          ! local
          INTEGER :: I
          INTEGER :: ILIMIT

          ILIMIT =  CEILING( SQRT( REAL(B) ) )
          DO I = 1, ILIMIT
             IF ( I*I == B ) EXIT
          ENDDO
          PERFECT_SQUARE = I
       END FUNCTION PERFECT_SQUARE

    END SUBROUTINE FERMAT_RAZOR

!=======================================================================
!
!> messy scalapack initialization, calculate all required workspace
!>
!> the first routine
!>    SUBROUTINE INIT_scala(myCOMM, N )
!> sets up a global structure for scaLAPACK communication
!> this is sufficiently flexible for most parts of VASP
!>
!> However in the more complicated cases, such as GW, it might
!> be necessary to use different contextes for instance
!> subdividing the nodes into different communication universes
!> in that case
!>   SUBROUTINE INIT_scala_DESC(myCOMM, N , DESCA, GS)
!> should be used with local DESCA and GS (scalapack_des) descriptors
!
!=======================================================================

    SUBROUTINE INIT_scala(myCOMM, N )
      TYPE(communic) :: myCOMM
      INTEGER N

      CALL INIT_scala_DESC(myCOMM, N , DESCSTD, GSD)

    END SUBROUTINE INIT_scala

!
!> here is the context safe version (of #INIT_scala)
!

    SUBROUTINE INIT_scala_DESC(myCOMM, N , DESCA, GS)

      USE main_mpi            ! to get world communicator
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      INTEGER              :: DESCA( DLEN_ ) !< distributed matrix descriptor array
      TYPE (scalapack_des) :: GS             !< descriptor that should be initialized by the routine
      TYPE(communic) :: myCOMM               !< VASP MPI communicator
      INTEGER N                              !< dimension of matrices to be handled
! local
      INTEGER INFO,NEIG
      INTEGER NP0,NQ0,NN
      INTEGER ierror

      INTRINSIC :: MAX,MIN
      INTEGER,EXTERNAL :: NUMROC,ICEIL
#ifdef oldversion
      INTEGER :: FACTORS(4), PWR(4)
      DATA FACTORS /2, 3, 5, 7/   ! a few prime factors likely to be be used
      INTEGER :: NNPCOL, NNPROW, MR, NDUM, NMIN, I
#endif
      EXTERNAL BLACS_PINFO
      EXTERNAL DESCINIT
#ifdef ELPA
      INTEGER ISUCCESS
#endif

      ! if gs%magic is set and the communicator changes or
      ! the size N of the matrix changes, then reset the GS structure
      IF (GS%MAGIC == 1789231) THEN
         IF (myCOMM%MPI_COMM /= GS%MPI_COMM .OR. GS%N /=N) THEN
            ! exit the BLACS context first
            CALL BLACS_GRIDEXIT(GS%ICTXT)
            ! require re-initialization
            GS%MAGIC=0
         ENDIF
      ENDIF

! if GS%MAGIC is not set to the right value reinitialize the GS structure
      IF (GS%MAGIC /= 1789231) THEN

      GS%MAGIC=1789231
      GS%N=N
      GS%MPI_COMM=myCOMM%MPI_COMM

! number of nodes (for one image)
      GS%NPROCS = myCOMM%NCPU

# ifdef oldversion
! determine processor grid, according to naive guess
      PWR=0
      MR = GS%NPROCS
      GS%NPROW=1; GS%NPCOL=1
      DO I=1,4
         DO
           IF( MOD(MR,FACTORS(I)) == 0 ) THEN
              PWR(I) = PWR(I) + 1
              MR = MR/FACTORS(I)
           ELSE
              EXIT
           ENDIF
         ENDDO
         NNPROW = PWR(I)/2
         NNPCOL = PWR(I) - NNPROW
         IF(MOD(I,2)==0) THEN         ! to make the grid more uniformly distributed (hopefully)
           NDUM = NNPCOL
           NNPCOL = NNPROW
           NNPROW = NDUM
         ENDIF
         GS%NPROW = GS%NPROW*FACTORS(I)**NNPROW
         GS%NPCOL = GS%NPCOL*FACTORS(I)**NNPCOL
      ENDDO
      IF( GS%NPROW*GS%NPCOL < GS%NPROCS )THEN  ! final adjustment to the total number of procs
         NMIN = MIN(GS%NPROW,GS%NPCOL)
         NDUM = GS%NPROCS/GS%NPROW/GS%NPCOL
         IF(GS%NPROW==NMIN)THEN          ! again, hoping that in this way the grid is more uniform
           GS%NPROW = GS%NPROW * NDUM   ! you might think to a better algorithm ( I dont think it
         ELSE                             ! really matters )
           GS%NPCOL = GS%NPCOL * NDUM
         ENDIF
      ENDIF
#else
! determine precessor grid using Fermat Sieve
      CALL FERMAT_RAZOR( GS%NPROCS, GS%NPROW, GS%NPCOL )
#endif

! make processor grids
      CALL MPI_barrier( myCOMM%MPI_COMM, ierror )

      CALL PROCMAP( myCOMM, GS%ICTXT, 2, GS%NPROW, GS%NPCOL)

! calculate local size of matrices

      CALL BLACS_GRIDINFO( GS%ICTXT, GS%NPROW, GS%NPCOL, GS%MYROW, GS%MYCOL )

      GS%NP = NUMROC(N,NB,GS%MYROW,0,GS%NPROW)   ! get number of rows on proc
      GS%NQ = NUMROC(N,NB,GS%MYCOL,0,GS%NPCOL)   ! get number of cols on proc

      CALL DESCINIT(DESCA,N,N,NB,NB,0,0,GS%ICTXT,MAX( 1, GS%NP ),INFO ) ! setup descriptor
      IF (INFO.NE.0) THEN
        CALL vtutor%bug("INIT_SCALA: DESCA, DESCINIT, INFO: " // str(INFO), __FILE__, __LINE__)
      ENDIF

#ifdef ELPA
! if one of the ranks has no data we can not use ELPA
      ISUCCESS=0; IF (GS%NP==0.OR.GS%NQ==0) ISUCCESS=1
      CALLMPI( M_sum_i(myCOMM,ISUCCESS,1) )
      IF (ISUCCESS>0) THEN
         LELPA=.FALSE.
         IF (myCOMM%NODE_ME==myCOMM%IONODE) &
            WRITE(*,*) 'ELPA: falling back to scaLAPACK'
      ELSE
         ELPA_MPI_COMM_ALL=myCOMM%MPI_COMM
      ENDIF
#endif

! calculate scalapack workspace and allocate, pDSSYEX_ZHEEVX only

   ! iwork
      GS%LIWORK=6*MAX(N,GS%NPROW*GS%NPCOL+1,4)
   ! wwork
      NN = MAX( N, NB, 2 )
      NEIG = N                ! number of eigenvectors requested
      IF ((DESCA(MB_).NE.NB).OR.(DESCA(NB_).NE.NB) .OR. &
          (DESCA(RSRC_).NE.0) .OR.(DESCA( CSRC_).NE.0)) THEN
        CALL vtutor%bug("INIT_SCALA: pSSYEZS_ETC, ERROR", __FILE__, __LINE__)
      ENDIF
      NP0 = NUMROC( NN, NB, 0, 0, GS%NPROW )
      NQ0 = MAX( NUMROC( NEIG, NB, 0, 0, GS%NPCOL ), NB )
#ifdef gammareal
      GS%LRWORK=1
      GS%LWWORK=5*N+MAX(5*NN,NP0*NQ0)+ICEIL(NEIG,GS%NPROW*GS%NPCOL)*NN+2*NB*NB+ &
         NCLUST*N
#else
      GS%LRWORK=4*N+MAX( 5*NN, NP0 * NQ0 )+ICEIL( NEIG, GS%NPROW*GS%NPCOL)*NN+ &
         NCLUST*N
      GS%LWWORK=N + MAX((NP0 +NQ0 +NB ) * NB, 3)
#endif

      GS%MALLOCPQ=MAX(GS%NP*GS%NQ,1)

! gK: dump global variables to see what has been set
!      WRITE(*,*) GS%NPROCS         ,' number of PEs'
!      WRITE(*,*) GS%ICTXT          ,' context handle of grid'
!      WRITE(*,*) GS%NPROW,GS%NPCOL ,' processor grid dimensions'
!      WRITE(*,*) GS%MYROW,GS%MYCOL ,' processor coordinates in ps. grid'
!      WRITE(*,*) GS%NP             ,' number of rows on the processor'
!      WRITE(*,*) GS%NQ             ,' number of cols on the processor'
!      WRITE(*,*) GS%LWWORK,GS%LIWORK,GS%LRWORK,GS%MALLOCPQ ,' sizes of scalapack workarrays'

      ENDIF

    END SUBROUTINE INIT_scala_DESC


!=======================================================================
!
!> @details any routine should either call INIT_scala or CHECK_scala
!> CHECK_scala, gracefully aborts if the GSD%N is incompatible
!> to the supplied value N
!
!=======================================================================

    SUBROUTINE CHECK_scala (N, STRING)
      USE tutor, ONLY: vtutor
      INTEGER N
      CHARACTER (LEN=*) :: STRING
      IF (N /= GSD%N) THEN
         CALL vtutor%bug("CHECK_scala: BG_INIT_SCALA_BG was not called beforehand for " &
            // "this k-point \n CHECK_scala was called by " // STRING, __FILE__, __LINE__)
      ENDIF

    END SUBROUTINE CHECK_scala

!=======================================================================
!
!> Routine to querry the GSD%NP variable
!
!=======================================================================

    FUNCTION SCALA_NP()
      INTEGER SCALA_NP
      SCALA_NP=MAX(GSD%NP,1)
    END FUNCTION SCALA_NP

!=======================================================================
!
!> Routine to querry the GSD%NQ variable
!
!=======================================================================

    FUNCTION SCALA_NQ()
      INTEGER SCALA_NQ
      SCALA_NQ=MAX(GSD%NQ,1)
    END FUNCTION SCALA_NQ

!=======================================================================
!
!> Maps processors onto a grid
!>
!> This subroutine maps processors onto a grid
!> and allows to perform scaLAPACK operations on a subset of nodes
!> taken from BLACS example code (written by Clint Whaley 7/26/94)
!> modified by gD and gK
!>
!> PROCMAP maps NPROW*NPCOL onto an VASP-Communicator
!> It is a nice system independent implementation
!>
!> @param[inout] CONTEXT (output) This integer is used by the BLACS to indicate
!>               a context. A context is a universe where messages exist and do not
!>              interact with other contexts messages. The context includes
!>              the definition of a grid, and each processs coordinates in it.
!>
!> @param[inout] MAPPING (input) Way to map processes to grid. Choices are:
!>               - `1` : row-major natural ordering
!>               - `2` : column-major natural ordering
!>
!> @param[inout] NPROW (input) The number of process rows the created grid
!>               should have.
!>
!> @param[inout] NPCOL (input) The number of process columns the created grid
!>               should have.
!
!=======================================================================

    SUBROUTINE PROCMAP(VCOMM,CONTEXT, MAPPING, NPROW, NPCOL)
      USE main_mpi     ! to get MPI communication handle
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(communic) :: VCOMM

!     .. Scalar Arguments ..
      INTEGER :: CONTEXT, MAPPING, NPROW, NPCOL
      INTEGER :: IMAP(NPROW, NPCOL)
      INTEGER ierror
!     ..
!
      INTEGER,EXTERNAL :: BLACS_PNUM
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDMAP
      INTEGER TMPCONTXT, NPROCS, I, J, K
      INTEGER BLACS_IDS(VCOMM%NCPU), MY_BLACS_ID
      INTEGER NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP
!
!     See how many processes there are in the system
!     NPROCS here refers to all processors being used

      NPROCS=VCOMM%NCPU
!
!     Temporarily map all processes into 1 x NPROCS grid
!
      TMPCONTXT=VCOMM%MPI_COMM
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )
      CALL BLACS_GRIDINIT( TMPCONTXT, 'Row', 1, NPROCS )
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )

!
!     If we want a row-major natural ordering
!
      CALL BLACS_GRIDINFO( TMPCONTXT, NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP)

      MY_BLACS_ID= BLACS_PNUM(TMPCONTXT, MYROW_TMP, MYCOL_TMP)
      CALL MPI_ALLGATHER(MY_BLACS_ID,1,MPI_INTEGER, &
     &     BLACS_IDS(1),1,MPI_INTEGER,VCOMM%MPI_COMM,IERROR)
      IF (IERROR.NE.0) CALL vtutor%error('ERROR allgather in procmap' // str(ierror))

      K=1
      IF (MAPPING .EQ. 1) THEN
         DO I = 1, NPROW
            DO J = 1, NPCOL
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
!
!     If we want a column-major natural ordering
!
      ELSE IF (MAPPING .EQ. 2) THEN
         DO J = 1, NPCOL
            DO I = 1, NPROW
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
      ELSE
         CALL vtutor%error("INIT_SCALA, PROCMAP: unknown mapping, STOP")
      END IF
!
!     Free temporary context
!
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )
      CALL BLACS_GRIDEXIT(TMPCONTXT)
!      WRITE (*,*) 'MAPPING',NPCOL,NPROW
!
!     Apply the new mapping to form desired context
!
      CONTEXT=VCOMM%MPI_COMM
      CALL BLACS_GRIDMAP( CONTEXT, IMAP, NPROW, NPROW, NPCOL )

!      WRITE(*,'("local map",16I3)') ((IMAP(I,J),I=1,NPROW),J=1,NPCOL)

      RETURN
    END SUBROUTINE PROCMAP

!=======================================================================
!
! the following routines rely on the global GSD communicator
! they however re-initialize GSD as well, so should be safe
! against changes of the matrix dimension
!
!=======================================================================
!
!> this subroutine determines U^-1: A= U+ U and invert U
!>
!> 1. calls INIT_scala
!> 2. calls the setup of distributed the distributed matrix
!> 3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!> 4. calls reconstruction of distributed data into patched matrix
!>
!> @note The sum over processors of the patched matrix is not realised
!>       in this module
!
!=======================================================================


    SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN, NDIM, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           !< leading dimension of AMATIN
      INTEGER N              !< NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N) !< input/output matrix
! local variables
      INTEGER INFO

      GDEF, ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      PROFILING_START('ppotrf_trtri')

      CALL INIT_scala(COMM, N)

!-----------------------------------------------------------------------
! allocate workarray (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------
      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,DESCSTD)

        INFO=0
#ifdef gammareal
        CALL PDPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#else
        CALL PZPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#endif

        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_TRTRI, POTRF, INFO: " // str(INFO) // "\n STOP")
        ENDIF
#ifdef gammareal
         CALL PDTRTRI &
#else
         CALL PZTRTRI &
#endif
     &    ('U', 'N', N, A, 1, 1, DESCSTD, INFO)
        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_TRTRI, TRTRI, INFO: " // str(INFO))
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,DESCSTD)
   ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)

      PROFILING_STOP('ppotrf_trtri')

      RETURN
    END SUBROUTINE pPOTRF_TRTRI


!=======================================================================
!
!> this subroutine determines A^-1: A= U+ U and invert (U+ U)
!>
!> 1. calls INIT_scala
!> 2. calls the setup of distributed the distributed matrix
!> 3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!> 4. calls reconstruction of distributed data into patched matrix
!>
!> @note The sum over processors of the patched matrix is not realised
!>       in this module
!
!=======================================================================


    SUBROUTINE pPOTRF_POTRI (COMM, AMATIN, NDIM, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           ! leading dimension of AMATIN
      INTEGER N              ! NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N) ! input/output matrix
! local variables
      INTEGER INFO

      GDEF, ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      PROFILING_START('ppotrf_potri')

      CALL INIT_scala(COMM, N)
!-----------------------------------------------------------------------
! allocate workarray
!-----------------------------------------------------------------------
      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,DESCSTD)

        INFO=0
#ifdef gammareal
        CALL PDPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#else
        CALL PZPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#endif

        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_TRTRI, POTRF, INFO: " // str(INFO) // "\n STOP")
        ENDIF
#ifdef gammareal
         CALL PDPOTRI &
#else
         CALL PZPOTRI &
#endif
     &    ('U', N, A, 1, 1, DESCSTD, INFO)
        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_POTRI, POTRI, INFO: " // str(INFO))
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,DESCSTD)
   ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)

      PROFILING_STOP('ppotrf_potri')

      RETURN
    END SUBROUTINE pPOTRF_POTRI


!=======================================================================
!
!> Call to pZHEEVX respectively pDSYEVX
!>
!> i.e. diagonalization of a symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE pDSSYEX_ZHEEVX(COMM, AMATIN, W, NDIM, N, ABSTOL_, AMATIN_SCALA)
#ifdef ELPA
      USE iso_c_binding
#endif
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          !< NxN matrix to be distributed
      INTEGER NDIM                       !< leading dimension of matrix
      GDEF    AMATIN(NDIM,N)             !< input/output matrix
      REAL(q) W(N)                       !< eigenvalues
      REAL(q), OPTIONAL      ::  ABSTOL_ !< tolerance
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(:) !< alternative input/output matrix (scaLAPACK format)
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2
      REAL(q) :: ABSTOL

      GDEF, ALLOCATABLE, TARGET :: A(:) ! A matrix to be diagonalized (allocatable)
      GDEF, ALLOCATABLE, TARGET :: Z(:) ! Z eigenvector matrix
      GDEF, POINTER        ::     AP(:) ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array

      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE

#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS

      GDEF, POINTER :: EP(:,:), EZ(:,:)
#endif

      PROFILING_START('pdssyex_zheevx')

      CALL INIT_scala(COMM, N)

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         CALL vtutor%error("pDSSYEX_ZHEEVX: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " "&
             // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         IF (SIZE(AMATIN_SCALA)<GSD%MALLOCPQ) CALL vtutor%bug("pDSSYEX_ZHEEVX: AMATIN_SCALA too small: " &
               // str(SIZE(AMATIN_SCALA)) // " " // str(GSD%MALLOCPQ),__FILE__,__LINE__)
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
      ALLOCATE(Z(GSD%MALLOCPQ))
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(GSD%LWWORK))
      ALLOCATE(RWORK(GSD%LRWORK))
DOESI Z = 0; IWORK = 0; WWORK = 0; RWORK = 0.0_q ! Breaks BO_phon (intel), uninitialized values trigger SIGFPE in scaLAPACK routine.

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

         !     get parts of the matrix into the local arrays
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N, A, DESCSTD)
         ENDIF
         INFO=0

#ifdef ELPA
        ! call ELPA 2-stage solver
        IF (LELPA) THEN
           IF (ELPA_INIT(20171201) /= ELPA_OK) &
              CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

           CALL c_f_pointer(c_loc(AP),EP,[GSD%NP,GSD%NQ])
           CALL c_f_pointer(c_loc( Z),EZ,[GSD%NP,GSD%NQ])

           ET => ELPA_ALLOCATE()

           CALL ET%SET("na",  N, ISUCCESS)
           CALL ET%SET("nev", N, ISUCCESS)

           CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
           CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

           CALL ET%SET("nblk", NB, ISUCCESS)

           CALL ET%SET("mpi_comm_parent", COMM%MPI_COMM, ISUCCESS)
           CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
           CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

           ISUCCESS = ET%SETUP()

           !CALL ET%SET("nvidia-gpu", 1, ISUCCESS)
           CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)

           CALL ET%EIGENVECTORS(EP, W, EZ, ISUCCESS)

           CALL ELPA_DEALLOCATE(ET)
           CALL ELPA_UNINIT()

           NULLIFY(EP,EZ)
        ENDIF
#endif
        ! call scaLAPACK routine
        IF (.NOT. LELPA) THEN
#ifdef gammareal
         CALL PDSYEVX( 'V', 'A', 'U', N, AP(1), 1, 1, DESCSTD, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, WWORK, &
              GSD%LWWORK, IWORK, GSD%LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
         CALL PZHEEVX( 'V', 'A', 'U', N, AP(1), 1, 1, DESCSTD, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, &
              WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
              IFAIL, ICLUSTR, GAP, INFO )
#endif
         IF (M.NE.N) THEN
            CALL vtutor%error("GSD%LWWORK " // str(GSD%LWWORK) // " " // str(GSD%LRWORK) // " " // &
               str(GSD%LIWORK) // " " // str(DESCSTD(M_)) // "\nERROR in subspace rotation PDSYEVX/ &
               &PZHEEVX: not enough eigenvalues found " // str(M) // " " // str(N))
         ENDIF
         IF (NZ.NE.N) THEN
            CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough eigenvectors &
               &computed " // str(M) // " " // str(N))
         ENDIF
         DO I2=1,N
            IF (IFAIL(I2).NE.0) THEN
               CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: I2,IFAIL= " // str(I2) &
                  // " " // str(IFAIL(I2)))
            ENDIF
         ENDDO
         IF (INFO<0) THEN
            CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
         ELSE IF (INFO>0) THEN
            IF (MOD(INFO,2).NE.0) THEN
               CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ELSEIF (MOD(INFO/2,2).NE.0) THEN
               ! this error condition happens very often, but does not imply that
               ! the result is not usefull
               ! IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
            ELSE
               CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ENDIF
         ENDIF
        ENDIF

        IF (PRESENT(AMATIN_SCALA)) THEN
           AMATIN_SCALA(1:GSD%MALLOCPQ)=Z(1:GSD%MALLOCPQ)
        ELSE
           CALL RECON(AMATIN, NDIM, N, Z, DESCSTD)
        ENDIF

     ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF
      DEALLOCATE(Z)
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)

      PROFILING_STOP('pdssyex_zheevx')

      RETURN
    END SUBROUTINE pDSSYEX_ZHEEVX

!=======================================================================
!
!> Call to pZHEEVD respectively pDSYEVD
!>
!> i.e. diagonalization of a symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE pDSYEV_ZHEEVD(COMM, AMATIN, W, NDIM, N, AMATIN_SCALA)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          !< NxN matrix to be distributed
      INTEGER NDIM                       !< leading dimension of matrix
      GDEF    AMATIN(NDIM,N)             !< input/output matrix
      REAL(q) W(N)                       !< eigenvalues
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(:) !< alternative input/output matrix (scaLAPACK format)
! local variables
      GDEF, ALLOCATABLE,TARGET ::  A(:) ! A matrix to be diagonalized (allocatable)
      GDEF, POINTER  ::  AP(:)          ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      INTEGER NODE_ME,IONODE

      CALL INIT_scala(COMM, N)
      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         CALL vtutor%error("pDSYEV_ZHEEVD: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " "&
             // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

         !     get parts of the matrix into the local arrays (AP=>A)
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N, A(1), DESCSTD)
         ENDIF

         CALL PDSYEV_ZHEEVD_DESC(AP(:), W, N, DESCSTD, COMM)
! if desired this can be also replaced by the standard matrix diagonalization
!         CALL PDSYEV_ZHEEV_DESC(AP(:), W, N, DESCSTD, COMM)

         IF (PRESENT(AMATIN_SCALA)) THEN
            AMATIN_SCALA=AP
         ELSE
            CALL RECON(AMATIN, NDIM, N, A(1), DESCSTD)
         ENDIF

     ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF

      RETURN
    END SUBROUTINE

!=======================================================================
!
!> Clean the matrix elements along the diagonal
!>
!> Clean the matrix elements along the diagonal (force them to be real),
!> as it is done in the serial version
!
!=======================================================================
#define testslice
#ifndef testslice
   SUBROUTINE BG_CHANGE_DIAGONALE(N, A, IU0, NSTART, NSTOP)
      IMPLICIT NONE
      GDEF    A(*)                !< distributed matrix
      INTEGER N                   !< dimension of the matrix
      INTEGER IU0                 !< IO unit for error report
      INTEGER, OPTIONAL :: NSTART !< start with global column NSTART
      INTEGER, OPTIONAL :: NSTOP  !< stop at global column NSTOP
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER IROW, ICOL, ICOL1, ICOL2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO
#ifndef gammareal
      CALL CHECK_scala (N, 'BG_CHANGE_DIAGONALE' )

      PROFILING_START('bg_change_diagonale')

      CALL BLACS_GRIDINFO(DESCSTD(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCSTD(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCSTD(NB_),MYCOL,0,NPCOL)

      ICOL1 = 1 ; IF (PRESENT(NSTART)) ICOL1 = NSTART
      ICOL2 = N ; IF (PRESENT(NSTOP) ) ICOL2 = NSTOP

      JCOL = 0

      ! loop over local columns
      DO J1=1,NQ,DESCSTD(NB_)
        J1RES=MIN(DESCSTD(NB_),NQ-J1+1)
        DO J2=1,J1RES

          IROW_JCOL = JCOL

          ! global column index
          ICOL = DESCSTD(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (ICOL < ICOL1 .OR. ICOL > ICOL2) GOTO 100

          ! loop over local rows
          DO I1=1,NP,DESCSTD(MB_)
            I1RES=MIN(DESCSTD(MB_),NP-I1+1)
            DO I2=1,I1RES

              IROW_JCOL = IROW_JCOL+1

              ! global row index
              IROW = DESCSTD(MB_)*MYROW+NPROW*(I1-1)+I2

              ! are we on the diagonal?
              IF( IROW == ICOL ) THEN
                 ! check whether the matrix is hermitian
                 IF (ABS(AIMAG(A(IROW_JCOL)))>1E-2_q .AND. IU0>=0) THEN
                    WRITE(IU0,*)'BG_CHANGE_DIAGONALE: WARNING: sub-space matrix is not hermitian', &
                    &                 AIMAG(A(IROW_JCOL)),DESCSTD(NB_)*MYCOL+NPCOL*(J1-1)+J2
                 ENDIF
                 ! make sure the diagonal is strictly real
                 A(IROW_JCOL) = REAL( A(IROW_JCOL) ,KIND=q)

                 ! since we have found the right row we may now move
                 ! on to the next column
                 GOTO 100 
              ENDIF

            ENDDO
          ENDDO
100       JCOL = JCOL + DESCSTD(LLD_)
        ENDDO
      ENDDO

      PROFILING_STOP('bg_change_diagonale')
#endif
      RETURN
    END SUBROUTINE BG_CHANGE_DIAGONALE
#else
    SUBROUTINE BG_CHANGE_DIAGONALE(N, A, IU0, NSTART, NSTOP)

GPACC USE moffload_struct_def

      IMPLICIT NONE

      GDEF    A(*)                !< distributed matrix
      INTEGER N                   !< dimension of the matrix
      INTEGER IU0                 !< IO unit for error report
      INTEGER, OPTIONAL :: NSTART !< start with global column NSTART
      INTEGER, OPTIONAL :: NSTOP  !< stop at global column NSTOP

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST, K
      INTEGER ROW, COL, ICOL1, ICOL2

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

#ifndef gammareal
      PROFILING_START('bg_change_diagonale')

      CALL CHECK_scala (N, 'BG_CHANGE_DIAGONALE' )

      CALL BLACS_GRIDINFO(DESCSTD(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCSTD(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCSTD(NB_),MYCOL,0,NPCOL)

      MB = DESCSTD(MB_ )
      NB = DESCSTD(NB_ )
      LD = DESCSTD(LLD_)

      ICOL1 = 1 ; IF (PRESENT(NSTART)) ICOL1 = NSTART
      ICOL2 = N ; IF (PRESENT(NSTOP) ) ICOL2 = NSTOP

!$ACC PARALLEL LOOP PRESENT(A) PRIVATE(JBLCK,JBASE,JREST,COL,I,IBLCK,IBASE,IREST,ROW,K) __IF_ASYNC__
      cols: DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         IF (COL < ICOL1 .OR. COL > ICOL2) CYCLE cols

         rows: DO I = 1, NP
            IBLCK = (I - 1) / MB
            IBASE = IBLCK * MB
            IREST = I - IBASE

            ! global row index
            ROW = MB * MYROW + IBASE * NPROW + IREST

            IF (COL == ROW) THEN
               K = I + (J - 1) * LD
#ifndef ACC_OFFLOAD
               IF ((ABS(AIMAG(A(K))) > 1E-2_q) .AND. (IU0 >= 0)) THEN
                  WRITE(IU0,*) 'BG_CHANGE_DIAGONALE: WARNING: subspace matrix is not hermitian', AIMAG(A(K)),COL
               ENDIF
#endif
               A(K) = REAL(A(K), KIND=q)
               ! move on to the next column, since we have now found
               ! the element on the diagonal in the current one
               EXIT rows
            ENDIF

         ENDDO rows
      ENDDO cols

      PROFILING_STOP('bg_change_diagonale')
#endif
      RETURN
    END SUBROUTINE BG_CHANGE_DIAGONALE
#endif

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version is virtually identical to pDSSYEX_ZHEEVX
!> when called with the additional AMATIN_SCALA argument,
!> i.e., `BG_pDSSYEX_ZHEEVX(COMM, A, W, N)` could be replaced by
!> `pDSSYEX_ZHEEVX(COMM, ADUMMY, W, NDUMMY, N, AMATIN_SCALA=A)`.
!> A is a two dimension array here, however.
!
!=======================================================================

    SUBROUTINE BG_pDSSYEX_ZHEEVX(COMM, A, W, N)

GPACC USE moffload_struct_def

      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               !< NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))  !< input/output matrix
      REAL(q) W(N)            !< eigenvalues
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2

      GDEF    Z(MAX(GSD%NP,1),MAX(GSD%NQ,1))
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE
#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS
#endif

      PROFILING_START('bg_pdssyex_zheevx')

#ifdef CUSOLVERMP
      !$ACC ENTER DATA COPYIN(Z,W) __IF_ASYNC__
#endif
      !$ACC UPDATE SELF(A) __IF_ASYNC__
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

      CALL CHECK_scala (N, 'BG_pDSSYEX_ZHEEVX' )

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
        CALL vtutor%error("pDSSYEX_ZHEEVX: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " " &
           // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(MAX(3,GSD%LWWORK)))
      ALLOCATE(RWORK(MAX(3,GSD%LRWORK)))
DOESI IWORK = 0; WWORK = 0; RWORK = 0.0_q ! Breaks bulk_GaAs_ACFDT (intel), uninitialized values trigger SIGFPE in scaLAPACK routine.
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

      INFO=0

#ifdef ELPA
      ! call ELPA 2-stage solver
      IF (LELPA) THEN
           IF (ELPA_INIT(20171201) /= ELPA_OK) &
              CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

           ET => ELPA_ALLOCATE()

           CALL ET%SET("na",  N, ISUCCESS)
           CALL ET%SET("nev", N, ISUCCESS)

           CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
           CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

           CALL ET%SET("nblk", NB, ISUCCESS)

           CALL ET%SET("mpi_comm_parent", COMM%MPI_COMM, ISUCCESS)
           CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
           CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

           ISUCCESS = ET%SETUP()

GPACC      IF (OFFLOAD_ON) THEN
GPACC         CALL ET%SET("nvidia-gpu", 1, ISUCCESS)
GPACC      ELSE
              CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)
GPACC      ENDIF

           CALL ET%EIGENVECTORS(A, W, Z, ISUCCESS)

           CALL ELPA_DEALLOCATE(ET)
           CALL ELPA_UNINIT()
      ENDIF
#endif
      ! call scaLAPACK routine
      IF (.NOT. LELPA) THEN

#ifdef gammareal
        CALL WPDSYEVX( 'V', 'A', 'U', N, A, 1, 1, DESCSTD, ZERO, ZERO, 13,  &
                    -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, WWORK, &
                    GSD%LWWORK, IWORK, GSD%LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
        CALL WPZHEEVX( 'V', 'A', 'U', N, A, 1, 1, DESCSTD, ZERO, ZERO, 13,  &
                    -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, &
                    WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
                    IFAIL, ICLUSTR, GAP, INFO )
#endif
! a work around for a bug in cusolverMp
! should be fixed in nvhpc 24.11
#if PGI_BEFORE_24_11 && defined(CUSOLVERMP)
    IF((DESCSTD(M_) <= DESCSTD(MB_)) .OR. (DESCSTD(N_) <= DESCSTD(NB_))) THEN
       CALL M_bcast_d_from(COMM , W(1), SIZE(W), 1)
    ENDIF
#endif
#ifndef CUSOLVERMP
        IF (M.NE.N) THEN
          CALL vtutor%error("GSD%LWWORK " // str(GSD%LWWORK) // " " // str(GSD%LRWORK) // " " // &
             str(GSD%LIWORK) // " " // str(DESCSTD(M_)) // "\nERROR in subspace rotation PDSYEVX: not &
             &enough eigenvalues found " // str(M) // " " // str(N))
        ENDIF
        IF (NZ.NE.N) THEN
          CALL vtutor%error("ERROR in subspace rotation PDSYEVX: not enough eigenvectors computed " &
             // str(M) // " " // str(N))
        ENDIF
        DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
            CALL vtutor%error("ERROR in subspace rotation PDSYEVX: I2,IFAIL= " // str(I2) // " " // &
               str(IFAIL(I2)))
          ENDIF
        ENDDO
        IF (INFO<0) THEN
           CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
        ELSE IF (INFO>0) THEN
           IF (MOD(INFO,2).NE.0) THEN
              CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
           ELSEIF (MOD(INFO/2,2).NE.0) THEN
              ! this error condition happens very often, but does not imply that
              ! the result is not usefull
              ! IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
           ELSE
              CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
           ENDIF
        ENDIF
#endif
      ENDIF

#ifdef CUSOLVERMP
      !$ACC KERNELS PRESENT(A, Z) __IF_ASYNC__
#endif
      A=Z
#ifdef CUSOLVERMP
      !$ACC END KERNELS
      !$ACC EXIT DATA COPYOUT(W,Z) __IF_ASYNC__
#endif

#ifdef CUSOLVERMP
      !$ACC UPDATE SELF(A) __IF_ASYNC__
#else
      !$ACC UPDATE DEVICE(A) __IF_ASYNC__
#endif

   ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)

!$ACC UPDATE DEVICE(A) __IF_ASYNC__

      PROFILING_STOP('bg_pdssyex_zheevx')

      RETURN
    END SUBROUTINE BG_pDSSYEX_ZHEEVX

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version is virtually identical to pDSSYED_ZHEEVD
!> when called with the additional AMATIN_SCALA argument,
!> i.e., `BG_pDSSYED_ZHEEVD(COMM, A, W, N)` could be replaced by
!> `pDSSYED_ZHEEVD(COMM, ADUMMY, W, NDUMMY, N, AMATIN_SCALA=A)`
!> A is a two dimension array here, however.
!
!=======================================================================

    SUBROUTINE BG_pDSYEV_ZHEEVD(COMM, A, W, N)

GPACC USE moffload_struct_def

      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               !< NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))  !< input/output matrix
      REAL(q) W(N)            !< eigenvalues
! local variables
      INTEGER NODE_ME,IONODE

      PROFILING_START('bg_pdsyev_zheevd')

!$ACC UPDATE SELF(A) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

      CALL CHECK_scala (N, 'BG_pDSYEV_ZHEEVD' )

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
        CALL vtutor%error("pDSYEV_ZHEEVD: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " " &
           // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

        CALL PDSYEV_ZHEEVD_DESC2(A, W, N, DESCSTD, COMM)
      ENDIF calc

#ifdef CUSOLVERMP
      !$ACC UPDATE SELF(A) __IF_ASYNC__
#else
      !$ACC UPDATE DEVICE(A) __IF_ASYNC__
#endif

      PROFILING_STOP('bg_pdsyev_zheevd')
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      RETURN
    END SUBROUTINE

!=======================================================================
!
!> Matrix LU decompositions scaLAPACK aware version
!>
!> Determine A = U+ U (POTRF) and invert U (TRTRI)
!
!=======================================================================

    SUBROUTINE BG_pPOTRF_TRTRI(A, N, INFO)

GPACC USE moffload_struct_def

      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      INTEGER N                              !< NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1)) !< input/output matrix
      INTEGER INFO

#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS
#endif

      PROFILING_START('bg_ppotrf_trtri')


      INFO=0

      CALL CHECK_scala (N, 'BG_pPOTRF_TRTRI' )

#ifdef ELPA
      ! call ELPA cholesky and invert
      IF (LELPA) THEN
         IF (ELPA_INIT(20171201) /= ELPA_OK) &
            CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

         ET => ELPA_ALLOCATE()

         CALL ET%SET("na",  N, ISUCCESS)
         CALL ET%SET("nev", N, ISUCCESS)

         CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
         CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

         CALL ET%SET("nblk", NB, ISUCCESS)

         CALL ET%SET("mpi_comm_parent", ELPA_MPI_COMM_ALL, ISUCCESS)
         CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
         CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

         ISUCCESS = ET%SETUP()

GPACC    IF (OFFLOAD_ON) CALL ET%SET("nvidia-gpu", 1, ISUCCESS)

         CALL ET%CHOLESKY(A, ISUCCESS)
         CALL ET%INVERT_TRIANGULAR(A, ISUCCESS)

         CALL ELPA_DEALLOCATE(ET)
         CALL ELPA_UNINIT()
      ENDIF
#endif

#if ! defined(CUSOLVERMP)
!$ACC UPDATE SELF(A) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
#endif

      IF (.NOT. LELPA) THEN
#ifdef gammareal
         CALL WPDPOTRF &
#else
         CALL WPZPOTRF &
#endif
           & ('U',N, A(1,1),1,1,DESCSTD,INFO)

#if ! defined(CUSOLVERMP)
!$ACC UPDATE DEVICE(A) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q)
#endif
#if ! defined(CUBLASMP)
!$ACC UPDATE SELF(A) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q)
#endif

         IF (INFO/=0) RETURN

#ifdef gammareal
         CALL WPDTRTRI &
#else
         CALL WPZTRTRI &
#endif
             & ('U','N',N,A(1,1),1,1,DESCSTD,INFO)
      ENDIF
#if ! defined(CUBLASMP)
!$ACC WAIT(ACC_ASYNC_Q)
!$ACC UPDATE DEVICE(A) __IF_ASYNC__
#endif

      PROFILING_STOP('bg_ppotrf_trtri')

      END SUBROUTINE BG_pPOTRF_TRTRI

!=======================================================================
!
!> Matrix LU decompositions scaLAPACK aware version
!>
!> Determine A = U+ U (POTRF) and invert A = U+U (POTRI)
!
!=======================================================================

    SUBROUTINE BG_pPOTRF_POTRI(A, N, INFO)
      IMPLICIT NONE

      INTEGER N                              ! NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1)) ! input/output matrix
      INTEGER INFO

      PROFILING_START('bg_ppotrf_potri')

      CALL CHECK_scala (N, 'BG_pPOTRF_POTRI' )
#ifdef gammareal
         CALL PDPOTRF &
#else
         CALL PZPOTRF &
#endif
          & ('U',N, A(1,1),1,1,DESCSTD,INFO)
#ifdef gammareal
         CALL PDPOTRI &
#else
         CALL PZPOTRI &
#endif
          & ('U',N, A(1,1), 1, 1, DESCSTD, INFO)

      PROFILING_STOP('bg_ppotrf_potri')

     END SUBROUTINE BG_pPOTRF_POTRI



!=======================================================================
!
! The following routines are called with descriptors DESCA.
! The routines can be called from any context or subgroup
! of MPI processed as long as DESCA is properly set up.
! They also do not refer to GSD. In that sense they are much safer to use, but
! matter of fact the caller must "know" scalapack.
! For the time being the check however that the passed matrix dimensions are
! consistent with those set up in GSD by calling CHECK_scala. This should be
! removed after more extensive testing.
!
!=======================================================================

    SUBROUTINE DISTRI(AMATIN, NDIM, N, A, DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEF    AMATIN(NDIM,N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("DISTRI: leading dimension of matrix too small", __FILE__, __LINE__)
      END IF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI

    !> same as #DISTRI using a single precision input array

    SUBROUTINE DISTRI_SINGLE(AMATIN, NDIM, N, A,DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM,N)
      GDEF    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("internal error in scala.F: leading dimension of matrix too small", __FILE__, __LINE__)
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE

    !> Same as #DISTRI using a single precision input array

    SUBROUTINE DISTRI_SINGLE_SINGLE(AMATIN, NDIM, N, A,DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM,N)
      GDEFS   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("internal error in scala.F: leading dimension of matrix too small", __FILE__, __LINE__)
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE_SINGLE


!=======================================================================
!
!> Merge distributed matrix into one large matrix
!
!=======================================================================

    SUBROUTINE RECON(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEF    AMATIN(NDIM, N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0._q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON

    !> Same as #RECON using a single precision output array

    SUBROUTINE RECON_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM, N)
      GDEF    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE


    !> Same as #RECON using a single precision output array

    SUBROUTINE RECON_SINGLE_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM, N)
      GDEFS   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE_SINGLE

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version expects a proper scaLAPACK descriptor.
!> @note This routine does stop when eigenvectors cannot be
!> determined
!
!=======================================================================

     SUBROUTINE PDSYEVX_ZHEEVX_DESC(A, W, N, DESCA_, COMM, JOBIN, ABSTOL_)
#if defined ELPA
       USE iso_c_binding
#endif
GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF, TARGET :: A(:)    !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
       REAL(q), OPTIONAL :: ABSTOL_
! local variables
       REAL(q) :: ABSTOL
       INTEGER ICLUSTR(2*COMM%NCPU)
       INTEGER IFAIL(N)
       REAL(q) :: GAP(COMM%NCPU)
       REAL(q) :: ZERO=0.0_q
       INTEGER INFO
       INTEGER M,NZ
       INTEGER I2

       GDEF, ALLOCATABLE, TARGET ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

#ifdef ELPA
       CLASS(elpa_t), POINTER :: ET
       INTEGER :: ISUCCESS

       GDEF, POINTER :: EP(:,:), EZ(:,:)
       INTEGER NP,NQ
       INTEGER MYROW, MYCOL, NPROW, NPCOL
       INTEGER    NUMROC
       EXTERNAL   NUMROC,BLACS_GRIDINFO
#endif

       PROFILING_START('pdsyevx_zheevx_desc')

       JOB='V'  ! default is eigenvalues and eigenvectors

       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(3), WWORK(3), IWORK(1))
DOESI  RWORK = 0.0_q; WWORK = 0; IWORK = 0
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
DOESI  Z = 0
       !$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF

#if defined ELPA
       CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
       NP = NUMROC(DESCA_(M_),DESCA(MB_),MYROW,0,NPROW)
       NQ = NUMROC(DESCA_(N_),DESCA(NB_),MYCOL,0,NPCOL)

       ! call ELPA 2-stage solver
       IF (LELPA) THEN
         IF (ELPA_INIT(20171201) /= ELPA_OK) &
           CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

         CALL c_f_pointer(c_loc(A(1)),EP,[NP,NQ])
         CALL c_f_pointer(c_loc(Z),EZ,[NP,NQ])

         ET => ELPA_ALLOCATE(ISUCCESS)

         CALL ET%SET("na",  N, ISUCCESS)
         CALL ET%SET("nev", N, ISUCCESS)

         CALL ET%SET("local_nrows", NP, ISUCCESS)
         CALL ET%SET("local_ncols", NQ, ISUCCESS)

         CALL ET%SET("nblk", NB, ISUCCESS)

         CALL ET%SET("mpi_comm_parent", COMM%MPI_COMM, ISUCCESS)
         CALL ET%SET("process_row", MYROW, ISUCCESS)
         CALL ET%SET("process_col", MYCOL, ISUCCESS)

!         CALL ET%SET("gpu", 1, ISUCCESS)
!         CALL ET%SET("debug", 1, ISUCCESS)
!         CALL ET%SET("verbose", 1, ISUCCESS)

!         CALL ET%SET("gpu", 1, ISUCCESS)
!         call ET%SET("use_gpu_id", int(gpuID,kind=c_int), ISUCCESS)
!         CALL ET%SET("gpu_cholesky",1, ISUCCESS)
!         CALL ET%SET("nvidia-gpu", 1, ISUCCESS)
!         CALL ET%SET("gpu_hermitian_multiply",1, ISUCCESS)

         CALL ET%SET("debug", 1, ISUCCESS)
         CALL ET%SET("verbose", 1, ISUCCESS)

         ISUCCESS = ET%SETUP()

!        CALL ET%PRINT_SETTINGS(ISUCCESS)

         CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)
         CALL ET%EIGENVECTORS(EP, W, EZ, ISUCCESS)

         CALL ELPA_DEALLOCATE(ET)
         CALL ELPA_UNINIT()

         NULLIFY(EP,EZ)
       ENDIF
#else
       ! call scaLAPACK routine
       IF (.NOT. LELPA) THEN
#ifdef gammareal
         CALL WPDSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, WWORK, &
              LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )

         RWORK(1)=1
         WWORK(1)=WWORK(1)+NCLUST*N
#else
         CALL WPZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
              WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
              IFAIL, ICLUSTR, GAP, INFO )
         RWORK(1)=RWORK(1)+NCLUST*N
#endif
         LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
         LIWORK= MAX(1,IWORK(1))
         LRWORK= MAX(1,NINT(RWORK(1)))

         DEALLOCATE(RWORK, WWORK, IWORK)
         ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))
DOESI    RWORK = 0.0_q; WWORK = 0; IWORK = 0


#ifdef gammareal
         CALL WPDSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, WWORK, &
              LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
         CALL WPZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
              WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
              IFAIL, ICLUSTR, GAP, INFO )
#endif
! a work around for a bug in cusolverMp
! should be fixed in nvhpc 24.11
#if PGI_BEFORE_24_11 && defined(CUSOLVERMP)
    IF((DESCSTD(M_) <= DESCSTD(MB_)) .OR. (DESCSTD(N_) <= DESCSTD(NB_))) THEN
       CALL M_bcast_d_from(COMM , W(1), SIZE(W), 1)
    ENDIF
#endif
#ifndef ACC_OFFLOAD
         IF (M.NE.N) THEN
            CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // str(LIWORK) // &
               " " // str(DESCA(M_)) // "\nERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough &
               &eigenvalues found " // str(M) // " " // str(N))
         ENDIF
         IF (NZ.NE.N) THEN
            CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough eigenvectors &
               &computed " // str(M) // " " // str(N))
         ENDIF
         DO I2=1,N
            IF (IFAIL(I2).NE.0) THEN
               CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: I2,IFAIL= " // str(I2) &
                  // " " // str(IFAIL(I2)))
            ENDIF
         ENDDO
         IF (INFO<0) THEN
            CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
         ELSE IF (INFO>0) THEN
            IF (MOD(INFO,2).NE.0) THEN
               CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ELSEIF (MOD(INFO/2,2).NE.0) THEN
            ! this error condition happens very often, but does not imply that
            ! the result is not usefull
            !  IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
            ELSE
               CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ENDIF
         ENDIF
#endif
       ENDIF ! IF (.NOT. LELPA) THEN
#endif

      !$ACC WAIT IF(OFFLOAD_ON)
      !$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
             A=Z
      !$ACC END KERNELS
      !$ACC EXIT DATA DELETE(Z) __IF_ASYNC__

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       PROFILING_STOP('pdsyevx_zheevx_desc')

     END SUBROUTINE PDSYEVX_ZHEEVX_DESC

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version expects a proper scaLAPACK descriptor.
!> @note This routine does stop when eigenvectors cannot be
!> determined
!>
!
!=======================================================================

     SUBROUTINE PSSYEVX_CHEEVX_DESC(A, W, N, DESCA_, COMM, JOBIN, ABSTOL_)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEFS    A(:)           !< input/output matrix
       REAL(qs) W(N)           !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
       REAL(q), OPTIONAL :: ABSTOL_
! local variables
       REAL(qs) :: ABSTOL
       INTEGER ICLUSTR(2*COMM%NCPU)
       INTEGER IFAIL(N)
       REAL(qs) :: GAP(COMM%NCPU)
       REAL(qs) :: ZERO=0.0_q
       INTEGER INFO
       INTEGER M,NZ
       INTEGER I2

       GDEFS, ALLOCATABLE ::    Z(:)
       REAL(qs), ALLOCATABLE::  RWORK(:) ! work array
       GDEFS, ALLOCATABLE   ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors

       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(3), WWORK(3), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
      !$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPSSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)+NCLUST*N
#else
       CALL WPCHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=RWORK(1)+NCLUST*N
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL WPSSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
       CALL WPCHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
#endif
! a work around for a bug in cusolverMp
! should be fixed in nvhpc 24.11
#if PGI_BEFORE_24_11 && defined(CUSOLVERMP)
    IF((DESCSTD(M_) <= DESCSTD(MB_)) .OR. (DESCSTD(N_) <= DESCSTD(NB_))) THEN
       CALL M_bcast_s_from(COMM , W(1), SIZE(W), 1)
    ENDIF
#endif
#ifndef ACC_OFFLOAD
       IF (M.NE.N) THEN
          CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // str(LIWORK) // &
             " " // str(DESCA(M_)) // "\nERROR in subspace rotation PSSYEVX/ PCHEEVX: not enough &
             &eigenvalues found " // str(M) // " " // str(N))
       ENDIF
       IF (NZ.NE.N) THEN
          CALL vtutor%error("ERROR in subspace rotation PSSYEVX/ PCHEEVX: not enough eigenvectors &
             &computed " // str(M) // " " // str(N))
       ENDIF
       DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
             CALL vtutor%error("ERROR in subspace rotation PSSYEVX/ PCHEEVX: I2,IFAIL= " // str(I2) &
                // " " // str(IFAIL(I2)))
          ENDIF
       ENDDO
       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (MOD(INFO,2).NE.0) THEN
             CALL vtutor%error("ERROR eigenvector not converged PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
          ELSEIF (MOD(INFO/2,2).NE.0) THEN
          ! this error condition happens very often, but does not imply that
          ! the result is not usefull
!             IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PSSYEVX/ PCHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
          ELSE
             CALL vtutor%error("ERROR in diagonalization PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
          ENDIF
       ENDIF
#endif

       !$ACC WAIT IF(OFFLOAD_ON)
       !$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
       !$ACC END KERNELS
       !$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

     END SUBROUTINE PSSYEVX_CHEEVX_DESC


!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> This version expects a proper scaLAPACK descriptor.
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:)            !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ )  !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEF, ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PDSYEV_ZHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PDSYEV_ZHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
!$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPDSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))
#ifdef gammareal
       CALL WPDSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

!$ACC WAIT IF(OFFLOAD_ON)
!$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSYEV_ZHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PDSYEV_ZHEEVD_DESC' )
    END SUBROUTINE PDSYEV_ZHEEVD_DESC

!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> Same as #PDSYEV_ZHEEVD_DESC, but complex version regardless of gammareal
!
!=======================================================================

     SUBROUTINE PZHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       COMPLEX(q) A(:)            !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ )  !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PZHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PZHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PZHEEVD failed using fallback PZHEEV"
          CALL PZHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PZHEEVD_DESC' )
    END SUBROUTINE PZHEEVD_DESC

!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> This version uses a two-dimensional input array.
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEVD_DESC2(A, W, N, DESCA_, COMM, JOBIN)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:,:)          !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEF, ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

       !$ACC ENTER DATA COPYIN(Z,W) __IF_ASYNC__

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

#ifdef CUSOLVERMP
      !$ACC KERNELS PRESENT(A, Z) __IF_ASYNC__
#endif
      A=Z
#ifdef CUSOLVERMP
      !$ACC END KERNELS
      !$ACC EXIT DATA DELETE(Z) COPYOUT(W) __IF_ASYNC__
#endif
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSYEV_ZHEEV_DESC2(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored in
! single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))


    END SUBROUTINE PDSYEV_ZHEEVD_DESC2

!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> Single precision version.
!
!=======================================================================

     SUBROUTINE PSSYEV_CHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
GPACC  USE moffload_struct_def
       IMPLICIT NONE

       INTEGER N                !< dimension of matrix
       GDEFS    A(:)            !< input/output matrix
       REAL(qs) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ )  !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEFS, ALLOCATABLE ::     Z(:)
       REAL(qs), ALLOCATABLE ::  RWORK(:) ! work array
       GDEFS, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PSSYEV_CHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PSSYEV_CHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
!$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPSSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL WPCHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       ! allocate optimal work space
       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL WPSSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL WPCHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

!$ACC WAIT IF(OFFLOAD_ON)
!$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEVD/ PCHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PSSYEVD/ PCHEEVD failed using fallback PSSYEV/ PCHEEV"
          CALL PSSYEV_CHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PSSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_s_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PSSYEV_CHEEVD_DESC' )
    END SUBROUTINE PSSYEV_CHEEVD_DESC

!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> Usually not very fast but robust
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:)            !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEF, ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START('pdsyev_zheev_desc')

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
!$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPDSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL WPZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL WPDSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL WPZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= " // str(INFO))
       ENDIF

!$ACC WAIT IF(OFFLOAD_ON)
!$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

       PROFILING_STOP('pdsyev_zheev_desc')

    END SUBROUTINE PDSYEV_ZHEEV_DESC

!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> Complex version regardless of gammareal
!
!=======================================================================

    SUBROUTINE PZHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       COMPLEX(q)    A(:)      !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PZHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PZHEEV_DESC

!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> This version uses a two-dimensional input array.
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEV_DESC2(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:,:)          !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEF, ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PDSYEV_ZHEEV_DESC2


!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> Single precision version
!
!=======================================================================

     SUBROUTINE PSSYEV_CHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEFS   A(:)            !< input/output matrix
       REAL(qs) W(N)           !< eigenvalues
       INTEGER :: DESCA_( DLEN_ )
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEFS, ALLOCATABLE ::     Z(:)
       REAL(qs), ALLOCATABLE :: RWORK(:) ! work array
       GDEFS, ALLOCATABLE   ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

#ifdef gammareal
       CALL PSSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL PCHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL PSSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL PCHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEV/ PCHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PSSYEV_CHEEV_DESC


!=======================================================================
!
!> Makes an (approximately unitary) transformation matrix
!> exactly unitary (or orthogonal)
!>
!> If the orbitals (observing `<n|m>= delta_mn` )
!> are transformed by a matrix A, they will no longer be
!> exactly orthonormal
!> ~~~
!> |n'> = |n> A_nn'  (see LINCOM)
!> ~~~
!> but rather observe
!> ~~~
!> <n'| m'> =  A^+_n'm A_mm' = (A+ A)_n'm'
!> ~~~
!> The routine calculate `S=  A+ A` then determines
!>    `S= U+ U`  (U upper triangular matrix)
!> inverts U and updates the transformation matrix to
!>    `A <-  A U^-1`
!> one can easily show that
!> ~~~
!>(A U^-1)^+ (A U^-1) = 1
!> ~~~
!
!=======================================================================

     SUBROUTINE MAKE_UNITARY_DESC(A, N, DESCA )
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:)            !< input/output matrix
       INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
       INTEGER :: INFO
! local variables
       GDEF, ALLOCATABLE ::     Z(:)

       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       ! calculate S= A+ A
       CALL PGGEMM( trans, 'N', N, N, N, one,  &
          A(1), 1, 1, DESCA, &
          A(1), 1, 1, DESCA, zero,  &
          Z(1), 1, 1, DESCA )

       ! determine S= U+ U  (U upper triangular matrix)
#ifdef gammareal
       CALL PDPOTRF &
#else
       CALL PZPOTRF &
#endif
          & ('U', N, Z(1), 1, 1, DESCA, INFO)

       IF (INFO/=0) THEN
          CALL vtutor%bug("internal error in VASP: MAKE_UNITARY_DESC PZPOTRF failed " // str(INFO), __FILE__, __LINE__)
       ENDIF
       ! invert U <- U^-1
#ifdef gammareal
       CALL PDTRTRI &
#else
       CALL PZTRTRI &
#endif
          & ('U','N', N, Z(1), 1, 1, DESCA, INFO)
       IF (INFO/=0) THEN
          CALL vtutor%bug("internal error in VASP: MAKE_UNITARY_DESC PZTRTRI failed " // str(INFO), __FILE__, __LINE__)
       ENDIF
! triangular update of matrix A <- A U^-1
! note that the lower part of Z is set so DGEMM does not work here
#ifdef gammareal
       CALL PDTRMM &
#else
       CALL PZTRMM &
#endif
     &                ('R', 'U', 'N', 'N' , N, N, one, &
     &                 Z, 1, 1, DESCA, A, 1, 1, DESCA)

       DEALLOCATE(Z)

     END SUBROUTINE MAKE_UNITARY_DESC

!=======================================================================
!
!> Add diagonal elements AD (usually eigenvalues) to the supplied
!> distributed matrix
!
!=======================================================================


    SUBROUTINE ADD_TO_DIAGONALE(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(:)        !< distributed matrix
      COMPLEX(q) AD(:)       !< diagonal elements to be added to A
      INTEGER N              !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE

!
!> Similary to #SET_DIAGONALE for a real valued A
!> Also, A can be of any shape as long as compatible with DESCA
!

    SUBROUTINE SET_DIAGONALE_D(N, A, AD, DESCA )
GPACC USE moffload_struct_def
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      REAL(q) A(*)              !< distributed matrix
      REAL(q) AD(*)          !< diagonal elements to be added to A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST, K
      INTEGER MB, NB, LD, ROW, COL

      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      MB = DESCSTD(MB_ )
      NB = DESCSTD(NB_ )
      LD = DESCSTD(LLD_)

!$ACC PARALLEL LOOP PRESENT(A,AD) PRIVATE(JBLCK,JBASE,JREST,COL,I,IBLCK,IBASE,IREST,ROW,K) __IF_ASYNC__
      cols: DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         rows: DO I = 1, NP
            IBLCK = (I - 1) / MB
            IBASE = IBLCK * MB
            IREST = I - IBASE

            ! global row index
            ROW = MB * MYROW + IBASE * NPROW + IREST

            IF (COL == ROW) THEN
               K = I + (J - 1) * LD
               A(K) = AD(ROW)
               ! move on to the next column, since we have now found
               ! the element on the diagonal in the current one
               EXIT rows
            ENDIF

         ENDDO rows
      ENDDO cols

      RETURN
    END SUBROUTINE SET_DIAGONALE_D

!
!> Similary to #SET_DIAGONALE for a complex valued A
!> Also, A can be of any shape as long as compatible with DESCA
!

    SUBROUTINE SET_DIAGONALE_Z(N, A, AD, DESCA )
GPACC USE moffload_struct_def
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(*)              !< distributed matrix
      REAL(q) AD(*)          !< diagonal elements to be added to A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST, K
      INTEGER MB, NB, LD, ROW, COL

      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      MB = DESCSTD(MB_ )
      NB = DESCSTD(NB_ )
      LD = DESCSTD(LLD_)

!$ACC PARALLEL LOOP PRESENT(A,AD) PRIVATE(JBLCK,JBASE,JREST,COL,I,IBLCK,IBASE,IREST,ROW,K) __IF_ASYNC__
      cols: DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         rows: DO I = 1, NP
            IBLCK = (I - 1) / MB
            IBASE = IBLCK * MB
            IREST = I - IBASE

            ! global row index
            ROW = MB * MYROW + IBASE * NPROW + IREST

            IF (COL == ROW) THEN
               K = I + (J - 1) * LD
               A(K) = AD(ROW)
               ! move on to the next column, since we have now found
               ! the element on the diagonal in the current one
               EXIT rows
            ENDIF

         ENDDO rows
      ENDDO cols

      RETURN
    END SUBROUTINE SET_DIAGONALE_Z

!
!> Same as #ADD_TO_DIAGONALE for a real valued function AD
!

    SUBROUTINE ADD_TO_DIAGONALE_REAL(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)       !< distributed matrix
      REAL(q) AD(:)      !< diagonal elements to be added to A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE_REAL

!
!> Similary to #ADD_TO_DIAGONALE for a GDEF valued function AD
!> Also, A can be of any shape as long as compatible with DESCA
!
       SUBROUTINE ADD_TO_DIAGONALE_GDEF(N, A, AD, DESCA )
         IMPLICIT NONE
         INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
         GDEF A(*)              !< distributed matrix
         REAL(q) AD(*)          !< diagonal elements to be added to A
         INTEGER N              !< dimension of the matrix
         ! local
         INTEGER NP,NQ
         INTEGER I1RES,J1RES,JCOL
         INTEGER MYROW, MYCOL, NPROW, NPCOL
         INTEGER I1,I2,J1,J2
   
         INTRINSIC  MIN
         INTEGER    NUMROC,IROW_JCOL
         EXTERNAL   NUMROC,BLACS_GRIDINFO
   
         CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
         NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
         NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
   
         JCOL = 0
         ! loop over local columns
         DO J1=1,NQ,DESCA(NB_)
           J1RES=MIN(DESCA(NB_),NQ-J1+1)
           DO J2=1,J1RES
             IROW_JCOL = JCOL
   
             ! loop over local rows
             DO I1=1,NP,DESCA(MB_)
               I1RES=MIN(DESCA(MB_),NP-I1+1)
               DO I2=1,I1RES
                 IROW_JCOL = IROW_JCOL+1
                 IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                      (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                   A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                   GOTO 100 ! next column, since we have found the right row
                 ENDIF
               ENDDO
             ENDDO
100          JCOL = JCOL + DESCA(LLD_)
           ENDDO
         ENDDO
   
         RETURN
       END SUBROUTINE ADD_TO_DIAGONALE_GDEF

!
!> Same as #ADD_TO_DIAGONALE_REAL, but set diagonal part
!

    SUBROUTINE SET_DIAGONALE_REAL(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)       !< distributed matrix
      REAL(q) AD(:)      !< diagonal elements to be subtracted from A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                IF (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2 <= SIZE(AD) ) THEN
                   A(IROW_JCOL) = AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                ENDIF
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_DIAGONALE_REAL


!=======================================================================
!
!> Determine diagonal elements AD of the supplied distributed matrix A
!>
!> @attention Global sum is required after the call!
!
!=======================================================================


    SUBROUTINE DETERMINE_DIAGONALE(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(:)    !< distributed matrix
      COMPLEX(q) AD(:)   !< on exit: diagonal elements
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DETERMINE_DIAGONALE

!
!> Same as #DETERMINE_DIAGONALE but for a GDEF input array (`REAL(q)` for gammaonly)
!

    SUBROUTINE DETERMINE_DIAGONALE_GDEF(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)       ! distributed matrix
      GDEF   AD(:)       ! on exit: diagonal elements
      INTEGER N          ! dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO
      PROFILING_START( 'determine_diagonale_gdef' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      PROFILING_STOP( 'determine_diagonale_gdef' )
    END SUBROUTINE DETERMINE_DIAGONALE_GDEF

!
!> Similary to #DETERMINE_DIAGONALE but for a GDEF input array and real output array AD
!> should be called for hermitian matrices A only
!> also A can have any shape that is compatible with DESCA 
!
    SUBROUTINE DETERMINE_DIAGONALE_GDEF_REAL(N, A, AD, DESCA)
      IMPLICIT NONE
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF A(*)    !< distributed matrix
      REAL(q) AD(:)   !< on exit: diagonal elements
      INTEGER N          !< dimension of the matrix
      !local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
   
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO
   
      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
   
      AD=0
   
      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
   
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO
   
    END SUBROUTINE DETERMINE_DIAGONALE_GDEF_REAL


!=========================================================================
!> Computes a Taylor expansion
!>
!> Computes the Taylor expansion of
!> ~~~
!> ( f_r - f_s ) / ( e_r - e_s )
!> ~~~
!> by expanding `f_s` around `e_r`, that is
!> ~~~
!> f_r ~ f_s + f'_s( e_r - e_s ) + f''_s( e_r - e_s)^2 / 2
!> ~~~
!=========================================================================
    FUNCTION TAYLOR_FR_FS_OVER_ER_ES( F_R, F_S, E_R, E_S, BETA )
        REAL(q),INTENT(IN) :: F_R
        REAL(q),INTENT(IN) :: F_S
        REAL(q),INTENT(IN) :: E_R
        REAL(q),INTENT(IN) :: E_S
        REAL(q),INTENT(IN) :: BETA
        REAL(q) :: TAYLOR_FR_FS_OVER_ER_ES
        ! local

        TAYLOR_FR_FS_OVER_ER_ES = BETA*F_S*(1-F_S)*&
          (BETA*(E_R-E_S)/2*(1-2*F_S)-1)

    END FUNCTION TAYLOR_FR_FS_OVER_ER_ES

!=========================================================================
!> Computes a Taylor expansion
!>
!> Computes the Taylor expansion of
!> ~~~
!> ( f_r*e_r - f_s*e_s ) / ( e_r - e_s )
!> ~~~
!> by expanding `f_s` around `e_r`, that is
!> ~~~
!> f_r ~ g_s - beta*f_s( 1 - f_s )*( e_r = e_m * beta /2 ( e_m - e_n ) * (1-2*f_s)
!> ~~~
!=========================================================================
    FUNCTION TAYLOR_FR_ER_FS_ES_OVER_ER_ES( F_R, F_S, E_R, E_S, BETA )
        REAL(q),INTENT(IN) :: F_R
        REAL(q),INTENT(IN) :: F_S
        REAL(q),INTENT(IN) :: E_R
        REAL(q),INTENT(IN) :: E_S
        REAL(q),INTENT(IN) :: BETA
        REAL(q) :: TAYLOR_FR_ER_FS_ES_OVER_ER_ES
        ! local

        TAYLOR_FR_ER_FS_ES_OVER_ER_ES = F_S - BETA * F_S * (1-F_S)*&
           ( E_R + BETA*E_S * ( E_S - E_R ) * ( 1 - 2*F_S)/2 )
    END FUNCTION TAYLOR_FR_ER_FS_ES_OVER_ER_ES

!=======================================================================
!
!> Divide distributed matrix by eigenvalue difference
!>
!> ~~~
!>   (       |       )         (                |                )
!>   (       |       )         ( A_ij*(f_i-f_j) | A_ia*(f_i-f_a) )
!>   (  A_ij | A_ia  )         ( -------------- | -------------- )
!>   (       |       )         (    e_i -e_j    |    e_i -e_a    )
!>   (_______|_______)  ----\  (________________|________________)
!>   (       |       )  ----/  (                |                )
!>   (       |       )         ( A_ai*(f_a-f_i) | A_ab*(f_a-f_b) )
!>   ( A_ai  | A_ab  )         ( -------------- | -------------- )
!>   (       |       )         (    e_a -e_i    |    e_a -e_b    )
!>   (       |       )         (                |                )
!> ~~~
!>
!> This is usually required in standard perturbation theory.
!> The operation is only done in the occupied-unoccupied part
!> of the matrix.
!> The first version is used to construct the density matrix
!> from a given perturbation.
!> The second version is used to construct an approximately
!> unitary rotation matrix.
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF(N, A, CELTOT, FERTOT, DESCA, &
      LFINITE_T, AUXTOT )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )!< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      COMPLEX(q) CELTOT(:)  !< eigenvalues
      REAL(q)    FERTOT(:)  !< occupancies
      INTEGER N             !< dimension of the matrix
      LOGICAL :: LFINITE_T  !< finite temperature mode
      REAL(q)    AUXTOT(:)  !< derivatives of occupancies
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW_JCOL = IROW_JCOL+1

               IF ( LFINITE_T ) THEN
                  CALL SET_A_RS_AT_T( A(IROW_JCOL),        &
                     FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),        &
                     REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),q),&
                     REAL(CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q),&
                     AUXTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     AUXTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
               ELSE
                  CALL SET_A_RS_AT_T0( A(IROW_JCOL),       &
                     FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),        &
                     REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),q),&
                     REAL(CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q) )

               ENDIF

            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO
      RETURN

      CONTAINS

      ! ------------------------------------------------------------------------
      !> Gives matrix element at T=0
      !>
      !> F_R and F_S has integer values 0 and 1 only
      ! ------------------------------------------------------------------------
      SUBROUTINE SET_A_RS_AT_T0( A, FI, FJ, EI, EJ )
         GDEF, INTENT(INOUT) :: A
         REAL(q),INTENT(IN)  :: FI
         REAL(q),INTENT(IN)  :: FJ
         REAL(q),INTENT(IN)  :: EI
         REAL(q),INTENT(IN)  :: EJ
         ! local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q
         REAL(q) :: DIFF

         ! identical occupancies, clear rotation matrix
         DIFF=( FI- FJ )

         IF (ABS(DIFF)<DIFMAX) THEN
            A = 0
         ELSE
            ! (e_i - e_j)/(f_i-f_j)
            DIFF = ( EI - EJ ) / DIFF

            ! huge rotation element, should never occur except for tiny
            ! smearing parameters of the order of DIFMAX
            IF (ABS(DIFF)<DIFMAX) THEN
               A = 0
            ELSE
               A = A /DIFF
            ENDIF
         END IF

      END SUBROUTINE SET_A_RS_AT_T0

      ! ------------------------------------------------------------------------
      !> Gives matrix element at T>0
      !>
      !> FI and FJ between 0 and 1
      ! ------------------------------------------------------------------------
      SUBROUTINE SET_A_RS_AT_T( A, FI, FJ, EI, EJ, DFI, DFJ )
         GDEF, INTENT(INOUT) :: A
         REAL(q),INTENT(IN)  :: FI
         REAL(q),INTENT(IN)  :: FJ
         REAL(q),INTENT(IN)  :: EI
         REAL(q),INTENT(IN)  :: EJ
         REAL(q),INTENT(IN)  :: DFI
         REAL(q),INTENT(IN)  :: DFJ
         ! local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q

         IF (ABS( EI - EJ )<DIFMAX) THEN
            ! this was tested for small transition energies against equation below
            A = A * (DFI + DFJ)/2
         ELSE
            A =  A* ( FI- FJ ) / ( EI - EJ )
         END IF

      END SUBROUTINE SET_A_RS_AT_T

    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF


!=======================================================================
!
!> Divide matrix A by eigenvalue difference times eigenvalue (uocc_occ block
!> version)
!>
!> Divide matrix A by eigenvalue difference times eigenvalue in the
!> uocc_occ block ( i,j \in occ ; a,b \in uocc), remove uocc_uocc block
!> ~~~
!>   (       |       )         (          |           )
!>   (       |       )         (          |  e_i*A_ia )
!>   (  A_ij | A_ia  )         (   A_ij   |  -------- )
!>   (       |       )         (          |  e_i-e_a  )
!>   (_______|_______)  ----\  (__________|___________)  ---    A'
!>   (       |       )  ----/  (          |           )  ---
!>   (       |       )         ( e_i*A_ai |     0     )
!>   ( A_ai  | A_ab  )         ( -------- |    0 0    )
!>   (       |       )         ( e_i-e_a  |     0     )
!>   (       |       )         (          |           )
!> ~~~
!> In index notation:
!> ~~~
!>       A' = A_ij * ( e_i f_i - e_j f_j) / (e_i - e_j)
!>          = A_ij/2 [ f_i+f_j + (e_i+e_j) (f_i-f_j) /(e_i - e_j)]
!> ~~~
!>
!> @note At finite temperature one has 0 < f_i, f_j < 1.
!> So careful approximation of the fraction has to be done.
!> Here f is assumed to be a Fermi function.
!>
!> And add it to matrix B:
!> ~~~
!> B = B + A'
!> ~~~
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL(N, A, B, CELTOT, FERTOT, DESCA,&
      LFINITE_T, AUXTOT)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed Hamilton matrix
      GDEF    B(:)          !< distributed result matrix
      COMPLEX(q) CELTOT(:)  !< eigenvalues
      REAL(q)    FERTOT(:)  !< occupancies
      REAL(q)    AUXTOT(:)  !< occupancies
      INTEGER N             !< dimension of the matrix
      LOGICAL LFINITE_T
      REAL(q) :: EFERMI
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTRINSIC  MIN
      INTEGER    I, J
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW_JCOL = IROW_JCOL+1
               ! global column index is DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
               ! global row index DESCA(MB_)*MYROW+NPROW*(I1-1)+I2

               J = DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
               I = DESCA(MB_)*MYROW+NPROW*(I1-1)+I2

               ! compute matrix elements
               IF ( LFINITE_T ) THEN
                  CALL MATRIX_ELEMENT_AT_T( B(IROW_JCOL), A(IROW_JCOL),&
                      FERTOT(J), FERTOT(I), REAL(CELTOT(J),q), REAL(CELTOT(I),q), &
                      AUXTOT(I), AUXTOT(J))
               ELSE
                  CALL MATRIX_ELEMENT_AT_T0( B(IROW_JCOL), A(IROW_JCOL),&
                       FERTOT(J), FERTOT(I), REAL(CELTOT(J),q), REAL(CELTOT(I),q)  )
               ENDIF

            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      CONTAINS
      ! ------------------------------------------------------------------------
      !> Gives matrix element at T=0
      !>
      !> FI and FJ have integer values 0 and 1 only
      ! ------------------------------------------------------------------------
      SUBROUTINE MATRIX_ELEMENT_AT_T0( B, A, FJ, FI, EJ, EI )
         GDEF,INTENT(INOUT) :: B
         GDEF,INTENT(IN) :: A
         REAL(q), INTENT(IN) :: FJ, FI
         REAL(q), INTENT(IN) :: EJ, EI
         !local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q

         IF (ABS(EJ-EI)<DIFMAX) THEN
            ! not the exact solution but a sensible simple approximation at T=0
            B = B +  A/2 * ( FI +FJ )
         ELSE
            B = B + A * ( FJ*EJ - FI*EI ) / (EJ-EI)
         END IF
      END SUBROUTINE MATRIX_ELEMENT_AT_T0

      ! ------------------------------------------------------------------------
      !> Gives matrix element at T>0
      !>
      !> The matrix element at T>0 is obtained from
      !> ~~~
      !> sum_n i w_n g( iw_n) g( iw_n ) e^{-iw_n eta} = 1/2 ( sign( eta) -1)
      !> + (f_n e_n - f_m e_m ) / ( e_n - e_m )
      !> ~~~
      !> the actual Matsubara sum is
      !> ~~~
      !> 1/beta \sum_n ( iw_n ) G(iw_n + mu ) F G( iw_n + mu) e^{-i w_n eta)
      !>    fj ej - fi ei
      !> = ----------------
      !>       ej - ei
      !> ~~~
      !> This result follows from the evaluation of the series for ~ -w_n^2
      !> and the series for ~ -I w_n^3 Exp[-iw_n eta] (eta>0)
      !> The term proportional to -w_n^2 gives the result -1/2, while the term
      !> proportional to -I w_n^3 Exp[-iw_n 0-] yiels +1/2.
      !>
      !> For terms on the diagonal ei ~ ej
      !> make a Taylor exansion of fj ej  around ei
      !> ~~~
      !> fj ej ~ [fi + fi' (ej-ei)  ][ ei + (ej-ei) ]
      !>       ~ fi ei + fi' (ej -ei ) ei + fi (ej-ei)
      !> ~~~
      !> and insert this in the first term of the result giving
      !> ~~~
      !>    fi ei + fi' (ej -ei ) ei + fi (ej-ei) - fi ei
      !> = -----------------------------------------------  = fi' ei + fi
      !>                        ej - ei
      !> ~~~
      !> symmetrizing yields
      !> ~~~
      !>    fi + fj    fi' ei + fj' ej
      !> =  ------- +  ---------------
      !>       2              2
      !> ~~~
      !> to include off-diagonal terms the seond term is symmetrized as follows
      !> not this is valid since ei ~ ej holds true
      !> ~~~
      !>    fi + fj    fi' ei + fj' ej + fi' ej + fj' ei
      !> =  ------- +  ---------------------------------
      !>       2                       4
      !> ~~~
      ! ------------------------------------------------------------------------
      SUBROUTINE MATRIX_ELEMENT_AT_T( B, A, FJ, FI, EJ, EI, DFJ, DFI )
         GDEF,INTENT(INOUT) :: B
         GDEF,INTENT(IN) :: A
         REAL(q), INTENT(IN) :: FJ, FI
         REAL(q), INTENT(IN) :: EJ, EI
         REAL(q), INTENT(IN) :: DFJ, DFI
         !local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q
         ! energies close to each other
         IF (ABS(EJ-EI)<DIFMAX) THEN
            !B = B +  A* ( FI+FJ + EI * DFI + EJ*DFJ ) / 2
            ! taylor expansion if occupancies are close (tested against below)
            B = B +  A* ( FI +FJ + (EI +EJ) * (DFI + DFJ)/2) / 2
         ELSE
            !B = B + A*( FJ*EJ - FI*EI )/(EJ-EI)
            B = B + A*(( FJ*EJ - FI*EI )/(EJ-EI))
         ENDIF
       END SUBROUTINE MATRIX_ELEMENT_AT_T

    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL

!=======================================================================
!
!> Divide distribute matrix by eigenvalue difference
!>
!> This is usually required in standard perturbation theory.
!> The operation is only done in the occupied-unoccupied part
!> of the matrix.
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF_SIGNED(N, A, CELTOT, DESCA, LFINITE_T )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      GDEF  CELTOT(:)       !< eigenvalues
      INTEGER N             !< dimension of the matrix
      LOGICAL :: LFINITE_T
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1

              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- &
                           CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)<DIFMAX) THEN
                 A(IROW_JCOL) = 0
              ELSE
                 ! the equations above also involve
                 ! the difference between the occupancies, for the density matrix they would drop out
                 ! in other words: rotations in the (un)occupied-(un)occupied  block are irrelevant
                 ! for energies
                 ! tests indicate that inclusion of DIFCELL makes indeed little difference for energies
                 A(IROW_JCOL) = A(IROW_JCOL)/DIFCELL
                 IF (ABS(A(IROW_JCOL))>0.4_q) THEN
                    ! damp rotation matrix as in Loewdin case
                    A(IROW_JCOL)=A(IROW_JCOL)/ABS(ABS(A(IROW_JCOL)))*0.4_q
                 ENDIF
              ENDIF
              ! set diagonal to 1
              IF (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2==DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) THEN
                 A(IROW_JCOL)=1
              ENDIF


            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF_SIGNED

!=======================================================================
!
!> calculate |H_ai|^2 / (epsilon_a-epsilon_i) *(f_a-f_i)
!>
!> Both the matrix H_ai and the eigenvalues, and the occupancies are passed
!> down.
!> The calculation is only performed for |epsilon_a - epsilon_i| >DIFFMAX.
!
!=======================================================================

    SUBROUTINE SUM_SECOND_ORDER(N, A, CELTOT, FERTOT, CSUM, DESCA, LFINITE_T )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      REAL(q) CELTOT(:)     !< eigenvalues epsilon
      REAL(q) FERTOT(:)     !< occupancies
      COMPLEX(q) CSUM       !< final result |H_ai|^2 / (epsilon_a-epsilon_i)
      INTEGER N             !< dimension of the matrix
      LOGICAL :: LFINITE_T  !< finite temperature mode
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CSUM=0
      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- &
                          CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)>=DIFMAX) THEN
                 CSUM=CSUM+A(IROW_JCOL)*GCONJG(A(IROW_JCOL))/DIFCELL* &
                 (FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)-FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SUM_SECOND_ORDER

!=======================================================================
!
!> Clear occupied-occupied and unoccupied-unoccupied block of
!> a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
                     ! do nothing
               ELSE
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC

!=======================================================================
!
!> Clear occupied-unoccupied and unoccupied-occupied blocks of
!> a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC


!=======================================================================
!
!> Clear occupied-unoccupied and unoccupied-occupied blocks of
!> a square matrix using supplied fermi-weights
!>
!> Specifically the supplied matrix is multiplied by
!> ~~~
!> ferwe(ncol)*ferwe(nrow) + (1-ferwe(ncol))*(1-ferwe(nrow))
!> ~~~
!
!=======================================================================

    SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC_FER(N, A, FERTOT, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(:)        !< distributed matrix
      INTEGER N              !< dimension of the matrix
      REAL(q)    FERTOT(:)   !< occupancies
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              A(IROW_JCOL)=A(IROW_JCOL)* &
                 (FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)* FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)+ &
                  (1-FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2))*(1-FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)))
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC_FER

!=======================================================================
!
!> Set occupied-occupied block of a square matrix to identity matrix
!
!=======================================================================

    SUBROUTINE SET_OCCOCC_IDENTITY(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<= LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0._q
                  IF ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) .EQ. (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)) THEN
                     A(IROW_JCOL)=1._q
                  ENDIF
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_OCCOCC_IDENTITY

!=======================================================================
!
!> Set occupied-occupied block to a square matrix
!>
!> ~~~
!> n+1        0             0           0
!> 0          n             0           0
!> 0          0             n-1         0
!> 0          0             0           1
!> ~~~
!> This will preserve the order of the states during diagonalization
!> of the density matrix.
!
!=======================================================================

    SUBROUTINE SET_OCCOCC_INCREASING(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<= LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0._q
                  IF ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) .EQ. (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)) THEN
                     A(IROW_JCOL)=LAST_OCCUPIED+1-(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
                  ENDIF
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_OCCOCC_INCREASING

!=======================================================================
!
!> Conjugate the upper triangle of a matrix
!>
!> This routine constructs from the change of the density matrix
!> the unitary transformation for the orbitals.
!> One can also simply take the lower triangle and change the sign
!> in the upper.
!
!=======================================================================

    SUBROUTINE MINUS_UPPER(N, A, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) < &
                  (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                   A(IROW_JCOL)=-(A(IROW_JCOL))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MINUS_UPPER


!=======================================================================
!
!> Multiply distribute matrix by eigenvalue difference
!>
!> ~~~
!>   1/(iw - epsilon_i) and  1/(iw - epsilon_a)
!> ~~~
!> from the left and right
!
!=======================================================================


    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) CELTOTINV(:) !< multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE

!
!> Same as #MULTIPLY_BY_IVEIGENVALUE, but use only real part of transformation
!
    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) CELTOTINV(:) !< multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*REAL(CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL

!
!> Similar to #MULTIPLY_BY_IVEIGENVALUE, but use GDEF of A and real CELTOTINV
!> also A can have any shape that is compatible with DESCA 
!
    SUBROUTINE MULTIPLY_A_BY_CELTOTINV_GDEF(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      GDEF A(*)               !< distributed matrix
      REAL(q) CELTOTINV(:)    !< multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

    END SUBROUTINE MULTIPLY_A_BY_CELTOTINV_GDEF

!=======================================================================
!
!> Right multiply matrix by a vector
!>
!> (most likely a BLACS routines exists to do this)
!> ~~~
!> A <- A R    or  A_ij = A_ij R_j
!> ~~~
!
!=======================================================================
!
    SUBROUTINE RIGHT_MULTIPLY_BY_R(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      GDEF A(:)               !< distributed matrix
      REAL(q) R(:)            !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_R


    SUBROUTINE RIGHT_MULTIPLY_BY_C(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) R(:)         !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_C

    SUBROUTINE LEFT_MULTIPLY_BY_R(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      GDEF A(:)               !< distributed matrix
      REAL(q) R(:)            !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE LEFT_MULTIPLY_BY_R

    SUBROUTINE LEFT_MULTIPLY_BY_C(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) R(:)         !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE LEFT_MULTIPLY_BY_C

!=======================================================================
!
!> Distribute slice of the global matrix AMATIN to a distributed matrix A
!>
!> As input the routine expects several columns of the matrix
!> to be distributed.
!> The columns of the input matrix are COLUMN_LOW up to COLUMN_HIGH.
!> All rows are passed to the routine.
!
!=======================================================================
#ifndef testslice
    SUBROUTINE DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array
    ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('distri_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     WRITE(*,*) 'distri',NDIM, N

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN

             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   A(IROW+(JCOL-1)*DESCA(LLD_))= &
                   AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      PROFILING_STOP('distri_slice')

      RETURN
    END SUBROUTINE DISTRI_SLICE
#else
    SUBROUTINE DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

GPACC USE moffload_struct_def

      IMPLICIT NONE

      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST
      INTEGER ROW, COL, COL_

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

      PROFILING_START('distri_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

      MB = DESCA(MB_ )
      NB = DESCA(NB_ )
      LD = DESCA(LLD_)

!$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) __IF_ASYNC__
      DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         ! storage column index
         COL_= COL - COLUMN_LOW + 1

         ! does this column match any of the columns of amatin,
         ! the current stripe of the global matrix
         IF (COL >= COLUMN_LOW .AND. COL <= COLUMN_HIGH) THEN

            DO I = 1, NP
               IBLCK = (I - 1) / MB
               IBASE = IBLCK * MB
               IREST = I - IBASE

               ! global row index
               ROW = MB * MYROW + IBASE * NPROW + IREST

               A(I + (J - 1) * LD) = AMATIN(ROW, COL_)
            ENDDO

         ENDIF
      ENDDO

      PROFILING_STOP('distri_slice')

      RETURN
    END SUBROUTINE DISTRI_SLICE
#endif

!
!> Distribute matrix and add hermitian conjugated elements
!
!
#ifndef testslice
    SUBROUTINE DISTRI_SLICE_HERM(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('distri_slice_herm')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part first
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                       A(IROW+(JCOL-1)*DESCA(LLD_))=GCONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                     A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      PROFILING_STOP('distri_slice_herm')

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM
#else
    SUBROUTINE DISTRI_SLICE_HERM(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

GPACC USE moffload_struct_def

      IMPLICIT NONE

      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST
      INTEGER ROW, COL, ROW_, COL_

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

      PROFILING_START('distri_slice_herm')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

      MB = DESCA(MB_ )
      NB = DESCA(NB_ )
      LD = DESCA(LLD_)

      ! conjugated part
!$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(IBLCK,IBASE,IREST,ROW,ROW_,J,JBLCK,JBASE,JREST,COL) __IF_ASYNC__
      DO I = 1, NP
         IBLCK = (I - 1) / MB
         IBASE = IBLCK * MB
         IREST = I - IBASE

         ! global row index
         ROW = MB * MYROW + IBASE * NPROW + IREST

         ! storage row index
         ROW_= ROW - COLUMN_LOW + 1

         ! does this row match any of the columns of amatin,
         ! the current stripe (amatin) of the global matrix
         IF (ROW >= COLUMN_LOW .AND. ROW <= COLUMN_HIGH) THEN

            DO J = 1, NQ
               JBLCK = (J - 1) / NB
               JBASE = JBLCK * NB
               JREST = J - JBASE

               ! global column index
               COL = NB * MYCOL + JBASE * NPCOL + JREST

               A(I + (J - 1) * LD) = GCONJG( AMATIN(COL, ROW_) )
            ENDDO

         ENDIF
      ENDDO

      ! non-conjugated part
!$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) __IF_ASYNC__
      DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         ! storage column index
         COL_= COL - COLUMN_LOW + 1

         ! does this column match any of the columns of amatin,
         ! the current stripe of the global matrix
         IF (COL >= COLUMN_LOW .AND. COL <= COLUMN_HIGH) THEN

            DO I = 1, NP
               IBLCK = (I -1) / MB
               IBASE = IBLCK * MB
               IREST = I - IBASE

               ! global row index
               ROW = MB * MYROW + IBASE * NPROW + IREST

               A(I + (J - 1) * LD) = AMATIN(ROW, COL_)
            ENDDO

         ENDIF
      ENDDO

      PROFILING_STOP('distri_slice_herm')

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM
#endif

!
!> Distribute matrix and add hermitian conjugated elements
!>
!> Similar to #DISTRI_SLICE_HERM. In this version it is assumed that the lower
!> triangle of the input matrix is set. This requires some extra precaution that
!> previously properly set elements are not overwritten. Only the lower part
!> of the input matrix is used as source.
!
    SUBROUTINE DISTRI_SLICE_HERM_L(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (COL>=ROW) A(IROW+(JCOL-1)*DESCA(LLD_))=GCONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW>=COL) A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM_L

!=======================================================================
!
!> Reconstruct a slice of the originally distributed matrix AMATIN
!
!=======================================================================
#ifndef testslice
    SUBROUTINE RECON_SLICE(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      GDEF    A(*)                       !< distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('recon_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      PROFILING_STOP('recon_slice')

      RETURN
    END SUBROUTINE RECON_SLICE
#else
    SUBROUTINE RECON_SLICE(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

GPACC USE moffload_struct_def

      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      GDEF    A(*)                       !< distributed matrix

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST
      INTEGER ROW, COL, COL_

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

      PROFILING_START('recon_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      MB = DESCA(MB_ )
      NB = DESCA(NB_ )
      LD = DESCA(LLD_)

!$ACC KERNELS PRESENT(AMATIN) __IF_ASYNC__
      AMATIN=0.0_q
!$ACC END KERNELS

      ! fill stripe of global matrix
!$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) __IF_ASYNC__
      DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         ! storage column index
         COL_= COL - COLUMN_LOW + 1

         ! does this column match any of the columns of amatin,
         ! the current stripe of the global matrix
         IF (COL >= COLUMN_LOW .AND. COL <= COLUMN_HIGH) THEN

            DO I = 1, NP
               IBLCK = (I - 1) / MB
               IBASE = IBLCK * MB
               IREST = I - IBASE

               ! global row index
               ROW = MB * MYROW + IBASE * NPROW + IREST

               AMATIN(ROW,COL_) = A(I + (J - 1) * LD)
            ENDDO

         ENDIF
      ENDDO

      PROFILING_STOP('recon_slice')

      RETURN
    END SUBROUTINE RECON_SLICE
#endif

!
!> Complex version of #RECON_SLICE
!

    SUBROUTINE RECON_SLICE_C(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      COMPLEX(q) AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      COMPLEX(q) A(*)                    !< distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_C

!
!> Similar to #RECON_SLICE_C but add complex conjugated elements
!> in the lower triangle in AMATIN
!
    SUBROUTINE RECON_SLICE_HERM(AMATIN,NDIM,N,A,DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_)              !< distributed matrix descriptor array
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
      GDEF    A(*)                       !< distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'RECON_SLICE_HERM' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     set conjugated part first restrict to lower triangle
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (ROW<COL) &
                     AMATIN(COL, ROW-COLUMN_LOW+1)= &
                        GCONJG(A(IROW+(JCOL-1)*DESCA(LLD_)))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     now the non conjugated part but restrict to upper triangle
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW<=COL) &
                      AMATIN(ROW,COL-COLUMN_LOW+1)= &
                            A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_HERM

!
!> Collect specific columns of the distribted matrix
!
    SUBROUTINE RECON_SLICE_REORDER(AMATIN, NDIM, N, A, DESCA, ICOL, NCOL)

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,NCOL)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      GDEF    A(*)                       !< distributed matrix
      INTEGER NCOL                       !< number of columns to collect
      INTEGER ICOL(NCOL)                 !< global indices of the columns to be collected
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2,IC
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('recon_slice_reorder')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          DO IC=1,NCOL
            IF (ICOL(IC)==COL) THEN
               IROW=0
               DO I1=1,NP,DESCA(MB_)
                  I1RES=MIN(DESCA(MB_),NP-I1+1)
                  DO I2=1,I1RES
                     IROW=IROW+1
                     ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                     AMATIN(ROW,IC)= A(IROW+(JCOL-1)*DESCA(LLD_))
                  ENDDO
               ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      PROFILING_STOP('recon_slice_reorder')

      RETURN
    END SUBROUTINE RECON_SLICE_REORDER


!=======================================================================
!
!> Shift clear off diagonal components below NLOW
!
!=======================================================================

    SUBROUTINE SCA_SHIFT_BANDS_LOW(N, NLOW, A, SHIFT, DESCA )
      IMPLICIT NONE

      INTEGER N              !< actual dimension of the matrix
      INTEGER NLOW           !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array

 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2, JCOL_GLOBAL, IROW_GLOBAL

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      ! loop over local columns
      JCOL = 0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL = JCOL + 1
          JCOL_GLOBAL=(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)

          ! loop over local rows
          IROW = 0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW = IROW + 1
               IROW_GLOBAL=(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)

               IF (IROW_GLOBAL < NLOW .OR. JCOL_GLOBAL < NLOW ) THEN
                  IF (IROW_GLOBAL== JCOL_GLOBAL) THEN
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=A(IROW+ (JCOL-1)*DESCA(LLD_) )+SHIFT
                  ELSE
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=0.0_q
                  ENDIF
               ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SCA_SHIFT_BANDS_LOW


    SUBROUTINE SCA_SHIFT_BANDS_HIGH(N, NHIGH, A, SHIFT, DESCA )
      IMPLICIT NONE

      INTEGER N              !< dimension of the matrix
      INTEGER NHIGH          !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array

 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2, JCOL_GLOBAL, IROW_GLOBAL

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      ! loop over local columns
      JCOL = 0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL = JCOL + 1
          JCOL_GLOBAL=(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)

          ! loop over local rows
          IROW = 0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW = IROW + 1
               IROW_GLOBAL=(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)

               IF (IROW_GLOBAL > NHIGH .OR. JCOL_GLOBAL > NHIGH ) THEN
                  IF (IROW_GLOBAL== JCOL_GLOBAL) THEN
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=A(IROW+ (JCOL-1)*DESCA(LLD_) )+SHIFT
                  ELSE
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=0.0_q
                  ENDIF
               ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SCA_SHIFT_BANDS_HIGH

!=======================================================================
!
!> Performs a pseudo inverse of CHI1D
!>
!> CHI1D is passed as pointer to 1D array.
!>
!> The following lines can be used to reshape a 2D array into 1D form:
!> ~~~
!> USE, INTRINSIC :: iso_c_binding
!> GDEFCPTR, ALLOCATABLE   :: CHI1D(:)
!> CHI1D = RESHAPE( CHI_INV, (/MY_NROWS*MY_NCOLS/) )
!> ~~~
!
!=======================================================================
   SUBROUTINE PSEUDO_INVERSE_DESC( CHI1D, NP, DESC, THRESHHOLD, IU6 )
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF, POINTER        :: CHI1D(:)
      INTEGER              :: NP             !< dimension of submatrix
      INTEGER              :: DESC( DLEN_ )  !< distributed matrix descriptor array
      REAL(q)              :: THRESHHOLD
      INTEGER              :: IU6
      ! local
      INTEGER              :: J
      INTEGER              :: PINFO                !scaLAPACK routine info variable
      INTEGER              :: LWORK, LRWORK        !for scaLAPACK routine
      REAL(q)              :: E(NP)
      GDEF   ,ALLOCATABLE  :: EV(:)
      GDEF   ,ALLOCATABLE  :: V(:)
      GDEF, ALLOCATABLE    :: WORK(:)
#ifndef gammareal
      COMPLEX(q), ALLOCATABLE :: CWORK(:)             !needed for complex computation
#endif
      ALLOCATE(EV(SIZE(CHI1D)))
      ALLOCATE(V(SIZE(CHI1D)))
      ! initialize eigenvalues and eigenvectors
      E = zero
      EV = zero
      !
      !first get the optimal working arrays
      !
      ALLOCATE(WORK(1))
#ifdef gammareal
      !for real values only
      CALL PDSYEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), -1, PINFO)
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP=" // str(NP) // " 1st call of &
            &PDSYEV returns " // str(PINFO))
      ENDIF
#else
      !for a hermitian matrix
      ALLOCATE(CWORK(1))
      CALL PZHEEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1), &
           EV(1), 1, 1, DESC, WORK(1), -1, CWORK(1), -1, PINFO )
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP=" // str(NP) // " 1st call of &
            &PZHEEV returns " // str(PINFO))
      ENDIF
#endif

      !LWORK is now the optimal size of WORK
      LWORK=WORK(1)
      IF ( WORK(1)  /= 0 ) THEN   !if it's not zero reallocate WORK array
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
      ENDIF

      !diagonalize
#ifdef gammareal
      !real version
      CALL PDSYEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), LWORK , PINFO)
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP " // str(NP) // " 2nd call of &
            &PDSYEV returns " // str(PINFO))
      ENDIF
#else
      !in the complex version do the same with the CWORK array
      LRWORK=CWORK(1)
      IF ( LRWORK  /= 0 ) THEN   !if it's not zero reallocate WORK array
        DEALLOCATE(CWORK)
        ALLOCATE(CWORK(LRWORK))
      ENDIF
      !hermitian version
      CALL PZHEEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), LWORK, CWORK(1), LRWORK, PINFO )
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP " // str(NP) // " 2nd call of &
            &PZHEEV returns " // str(PINFO))
      ENDIF
      DEALLOCATE( CWORK )
#endif

      IF ( PINFO > 0 ) THEN
         WRITE(*,'(" Eigenvector #",I4," may not be converged" )') PINFO
      ENDIF

      IF (IU6>=0) WRITE(IU6,'(A)') ' eigenvalues of the response function'
      IF (IU6>=0) WRITE(IU6,'(2X,10F12.6)') E(1:NP)

      ! invert eigenvalues based on passed threshhold
      DO J=1,NP
         IF (ABS(E(J))<THRESHHOLD) THEN
            E(J)=0
         ELSE
            ! uncomment the following line to represent original matrix
            E(J)=1/E(J)
         ENDIF
      ENDDO

      ! form diagonal matrix with eigenvalues
      V = zero
      CALL ADD_TO_DIAGONALE_REAL(NP, V, E,  DESC )
      ! re-shape back to 2D array
       ! EIG * V,  multiply with rotation matrix
      CALL PGGEMM( 'N', 'N', NP, NP, NP, one,  &
         EV(1), 1, 1, DESC,              &
         V(1), 1, 1, DESC,         &
         zero, CHI1D(1), 1, 1, DESC )
      V = CHI1D
       ! 1/EIG * V,  multiply with rotation matrix
      CALL PGGEMM( 'N', trans, NP, NP, NP, one,  &
         EV(1), 1, 1, DESC,              &
         V(1), 1, 1, DESC,         &
         zero, CHI1D(1), 1, 1, DESC )

      DEALLOCATE( WORK )
      DEALLOCATE( EV )
      DEALLOCATE( V )
   END SUBROUTINE PSEUDO_INVERSE_DESC

!=======================================================================
!
!> Performs a Cholesky decomposition using the scalapack Descriptor DESC
!>
!=======================================================================
   SUBROUTINE PDPOTRF_PZPOTRF_DESC( CHI, NP,  DESC, JOBIN)
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF                 :: CHI(*)
      INTEGER              :: NP               !< dimension of submatrix
      INTEGER              :: DESC( DLEN_ )  !< distributed matrix descriptor array
      CHARACTER (LEN=1), OPTIONAL :: JOBIN
      ! local
      INTEGER              :: PINFO 
      CHARACTER (LEN=1)    :: JOB

      ! default job is upper triangular part 
      JOB = 'U'
      ! optionally pass lower triangular part 
      IF ( PRESENT( JOBIN) ) JOB = JOBIN
      !call Cholesky for gamma or complex version 
      CALL PGPOTRF( JOB, NP, CHI, 1, 1, DESC, PINFO)

      IF ( PINFO /=0 ) THEN
         CALL vtutor%error(" PGPOTRF failed for NP=" // str(NP) // " with error "// str(PINFO) )
      ENDIF

   END SUBROUTINE PDPOTRF_PZPOTRF_DESC

!=======================================================================
!
!> Determines the Frobenius Norm of a distributed matrix of dimension 
!> N x N and returns the result in R 
!> requires a MPI sum afterwards !!
!=======================================================================
   SUBROUTINE MATRIX_NORM_DESC( A, N,  DESC, R, NORMIN ) 
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF                   :: A(*)           !< distributed matrix
      INTEGER, INTENT(IN)    :: N              !< dimension of submatrix
      INTEGER                :: DESC( DLEN_ )  !< distributed matrix descriptor array
      REAL(q), INTENT(INOUT) :: R
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORMIN
      ! local
      CHARACTER(LEN=1)       :: NORM
      CHARACTER(LEN=1)       :: UPLO
      REAL(q)                :: WORK(1)
      REAL(q), EXTERNAL      :: PGLANXX

      NORM = 'F'
      IF( PRESENT( NORMIN) ) THEN
         IF( NORMIN == 'e' .OR. NORMIN == 'E' .OR. NORMIN == 'f' .OR. NORMIN == 'F' ) THEN
            NORM = NORMIN
         ELSE
            CALL VTUTOR%BUG( "requested norm "//NORMIN//" not implemented, only NORM = e,f supported", __FILE__, __LINE__ )
         ENDIF
      ENDIF

      UPLO='U'
      !frobenius norm of symmetric/hermitian real-valued matrix 
      R = PGLANXX( NORM, UPLO, N, A, 1, 1, DESC, WORK ) 
      R = R**2
   END SUBROUTINE MATRIX_NORM_DESC

!=======================================================================
!
!> Performs a matrix inversion of CHI
!>
!=======================================================================
   SUBROUTINE PDGETRI_PZGETRI_DESC( INV_EPS, NP,  DESC)
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF                 :: INV_EPS(*)
      INTEGER, INTENT(IN)  :: NP               !< dimension of submatrix
      INTEGER, INTENT(IN)  :: DESC( DLEN_ )  !< distributed matrix descriptor array
      ! local
      INTEGER              :: PINFO 
      GDEF, ALLOCATABLE    :: WORK(:)              !working array and eigenvectors 
      GDEF, ALLOCATABLE    :: RWORK
      INTEGER,ALLOCATABLE  :: IPIV(:),IWORK(:)     !for scaLAPACK
      INTEGER              :: LWORK, LIWORK        !for scaLAPACK routine 
      INTEGER              :: NPROW, NPCOL, MYROW, MYCOL 
      INTEGER              :: MY_NROWS 
      INTEGER              :: MY_NCOLS
      INTEGER, EXTERNAL    :: NUMROC 

      PROFILING_START( 'pdgetri_pzgetri_desc') 

      CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      !obtain local matrix sizes
      MY_NROWS = NUMROC(DESC(M_),DESC(MB_),MYROW,0,NPROW)
      MY_NCOLS = NUMROC(DESC(N_),DESC(NB_),MYCOL,0,NPCOL) 

      !allocate IPIV
      ALLOCATE( IPIV(MY_NROWS+DESC(MB_)))
      
      !LU decomposition of matrix 
      CALL PGGETRF( NP, NP, INV_EPS, 1, 1, DESC, IPIV, PINFO)
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("first call to PGGETRF returns " // str(PINFO) )
      ENDIF

      !obtain optimal size of working arrays
      ALLOCATE(WORK(1),IWORK(1))
      CALL PGGETRI( NP, INV_EPS, 1, 1, DESC, IPIV(1), WORK(1), -1 ,&
         IWORK(1), -1, PINFO)
   
      !allocate working arrays with optimum size 
      RWORK = WORK(1)
      LWORK=WORK(1)
      LIWORK=IWORK(1)
      IF ( LWORK == 0 .OR. LIWORK ==0 ) THEN
         CALL vtutor%error("LWORK or LIWORK estimated from "//str([RWORK])//" is zero:" // str([LWORK, LIWORK]) )
      ENDIF

      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      ALLOCATE(WORK(LWORK))
      ALLOCATE(IWORK(LIWORK))
      
      !invert now
      CALL PGGETRI( NP, INV_EPS, 1, 1, DESC, IPIV(1), WORK(1), LWORK ,&
         IWORK(1), LIWORK, PINFO)
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error( "PGGETRI returns "//str(PINFO) )
      ENDIF

      !deallocate scaLAPACK's  axuillary arrays
      IF ( ALLOCATED( IPIV  )) DEALLOCATE( IPIV  )
      IF ( ALLOCATED( WORK  )) DEALLOCATE( WORK  )
      IF ( ALLOCATED( IWORK )) DEALLOCATE( IWORK ) 

      PROFILING_STOP( 'pdgetri_pzgetri_desc') 
   END SUBROUTINE PDGETRI_PZGETRI_DESC

!=======================================================================
!
!> allocate an array for the Green function stored in the orbital
!> basis
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis
!> @param[in]  NKPTS_IRZ number of k-points
!> @param[in]  ISPIN     number of spin degrees of freedom
!> @param[in]  NOMEGA    number of frequency points (or 1)
!
!=======================================================================

  SUBROUTINE ALLOCATE_GREENS_MAT(GDES_MAT, CHAM_MAT,  NKPTS_IRZ, ISPIN, NOMEGA)
     USE ini
     USE tutor, ONLY: vtutor
     USE string, ONLY: str
     TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
     COMPLEX(q), POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis
     INTEGER             :: NKPTS_IRZ           ! in: number of k-points
     INTEGER             :: ISPIN               ! in: number of spin degrees of freedom
     INTEGER             :: NOMEGA              ! in: number of frequency points (or 1)
   ! local
     INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
     INTEGER,EXTERNAL :: NUMROC
     INTEGER          :: ISTAT 
 
 ! calculate local size of matrices
 
     CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )
 
     NP = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(MB_), MYROW, 0, NPROW)
     NQ = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(NB_), MYCOL, 0, NPCOL)
 
     ALLOCATE(CHAM_MAT(NP*NQ, NKPTS_IRZ, ISPIN, NOMEGA), STAT = ISTAT ) 
     IF ( ISTAT/=0 ) THEN
        CALL vtutor%error( "ALLOCATE_GREENS_MAT is not able to allocate "//&
          str( 16._q*NP*NQ*NKPTS_IRZ*ISPIN*NOMEGA/1024 ) //&
          " kB of data on MPI rank 0." )
     ENDIF
     CALL REGISTER_ALLOCATE(16._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
    
  END SUBROUTINE ALLOCATE_GREENS_MAT

!=======================================================================
!> deallocate an array for the Green function stored in the orbital
!> basis
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis, that is deallocated
!=======================================================================

  SUBROUTINE DEALLOCATE_GREENS_MAT(GDES_MAT, CHAM_MAT)
    USE ini
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    COMPLEX(q), POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis

    IF (ASSOCIATED(CHAM_MAT)) THEN
       CALL DEREGISTER_ALLOCATE(16._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
       DEALLOCATE(CHAM_MAT)
       NULLIFY(CHAM_MAT)
    ELSE
       WRITE(*,*) 'internal warning in DEALLOCATE_GREENS_MAT: matrix is not allocated'
    ENDIF
    
  END SUBROUTINE DEALLOCATE_GREENS_MAT

!=======================================================================
!
!> same routines as ALLOCATE_GREENS_MAT
!> in imaginary time the matrix is real valued (and symmetric) at the Gamma point
!> hence once can use GDEF instead of COMPLEX(q)
!
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis
!> @param[in]  NKPTS_IRZ number of k-points
!> @param[in]  ISPIN     number of spin degrees of freedom
!> @param[in]  NOMEGA    number of frequency points (or 1)
!=======================================================================
  SUBROUTINE ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT,  NKPTS_IRZ, ISPIN, NOMEGA)
    USE ini
    USE tutor, ONLY: vtutor
    USE string, ONLY: str 
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    GDEF, POINTER       :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis
    INTEGER             :: NKPTS_IRZ           ! in: number of k-points
    INTEGER             :: ISPIN               ! in: number of spin degrees of freedom
    INTEGER             :: NOMEGA              ! in: number of frequency points (or 1)
  ! local
    INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC
    INTEGER          :: ISTAT 
    
! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    NP = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(MB_), MYROW, 0, NPROW)
    NQ = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(NB_), MYCOL, 0, NPCOL)

    ALLOCATE(CHAM_MAT(NP*NQ, NKPTS_IRZ, ISPIN, NOMEGA), STAT = ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_GREENS_MAT_GDEF is not able to allocate "//&
         str( ndata*8._q*NP*NQ*NKPTS_IRZ*ISPIN*NOMEGA/1024 ) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CHAM_MAT=0
    CALL REGISTER_ALLOCATE(ndata*8._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
    
  END SUBROUTINE ALLOCATE_GREENS_MAT_GDEF


!=======================================================================
!> same routines as DEALLOCATE_GREENS_MAT
!> in imaginary time the matrix is real valued (and symmetric) at the Gamma point
!> hence once can use GDEF instead of COMPLEX(q)
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis, that is deallocated
!=======================================================================
  SUBROUTINE DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT)
    USE ini
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    GDEF, POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis

    IF (ASSOCIATED(CHAM_MAT)) THEN
       CALL DEREGISTER_ALLOCATE(ndata*8._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
       DEALLOCATE(CHAM_MAT)
       NULLIFY(CHAM_MAT)
    ELSE
       WRITE(*,*) 'internal warning in DEALLOCATE_GREENS_MAT: matrix is not allocated'
    ENDIF
    
  END SUBROUTINE DEALLOCATE_GREENS_MAT_GDEF



!
!> destroys distributed response descriptor
!
   SUBROUTINE DESTROY_GDES_MAT( RDES )
      TYPE(greens_mat_des), POINTER :: RDES 
 
      IF ( .NOT. ASSOCIATED( RDES ) ) RETURN

      CALL BLACS_GRIDEXIT( RDES%ICTXT ) 
      DEALLOCATE( RDES ) 
      NULLIFY( RDES )
   END SUBROUTINE DESTROY_GDES_MAT

!
!>  initalize column major distrubuted scalapack matrix descriptor
!
   SUBROUTINE RESPONSE_DES_INIT_COLMAJ( RDES, LOOPD, GDES )
      USE chi_glb
      USE minimax_struct
      USE tutor, ONLY: vtutor
      USE string, ONLY: str 
      TYPE(greens_mat_des), POINTER :: RDES ! distributed response
      TYPE(loop_des)   :: LOOPD
      TYPE(greensfdes) :: GDES
      ! local 
      INTEGER        :: NPROW0 
      INTEGER        :: NPCOL0
      INTEGER        :: MYROW0
      INTEGER        :: MYCOL0
      INTEGER        :: LLD
      INTEGER        :: INFO 
      INTEGER        :: NROWS
      INTEGER        :: NCOLS
      INTEGER, EXTERNAL :: NUMROC

      ! mark descriptor as initialized 
      ALLOCATE( RDES )

      ! set communicators
      RDES%COMM_INTRA=LOOPD%COMM_IN_GROUP
      RDES%COMM_INTER=LOOPD%COMM_BETWEEN_GROUPS
      RDES%COMM=LOOPD%COMM

      !---------------------------------------------------------------------------------------
      !initialzie 1D column-major processor grid
      ! row major ordering (but column-major ordering is identical if there is only 1 column)
      NPROW0=1
      NPCOL0=LOOPD%COMM_IN_GROUP%NCPU
      ! actual size of local matrix
      NROWS  = GDES%RES_NRPLWV_ROW_DATA_POINTS
      !number of cols ( in parallel blocking factor )
      NCOLS  = GDES%RES_NRPLWV_COL_DATA_POINTS
      ! 
      RDES%NROWS = NROWS

      CALL PROCMAP( RDES%COMM_INTRA, RDES%ICTXT, 2, NPROW0, NPCOL0)
      CALL BLACS_GRIDINFO( RDES%ICTXT, NPROW0, NPCOL0, MYROW0, MYCOL0) !get the ids for scaLAPACK
      !define corresponding array descritpor DESCA
      !since CHI_WORK is distributed blockwisely,  
      !the leading dimension of rows is simply NROWS, 
      !check this 
      LLD = MAX( 1, NUMROC(NROWS, NCOLS, MYROW0, 0, NPROW0) )
      IF ( LLD /= NROWS ) THEN 
         CALL vtutor%bug("RESPONSE_DES_INIT_COLMAJ LLD has wrong size for NODE " // str(RDES%COMM%NODE_ME) // " " // str(LLD) // " " // str(NROWS), __FILE__, __LINE__)
      ENDIF 
      CALL DESCINIT(RDES%DESC, NROWS , NROWS, NCOLS, NCOLS, 0, 0, RDES%ICTXT, &
         LLD, INFO )
      !Eventually take care for errors 
      IF ( INFO /= 0 ) THEN 
         CALL vtutor%bug("Error in 1D col-major initialization " // str(INFO), __FILE__,__LINE__)
      ENDIF

      !check if CHI_WORK has proper size
      RDES%MY_NROWS = NUMROC(RDES%DESC(M_),RDES%DESC(MB_),MYROW0,0,NPROW0)
      RDES%MY_NCOLS = NUMROC(RDES%DESC(N_),RDES%DESC(NB_),MYCOL0,0,NPCOL0) 

   END SUBROUTINE RESPONSE_DES_INIT_COLMAJ

!
!>  initalize block cyclic distrubuted scalapack matrix descriptor
!
   SUBROUTINE RESPONSE_DES_INIT_BC( RDES, LOOPD, GDES )
      USE chi_glb
      USE minimax_struct
      USE tutor, ONLY: vtutor
      USE string, ONLY: str 
      TYPE(greens_mat_des), POINTER :: RDES ! distributed response
      TYPE(loop_des)   :: LOOPD
      TYPE(greensfdes) :: GDES
      ! local 
      INTEGER        :: LLD
      INTEGER        :: INFO 
      INTEGER, EXTERNAL :: NUMROC
      INTEGER :: INU=-1
      INTEGER, PARAMETER :: ITRY = 6
      INTEGER :: I 
      ! mark descriptor as initialized 
      ALLOCATE( RDES )

      ! set communicators
      RDES%COMM_INTRA=LOOPD%COMM_IN_GROUP
      RDES%COMM_INTER=LOOPD%COMM_BETWEEN_GROUPS
      RDES%COMM=LOOPD%COMM

      !---------------------------------------------------------------------------------------
      !initialzie block-cyclic 2D processor grid
      CALL FERMAT_RAZOR( RDES%COMM_INTRA%NCPU, RDES%NPROW, RDES%NPCOL )
      IF ( INU == 1 ) THEN
         WRITE(*,'(" Using a ",I3,"x",I3," processor grid for a ",I5," x ",I3,"*",I5," matrix")')&
            RDES%NPROW, RDES%NPCOL, GDES%RES_NRPLWV_ROW_DATA_POINTS,  RDES%COMM_INTRA%NCPU, GDES%RES_NRPLWV_COL_DATA_POINTS
      ENDIF
   
      ! update blocksize and processor grid, might be to large
      IF ( INU >= 0 ) THEN
         WRITE(*,'(" distribution is done in blocks of ",I3," x ", I3)')NB, NB
      ENDIF

      ! actual size of local matrix
      RDES%NROWS = GDES%RES_NRPLWV_ROW_DATA_POINTS

      CALL PROCMAP( RDES%COMM_INTRA, RDES%ICTXT, 2, RDES%NPROW, RDES%NPCOL)
      CALL BLACS_GRIDINFO( RDES%ICTXT, RDES%NPROW, RDES%NPCOL, RDES%MYROW, RDES%MYCOL) !get the ids for scaLAPACK
       
      !define corresponding array descritpor DESCA
      !since CHI_WORK is distributed blockwisely,  
      !the leading dimension of rows is simply NROWS, 
      !check this 
      !leading dimension of local array for block cyclic distribution
      LLD = MAX( 1 , NUMROC(RDES%NROWS, NB, RDES%MYROW, 0, RDES%NPROW) )
      CALL DESCINIT(RDES%DESC, RDES%NROWS , RDES%NROWS, NB, NB, 0, 0, RDES%ICTXT, LLD, INFO )
      !Eventually take care for errors 
      IF ( INFO /= 0 ) THEN 
         CALL vtutor%error("Error in 2D block-cyclic initialization " // str(INFO))
      ENDIF

      !check if CHI_WORK has proper size
      RDES%MY_NROWS = NUMROC(RDES%DESC(M_),RDES%DESC(MB_),RDES%MYROW,0,RDES%NPROW)
      RDES%MY_NCOLS = NUMROC(RDES%DESC(N_),RDES%DESC(NB_),RDES%MYCOL,0,RDES%NPCOL) 
   END SUBROUTINE RESPONSE_DES_INIT_BC

!*********************************************************************
!
!> initialize Green's function descriptor in orbital representation
!> @param[out] GDES_MAT  pointer to generated descriptor
!> @param[in] IMAG_GRIDS time and freuency grid handle
!> @param[in] NB_TOT     maximum number of orbitals
!
!*********************************************************************

  SUBROUTINE GREENS_MAT_DES_INIT( GDES_MAT, IGRID_DES, NB_TOT)
    USE chi_glb
    USE minimax_struct
    USE tutor, ONLY: vtutor
    USE string, ONLY: str 

    TYPE (greens_mat_des), POINTER :: GDES_MAT   ! out: pointer to generated descriptor
    TYPE (loop_des)                :: IGRID_DES
    INTEGER                        :: NB_TOT     ! in: maximum number of orbitals
  ! local
    INTEGER :: NPROCS, NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC
    INTEGER :: IERROR
    INTEGER :: NB_SCALA
   
    ALLOCATE(GDES_MAT)
    
    GDES_MAT%COMM=IGRID_DES%COMM
    GDES_MAT%COMM_INTRA=IGRID_DES%COMM_IN_GROUP
    GDES_MAT%COMM_INTER=IGRID_DES%COMM_BETWEEN_GROUPS

    ! few sanity tests
    IF (.NOT. LscaLAPACK) THEN
       CALL vtutor%bug("GREENS_MAT_DES_INIT requires scaLAPACK", __FILE__, __LINE__)
    ENDIF

    IF (GDES_MAT%COMM%NCPU /= GDES_MAT%COMM_INTRA%NCPU * GDES_MAT%COMM_INTER%NCPU) THEN
       CALL vtutor%bug("core number incompatible " // str(GDES_MAT%COMM%NCPU) // " " // str(GDES_MAT%COMM_INTRA%NCPU) // " " // str(GDES_MAT%COMM_INTER%NCPU), __FILE__, __LINE__)
    ENDIF

    GDES_MAT%NROWS=NB_TOT
    GDES_MAT%ICTXT=0
    NPROCS=GDES_MAT%COMM_INTRA%NCPU
!#define distribute_columns
#ifdef distribute_columns    
    ! we chose to distribute the columns in a round robin fashion
    ! at least for the time being
    NPROW=1
    NPCOL=NPROCS
    
    CALL MPI_barrier( GDES_MAT%COMM_INTRA%MPI_COMM, IERROR )
    CALL PROCMAP( GDES_MAT%COMM_INTRA, GDES_MAT%ICTXT, 2, NPROW, NPCOL)

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    NP = NUMROC(NB_TOT,NB_TOT,MYROW,0,NPROW)   ! get number of rows on proc
    NQ = NUMROC(NB_TOT,1, MYCOL,0, NPCOL)      ! get number of cols on proc

    CALL DESCINIT(GDES_MAT%DESC, NB_TOT,NB_TOT, NB_TOT,1, 0,0, GDES_MAT%ICTXT,MAX( 1, NP ), IERROR ) ! setup descriptor
    IF (IERROR.NE.0) THEN
      CALL vtutor%bug("DESC, DESCINIT, INFO: " // str(IERROR), __FILE__, __LINE__)
    ENDIF
#else
! determine processor grid, according to naive guess
    CALL FERMAT_RAZOR( NPROCS, NPROW, NPCOL )

    CALL MPI_barrier( GDES_MAT%COMM_INTRA%MPI_COMM, IERROR )
    CALL PROCMAP( GDES_MAT%COMM_INTRA, GDES_MAT%ICTXT, 2, NPROW, NPCOL)

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    ! set blocking factor compatible to scala.F
    ! this allows to call INIT_scala and use simple functions
    ! defined in scala.F
    CALL QUERRY_SCALA_NB( NB_SCALA )

    NP = NUMROC(NB_TOT,NB_SCALA, MYROW,0,NPROW)   ! get number of rows on proc
    NQ = NUMROC(NB_TOT,NB_SCALA, MYCOL,0,NPCOL)   ! get number of cols on proc

    CALL DESCINIT(GDES_MAT%DESC, NB_TOT,NB_TOT, NB_SCALA,NB_SCALA, 0,0, GDES_MAT%ICTXT, MAX( 1, NP ), IERROR ) ! setup descriptor
    IF (IERROR.NE.0) THEN
      CALL vtutor%bug("DESC, DESCINIT, INFO: " // str(IERROR) // " " // str(NB_TOT) // " " // str(NB_SCALA) // " " // str(NP),__FILE__, __LINE__)
    ENDIF
#endif
    GDES_MAT%MY_NROWS = NP
    GDES_MAT%MY_NCOLS = NQ
  END SUBROUTINE GREENS_MAT_DES_INIT

!=======================================================================
!
!> small routine to dump a distributed matrix 
!> the descriptor in DESCA must properly describe the matrix
!> since RECON_SLICE calls check (GDEF version)
!
!> @param[in]   STRING          string to dump
!> @param[in]   CHAM_DISTRI     distributed matrix
!> @param[in]   GDES_MAT        descriptor
!> @param[in]   NB_TOT          leading dimension of distributed matrix 
!> @param[in]   IU_             unit to write to (not dump for IU<0)
!> @param[in]   IUNIT           unit to write to (not dump for IU<0)
!=======================================================================
    
  SUBROUTINE DUMP_GREENS_HAM_GDEF( STRING, CHAM_DISTRI, NB_TOT, GDES_MAT, IU_,IUNIT)
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    GDEF              :: CHAM_DISTRI(*)    ! distributed matrix
    TYPE (greens_mat_des):: GDES_MAT       ! in: descriptor
    INTEGER           :: NB_TOT
    INTEGER           :: IU_                ! unit to write to (not dump for IU<0)
    INTEGER,OPTIONAL  :: IUNIT             ! unit to write to (not dump for IU<0)
  ! local
    INTEGER :: NDUMP=16
    GDEF,ALLOCATABLE   :: CHAM(:,:)
    INTEGER IU
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW
    LOGICAL :: LGAMMA_PRINT = .FALSE.

    IF ( PRESENT( IUNIT ) ) THEN 
       IU = IUNIT
       NDUMP =  NB_TOT
    ELSE
       IU = IU_
    ENDIF
    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  GDES_MAT%DESC, COLUMN_LOW, COLUMN_HIGH)

    CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, CHAM(1,1), SIZE(CHAM)))


    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
    IF ( PRESENT( IUNIT ) ) THEN 
       DO N2 = 1, NDUMP
          WRITE(IU,10)N1+COLUMN_LOW-1, N2+COLUMN_LOW-1, CHAM(N1+COLUMN_LOW-1,N2)
       ENDDO
    ELSE
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1+COLUMN_LOW-1,N2) ,KIND=q) ,N2=1,NDUMP)
    ENDIF
    ENDDO
    WRITE(IU,*)
#ifndef gammareal
IF ( LGAMMA_PRINT ) THEN
    DO N1=1,NDUMP
    IF ( .NOT. PRESENT( IUNIT ) ) THEN 
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1+COLUMN_LOW-1,N2)),N2=1,NDUMP)
    ENDIF
    ENDDO
    WRITE(IU,*)
ENDIF
#endif
    ENDIF

    DEALLOCATE(CHAM)
1   FORMAT(1I3,3X,40F10.4)
2   FORMAT(1I3,3X,40F10.4)
10   FORMAT(2I4,2F20.10)

  END SUBROUTINE DUMP_GREENS_HAM_GDEF

!> small routine to dump a distributed matrix
!> the descriptor in DESCA must properly describe the matrix
!> since RECON_SLICE calls check 
!> @param[in]   STRING          string to dump
!> @param[in]   CHAM_DISTRI     distributed matrix
!> @param[in]   GDES_MAT        descriptor
!> @param[in]   NB_TOT          leading dimension of distributed matrix 
!> @param[in]   IU_             unit to write to (not dump for IU<0)
!> @param[in]   IUNIT           unit to write to (not dump for IU<0)
  SUBROUTINE DUMP_GREENS_HAM( STRING, CHAM_DISTRI, NB_TOT, GDES_MAT, IU_,IUNIT)
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    COMPLEX(q)        :: CHAM_DISTRI(*)    ! distributed matrix
    TYPE (greens_mat_des):: GDES_MAT       ! in: descriptor
    INTEGER           :: NB_TOT
    INTEGER           :: IU_               ! unit to write to (not dump for IU<0)
    INTEGER,OPTIONAL  :: IUNIT             ! unit to write to (not dump for IU<0)
  ! local
    INTEGER           :: IU                ! unit to write to (not dump for IU<0)
    INTEGER, PARAMETER :: NDUMP=8
    COMPLEX(q),ALLOCATABLE  :: CHAM(:,:)
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW

    IF ( PRESENT( IUNIT ) ) THEN 
       IU = IUNIT
    ELSE
       IU = IU_
    ENDIF
    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE_C(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  GDES_MAT%DESC, COLUMN_LOW, COLUMN_HIGH)

    CALLMPI( M_sum_z(GDES_MAT%COMM_INTRA, CHAM(1,1), SIZE(CHAM)))

    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1,N2+COLUMN_LOW-1) ,KIND=q) ,N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#ifndef gammareal
    DO N1=1,NDUMP
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1,N2+COLUMN_LOW-1)),N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#endif
    ENDIF

    DEALLOCATE(CHAM)
1   FORMAT(1I2,3X,40F16.8)
2   FORMAT(1I2,3X,40F16.8)

  END SUBROUTINE DUMP_GREENS_HAM

!=======================================================================
!
!> Determines the Frobenius Norm of a distributed matrix of dimension 
!> N x N and returns the result in R 
!> requires a MPI sum afterwards !!
!> @param[in]   A               distributed matrix
!> @param[in]   N               dimension of sub matrix of A for which norm is calculated
!> @param[in]   RDES2D          descriptor
!> @param[out]  R               squared Frobenius norm
!> @param[in]   NORMIN          optional argument that specifies which norm is calculated
!=======================================================================
   SUBROUTINE MATRIX_NORM_DESC_OLD( A, N,  RDES2D, R, NORMIN ) 
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF, INTENT(IN)       :: A(:)
      INTEGER, INTENT(IN)    :: N               !< dimension of submatrix
      TYPE( greens_mat_des ) :: RDES2D
      REAL(q), INTENT(INOUT) :: R
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORMIN
      INTEGER               :: I, J 
      INTEGER               :: NCOLS, NROWS 
      INTEGER, EXTERNAL     :: NUMROC 
      ! local
      ! padestrian computation of squared Frobenius norm of A: 
      !CALL BLACS_GRIDINFO( RDES2D%ICTXT, RDES2D%NPROW, RDES2D%NPCOL, RDES2D%MYROW, RDES2D%MYCOL )
      NROWS = NUMROC(N,RDES2D%DESC(MB_), RDES2D%MYROW, 0, RDES2D%NPROW )
      NCOLS = NUMROC(N,RDES2D%DESC(NB_), RDES2D%MYCOL, 0, RDES2D%NPCOL )
       
      R = 0 
      DO J = 1, NCOLS
         DO I = 1, NROWS
            R = R + A( I + (J-1)*RDES2D%MY_NROWS )*GCONJG( A( I + (J-1)*RDES2D%MY_NROWS ) )
         ENDDO
      ENDDO
      ! accumulate results from each rank
      CALLMPI( M_sum_d(RDES2D%COMM_INTRA, R , 1 ))
   END SUBROUTINE MATRIX_NORM_DESC_OLD

END MODULE


#endif


#ifdef scaLAPACK

!========================== ORTH1_DISTRI ===============================
!
!> ORTH1_DISTRI is essentially identical to #ORTH1
!>
!> It calculates a slice or stripe of the matrix
!> ~~~
!> O(I,J) =  <CPTWFP(I) | 1+ Q | CFW(J) >
!> ~~~
!> and **then distributes the slice** in a scaLAPACK compatible manner
!> using #DISTRI_SLICE.
!> #ORTH1 requires the storage of a global NB_TOT time NB_TOT matrix
!> whereas ORTH1_DISTRI only the storage of the distributed matrix.
!>
!> Arguments are identical to #ORTH1, except for additional communicator
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)

GPACC   USE moffload_struct_def

        USE prec
        USE scala
        USE tutor, ONLY: vtutor
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        PROFILING_START('orth1_distri')

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI' )

!$ACC ENTER DATA CREATE(COVL) __IF_ASYNC__
!$ACC KERNELS PRESENT(COVL) __IF_ASYNC__
        COVL=0
!$ACC END KERNELS
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        PROFILING_START('m_sum_g')
        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))
        PROFILING_STOP('m_sum_g')

#ifdef ELPA
        IF (LELPA) THEN
           IF (CSEL=='U') THEN
              CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
           ELSE
              CALL vtutor%bug("internal error in VASP: ORTH1_DISTRI is only tested for CSEL=U try " &
                 // "DISTRI_SLICE_HERM_L", __FILE__, __LINE__)
           ENDIF
        ENDIF
#endif
        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF

!$ACC EXIT DATA DELETE(COVL) __IF_ASYNC__

        PROFILING_STOP('orth1_distri')

      END SUBROUTINE ORTH1_DISTRI

!
!> Same as #ORTH1_DISTRI_HERM, this version always adds the Hermitian conjugated elements
!> to the matrix
!>
!> Costs little extra time, and is more secure, in the sense that
!> we can use the constructed matrix in places where the lower triangle is required as well
!

      SUBROUTINE ORTH1_DISTRI_HERM(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_HERM' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF
      END SUBROUTINE ORTH1_DISTRI_HERM


!========================== ORTH1_DISTRI_HERM_DESC =====================
!
!> Same as #ORTH1_DISTRI_HERM, this version requires an additional argument, the descriptor DESCA
!> of the matrix
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI_HERM_DESC(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK, DESCA)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL
        INTEGER DESCA(*)

        PROFILING_START('orth1_distri_herm_desc')

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF

        PROFILING_STOP('orth1_distri_herm_desc')

      END SUBROUTINE ORTH1_DISTRI_HERM_DESC


!========================== ORTH1_DISTRI_DAVIDSON ======================
!
!> ORTH1_DISTRI_DAVIDSON is identical to #ORTH1_DISTRI
!> but also changes the intermediate COVL matrix as required in
!> the davidson.F routine (requires vector R)
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI_DAVIDSON(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,R,MY_COMM, NBANDSK )
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK, I
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        REAL(q) :: R(NSTRIP)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_DAVIDSON' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        COVL(NPOS:NPOS+NSTRIP-1,1:NSTRIP)=0
        IF (MY_COMM%NODE_ME==MY_COMM%IONODE) THEN
           DO I=1,NSTRIP
              COVL(NPOS-1+I,I)=R(I)
           ENDDO
        ENDIF

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

#ifdef ELPA
        IF (LELPA) THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF
#endif
        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF

      END SUBROUTINE ORTH1_DISTRI_DAVIDSON

!========================== LINCOM_DISTRI ==============================
!
!> Variant of dfast::LINCOM accepting distributed matrix
!>
!> LINCOM_DISTRI does essentially the same as dfast::LINCOM but
!> accepts a distributed matrix COVL2 instead of a global
!> NB_TOT times NB_TOT matrix.
!> The required matrix COVL is reconstructed on the fly in memory
!> slice by slice.
!> The disadvantage is that the transformed wavefunctions need
!> to be stored in a temporary array.
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN,MY_COMM,NBLK)

GPACC   USE moffload_struct_def

        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q)  CF(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN)
        GDEF        COVL(*)

      ! local variables
        COMPLEX(q)  CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF_RESULT(NPROD,NIN)
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

DOESI   CF_RESULT=(0.0_q, 0.0_q); CPROF_RESULT=0; COVL_GLOBAL=0

        PROFILING_START('lincom_distri')

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI' )

!$ACC ENTER DATA CREATE(CF_RESULT,CPROF_RESULT,COVL_GLOBAL) __IF_ASYNC__
      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK

          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

          PROFILING_START('m_sum_g')
          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))
          PROFILING_STOP('m_sum_g')

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

        PROFILING_START('m_sum_g')
        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))
        PROFILING_STOP('m_sum_g')

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

!$ACC KERNELS PRESENT(CF,CPROF,CF_RESULT,CPROF_RESULT) __IF_ASYNC__
        CF   =CF_RESULT
        CPROF=CPROF_RESULT
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(CF_RESULT,CPROF_RESULT,COVL_GLOBAL) __IF_ASYNC__

        PROFILING_STOP('lincom_distri')

      END SUBROUTINE LINCOM_DISTRI

!
!> similar to #LINCOM_DISTRI but add the Hermitian conjugated elements
!> in the lower triangle
!>
!> Can be used for "full" mode (MODE="F") if the transformation
!> is Hermitian and only upper triangle is set (ZPOTRI does that).
!
      SUBROUTINE LINCOM_DISTRI_HERM(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(*)
      ! local
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

        PROFILING_START('lincom_distri_herm')

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI_HERM' )

      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

        PROFILING_STOP('lincom_distri_herm')

      END SUBROUTINE LINCOM_DISTRI_HERM


!========================== LINCOM_DISTRI_DESC =========================
!
!> LINCOM_DISTRI_DESC essentially the same as #LINCOM_DISTRI
!> but additionally requires the distributed matrix descriptor
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI_DESC(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK, DESCA)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(*)
        INTEGER DESCA(*)                               !< distributed matrix descriptor array
      ! local
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)

          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)

        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

      END SUBROUTINE LINCOM_DISTRI_DESC

!=======================================================================
!
!> Routine to call #DISTRI_SLICE without explicitly checking
!> the interface
!
!=======================================================================
      SUBROUTINE DISTRI_SLICE_NOINT(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
        USE prec
        USE scala
        IMPLICIT NONE

        INTEGER NDIM                                   !< first dimension of matrix
        INTEGER N                                      !< rank of matrix
        INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
        GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
        GDEF    A(*)                                   !< output distributed matrix
        INTEGER DESCA(*)                               !< distributed matrix descriptor array

        CALL DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

      END SUBROUTINE DISTRI_SLICE_NOINT


      SUBROUTINE BG_pDSSYEX_ZHEEVX_NOINT(COMM,A,W,N)
      USE prec
      USE scala
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               !< NxN matrix to be distributed
      GDEF    A(*)            !< input/output matrix
      REAL(q) W(N)            !< eigenvalues

      CALL BG_pDSSYEX_ZHEEVX(COMM,A(1),W,N)

      END SUBROUTINE BG_pDSSYEX_ZHEEVX_NOINT


      SUBROUTINE SCA_SHIFT_BANDS_LOW_NOINT(N, NLOW, A, SHIFT, DESCA )
      USE scala
      IMPLICIT NONE

      INTEGER N              !< actual dimension of the matrix
      INTEGER NLOW           !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA(*)       !< distributed matrix descriptor array

      CALL SCA_SHIFT_BANDS_LOW(N, NLOW, A, SHIFT, DESCA )
      END SUBROUTINE

      SUBROUTINE SCA_SHIFT_BANDS_HIGH_NOINT(N, NHIGH, A, SHIFT, DESCA )
      USE scala
      IMPLICIT NONE

      INTEGER N              !< actual dimension of the matrix
      INTEGER NHIGH          !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA(*)       !< distributed matrix descriptor array

      CALL SCA_SHIFT_BANDS_HIGH(N, NHIGH, A, SHIFT, DESCA )
      END SUBROUTINE


!=======================================================================
!
!> Dummy routines if scaLAPACK is not compiled in
!
!=======================================================================
#else
      SUBROUTINE ORTH1_DISTRI(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM, NBANDSK)
        USE prec
        USE scala

        IMPLICIT NONE

        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        TYPE(communic) :: MY_COMM
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI is not available"

      END SUBROUTINE ORTH1_DISTRI

      SUBROUTINE ORTH1_DISTRI_DAVIDSON(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,R,MY_COMM, NBANDSK)
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        REAL(q) :: R(NSTRIP)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI_DAVIDSON is not available"

      END SUBROUTINE

      SUBROUTINE ORTH1_DISTRI_DAVMAT(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2, &
           CHAM_CORRECTION,MY_COMM, NBANDSK )
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        GDEF      CHAM_CORRECTION(NSTRIP,NSTRIP)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI_DAVIDSON is not available"

      END SUBROUTINE ORTH1_DISTRI_DAVMAT

      SUBROUTINE LINCOM_DISTRI(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(1,1)

        WRITE(0,*) "internal error in scala.F:  LINCOM_DISTRI is not available"

      END SUBROUTINE LINCOM_DISTRI


      SUBROUTINE LINCOM_DISTRI_HERM(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(1,1)

        WRITE(0,*) "internal error in scala.F:  LINCOM_DISTRI is not available"

      END SUBROUTINE LINCOM_DISTRI_HERM
#endif
