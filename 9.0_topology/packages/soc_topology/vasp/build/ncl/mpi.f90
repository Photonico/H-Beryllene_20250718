# 1 "mpi.F"
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


# 2 "mpi.F" 2 




!=======================================================================
!
! in most 1 implementations the collective communcation
! is slow, hence by default we avoid them
! if you want to use them define 1 (here or in the makefile)
!
! 1 use MPI_alltoall
! avoid_async    tries to post syncronised send and read operations
!                such that collisions are avoided, this is usually slower
! PROC_GROUP     does communication is block wise in a group of
!                roughly PROC_GROUP processors
!                this reduces collisions
!
! in addition our own implementation of the collective communication
! routines allow
!=======================================================================

!#define 1
!#define avoid_async





# 40

!
! OPENMPI seems to have a bcast bug
! whenever a node initiates a bcast the master node seems
! to return immediately, whereas all the other nodes seem
! to wait until the master node initiates ANOTHER 1 call
! a simple work around is to initiate a barrier command after
! every single bcast
! to do this set MPI_bcast_with_barrier the makefile
! or undocument this line

!#define MPI_avoid_bcast
!
! alternatively the main VASP code can initiate a barrier after
! a block of bcast calls
! this might be more efficient since less barrier are initiated
! this requires to define MPI_barrier_after_bcast
! in the makefile.
! The corresponding barriers have been inserted in
! wave_high.F and wave_mpi.F

!
! In CUDA-aware mode MPI_Ireduce is not supported yet


!
! To switch to blocking bcast in M_ibcast_z/d_from




!
! To switch to blocking reduce in M_ireduce_z/d_to




!=======================================================================
!
!> 1 communication routines for VASP
!>
!> all communication should be 1._q using this interface to allow
!> adaption of other communication routines
!> routines were entirely rewritten by Kresse Georg,
!> but functionallity is similar to a module written by Peter
!> Lockey at Daresbury
!
!======================================================================
      MODULE mpimy
# 95

      USE prec
!
! Nbranch is the number of branches at each node in the gsum routines
! two should be always fine
!
      INTEGER Nbranch
      PARAMETER( Nbranch=2 )

      TYPE communic
!only COMM
        INTEGER MPI_COMM          !< MPI_Communicator
        INTEGER NODE_ME           !< node id starting from 1 ... NCPU
        INTEGER IONODE            !< node which has to do IO (set to 0 for no IO)
        INTEGER NCPU              !< total number of proc in this communicator
# 114

      END TYPE

! Standard 1 include file.
! I would like to have everything in the header but "freaking" SGI
! compiler can not handle this, thus I have to use an include file
      INCLUDE "mpif.h"

# 125

!> create operation for collective summation of quadruples
      INTEGER :: M_sum_qd_op

# 135

! There are no global local sum routines in 1, thus some workspace
! is required to store the results of the global sum
      INTEGER,PARAMETER ::  NZTMP=8000/2, NDTMP=8000, NITMP=8000, NLTMP=8000
! workspace for integer, complex, and real
      COMPLEX(q),SAVE :: ZTMP_m(NZTMP)
      REAL(q),SAVE    :: DTMP_m(NDTMP)
      INTEGER,SAVE    :: ITMP_m(NITMP)
      LOGICAL,SAVE    :: LTMP_m(NLTMP)


      EQUIVALENCE (DTMP_m,ZTMP_m)
      EQUIVALENCE (ITMP_m,ZTMP_m)
      EQUIVALENCE (LTMP_m,ZTMP_m)



      CONTAINS

!
!> Common error handling for all 1 calls.
!>
!> The macro wraps the subroutine with the same name, to make it work with the
!> preprocessor, we choose a different capitalization for the subroutine than
!> for the macro.
!

      SUBROUTINE check_mpi_error(IERROR, VASP_FUNC, MPI_FUNC, FILENAME, LINE)
         USE string, ONLY: str
         USE tutor, ONLY: vtutor_bug
         INTEGER, INTENT(IN) :: IERROR, LINE
         CHARACTER(LEN=*), INTENT(IN) :: VASP_FUNC, MPI_FUNC, FILENAME
         IF (IERROR /= MPI_success) &
            CALL vtutor_bug(VASP_FUNC // ": " // MPI_FUNC // " returns: " // str(IERROR), FILENAME, LINE)
      END SUBROUTINE check_mpi_error

!
!> define an 1 operation for emulated quad type with the same
!> effect as the intrinsically defined MPI_SUM
!
      SUBROUTINE M_sum_qd_op_function( INVEC, INOUTVEC, LEN, DATATYPE )
# 178

         REAL(qd) :: INVEC( LEN ), INOUTVEC( LEN )
         INTEGER :: LEN
         INTEGER :: DATATYPE
! local
         INTEGER :: I
# 196

         DO I = 1, LEN
            INOUTVEC( I ) = INVEC( I ) + INOUTVEC( I )
         ENDDO

      END SUBROUTINE M_sum_qd_op_function

!----------------------------------------------------------------------
!
!> initialise the basic communications (number of nodes, determine ionode)
!
!----------------------------------------------------------------------

      SUBROUTINE M_init( COMM )
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror


      CALL MPI_init( ierror )
      CALL check_mpi_error(ierror, 'M_init', 'MPI_init', "mpi.F", 219)
# 227

!
! initial communicator is world wide
! set only NCPU, NODE_ME and IONODE
! no internal setup 1._q at this point
!
      COMM%MPI_COMM= MPI_comm_world

      CALL MPI_comm_rank( COMM%MPI_COMM, COMM%NODE_ME, ierror )
      CALL check_mpi_error(ierror, 'M_init', 'MPI_comm_rank', "mpi.F", 236)
      COMM%NODE_ME= COMM%NODE_ME+1

      CALL MPI_comm_size( COMM%MPI_COMM, COMM%NCPU , ierror )
      CALL check_mpi_error(ierror, 'M_init', 'MPI_comm_size', "mpi.F", 240)

      COMM%IONODE = 1

# 256

! my own 1 summation of quadruple precision variables
      CALL MPI_op_create( M_sum_qd_op_function, .TRUE., M_sum_qd_op, ierror )
      CALL check_mpi_error(ierror, 'M_init', 'MPI_op_create', "mpi.F", 259)

      END SUBROUTINE M_init


!----------------------------------------------------------------------
!
!> creates a 2 dimensional cartesian topology and a communicator along
!> rows and columns of the process matrix
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide( COMM, NPAR, COMM_INTER, COMM_INB, reorder)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM, COMM_INTER, COMM_INB, COMM_CART
      INTEGER NPAR,NPAR_2
      INTEGER, PARAMETER :: ndims=2
      INTEGER :: dims(ndims)
      LOGICAL :: periods(ndims), reorder, remain_dims(ndims)
      INTEGER :: ierror

      IF (NPAR >= COMM%NCPU) NPAR=COMM%NCPU
      dims(1)       = NPAR
      dims(2)       = COMM%NCPU/ NPAR
      IF (dims(1)*dims(2) /= COMM%NCPU ) THEN
         WRITE(0,*) 'M_divide: can not subdivide ',COMM%NCPU,'nodes by',NPAR
      ENDIF

      periods(ndims)=.FALSE.

      CALL MPI_cart_create( COMM%MPI_COMM , ndims, dims, periods, reorder, &
                COMM_CART%MPI_COMM , ierror)
      CALL check_mpi_error(ierror, 'M_divide', 'MPI_cart_create', "mpi.F", 293)
! create the in-band communicator
      remain_dims(1)= .FALSE.
      remain_dims(2)= .TRUE.

      CALL MPI_cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INB%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_divide', 'MPI_cart_sub (1)', "mpi.F", 299)

! create the inter-band communicator
      remain_dims(1)= .TRUE.
      remain_dims(2)= .FALSE.

      CALL MPI_cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INTER%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_divide', 'MPI_cart_sub (2)', "mpi.F", 306)
! overwrite initial communicator by new (1._q,0._q)
      COMM=COMM_CART

      END SUBROUTINE M_divide


!----------------------------------------------------------------------
!
!> subdivides the communicator into two "images" (groups)
!>
!> the first (1._q,0._q) includes cores 1...NCORE, the second (1._q,0._q) all the remaining cores
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide2( COMM, NCORE, COMM_INTER, COMM_INB)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM, COMM_INTER, COMM_INB
      INTEGER NCORE
      INTEGER :: color, key
      INTEGER, PARAMETER :: ndims=2
      LOGICAL :: remain_dims(ndims)
      INTEGER :: i,ierror

      IF (NCORE >= COMM%NCPU) NCORE=COMM%NCPU
! first create the communicator COMM_INB which communicates inside (1._q,0._q) image
      IF (COMM%NODE_ME<=NCORE) THEN
         color=1
         key=COMM%NODE_ME
      ELSE
         color=2
         key=COMM%NODE_ME-NCORE
      ENDIF

      CALL MPI_comm_split(COMM%MPI_COMM, color, key, COMM_INB%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_divide2', 'MPI_comm_split (1)', "mpi.F", 343)

      CALL M_initc( COMM_INB)
!
! create the  communicator COMM_INTER
! this is only defined for the first node in each image
      color=MPI_UNDEFINED
      key=COMM%NODE_ME
      IF (COMM_INB%NODE_ME==1) THEN
         color=1
      ENDIF

      CALL MPI_comm_split(COMM%MPI_COMM, color, key, COMM_INTER%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_divide2', 'MPI_comm_split (2)', "mpi.F", 356)

      IF (COMM_INTER%MPI_COMM /= MPI_COMM_NULL) THEN
         CALL M_initc( COMM_INTER)
      ELSE
         COMM_INTER%NCPU=0
         COMM_INTER%IONODE=1
! for VASP internal reasons it is important that
! COMM%NODE_ME is equivalent to the "image"
         IF (COMM%NODE_ME<=NCORE) THEN
            COMM_INTER%NODE_ME=1
         ELSE
            COMM_INTER%NODE_ME=2
         ENDIF
      ENDIF

      END SUBROUTINE M_divide2


!----------------------------------------------------------------------
!> splits the communicator COMM into NTAUPAR groups
!>
!> splitting is cartesian and INTER is a proper communicator
!> that allows for collective communication among groups
!----------------------------------------------------------------------

      SUBROUTINE M_divide_general( COMM, NTAUPAR, COMM_INTER, COMM_INTRA, MAP )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      INTEGER NTAUPAR
      INTEGER MAP(NTAUPAR)  !< distribution matrix
      TYPE(communic) COMM, COMM_INTRA, COMM_INTER, COMM_CART
      INTEGER :: NDONE , NREMAIN, NGROUPS
      INTEGER :: ID, I
      INTEGER :: ierr

!----------------------------------------------------------------------------------------
!set up node IDs between groups (colums of processor grid)
!note that the # of tau groups, i.e. COMM_INTER%NCPU is given
!by NTAUPAR.
      COMM_INTER%NCPU = NTAUPAR
      COMM_INTER%NODE_ME = 0

      IF ( MOD( COMM%NCPU , NTAUPAR ) == 0 ) THEN
!if NTAUPAR is a divisor of NCPU, there are COMM%NCPU / NTAUPAR cores in each group
!hence the gropu id is basically NODE_ME/ (COMM%NCPU / NTAUPAR) + 1

!care must be taken if COMM%NCPU / NTAUPAR  is a divisor of the ID.
!in this case the group id is simply COMM%NODE_ME / ( COMM%NCPU / NTAUPAR )
         IF ( MOD( COMM%NODE_ME , COMM%NCPU / NTAUPAR ) == 0 ) THEN
            COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR )
         ELSE
            COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR )  + 1
         ENDIF

      ELSE
         IF ( COMM%NODE_ME <=  MOD( COMM%NCPU , NTAUPAR )*(COMM%NCPU / NTAUPAR + 1 )) THEN
!in the case when NTAUPAR is not a divisor of NCPU the distribution is such that
!the first MOD( NCPU , NTAUPAR ) groups have NCPU / NTAUPAR + 1 cores
!so the first MOD( NCPU , NTAUPAR )*( NCPU / NTAUPAR + 1 ) cores have
!following group IDs
            IF ( MOD( COMM%NODE_ME , COMM%NCPU / NTAUPAR + 1 ) == 0 ) THEN
               COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR + 1 )
            ELSE
               COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR + 1 )  + 1
            ENDIF

         ELSE
!the remaining NCPU - MOD( NCPU , NTAUPAR )*( NCPU / NTAUPAR + 1 ) cores
!are distributed in groups of NCPU / NTAUPAR cores per group
!so the groups MOD( NCPU , NTAUPAR ) + 1 ... NTAUPAR contain
!NCPU / NTAUPAR cores each

!NDONE CPUs 1._q
            NDONE = MOD( COMM%NCPU , NTAUPAR )*( COMM%NCPU / NTAUPAR + 1 )

!NREMAIN remaining CPUs
            NREMAIN = COMM%NCPU - NDONE

!they need to be grouped into NGROUPS
            NGROUPS = NTAUPAR - MOD( COMM%NCPU , NTAUPAR )

!where each of the NGROUPS groups contain NREMAIN / NGROUPS CPUs
!check if any cpu is left out
            IF ( MOD( NREMAIN , NGROUPS ) /= 0 ) THEN
               CALL vtutor%bug("M_divide_general: some nodes are missing " // str(NREMAIN) // &
                  " " // str(NGROUPS) // " " // str(NREMAIN / NGROUPS), "mpi.F", 444)
            ENDIF

            IF ( .NOT.( (1 <= COMM%NODE_ME - NDONE).AND.(COMM%NODE_ME - NDONE <=  NREMAIN) ) ) THEN
               CALL vtutor%error("M_divide_general: NODE " // str(COMM%NODE_ME) // "has no group id, choose another NTAUPAR")
            ENDIF

!distribute them according
             IF ( MOD( COMM%NODE_ME - NDONE , NREMAIN / NGROUPS ) == 0 ) THEN
                COMM_INTER%NODE_ME = ( COMM%NODE_ME - NDONE ) / ( NREMAIN / NGROUPS ) + &
                   MOD( COMM%NCPU , NTAUPAR )
             ELSE
                COMM_INTER%NODE_ME = ( COMM%NODE_ME - NDONE ) / ( NREMAIN / NGROUPS ) + 1 +&
                   MOD( COMM%NCPU , NTAUPAR )
             ENDIF
         ENDIF
      ENDIF

!----------------------------------------------------------------------------------------
!set up node IDs in groups (rows of processor grid)

!the distribution matrix MAP(I)
!gives us the number of nodes in the Ith tau group
!as determined above there are COMM_INTER%NCPU tau groups
      COMM_INTRA%NCPU = MAP(COMM_INTER%NODE_ME)

!the IDs inside each tau group can be set elegantly
!using the distribution matrix
      IF ( COMM_INTER%NODE_ME == 1 ) THEN
!in the first group the global ids conincide with the intau-ids
         COMM_INTRA%NODE_ME = COMM%NODE_ME
      ELSE
!in the Ith group the total number of CPU from group 1 to I-1
!are subtracted from the global node
         COMM_INTRA%NODE_ME = COMM%NODE_ME - SUM( MAP( 1 : COMM_INTER%NODE_ME - 1 ) )
      ENDIF

!----------------------------------------------------------------------------------------
!divide the communicator into columns, i.e. each tau group shares a separate communicator
      CALL MPI_comm_split( COMM%MPI_COMM , COMM_INTER%NODE_ME , COMM_INTRA%NODE_ME , &
         COMM_INTRA%MPI_COMM, ierr )
!check if everything went ok
      IF ( ierr /= 0 ) THEN
         CALL vtutor%bug("M_divide_general: initialization of communicator in tau groups " &
            // "failed " // str(ierr), "mpi.F", 488)
      ENDIF

!----------------------------------------------------------------------------------------
!divide the communicator into rows, i.e. between tau groups
!this means all nodes with the same COMM_INTRA%NODE_ME are in the new
!communicator COMM_BEWWEENTAU%MPI_COMM
      CALL MPI_comm_split( COMM%MPI_COMM , COMM_INTRA%NODE_ME , COMM_INTER%NODE_ME , &
         COMM_INTER%MPI_COMM, ierr )
!check if everything went ok
      IF ( ierr /= 0 ) THEN
         CALL vtutor%bug("M_divide_general: initialization of communicator between tau groups " &
            // "failed " // str(ierr), "mpi.F", 500)
      ENDIF

!initialize the communicators, just to be sure
!(although this should have happen already in M_DIVIDE_GENERAL)
      CALL M_initc( COMM_INTER )
      CALL M_initc( COMM_INTRA )

      END SUBROUTINE M_divide_general


!----------------------------------------------------------------------
!> splits communicator into NGROUPS groups
!>
!> if splitting is not cartesian so COMM_INTER is not a proper communicator
!> collective communication among groups using COMM_INTER will crash!
!> cartesian splitting satisfies MOD( COMM%NCPU , NGROUPS ) == 0
!----------------------------------------------------------------------

      SUBROUTINE M_divide_into_groups( COMM, COMM_INTER, COMM_INTRA, NGROUPS )
      IMPLICIT NONE
      TYPE(communic) :: COMM
      TYPE(communic) :: COMM_INTER
      TYPE(communic) :: COMM_INTRA
      INTEGER        :: NGROUPS
! local
      INTEGER :: ierror
      INTEGER :: MAP(COMM%NCPU,2)
      INTEGER :: MAP_(COMM%NCPU)
      INTEGER :: LMAP(NGROUPS)
      INTEGER :: ICOLOR, I, NCPU

      COMM_INTER%NCPU = 0
      COMM_INTRA%NCPU = 0

! map contains layout how COMM is subdivided into groups
! because we allow for non-cartesian splitting,
! there is no proper inter-communicator COMM_INTER
! which would allow for collective inter-group communications
! since (1._q,0._q) or more groups may contain more ranks than the other

! how may ranks are in each group?
      LMAP = 0
      IF ( MOD( COMM%NCPU , NGROUPS ) == 0 ) THEN
         LMAP(1:NGROUPS) = COMM%NCPU/NGROUPS
      ELSE
         LMAP( 1:MOD(COMM%NCPU,NGROUPS) ) = COMM%NCPU/NGROUPS +1
         LMAP( MOD(COMM%NCPU,NGROUPS)+1:NGROUPS ) = COMM%NCPU/NGROUPS
      ENDIF

! MAP stores the mapping of global communicator into groups
      MAP = 0

! first row of MAP tells us the new group id of each COMM%NODE_ME
      MAP_ = 0
      IF ( LMAP(1)+1 > COMM%NODE_ME ) MAP_( COMM%NODE_ME ) = 1
      DO ICOLOR = 1, NGROUPS-1
         IF ( SUM( LMAP( 1:ICOLOR ) ) < COMM%NODE_ME .AND.&
              SUM( LMAP( 1:ICOLOR+1 ))+1 > COMM%NODE_ME ) THEN
            MAP_( COMM%NODE_ME ) = ICOLOR+1
         ENDIF
      ENDDO
      CALL M_sum_i( COMM, MAP_, SIZE( MAP_ ) )

! second row of MAP tells us the new local id of each COMM%NODE_ME
      MAP( COMM%NODE_ME, 2 ) = SUM( LMAP(1:MAP_(COMM%NODE_ME)) )-COMM%NODE_ME+1
      CALL M_sum_i( COMM, MAP, SIZE( MAP ) )
      MAP( :, 1 ) = MAP_(:)

! initialize colors
      COMM_INTER%NODE_ME = MAP( COMM%NODE_ME, 1)
! initialize node ids in each color
      COMM_INTRA%NODE_ME = MAP( COMM%NODE_ME, 2)

      CALL MPI_comm_split( COMM%MPI_COMM, COMM_INTER%NODE_ME, &
         COMM_INTRA%NODE_ME, COMM_INTRA%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_divide_into_groups', 'MPI_comm_split (1)', "mpi.F", 576)

      CALL M_initc( COMM_INTRA )
      COMM_INTRA%NCPU = LMAP( COMM_INTER%NODE_ME )

! in case splitting is cartesian initialize inter communicator as well
      IF ( MOD( COMM%NCPU, NGROUPS ) == 0 ) THEN
         ICOLOR=COMM_INTRA%NODE_ME
      ELSE
! this is only defined for the first nodes in each group
         ICOLOR=MPI_UNDEFINED
         NCPU=LMAP( 1 )
         DO I = 1, NGROUPS-1
            NCPU = MIN( NCPU, LMAP( I+1 ) )
         ENDDO
         IF ( COMM_INTRA%NODE_ME <= NCPU ) THEN
            ICOLOR=COMM_INTRA%NODE_ME
         ENDIF
      ENDIF
      CALL MPI_comm_split( COMM%MPI_COMM, ICOLOR, &
         COMM_INTER%NODE_ME, COMM_INTER%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_divide_into_groups', 'MPI_comm_split (2)', "mpi.F", 597)

      IF (COMM_INTER%MPI_COMM /= MPI_COMM_NULL) THEN
         CALL M_initc( COMM_INTER )
      ELSE
         COMM_INTER%NCPU=0
      ENDIF

      END SUBROUTINE M_divide_into_groups


!----------------------------------------------------------------------
!
! M_divide_shmem:
!
!----------------------------------------------------------------------
# 710


!----------------------------------------------------------------------
!
! M_divide_intra_inter_node:
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide_intra_inter_node(COMM_WORLD,COMM,COMM_intra,COMM_inter)
      USE c2f_interface
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"
      TYPE (communic) COMM_WORLD,COMM,COMM_intra,COMM_inter,COMM_test
! local variables
      CHARACTER*(MPI_MAX_PROCESSOR_NAME) myname,pname
      INTEGER :: COMM_group,COMM_intra_group,COMM_inter_group
      INTEGER :: group(0:COMM%NCPU-1)
      INTEGER :: resultlen,myid
      INTEGER :: ierror
      INTEGER :: I,IGRP
! shmem variables for sanity check
      INTEGER(c_int)    :: shmid
      TYPE(c_ptr)       :: address
      INTEGER(c_size_t) :: k

      INTEGER*4, POINTER :: IDSHMEM(:)
      INTEGER :: ID

! attempt automatic division of COMM

# 760

      CALL MPI_comm_split_type(COMM%MPI_COMM,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,COMM_intra%MPI_COMM,ierror)
      CALL MPI_comm_group(COMM_intra%MPI_COMM,COMM_intra_group,ierror)
      CALL MPI_comm_group(COMM%MPI_COMM,COMM_group,ierror)

      CALL M_initc(COMM_intra)

# 791


      IGRP=0
      CALL MPI_group_rank(COMM_intra_group,myid,ierror)
      DO I=1,COMM%NCPU
         IF (I==COMM%NODE_ME) ID=myid
         CALL MPI_bcast(ID,1,MPI_INTEGER,I-1,COMM%MPI_COMM,ierror)
         CALL check_mpi_error(ierror, 'M_divide_intra_inter_node', 'MPI_bcast (2)', "mpi.F", 798)

         CALL MPI_barrier(COMM%MPI_COMM,ierror)

         IF (ID==myid) THEN
            group(IGRP)=I-1
            IGRP=IGRP+1
         ENDIF
      ENDDO

      CALL MPI_group_incl(COMM_group,IGRP,group,COMM_inter_group,ierror)
      CALL MPI_comm_create(COMM%MPI_COMM,COMM_inter_group,COMM_inter%MPI_COMM,ierror)

      CALL M_initc(COMM_inter)

!      DO I=1,COMM_WORLD%NCPU
!         IF (COMM_WORLD%NODE_ME==I) THEN
!            WRITE(*,'(A,I4,X,A,I4,X,A,I4)') 'global id:',comm_world%node_me,'intra node id:',comm_intra%node_me,'inter node id:',comm_inter%node_me
!         ENDIF
!         CALL MPI_barrier(COMM_WORLD%MPI_COMM,ierror)
!      ENDDO

      END SUBROUTINE M_divide_intra_inter_node


!----------------------------------------------------------------------
!
!> initialise a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_initc( COMM)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror
      INTEGER id_in_group

      CALL MPI_comm_rank( COMM%MPI_COMM, id_in_group, ierror )
      CALL check_mpi_error(ierror, 'M_initc', 'MPI_comm_rank', "mpi.F", 838)

      CALL MPI_comm_size( COMM%MPI_COMM, COMM%NCPU, ierror )
      CALL check_mpi_error(ierror, 'M_initc', 'MPI_comm_size', "mpi.F", 841)

      COMM%NODE_ME = id_in_group + 1

      CALL MPI_barrier( COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_initc', 'MPI_barrier', "mpi.F", 846)
      COMM%IONODE =1

      END SUBROUTINE M_initc


!----------------------------------------------------------------------
!
!> frees a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_freec( COMM )
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror
      INTEGER id_in_group

      CALL MPI_comm_free( COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_freec', 'MPI_comm_free', "mpi.F", 867)

      COMM%NODE_ME = 0
      COMM%IONODE = 0
      COMM%NCPU = 0

      END SUBROUTINE M_freec


!----------------------------------------------------------------------
!
!> initialise a NCCL_COMM communicator
!
!----------------------------------------------------------------------
# 917


      END MODULE mpimy


!======================================================================
!
! all other routines are often called with either
! real or complex arrays, vectors or scalars, so I can not put
! them into the F90 module
!
!======================================================================
!----------------------------------------------------------------------
!
!> exit 1, very simple just exit 1
!
!----------------------------------------------------------------------

      SUBROUTINE M_exit()
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      INTEGER ierror

!      CALL MPI_barrier(MPI_comm_world, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('M_exit: MPI_barrier returns: ',ierror)

! freq quadruple summation operator
      CALL MPI_op_free( M_sum_qd_op, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('M_exit: MPI_op_free returns: ',ierror)

      CALL MPI_finalize( ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('M_exit: MPI_finalize returns: ',ierror)
      STOP

      END SUBROUTINE M_exit


!----------------------------------------------------------------------
!
!> exits 1 and program because of error, a message is printed
!
!----------------------------------------------------------------------

      SUBROUTINE M_stop(message)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER (LEN=*) message
      INTEGER request, status, ierror
      LOGICAL flag
      INTEGER i

      CALL MPI_ibarrier(MPI_comm_world, request, ierror )

      WRITE (*,*) message

      DO i=1,30
         CALL MPI_Test(request, flag, status, ierror )
         IF (flag) CALL M_exit()
         CALL SLEEP(1)
      ENDDO

      CALL MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      END SUBROUTINE M_stop


      SUBROUTINE M_stop_ierr(message, ierror)
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER (LEN=*) message
      INTEGER ierror

      WRITE (*,*) message, ierror

      CALL MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      END SUBROUTINE M_stop_ierr


!======================================================================
!
! Send and Receive routines, map directly onto 1
!
!======================================================================

!----------------------------------------------------------------------
!
!> send n integers stored in ivec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_i (COMM, node, ivec, n)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)

      INTEGER status(MPI_status_size), ierror

      CALL MPI_send( ivec(1), n, MPI_integer, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_send_i', 'MPI_send', "mpi.F", 1033)

      END SUBROUTINE M_send_i


!----------------------------------------------------------------------
!
!> receive n integers into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_i(COMM, node, ivec, n )
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)
      INTEGER status(MPI_status_size), ierror

      CALL MPI_recv( ivec(1), n, MPI_integer , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      CALL check_mpi_error(ierror, 'M_recv_i', 'MPI_recv', "mpi.F", 1057)

      END SUBROUTINE M_recv_i


!----------------------------------------------------------------------
!
!> send n double complex stored in zvec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_z (COMM, node, zvec, n)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      COMPLEX(q) :: zvec(n)

      INTEGER status(MPI_status_size), ierror
# 1083


# 1105


# 1116


      CALL MPI_send( zvec(1), n, MPI_double_complex, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_send_z', 'MPI_send', "mpi.F", 1120)

      END SUBROUTINE M_send_z


!----------------------------------------------------------------------
!
!> receive n double complex into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_z(COMM, node, zvec, n )

!! USE moffload_struct_def

      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      COMPLEX(q) :: zvec(n)
      INTEGER status(MPI_status_size), ierror
# 1146


# 1168


# 1179


      CALL MPI_recv( zvec(1), n, MPI_double_complex , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      CALL check_mpi_error(ierror, 'M_recv_z', 'MPI_recv', "mpi.F", 1183)

      END SUBROUTINE M_recv_z


!----------------------------------------------------------------------
!
!> send n double stored in ivec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_d (COMM, node, dvec, n)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      REAL(q) :: dvec(n)

      INTEGER status(MPI_status_size), ierror
# 1209


# 1231


# 1242


      CALL MPI_send( dvec(1), n, MPI_double_precision, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_send_d', 'MPI_send', "mpi.F", 1246)

      END SUBROUTINE M_send_d


!----------------------------------------------------------------------
!
!> receive n double  into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_d(COMM, node, dvec, n )

!! USE moffload_struct_def

      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      REAL(q) :: dvec(n)
      INTEGER status(MPI_status_size), ierror
# 1272


# 1294


# 1305


      CALL MPI_recv( dvec(1), n, MPI_double_precision , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      CALL check_mpi_error(ierror, 'M_recv_d', 'MPI_recv', "mpi.F", 1309)

      END SUBROUTINE M_recv_d


!======================================================================
!
! global sum and maximum routines
!
!======================================================================

# 1955


! split the arrays and copy for MPI_allreduce

!----------------------------------------------------------------------
!
!> performs a global and on n logicals in vector lvec
!
!----------------------------------------------------------------------

      SUBROUTINE M_and(COMM, lvec, n)
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      LOGICAL lvec(n)
      INTEGER j,k

      INTEGER ierror, ichunk

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_and: invalid vector size n " // str(n), "mpi.F", 1991)
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NLTMP
         ichunk = MIN( n-j+1 , NLTMP)

         CALL MPI_allreduce(lvec(j), LTMP_m(1), ichunk, MPI_logical, &
     &                      MPI_land, COMM%MPI_COMM, ierror)
         CALL check_mpi_error(ierror, 'M_and', 'MPI_allreduce', "mpi.F", 2002)

         DO k = 0, ichunk-1
            lvec(j+k) = LTMP_m(k+1)
         ENDDO

      ENDDO

      END SUBROUTINE M_and


!----------------------------------------------------------------------
!
!> performs a global sum on n integers in vector ivec
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_i(COMM, ivec, n )
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sum_i: invalid vector size n " // str(n), "mpi.F", 2045)
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
     &                       MPI_sum, COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_sum_i', 'MPI_allreduce', "mpi.F", 2056)

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO

      ENDDO

      END SUBROUTINE M_sum_i


!----------------------------------------------------------------------
!
!> performs a global min on n integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_i(COMM, ivec, n )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk
# 2091


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_i: invalid vector size n " // str(n), "mpi.F", 2100)
      END IF

      

# 2138


# 2161


!  there is no inplace global min in 1, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
                             MPI_min, COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_min_i', 'MPI_allreduce', "mpi.F", 2171)

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO
      ENDDO

      

      END SUBROUTINE M_min_i


!----------------------------------------------------------------------
!
! M_max_i: performs a global max on n integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_i(COMM, ivec, n )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk
# 2207


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_i: invalid vector size n " // str(n), "mpi.F", 2216)
      END IF

      

# 2254


# 2277


!  there is no inplace global max in 1, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
                             MPI_max, COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_max_i', 'MPI_allreduce', "mpi.F", 2287)

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO
      ENDDO

      

      END SUBROUTINE M_max_i


!----------------------------------------------------------------------
!
!> performs a global max search on n doubles in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_d(COMM, vec, n )
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_d: invalid vector size n " // str(n), "mpi.F", 2326)
      END IF

!  there is no inplace global max in 1, thus we have to use
!  a work array
      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_max, &
                             COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_max_d', 'MPI_allreduce', "mpi.F", 2337)

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_max_d


!----------------------------------------------------------------------
!
!> performs a global min search on n doubles in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_d(COMM, vec, n )
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_d: invalid vector size n " // str(n), "mpi.F", 2372)
      END IF

!  there is no inplace global max in 1, thus we have to use
!  a work array
      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_min, &
                             COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_min_d', 'MPI_allreduce', "mpi.F", 2383)

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_min_d


!----------------------------------------------------------------------
!
!> performs a global sum on n doubles in vector vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!>  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_d(COMM, vec, n )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_d: invalid vector size n " // str(n), "mpi.F", 2424)
      END IF

      

# 2447


# 2466


!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_sumb_d', 'MPI_allreduce', "mpi.F", 2477)

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      

      END SUBROUTINE M_sumb_d


!----------------------------------------------------------------------
!
!> performs a global sum on n singles in vector vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!>  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_s(COMM, vec, n )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_s: invalid vector size n " // str(n), "mpi.F", 2520)
      END IF

      

# 2543


# 2561


!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_REAL, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_sumb_s', 'MPI_allreduce', "mpi.F", 2572)

         CALL SCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      

      END SUBROUTINE M_sumb_s


!----------------------------------------------------------------------
!
!> performs a global sum on n complex items in vector vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!>  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_z(COMM, vec, n )
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER j

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_z: invalid vector size n " // str(n), "mpi.F", 2612)
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array
      DO j = 1, n, NZTMP
         ichunk = MIN( n-j+1 , NZTMP)

         CALL MPI_allreduce( vec(j), ZTMP_m(1), ichunk, &
                             MPI_double_complex, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_sumb_z', 'MPI_allreduce', "mpi.F", 2623)

         CALL ZCOPY(ichunk , ZTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_sumb_z


!----------------------------------------------------------------------
!
!> performs a global product on n complex numbers in vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!
!----------------------------------------------------------------------

      SUBROUTINE M_prodb_z(COMM, vec, n )
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_prodb_z: invalid vector size n " // str(n), "mpi.F", 2660)
      END IF

!  there is no inplace global sum in 1, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), ZTMP_m(1), ichunk, &
                             MPI_double_complex, MPI_prod, &
                             COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_prodb_z', 'MPI_allreduce', "mpi.F", 2672)

         CALL ZCOPY(ichunk , ZTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_prodb_z



!----------------------------------------------------------------------
!
!> to make live easier, a global sum for 2 scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_2(COMM, v1, v2)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      REAL(q) vec(2),v1,v2

      vec(1)=v1
      vec(2)=v2

      CALL M_sumb_d(COMM, vec, 2)

      v1=vec(1)
      v2=vec(2)

      END SUBROUTINE M_sum_2


!----------------------------------------------------------------------
!
!> a global sum for 3 scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_3(COMM, v1, v2, v3)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      REAL(q) vec(3),v1,v2,v3

      vec(1)=v1
      vec(2)=v2
      vec(3)=v3

      CALL M_sumb_d(COMM, vec, 3)

      v1=vec(1)
      v2=vec(2)
      v3=vec(3)

      END SUBROUTINE M_sum_3


!----------------------------------------------------------------------
!
!> a global sum for 4 scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_4(COMM, v1, v2, v3, v4)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      REAL(q) vec(4),v1,v2,v3,v4

      vec(1)=v1
      vec(2)=v2
      vec(3)=v3
      vec(4)=v4

      CALL M_sumb_d(COMM, vec, 4)

      v1=vec(1)
      v2=vec(2)
      v3=vec(3)
      v4=vec(4)

      END SUBROUTINE M_sum_4


!======================================================================
!
!> Global barrier routine
!
!======================================================================

      SUBROUTINE M_barrier(COMM )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      CALL MPI_barrier( COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_barrier', 'MPI_barrier', "mpi.F", 2777)

      END SUBROUTINE M_barrier

!======================================================================
!
! Global Copy Routines
!
!======================================================================

!----------------------------------------------------------------------
!
!> copy n integers from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i(COMM, vec, n )
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER vec(n)

      INTEGER ierror
# 2807

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_i: invalid vector size n " // str(n), "mpi.F", 2815)
      END IF
# 2823

      CALL MPI_bcast( vec(1), n, MPI_integer, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_i', 'MPI_bcast', "mpi.F", 2826)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_i


!----------------------------------------------------------------------
!
!> copy n logical from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_l(COMM, vec, n )
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      LOGICAL vec(n)

      INTEGER ierror
# 2855

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_l: invalid vector size n " // str(n), "mpi.F", 2863)
      END IF
# 2871

      CALL MPI_bcast( vec(1), n, MPI_logical, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_l', 'MPI_bcast', "mpi.F", 2874)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_l


!----------------------------------------------------------------------
!
!> copy n integers from node inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i_from(COMM, vec, n , inode)
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER inode
      INTEGER vec(n)

      INTEGER ierror
# 2904

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_i_from: invalid vector size n " // str(n), "mpi.F", 2912)
      END IF
# 2920

      CALL MPI_bcast( vec(1), n, MPI_integer, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_i_from', 'MPI_bcast', "mpi.F", 2923)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_i_from


!----------------------------------------------------------------------
!
!> copy n double precision from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d(COMM, vec, n )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER ierror
# 2955

# 2958


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_d: invalid vector size n " // str(n), "mpi.F", 2967)
      END IF

      

# 2998


# 3013


# 3021

      CALL MPI_bcast( vec(1), n,  MPI_double_precision, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_d', 'MPI_bcast', "mpi.F", 3024)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )



      

      END SUBROUTINE M_bcast_d


!----------------------------------------------------------------------
!
!> copy n single precision from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_s(COMM, vec, n )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)

      INTEGER ierror
# 3059

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_s: invalid vector size n " // str(n), "mpi.F", 3067)
      END IF

      

# 3086


# 3101


# 3109

      CALL MPI_bcast( vec(1), n,  MPI_real, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_s', 'MPI_bcast', "mpi.F", 3112)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      

      END SUBROUTINE M_bcast_s


!----------------------------------------------------------------------
!
!> copy n double precision complex from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z(COMM, vec, n )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER ierror
# 3146

# 3149


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_z: invalid vector size n " // str(n), "mpi.F", 3158)
      END IF

      

# 3190


# 3205


# 3213

      CALL MPI_bcast( vec(1), n,  MPI_double_complex, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_z', 'MPI_bcast', "mpi.F", 3216)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      

      END SUBROUTINE M_bcast_z


!----------------------------------------------------------------------
!
!> copy n single precision complex from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_c_from(COMM, vec, n, inode )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(qs) vec(n)

      INTEGER inode, ierror
# 3250

# 3253

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_z_from: invalid vector size n " // str(n), "mpi.F", 3261)
      END IF

# 3286


# 3300


# 3308

      CALL MPI_bcast( vec(1), n,  MPI_complex, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_z_from', 'MPI_bcast', "mpi.F", 3311)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_c_from


!----------------------------------------------------------------------
!
!> copy n double precision complex from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z_from(COMM, vec, n, inode )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER inode, ierror
# 3343

# 3346

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_z_from: invalid vector size n " // str(n), "mpi.F", 3354)
      END IF

# 3379


# 3393


# 3401

      CALL MPI_bcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_z_from', 'MPI_bcast', "mpi.F", 3404)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_z_from


!----------------------------------------------------------------------
!
!> copy n double precision complex from inode to all nodes the number
!> of items is INTEGER (KIND=8)
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z8_from(COMM, vec, n, inode )

!! USE moffload_struct_def

      USE mpimy
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER (KIND=selected_int_kind(15)) :: n
      COMPLEX(q) vec(n)
! local
      INTEGER (KIND=selected_int_kind(15)) :: npos, nstep=2**24 ! nmax about 268 Mbyte bcast
      INTEGER :: nsend

      INTEGER inode, ierror
# 3439

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_z8_from: invalid vector size", "mpi.F", 3447)
      END IF

# 3468


# 3487


# 3500

      DO npos=1,n,nstep

         nsend=MIN(n-npos+1,nstep)

         CALL MPI_bcast( vec(npos), nsend,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
              &        ierror )
         CALL check_mpi_error(ierror, 'M_bcast_z8_from', 'MPI_bcast', "mpi.F", 3507)
      ENDDO

      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_z8_from


!----------------------------------------------------------------------
!
!> copy n single precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_s_from(COMM, vec, n, inode )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)

      INTEGER inode, ierror
# 3539

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_s_from: invalid vector size n " // str(n), "mpi.F", 3547)
      END IF

# 3562


# 3575


# 3583

      CALL MPI_bcast( vec(1), n,  MPI_real, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_s_from', 'MPI_bcast', "mpi.F", 3586)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_s_from


!----------------------------------------------------------------------
!
!> copy n double precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d_from(COMM, vec, n, inode )

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER inode, ierror
# 3618

# 3621

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_d_from: invalid vector size n " // str(n), "mpi.F", 3629)
      END IF

# 3654


# 3667


# 3675

      CALL MPI_bcast( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CALL check_mpi_error(ierror, 'M_bcast_d_from', 'MPI_bcast', "mpi.F", 3678)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_d_from


!-----------------------------------------------------------------------
!
!> copy n double precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d8_from(COMM, vec, n, inode )

!! USE moffload_struct_def

      USE mpimy
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER (KIND=selected_int_kind(15)) :: n
      REAL(q) vec(n)
! local
      INTEGER (KIND=selected_int_kind(15)) :: npos, nstep=2**24 ! nmax about 268 Mbyte bcast
      INTEGER :: nsend

      INTEGER inode, ierror
# 3712

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_d8_from: invalid vector size", "mpi.F", 3720)
      END IF

# 3741


# 3760


# 3773

      DO npos=1,n,nstep

         nsend=MIN(n-npos+1,nstep)

         CALL MPI_bcast( vec(npos), nsend,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
              &                ierror )
         CALL check_mpi_error(ierror, 'M_bcast_d8_from', 'MPI_bcast', "mpi.F", 3780)
      ENDDO

      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      END SUBROUTINE M_bcast_d8_from


!----------------------------------------------------------------------
!
!> copy n double precision complex from inode to all nodes (non-blocking)
!
!----------------------------------------------------------------------

      SUBROUTINE M_ibcast_z_from(COMM, vec, n, inode, request)

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER inode
      INTEGER request

      INTEGER ierror
# 3815


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_ibcast_z_from: invalid vector size n " // str(n), "mpi.F", 3824)
      END IF

      

# 3857


# 3875


# 3880

      CALL MPI_bcast ( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                          ierror )
      request = MPI_REQUEST_NULL

      CALL check_mpi_error(ierror, 'M_ibcast_z_from', 'MPI_ibcast', "mpi.F", 3885)

      

      END SUBROUTINE M_ibcast_z_from


!----------------------------------------------------------------------
!
!> copy n double precision from inode to all nodes (non-blocking)
!
!----------------------------------------------------------------------

      SUBROUTINE M_ibcast_d_from(COMM, vec, n, inode, request)

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER inode
      INTEGER request

      INTEGER ierror
# 3918


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_ibcast_d_from: invalid vector size n " // str(n), "mpi.F", 3927)
      END IF

      

# 3960


# 3978


# 3983

      CALL MPI_bcast ( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
     &                          ierror )
      request = MPI_REQUEST_NULL

      CALL check_mpi_error(ierror, 'M_ibcast_d_from', 'MPI_ibcast', "mpi.F", 3988)

      

      END SUBROUTINE M_ibcast_d_from


!-----------------------------------------------------------------------
!
!> copy n quad precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_qd_from(COMM, vec, n, inode )
      USE mpimy
      USE prec
# 4006

      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n, inode
      REAL(qd) vec(n)

      INTEGER i, ierror
# 4019

      
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_qd_from: invalid vector size n " // str(n), "mpi.F", 4028)
      END IF
# 4050


      CALL MPI_bcast( vec(1), n, MPI_real16, inode-1, COMM%MPI_COMM, &
     &                ierror )
# 4063

      CALL check_mpi_error(ierror, 'M_bcast_qd_from', 'MPI_bcast', "mpi.F", 4064)


      CALL MPI_barrier( COMM%MPI_COMM, ierror )


      

      END SUBROUTINE M_bcast_qd_from

!======================================================================
!
! Global Exchange Routine
!
!======================================================================

!----------------------------------------------------------------------
!
!> complex global exchange routine
!>
!> This maps directly onto MPI_alltoallv which is usually very slow;
!> therefore an alternative implementation optimised for clusters exists
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_z(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      COMPLEX(q) xsnd(*)         !< send buffer
      COMPLEX(q) xrcv(*)         !< receive buffer

      INTEGER psnd(COMM%NCPU+1) !< location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) !< location of data in recv buffer (0 based)
      INTEGER nsnd(COMM%NCPU+1) !< number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) !< number of data recv from each node

!----------------------------------------------------------------------

!----------------------------------------------------------------------
      INTEGER ierror
# 4113


      

# 4144


# 4157


      CALL MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_double_complex, &
                          xrcv(1), nrcv(1), prcv(1), MPI_double_complex, &
                          COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_alltoallv_z', 'MPI_alltoallv', "mpi.F", 4162)

!----------------------------------------------------------------------
# 4234

!----------------------------------------------------------------------

      

      END SUBROUTINE M_alltoallv_z


!----------------------------------------------------------------------
!
!> real cyclic exchange routine which maps directly onto
!
!----------------------------------------------------------------------

      SUBROUTINE M_cycle_d(COMM, xsnd, nsnd)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"


      TYPE(communic) COMM
      REAL(q) xsnd(nsnd)         !< send/ receive buffer
      INTEGER nsnd

      INTEGER ierror, i, j, ichunk
      INTEGER :: tag=202
      INTEGER :: request(2)
      INTEGER :: A_OF_STATUSES(2)
      INTEGER :: NDTMP_=3

      DO j = 1, nsnd, NDTMP_
         ichunk = MIN( nsnd-j+1 , NDTMP_)

! initiate the receive from node-1 (note (0._q,0._q) based)
         i = MOD(COMM%NODE_ME-1-1+COMM%NCPU , COMM%NCPU)
         CALL MPI_irecv( DTMP_m(1), ichunk, MPI_double_precision, &
              i, tag,  COMM%MPI_COMM, request(1), ierror )
         CALL check_mpi_error(ierror, 'M_cycle_d', 'MPI_irecv', "mpi.F", 4271)

! initiate the send
         i = MOD(COMM%NODE_ME-1+1 , COMM%NCPU)  ! i (0._q,0._q) based
         CALL MPI_isend( xsnd(j), ichunk, MPI_double_precision, &
              i, tag,  COMM%MPI_COMM, request(2), ierror )
         CALL check_mpi_error(ierror, 'M_cycle_d', 'MPI_isend', "mpi.F", 4277)

! wait for send and receive to finish
         CALL MPI_waitall(2 , request, A_OF_STATUSES, ierror)
         CALL check_mpi_error(ierror, 'M_cycle_d', 'MPI_waitall', "mpi.F", 4281)

! copy result back
         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  xsnd(j) , 1)
      END DO

      END SUBROUTINE M_cycle_d


!----------------------------------------------------------------------
!
!> integer routine which maps directly onto MPI_alltoallv
!>
!> this is used only once by VASP
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_i(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER xsnd(*)           !< send buffer
      INTEGER xrcv(*)           !< receive buffer

      INTEGER psnd(COMM%NCPU+1) !< location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) !< location of data in recv buffer (0 based)
      INTEGER nsnd(COMM%NCPU+1) !< number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) !< number of data recv from each node

! local data
      INTEGER ierror

      CALL MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_integer, &
                          xrcv(1), nrcv(1), prcv(1), MPI_integer, &
                          COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_alltoall_i', 'MPI_alltoallv', "mpi.F", 4318)

      END SUBROUTINE M_alltoall_i


!----------------------------------------------------------------------
!
!> requires as only input the number of data nsnd send from each node to each other node
!>
!> it assumes a continous data arrangement on sender and receiver
!> and sets up the arrays which are required for MPI_alltoall
!> nrcv, psnd, prcv
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_simple(COMM, nsnd, nrcv, psnd, prcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
! input
      INTEGER nsnd(COMM%NCPU)   !< number of data send to each node
! output
      INTEGER psnd(COMM%NCPU+1) !< location of data in send buffer (0 based)
      INTEGER nrcv(COMM%NCPU)   !< number of data recv from each node
      INTEGER prcv(COMM%NCPU+1) !< location of data in recv buffer (0 based)
! local variable
      INTEGER ierror,i

! only (1._q,0._q) simple MPI_alltoall is required
! to find number of received data on each node
      CALL MPI_alltoall( nsnd(1),  1, MPI_integer, &
                         nrcv(1),  1, MPI_integer, &
                         COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_alltoallv_simple', 'MPI_alltoall', "mpi.F", 4353)

! now set the locations assuming linear arrangement of data
      psnd(1)=0
      prcv(1)=0

      DO i=1,COMM%NCPU
        psnd(i+1)=psnd(i)+nsnd(i)
        prcv(i+1)=prcv(i)+nrcv(i)
      ENDDO

      END SUBROUTINE M_alltoallv_simple


!----------------------------------------------------------------------
!
!> real global exchange routine
!>
!> redistributes an array from distribution over bands to
!> distribution over coefficient (or vice versa)
!>      original distribution             final distribution
!>     |  1  |  2  |  3  |  4  |       |  1  |  1  |  1  |  1  |
!>     |  1  |  2  |  3  |  4  |       |  2  |  2  |  2  |  2  |
!>     |  1  |  2  |  3  |  4  |  <->  |  3  |  3  |  3  |  3  |
!>     |  1  |  2  |  3  |  4  |       |  4  |  4  |  4  |  4  |
!>
!> mind that only (n/NCPU) *NCPU data are exchanged
!> it is the responsability of the user to guarantee that n is
!> correct
!>
!> \param xsnd the array to be redistributed (having n elements)
!> \param xrcv the result array with n/NCPU elements received from each processor
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d(COMM, n, xsnd, xrcv )

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i

      INTEGER, SAVE :: tag=201
      INTEGER       :: in
      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=8000
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests
# 4413


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

!----------------------------------------------------------------------

!----------------------------------------------------------------------
# 4453


# 4466


      CALL MPI_alltoall( xsnd(1), sndcount, MPI_double_precision, &
     &                   xrcv(1), rcvcount, MPI_double_precision, &
     &                   COMM%MPI_COMM, ierror )

      CALL check_mpi_error(ierror, 'M_alltoall_d', 'MPI_alltoall', "mpi.F", 4472)

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CALL check_mpi_error(ierror, 'M_alltoall_d', 'MPI_barrier', "mpi.F", 4475)
!----------------------------------------------------------------------
# 4643

!----------------------------------------------------------------------

      

      END SUBROUTINE M_alltoall_d


      SUBROUTINE M_alltoall_d_omp(COMM, n, xsnd, xrcv )
      USE mpimy
      USE mopenmp_struct_def, ONLY : omp_nthreads_alltoall
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i

      INTEGER, SAVE :: tag=201
      INTEGER       :: in

      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
! test_
!     INTEGER, PARAMETER :: max_=8000
      INTEGER       :: max_
! test_
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests

!$    INTEGER       :: omp_id
!$    INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU
! test_
      max_=(sndcount+omp_nthreads_alltoall-1)/omp_nthreads_alltoall
! test_

!----------------------------------------------------------------------
!  this version  throws all data in (1._q,0._q) go onto all nodes
!  then waits for all sends and receives to finish
!----------------------------------------------------------------------
! initiate the receive on all nodes
! local copy has already been 1._q on each node handle remaining  NCPU-1 packages
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) NUM_THREADS(omp_nthreads_alltoall) &
!$OMP SHARED(sndcount,max_,COMM,xrcv,xsnd) PRIVATE(block,omp_id,sndcount_,p,in,i,request,ierror,A_OF_STATUSES)
      DO block = 0, sndcount-1, max_
!$    omp_id   = OMP_GET_THREAD_NUM()
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

      DO in = 1, COMM%NCPU-1
! send to node in + own node id
! such a construct should allow for an efficient use of the network

         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i (0._q,0._q) based

         CALL MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, 1,  COMM%MPI_COMM, request(in), ierror )
         CALL check_mpi_error(ierror, 'M_alltoall_d_omp', 'MPI_irecv', "mpi.F", 4709)
      ENDDO

! initiate the send on all nodes
      DO in = 1, COMM%NCPU-1
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i (0._q,0._q) based

         CALL MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, 1,  COMM%MPI_COMM, request(in+ COMM%NCPU-1), ierror )
         CALL check_mpi_error(ierror, 'M_alltoall_d_omp', 'MPI_isend', "mpi.F", 4718)
      ENDDO

! local memory copy for data kept on the local node
! overlaps with communication
      CALL DCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )

      CALL MPI_waitall((COMM%NCPU-1)*2, request, A_OF_STATUSES, ierror)
      CALL check_mpi_error(ierror, 'M_alltoall_d_omp', 'MPI_waitall', "mpi.F", 4726)
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE M_alltoall_d_omp


!----------------------------------------------------------------------
!
!> real global exchange routine (single precision)
!>
!> redistributes an array from distribution over bands to
!> distribution over coefficient (or vice versa)
!>      original distribution             final distribution
!>     |  1  |  2  |  3  |  4  |       |  1  |  1  |  1  |  1  |
!>     |  1  |  2  |  3  |  4  |       |  2  |  2  |  2  |  2  |
!>     |  1  |  2  |  3  |  4  |  <->  |  3  |  3  |  3  |  3  |
!>     |  1  |  2  |  3  |  4  |       |  4  |  4  |  4  |  4  |
!>
!> mind that only (n/NCPU) *NCPU data are exchanged
!> it is the responsability of the user to guarantee that n is
!> correct
!>
!> \param xsnd the array to be redistributed (having n elements)
!> \param xrcv the result array with n/NCPU elements received from each processor
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_s(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i

      INTEGER, SAVE :: tag=201
      INTEGER       :: in
      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=8000
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

!----------------------------------------------------------------------

!----------------------------------------------------------------------
      CALL MPI_alltoall( xsnd(1), sndcount, MPI_real, &
     &                   xrcv(1), rcvcount, MPI_real, &
     &                   COMM%MPI_COMM, ierror )
      CALL check_mpi_error(ierror, 'M_alltoall_s', 'MPI_alltoall', "mpi.F", 4786)

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CALL check_mpi_error(ierror, 'M_alltoall_s', 'MPI_barrier', "mpi.F", 4789)

!----------------------------------------------------------------------
# 4957

!----------------------------------------------------------------------

      END SUBROUTINE M_alltoall_s


!----------------------------------------------------------------------
!
!> complex global exchange routine
!>
!> uses #M_alltoall_d with twice as many elements
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_z(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) xsnd(n), xrcv(n)

      CALL M_alltoall_d(COMM, n*2, xsnd, xrcv)

      END SUBROUTINE M_alltoall_z


!----------------------------------------------------------------------
!
!> real global exchange routine (nonblocking)
!>
!> redistributes an array from distribution over bands to
!> distribution over coefficient (or vice versa)
!>      original distribution             final distribution
!>     |  1  |  2  |  3  |  4  |       |  1  |  1  |  1  |  1  |
!>     |  1  |  2  |  3  |  4  |       |  2  |  2  |  2  |  2  |
!>     |  1  |  2  |  3  |  4  |  <->  |  3  |  3  |  3  |  3  |
!>     |  1  |  2  |  3  |  4  |       |  4  |  4  |  4  |  4  |
!>
!> mind that only (n/NCPU) *NCPU data are exchanged
!> it is the responsability of the user to guarantee that n is
!> correct
!>
!> \param xsnd the array to be redistributed (having n elements)
!> \param xrcv the result array with n/NCPU elements received from each processor
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d_async(COMM, n, xsnd, xrcv, tag, srequest, rrequest )

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER tag
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror
      INTEGER i,j,in
# 5023


! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

! local memory copy for data kept on the local node
!! !$ACC PARALLEL LOOP PRESENT(xrcv,xsnd) 
!! !$OMP TARGET TEAMS DISTRIBUTE SIMD IF(OFFLOAD_ON)
      DO i = 1,sndcount
         xrcv((COMM%NODE_ME-1)*sndcount + i) = xsnd((COMM%NODE_ME-1)*sndcount + i)
      ENDDO

! initiate send and receive on all nodes
! local copy has already been 1._q each node send NCPU-1 packages
      j=1
# 5087


# 5110


      DO in = 0, COMM%NCPU-1
! send to node in + own node id
! such a construct should allow for an efficient use of the network

         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)
         IF ( COMM%NODE_ME-1 /= i) THEN
            CALL MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, srequest(j), ierror )
            CALL check_mpi_error(ierror, 'M_alltoall_d_async', 'MPI_isend', "mpi.F", 5120)
            CALL MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, rrequest(j), ierror )
            CALL check_mpi_error(ierror, 'M_alltoall_d_async', 'MPI_irecv', "mpi.F", 5123)
            j=j+1
         ENDIF
      ENDDO

      

      END SUBROUTINE M_alltoall_d_async


      SUBROUTINE M_alltoall_wait(COMM, srequest, rrequest )

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      INTEGER ierror
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,COMM%NCPU)

! wait for the NCPU-1 outstanding packages

# 5157

      CALL MPI_waitall(COMM%NCPU-1, srequest, A_OF_STATUSES, ierror)
      CALL check_mpi_error(ierror, 'M_alltoall_wait', 'MPI_waitall (1)', "mpi.F", 5159)
      CALL MPI_waitall(COMM%NCPU-1, rrequest, A_OF_STATUSES, ierror)
      CALL check_mpi_error(ierror, 'M_alltoall_wait', 'MPI_waitall (2)', "mpi.F", 5161)

      END SUBROUTINE M_alltoall_wait


!----------------------------------------------------------------------
!
!> performs a fast global sum on n doubles in vector vec
!  (algorithm by Kresse Georg)
!>
!> uses complete interchange algorithm
!  (my own invention, but I guess some people must know it)
!> exchange data between nodes, sum locally and
!> interchange back, this algorithm is faster than typical 1 based
!> algorithms (on 8 nodes under MPICH a factor 4)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_d(COMM, vec, n)
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n,ncount,nsummed,ndo,i,j, info, n_,mmax
      REAL(q) vec(n)

      REAL(q), ALLOCATABLE :: vec_inter(:)
! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=8000

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      mmax=MIN(n/COMM%NCPU,max_)
      ALLOCATE(vec_inter(mmax*COMM%NCPU))

      nsummed=0
      n_=n/COMM%NCPU

      DO ndo=0,n_-1,mmax
! forward exchange
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
! sum localy
         DO i=2, COMM%NCPU
           CALL DAXPY(ncount, 1.0_q, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
! backward exchange
         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL DCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_d: " // str(n_) // " " // str(nsummed), "mpi.F", 5230)
      ENDIF

      IF (n-nsummed /= 0 ) &
        CALL M_sumb_d(COMM, vec(nsummed+1), n-nsummed)

      DEALLOCATE(vec_inter)

      END SUBROUTINE M_sumf_d


!----------------------------------------------------------------------
!
!> performs a fast global sum on n doubles in vector vec
! (algorithm by Kresse Georg)
!>
!> This a special version for giant arrays exceeding an element count n
!> of 2^31-1 (the maximum size that can be handled with 32-bit INTEGER).
!> Mainly, the major difference to routine #M_sumf_d is that the last
!> argument ("n") is now declared as 64-bit INTEGER -- and that
!> in addition some intermediate integer operation need also to
!> be 1._q with 64-bit INTEGER temporary variables
!  (adapted by jF)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_d8(COMM, vec, n)

!! USE moffload_struct_def

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
! this must all declared as 64-bit INTEGER
      INTEGER(qi8) n,nsummed,ndo,i,j,n_,mmax
! this must stay to be 32-bit INTEGER and a new "n4" is needed for 1 calls
      INTEGER ncount,info,n4
      REAL(q) vec(n)

      REAL(q), ALLOCATABLE :: vec_inter(:)
! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=8000

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n/COMM%NCPU < 2147483648_qi8) THEN
! here the case that n/COMM%CPU still fits into the 32-bit INTEGER range ...
         n4=n/COMM%NCPU
         mmax=MIN(n4,max_)
      ELSE
! ... if not: (since max_ is of type 32-bit INTEGER) max_ *must* be the minimum
         mmax=max_
      ENDIF
      ALLOCATE(vec_inter(mmax*COMM%NCPU))
!$ACC ENTER DATA CREATE(vec_inter) 

      nsummed=0_qi8
      n_=n/COMM%NCPU

      DO ndo=0_qi8,n_-1_qi8,mmax
! forward exchange
! since the maximum value of mmax is finally determined by max_ (8000)
! which is a 32-bit INTEGER (and hence can never violate the allowed range
! of 32-bit INTEGER numbers) I could safely assume that "ncount" may be of
! type "INTEGER" -- this was essential since DAXPY and DCOPY as well as
! all 1 calls require 32-bit INTEGER arguments (for a standard BLAS);
! be warned that the situation might change if COMM%NCPU goes into the
! millions (since actually not ncount alone but ncount*COMM%NCPU is handed
! over as argument!) but then not only here but in ALL other subroutines we
! would run into trouble -- or needed BLAS/1 routines that can accept
! 64-bit INTEGERS (still rare but already existing ...); at the moment I
! could just recommend to set a small enough value for max_ (8000)
! which must be simply smaller than INT(2147483647/NCPU) ...
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
! sum localy
         DO i=2, COMM%NCPU
           CALL DAXPY(ncount, 1.0_q, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
! replicate data (will be send to each proc)
!$ACC PARALLEL LOOP PRESENT(vec,vec_inter) 
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
! backward exchange
         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL DCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_d8: " // str(n_) // " " // str(nsummed), "mpi.F", 5334)
      ENDIF

! Now the remaining few elements (we can use M_sumb_d for this task since the
! number of elements left cannot be larger than "max_" [being a 32-bit INTEGER])
      IF (n-nsummed /= 0_qi8 ) THEN
        n4=n-nsummed
        CALL M_sumb_d(COMM, vec(nsummed+1), n4)
      ENDIF

!$ACC EXIT DATA DELETE(vec_inter) 
      DEALLOCATE(vec_inter)

      END SUBROUTINE M_sumf_d8


!----------------------------------------------------------------------
!
!> performs a fast global sum on n singles in vector vec
!  (algorithm by Kresse Georg)
!>
!> uses complete interchange algorithm
!  (my own invention, but I guess some people must know it)
!> exchange data between nodes, sum locally and
!> interchange back, this algorithm is faster than typical 1 based
!> algorithms (on 8 nodes under MPICH a factor 4)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_s(COMM, vec, n)
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n,ncount,nsummed,ndo,i,j, info, n_,mmax
      REAL(qs) vec(n)
      REAL(qs), ALLOCATABLE :: vec_inter(:)
! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=8000

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      mmax=MIN(n/COMM%NCPU,max_)
      ALLOCATE(vec_inter(mmax*COMM%NCPU))

      nsummed=0
      n_=n/COMM%NCPU

      DO ndo=0,n_-1,mmax
! forward exchange
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
! sum localy
         DO i=2, COMM%NCPU
           CALL SAXPY(ncount, 1.0, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
! backward exchange
         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL SCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_s: " // str(n_) // " " // str(nsummed), "mpi.F", 5413)
      ENDIF

      IF (n-nsummed /= 0 ) &
        CALL M_sumb_s(COMM, vec(nsummed+1), n-nsummed)

      DEALLOCATE(vec_inter)
      END SUBROUTINE M_sumf_s


!----------------------------------------------------------------------
!
!> performs a fast global sum on n singles in vector vec
!  (algorithm by Kresse Georg)
!> This a special version for giant arrays exceeding an element count n of
!> 2^31-1 (the maximum size that can be handled with 32-bit INTEGER).
!> Mainly, the major difference to routine M_sumf_s is that the last
!> argument ("n") is now declared as 64-bit INTEGER -- and that
!> in addition some intermediate integer operation need also to
!> be 1._q with 64-bit INTEGER temporary variables
!  (adapted by jF)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_s8(COMM, vec, n)
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
! this must all declared as 64-bit INTEGER
      INTEGER(qi8) n,nsummed,ndo,i,j,n_,mmax
! this must stay to be 32-bit INTEGER and a new "n4" is needed for 1 calls
      INTEGER ncount,info,n4
      REAL(qs) vec(n)
      REAL(qs), ALLOCATABLE :: vec_inter(:)
! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=8000

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n/COMM%NCPU < 2147483648_qi8) THEN
! here the case that n/COMM%CPU still fits into the 32-bit INTEGER range ...
         n4=n/COMM%NCPU
         mmax=MIN(n4,max_)
      ELSE
! ... if not: (since max_ is of type 32-bit INTEGER) max_ *must* be the minimum
         mmax=max_
      ENDIF
      ALLOCATE(vec_inter(mmax*COMM%NCPU))

      nsummed=0_qi8
      n_=n/COMM%NCPU

      DO ndo=0_qi8,n_-1_qi8,mmax
! forward exchange
! same story as in M_sumf_d8: we can again safely keep ncount to be INTEGER
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
! sum localy
         DO i=2, COMM%NCPU
           CALL SAXPY(ncount, 1.0, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
! backward exchange
         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL SCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_s8: " // str(n_) // " " // str(nsummed), "mpi.F", 5498)
      ENDIF

! Now the remaining few elements (we can use M_sumb_s for this task since the
! number of elements left cannot be larger than "max_" [being a 32-bit INTEGER])
      IF (n-nsummed /= 0_qi8 ) THEN
        n4=n-nsummed
        CALL M_sumb_s(COMM, vec(nsummed+1), n4)
      ENDIF

      DEALLOCATE(vec_inter)
      END SUBROUTINE M_sumf_s8


!----------------------------------------------------------------------
!
!> performs a fast global sum on n (double) complex in vector 'vec' (see #M_sumf_d)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(2*n)
      CALL M_sumf_d(COMM, vec, 2*n)
      END SUBROUTINE M_sumf_z


!----------------------------------------------------------------------
!
! performs a fast global sum on n (single) complex in vector 'vec' (see #M_sumf_s)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_c(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(2*n)
      CALL M_sumf_s(COMM, vec, 2*n)
      END SUBROUTINE M_sumf_c


!----------------------------------------------------------------------
!
!> performs a sum on n double complex numbers; it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(2*n)


      CALL M_sumb_d(COMM, vec, 2*n)
# 5568

      END SUBROUTINE M_sum_z


!----------------------------------------------------------------------
!
!> performs a sum on n double complex numbers
!>
!> This is for giant arrays with giant dimensions (parameter 2*n potentially
!> exceeding the range of 32-bit INTEGER numbers); this uses the special
!> version M_sum_d8 of routine M_sum_d and just uses 2*n as the
!> corresponding number of doubles to be communicated (like in M_sumf_z)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_z8(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER(qi8) n
      REAL(q) vec(2_qi8*n)

      CALL M_sum_d8(COMM, vec, 2_qi8*n)

      END SUBROUTINE M_sum_z8


!----------------------------------------------------------------------
!
!> performs a sum on n single complex numbers; it uses either sumb_s or sumf_s
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_c(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(2*n)


      CALL M_sumb_s(COMM, vec, 2*n)
# 5618

      END SUBROUTINE M_sum_c


!----------------------------------------------------------------------
!
!> performs a sum on n quadruple complex numbers
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_qd(COMM, vec, n)
      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
# 5634

      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      REAL(qd)             :: vec(n)

      INTEGER  ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sum_qd: invalid vector size n " // str(n), "mpi.F", 5651)
      END IF

      


! invoke in-place version of MPI_allreduce, with own operation
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_real16, &
         &                M_sum_qd_op, COMM%MPI_COMM, ierror )
# 5664

      CALL check_mpi_error(ierror, 'M_sum_qd', 'MPI_allreduce', "mpi.F", 5665)

      

      END SUBROUTINE M_sum_qd


!----------------------------------------------------------------------
!
!> performs a sum on n double prec numbers; it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF


      CALL M_sumb_d(COMM, vec, n)
# 5699

      END SUBROUTINE M_sum_d


!----------------------------------------------------------------------
!
!> performs a sum on n double precision numbers
!>
!> This is for giant arrays with giant dimensions (parameter n potentially
!> exceeding the range of 32-bit INTEGER numbers); this mainly uses special
!> version M_sumf_d_giant of routine M_sumf_d (use of M_sumb_d not always
!> possible for "use_collective_sum" -- only second call always possible)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_d8(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER(qi8) n
      REAL(q) vec(n)
! we need additional local variables
      INTEGER(qi8) max_
      INTEGER n4

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF


! this is now a bit more complicated than in routine M_sum_d ...
      max_=2147483647_qi8          ! 2^31-1 = maximum possible 32-bit integer
      IF ( n>max_) THEN
! here we have no choice; we cannot use M_sumb_d (1 limitations to 32-bit
! integers!) and therefore we have to use M_sumf_d8 instead ... !!
         CALL M_sumf_d8(COMM, vec, n)
      ELSE
! only "n" smaller than 2^31 (limit for 32-bit integers) allows use of M_sumb_d
         n4=n
         CALL M_sumb_d(COMM, vec, n4)
      ENDIF
# 5754

      END SUBROUTINE M_sum_d8


!----------------------------------------------------------------------
!
!> performs a sum on n single prec numbers; it uses either sumb_s or sumf_s
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_single(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)


      CALL M_sumb_s(COMM, vec, n)
# 5780

      END SUBROUTINE M_sum_single


!----------------------------------------------------------------------
!
!> performs a sum on n single prec numbers
!>
!> This is for giant arrays with giant dimensions (parameter n potentially
!> exceeding the range of 32-bit INTEGER numbers); this mainly uses special
!> version M_sumf_s8 of routine M_sumf_s (use of M_sumb_s not always
!> possible for "use_collective_sum" -- only second call always possible);
!> the structure is basically the same as in M_sum_d8 ("same story")
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_single8(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER(qi8) n
      REAL(qs) vec(n)
! local
      INTEGER(qi8) max_
      INTEGER n4


      max_=2147483647_qi8
      IF ( n>max_) THEN
         CALL M_sumf_s8(COMM, vec, n)
      ELSE
         n4=n
         CALL M_sumb_s(COMM, vec, n4)
      ENDIF
# 5823

      END SUBROUTINE M_sum_single8


!----------------------------------------------------------------------
!
!> performs a sum on n double prec numbers to master
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_master_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER :: I,thisdata,ierror

      IF (COMM%NCPU == 1 ) RETURN

      DO i=1,n,8000
         thisdata=min(n-i+1,8000)

         CALL MPI_reduce( vec(i), DTMP_m(1), thisdata, &
                             MPI_double_precision, MPI_sum, &
                             0, COMM%MPI_COMM, ierror )
         CALL check_mpi_error(ierror, 'M_sum_master_d', 'MPI_reduce', "mpi.F", 5850)

         IF (COMM%NODE_ME==1) THEN
            vec(i:i+thisdata-1)=DTMP_m(1:thisdata)
         ENDIF
      END DO

      END SUBROUTINE M_sum_master_d


!----------------------------------------------------------------------
!
!> performs a sum on n double complex numbers to master
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_master_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      CALL M_sum_master_d(COMM,vec,2*n)

      END SUBROUTINE M_sum_master_z

!----------------------------------------------------------------------
!
!> copy nrcv double complex from node i=1,COMM%NCPU to x((i-1)*nrcv+1:i*nrcv)
!> on all nodes using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgather_z(COMM, nrcv, x)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv

      COMPLEX(q) x(*)

! local variables
      INTEGER ierror
# 5901


      IF (COMM%NCPU == 1 ) RETURN

# 5926


# 5938


      CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_double_complex, &
                         x(1), nrcv, MPI_double_complex, &
                         COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_allgather_z', 'MPI_allgather', "mpi.F", 5943)

      END SUBROUTINE M_allgather_z


!----------------------------------------------------------------------
!
!> copy nrcv double precision from node i=1,COMM%NCPU to x((i-1)*nrcv+1:i*nrcv)
!> on all nodes using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgather_d(COMM, nrcv, x)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv

      REAL(q) x(*)

! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

      CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_double_precision, &
                         x(1), nrcv, MPI_double_precision, &
                         COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_allgather_d', 'MPI_allgather', "mpi.F", 5972)

      END SUBROUTINE M_allgather_d

      SUBROUTINE M_allgather_i(COMM, nrcv, x)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv

      INTEGER x(*)

! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

      CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_integer, &
                         x(1), nrcv, MPI_integer, &
                         COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_allgather_i', 'MPI_allgather', "mpi.F", 5993)

      END SUBROUTINE M_allgather_i


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgather_z
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgather_z(COMM, nrcv, x, request)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv
      INTEGER request

      COMPLEX(q) x(*)

! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

      CALL MPI_iallgather(MPI_IN_PLACE, 0, MPI_double_complex, &
                          x(1), nrcv, MPI_double_complex, &
                          COMM%MPI_COMM, request, ierror)
      CALL check_mpi_error(ierror, 'M_iallgather_z', 'MPI_iallgather', "mpi.F", 6022)

      END SUBROUTINE M_iallgather_z


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgather_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgather_d(COMM, nrcv, x, request)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv
      INTEGER request

      REAL(q) x(*)

! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      CALL MPI_iallgather(MPI_IN_PLACE, 0, MPI_double_precision, &
                          x(1), nrcv, MPI_double_precision, &
                          COMM%MPI_COMM, request, ierror)
      CALL check_mpi_error(ierror, 'M_iallgather_d', 'MPI_iallgather', "mpi.F", 6051)

      END SUBROUTINE M_iallgather_d


!----------------------------------------------------------------------
!
!> non-blocking version of #M_gather_z , gathers x from ranks in comm
!> into y
!
!----------------------------------------------------------------------

      SUBROUTINE M_gathero_z(COMM, n, x, y)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      COMPLEX(q) x(*),y(*)

! local variables
      INTEGER ierror
# 6079


      IF (COMM%NCPU == 1) RETURN

      

# 6107


# 6119


      CALL MPI_gather( x(1), n, MPI_double_complex, &
                       y(1), n, MPI_double_complex, &
                       COMM%IONODE-1, COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_gathero_z', 'MPI_gather', "mpi.F", 6124)

      

      END SUBROUTINE M_gathero_z


!----------------------------------------------------------------------
!
!> non-blocking version of #M_gather_d , gathers x from ranks in comm
!> into y
!
!----------------------------------------------------------------------

      SUBROUTINE M_gathero_d(COMM, n, x, y)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      REAL(q) x(*),y(*)

! local variables
      INTEGER ierror
# 6154


      IF (COMM%NCPU == 1) RETURN

      

# 6182


# 6194


      CALL MPI_gather( x(1), n, MPI_double_precision, &
                       y(1), n, MPI_double_precision, &
                       COMM%IONODE-1, COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_gathero_d', 'MPI_gather', "mpi.F", 6199)

      

      END SUBROUTINE M_gathero_d


!----------------------------------------------------------------------
!
!> copy n double complex x(1:n) from nodes i=1,COMM%NCPU to y((i-1)*n+1:i*n)
!> on all nodes using out-of-place communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgathero_z(COMM, n, x, y)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      COMPLEX(q) x(*),y(*)

! local variables
      INTEGER ierror
# 6229


      IF (COMM%NCPU == 1) RETURN

      

# 6257


# 6269


      CALL MPI_allgather( x(1), n, MPI_double_complex, &
                          y(1), n, MPI_double_complex, &
                          COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_allgathero_z', 'MPI_allgather', "mpi.F", 6274)

      

      END SUBROUTINE M_allgathero_z


!----------------------------------------------------------------------
!
!> copy n double precision x(1:n) from nodes i=1,COMM%NCPU to y((i-1)*n+1:i*n)
!> on all nodes using out-of-place communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgathero_d(COMM, n, x, y)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      REAL(q) x(*),y(*)

! local variables
      INTEGER ierror
# 6304


      IF (COMM%NCPU == 1) RETURN

      

# 6332


# 6344


      CALL MPI_allgather( x(1), n, MPI_double_precision, &
                          y(1), n, MPI_double_precision, &
                          COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_allgathero_d', 'MPI_allgather', "mpi.F", 6349)

      

      END SUBROUTINE M_allgathero_d


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgathero_z
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgathero_z(COMM, n, x, y, request)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      INTEGER request

      COMPLEX(q) x(*),y(*)

! local variables
      INTEGER ierror
# 6379


      IF (COMM%NCPU == 1) RETURN

      

# 6407


# 6419


      CALL MPI_iallgather( x(1), n, MPI_double_complex, &
                           y(1), n, MPI_double_complex, &
                           COMM%MPI_COMM, request, ierror)
      CALL check_mpi_error(ierror, 'M_iallgathero_z', 'MPI_iallgather', "mpi.F", 6424)

      

      END SUBROUTINE M_iallgathero_z


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgathero_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgathero_d(COMM, n, x, y, request)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      INTEGER request

      REAL(q) x(*),y(*)

! local variables
      INTEGER ierror
# 6454


      IF (COMM%NCPU == 1) RETURN

      

# 6482


# 6494


      CALL MPI_iallgather( x(1), n, MPI_double_precision, &
                           y(1), n, MPI_double_precision, &
                           COMM%MPI_COMM, request, ierror)
      CALL check_mpi_error(ierror, 'M_iallgathero_d', 'MPI_iallgather', "mpi.F", 6499)

      

      END SUBROUTINE M_iallgathero_d


!----------------------------------------------------------------------
!
!> copy nrcv(i) n double complex from node i=1,COMM%NCPU to
!> c(prcv(i)+1:prcv(i)+nrcv(i)) on all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgatherv_z(COMM, x, c, nrcv, prcv)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv(COMM%NCPU)
      INTEGER prcv(COMM%NCPU)

      COMPLEX(q) x(*), c(*)

! local variables
      INTEGER ierror
# 6531


      IF (COMM%NCPU == 1 ) THEN
!! !$ACC KERNELS PRESENT(x,c) 
!! !$OMP TARGET IF(OFFLOAD_ON)
         c(prcv(1)+1:prcv(1)+nrcv(1)) = x(1:nrcv(1))
!! !$OMP END TARGET
!! !$ACC END KERNELS
         RETURN
      ENDIF

# 6566


# 6577


      CALL MPI_allgatherv(x(1), nrcv(COMM%NODE_ME), MPI_double_complex, &
                          c(1), nrcv(1), prcv(1), MPI_double_complex, &
                          COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_allgatherv_z', 'MPI_allgatherv', "mpi.F", 6582)

      END SUBROUTINE M_allgatherv_z


!----------------------------------------------------------------------
!
!> copy nrcv(i) n double precision from node i=1,COMM%NCPU to
!> c(prcv(i)+1:prcv(i)+nrcv(i)) on all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgatherv_d(COMM, x, c, nrcv, prcv)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv(COMM%NCPU)
      INTEGER prcv(COMM%NCPU)

      REAL(q) x(*), c(*)

! local variables
      INTEGER ierror
# 6612


      IF (COMM%NCPU == 1 ) THEN
!! !$ACC KERNELS PRESENT(x,c) 
!! !$OMP TARGET IF(OFFLOAD_ON)
         c(prcv(1)+1:prcv(1)+nrcv(1)) = x(1:nrcv(1))
!! !$OMP END TARGET
!! !$ACC END KERNELS
         RETURN
      ENDIF

# 6647


# 6658


      CALL MPI_allgatherv(x(1), nrcv(COMM%NODE_ME), MPI_double_precision, &
                           c(1), nrcv(1), prcv(1), MPI_double_precision, &
                           COMM%MPI_COMM, ierror)
      CALL check_mpi_error(ierror, 'M_allgatherv_d', 'MPI_allgatherv', "mpi.F", 6663)

      END SUBROUTINE M_allgatherv_d


!----------------------------------------------------------------------
!
!> reduce (MPI_sum) n complex doubles to node i from all nodes to node i
!> using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_reduce_z_to(COMM, vec, n, inode)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER inode

! local variables
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      IF (COMM%NODE_ME==inode) THEN
         CALL MPI_reduce( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ELSE
         CALL MPI_reduce( vec(1),       vec(1), n, MPI_double_complex, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ENDIF

      CALL check_mpi_error(ierror, 'M_reduce_z_to', 'MPI_reduce', "mpi.F", 6698)

      END SUBROUTINE M_reduce_z_to


!----------------------------------------------------------------------
!
!> reduce (MPI_sum) n doubles to node i from all nodes to node i using
!> inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_reduce_d_to(COMM, vec, n, inode)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER inode

! local variables
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      IF (COMM%NODE_ME==inode) THEN
         CALL MPI_reduce( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ELSE
         CALL MPI_reduce( vec(1),       vec(1), n, MPI_double_precision, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ENDIF

      CALL check_mpi_error(ierror, 'M_reduce_d_to', 'MPI_reduce', "mpi.F", 6733)

      END SUBROUTINE M_reduce_d_to


!----------------------------------------------------------------------
!
!> non-blocking reduce (MPI_sum) n complex doubles from all nodes to node i
!> using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_ireduce_z_to(COMM, vec, n, inode, request)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER inode
      INTEGER request

! local variables
      INTEGER ierror
# 6762


! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      

# 6809


      IF (COMM%NODE_ME==inode) THEN
# 6815

         CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)

      ELSE
# 6823

         CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_complex, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)

      ENDIF

      CALL check_mpi_error(ierror, 'M_ireduce_z_to', 'MPI_ireduce', "mpi.F", 6829)

      request = MPI_REQUEST_NULL

      

      END SUBROUTINE M_ireduce_z_to


!----------------------------------------------------------------------
!
!> non-blocking reduce (MPI_sum) n doubles to node i from all nodes to node i
!> using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_ireduce_d_to(COMM, vec, n, inode, request)

!! USE moffload_struct_def

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER inode
      INTEGER request

! local variables
      INTEGER ierror
# 6862


! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      

# 6909


      IF (COMM%NODE_ME==inode) THEN
# 6915

         CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)

      ELSE
# 6923

         CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_precision, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)

      ENDIF

      CALL check_mpi_error(ierror, 'M_ireduce_d_to', 'MPI_ireduce', "mpi.F", 6929)

      request = MPI_REQUEST_NULL

      

      END SUBROUTINE M_ireduce_d_to


!----------------------------------------------------------------------
!
!> just a wrapper of MPI_waitall
!
!----------------------------------------------------------------------

      SUBROUTINE M_waitall(n,requests)
      USE mpimy
      IMPLICIT NONE
      INTEGER n,requests(n)
! local variables
      INTEGER ierror

      CALL MPI_waitall(n, requests, MPI_STATUSES_IGNORE, ierror)

      CALL check_mpi_error(ierror, 'M_waitall', 'MPI_waitall', "mpi.F", 6953)

      END SUBROUTINE M_waitall

# 6981

