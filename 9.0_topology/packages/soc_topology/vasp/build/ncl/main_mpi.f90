# 1 "main_mpi.F"
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


# 2 "main_mpi.F" 2 
!***********************************************************************
! RCS:  $Id: main_mpi.F,v 1.3 2002/04/16 07:28:45 kresse Exp $
!
!> This module initializes the communication universe used in VASP
!
!***********************************************************************
      MODULE main_mpi
      USE prec
      USE base
      USE mpimy
      IMPLICIT NONE


      TYPE (communic),TARGET :: COMM_WORLD !< our world wide communicator
      TYPE (communic),TARGET :: COMM_CHAIN !< communication between images
      TYPE (communic),TARGET :: COMM       !< (1._q,0._q) image communicator
      TYPE (communic),TARGET :: COMM_KINTER !< between k-points communicator
      TYPE (communic),TARGET :: COMM_KIN   !< in k-point communicator
      TYPE (communic),TARGET :: COMM_INTER !< between-band communicator
      TYPE (communic),TARGET :: COMM_INB   !< in-band communicator
# 25

      TYPE (communic),TARGET :: COMM_intra_node_world,COMM_inter_node_world
      INTEGER :: IMAGES=0                  !< total number of images
      INTEGER :: KIMAGES=0                 !< distribution of k-points
!> if IMAGES is 2 it is possible to specify the number of cores that are in the first image
      INTEGER :: NCORE_IN_IMAGE1=0
      CHARACTER(LEN=10) ::  DIR_APP
      INTEGER           ::  DIR_LEN=0
      INTEGER           ::  KPAR,NCORE
# 36

      INTEGER, PARAMETER :: NCSHMEM=1


      CONTAINS
!***********************************************************************
!
!> initialize 1 and sub-divide communicators based on IMAGES, KPAR, and NCORE
!>
!> because redistribution of nodes might be 1._q (1._q,0._q) must do this as early as possible
!> @ref openmp :
!> forces main_mpi::ncshmem=1 and main_mpi::NCORE=1, and gives
!> additional output.
!
!***********************************************************************


      SUBROUTINE INIT_MPI(NPAR,IO)

!! USE moffload_struct_def

      USE command_line
      USE reader_tags
      USE incar_reader, ONLY: INCAR_FROM_FILE
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      USE vhdf5_base


      TYPE (in_struct) :: IO
      INTEGER          :: NPAR
! local
      INTEGER :: IERR, IDUM, UNIT_
      LOGICAL :: LVCAIMAGES, INCAR_SUBDIR
      REAL(q) :: VCAIMAGES

      INTEGER :: IH5ERR
# 76

!$    INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

# 82


      CALL M_init(COMM_WORLD)

      CALL PARSE_COMMAND_LINE(IO, COMM_WORLD%NODE_ME == COMM_WORLD%IONODE)

! in case a demo is compiled
      CALL VASP_DEMO( COMM_WORLD, IO%IU0 )


      IH5ERR = VH5_START()
!
! incar, vaspin
!
      INQUIRE(FILE=INCAR,EXIST=INCAR_FOUND)
      INQUIRE(FILE=VASPIN,EXIST=HDF5_FOUND)
      IF (.NOT.HDF5_FOUND) THEN
       IF (.NOT.INCAR_FOUND) CALL vtutor%error("No INCAR or " // VASPIN // " found, STOPPING")
      ELSE
        IH5ERR = VH5_FILE_OPEN_READ(VASPIN, IH5INFILEID)
        IH5ERR = VH5_GROUP_OPEN(IH5INFILEID, GRP_INPUT, IH5ININPUTGROUP_ID)
        INCAR_F%FROM_HDF5 = .TRUE.
        INCAR_F%ERROR = ""
        ALLOCATE(INCAR_F%TAGS(0))
      ENDIF
# 112

      IF (INCAR_FOUND) INCAR_F = INCAR_FROM_FILE(INCAR)

      VCAIMAGES=-1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'VCAIMAGES',VCAIMAGES,IERR)

      IF (VCAIMAGES==-1) THEN
         LVCAIMAGES=.FALSE.
      ELSE
         LVCAIMAGES=.TRUE.
      ENDIF

      NCORE_IN_IMAGE1=0
      IF (LVCAIMAGES) THEN
         IMAGES=2
         CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'NCORE_IN_IMAGE1',NCORE_IN_IMAGE1,IERR)
      ELSE
         IMAGES=0
         CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'IMAGES',IMAGES,IERR)
      ENDIF

      KIMAGES=0
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'KIMAGES',KIMAGES,IERR)
      IF (KIMAGES>0) THEN
         IDUM=0
         CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'FOURORBIT',IDUM,IERR)
         IF (IDUM/=1) THEN
            CALL vtutor%error("Distribution of k-points over KIMAGES only works \n in combination with &
               &FOURORBIT=1, sorry, stopping ...")
         ENDIF
      ENDIF
! KPAR division of kpoints,  default to unity in case only 1 k-point
      KPAR=1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'KPAR',KPAR,IERR)
      IF (KPAR>1.AND.KIMAGES>0) THEN
         CALL vtutor%error("Untested combination of FOURORBIT with KPAR \n (k-point parallelization),&
            & sorry, stopping...")
      END IF

# 153


!----------------------------------------------------------------------
!> # Image parallelization
!>
!> Creates a 2 dimensional cartesian topology for seperate images or
!> work groups. Each work group (image) will run VASP independently in
!> (1._q,0._q) sub directory (01-99) of the current directory. This mode is
!> required to support either independent calculations on parallel machines
!> or the nudged elastic band method.
!----------------------------------------------------------------------
      IF (ABS(IMAGES)>0) THEN
         IF (IMAGES==2 .AND. NCORE_IN_IMAGE1/=0) THEN
            CALL M_initc( COMM_WORLD)
! M_divide2 calls M_initc for COMM and COMM_CHAIN
            CALL M_divide2( COMM_WORLD, NCORE_IN_IMAGE1, COMM_CHAIN, COMM)
         ELSE
            CALL M_divide( COMM_WORLD, ABS(IMAGES), COMM_CHAIN, COMM, .TRUE. )
            CALL M_initc( COMM_WORLD)
            CALL M_initc( COMM)
            CALL M_initc( COMM_CHAIN)
         ENDIF
      ELSEIF (KIMAGES>0) THEN
         CALL M_divide( COMM_WORLD, KIMAGES, COMM_CHAIN, COMM, .TRUE. )
         CALL M_initc( COMM_WORLD)
         CALL M_initc( COMM)
         CALL M_initc( COMM_CHAIN)
      ELSE
         CALL M_initc( COMM_WORLD)
         COMM=COMM_WORLD
         COMM_CHAIN=COMM
      ENDIF

      CALL M_divide_intra_inter_node(COMM_WORLD,COMM_WORLD,COMM_intra_node_world,COMM_inter_node_world)

# 192

      IF ( COMM_WORLD%NODE_ME == 1 ) &
         WRITE(IO%IU0,'(" running ",I4," mpi-ranks, on ",I4," nodes")') &
             COMM_WORLD%NCPU,COMM_inter_node_world%NCPU

      IF ( COMM_WORLD%NODE_ME == 1 .AND. (ABS(IMAGES)>0.OR.KIMAGES>0) ) &
         WRITE(IO%IU0,'(" each image running on ",I4," cores")') COMM%NCPU

      IF (KPAR>=1) THEN
!----------------------------------------------------------------------
!> # k-point parallelization
!>
!> Creates a 2 dimensional cartesian topology within (1._q,0._q) work group (image).
!> This is required for simultaneous distribution over k-points and bands
!----------------------------------------------------------------------
         CALL M_divide( COMM, KPAR, COMM_KINTER, COMM_KIN, .FALSE.)
         CALL M_initc( COMM)   ! probably not required but who knows
         CALL M_initc( COMM_KINTER)
         CALL M_initc( COMM_KIN)
         IF ( COMM_WORLD%NODE_ME == 1 ) &
         WRITE(IO%IU0,'(" distrk:  each k-point on ",I4," cores, ",I4," groups")') &
                        COMM_KIN%NCPU,COMM_KINTER%NCPU
      ELSE
         COMM_KINTER = COMM
         COMM_KIN    = COMM
      ENDIF

! NCORE species onto how many cores a band is distributed
! often this values can be now set to the number of cores per node
! this is more handy than NPAR in most cases
      NCORE=1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'NCORES_PER_BAND',NCORE,IERR)
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'NCORE',NCORE,IERR)
      NCORE = MAX(MIN(COMM_KIN%NCPU, NCORE), 1)
      IF (MOD(COMM_KIN%NCPU, NCORE) /= 0) NCORE = 1

! NPAR number of bands distributed over processors, defaults to COMM_KIN%NCPU/NCORE
      NPAR = 0
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'NPAR',NPAR,IERR)
      IF ((NPAR < 1).OR.(NPAR > COMM_KIN%NCPU)) THEN
         NPAR = MAX(COMM_KIN%NCPU / NCORE, 1)
      ELSEIF (MOD(COMM_KIN%NCPU, NPAR) /= 0) THEN
         NPAR = MAX(COMM_KIN%NCPU / NCORE, 1)
      ELSE
         NCORE = COMM_KIN%NCPU / NPAR
      ENDIF

# 244

# 250


      IF (NPAR>=1) THEN
!----------------------------------------------------------------------
!> # Band and FFT parallelization
!>
!> Creates a 2 dimensional cartesian topology within (1._q,0._q) work group (image)
!> This is required for simultaneous distribution over bands and plane
!> wave coefficients. Communicators are created for inter-band in intra-band
!> communication. *NPAR* is the number of bands over which is parallelized.
!> The resulting layout will be the following (NPAR=4, NCPU=8):
!> (1 uses allways Row-major layout)
!>
!>     wave1        0 (0,0)        1 (0,1)
!>     wave2        2 (1,0)        3 (1,1)
!>     wave3        4 (2,0)        5 (2,1)
!>     wave4        6 (3,0)        7 (3,1)
!>     wave5        0 (0,0)        1 (0,1)
!>     etc.
!>
!> The sub-communicators #COMM_INB are (1._q,0._q) dimensional communicators
!> which allow communication within (1._q,0._q) row i.e. nodes are grouped to
!>
!>     0-1      2-3     4-5      6-7
!>
!> For shortness these groups will we called in-band-groups and
!> their communicators are called in-band-communicator.
!>
!> The sub-communicators #COMM_INTER are (1._q,0._q) dimensional communicators
!> which allow communication within (1._q,0._q) column i.e nodes are grouped
!>
!>     0-2-4-6      1-3-5-7
!>
!> These groups will we called inter-band-groups.
!>
!> The most complicated thing is the FFT of soft charge densities.
!> The following algorithm is used:
!> * The soft charge density is calculated in real space
!> * Charge from all bands is merged to processor 0 and 1 using #COMM_INTER
!> * A FFT involving all processors is 1._q (on processors 2-7 no components
!>   exist in real space)
!> * The final result in reciprocal space is defined on all processors
!>   (see #mgrid::SET_RL_GRID for more information)
!----------------------------------------------------------------------
         CALL M_divide( COMM_KIN, NPAR, COMM_INTER, COMM_INB, .FALSE.)
         CALL M_initc( COMM_KIN)   ! propably not required but who knows
         CALL M_initc( COMM_INTER)
         CALL M_initc( COMM_INB)
         IF ( COMM_WORLD%NODE_ME == 1 ) &
         WRITE(IO%IU0,'(" distr:  one band on ",I4," cores, ",I4," groups")') &
                        COMM_INB%NCPU,COMM_INTER%NCPU
      ELSE
         COMM_INTER = COMM
         COMM_INB   = COMM
      ENDIF

# 312


      IF (COMM%NODE_ME /= COMM%IONODE) IO%IU6 = -1
      IF (COMM%NODE_ME /= COMM%IONODE) IO%IU0 = -1

      IF (KIMAGES>0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU6 = -1
      IF (KIMAGES>0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU0 = -1

      IF (IMAGES<0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU6 = -1
      IF (IMAGES<0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU0 = -1

      IF (KPAR>1.AND.COMM%NODE_ME/=COMM%IONODE) IO%IU6 = -1
      IF (KPAR>1.AND.COMM%NODE_ME/=COMM%IONODE) IO%IU0 = -1

      IF ( IO%IU0/=-1 .AND. IO%IU6 == -1) THEN
         CALL vtutor%bug("internal ERROR: io-unit problem " // str(IO%IU0) // " " // str(IO%IU6), "main_mpi.F", 327)
      ENDIF

      CALL MAKE_DIR_APP(COMM_CHAIN%NODE_ME)
      IMAGES = ABS(IMAGES)

! if all nodes should write (giving complete mess) do not use the
! following line
      IF (COMM_WORLD%NODE_ME /= COMM_WORLD%IONODE .AND. IO%IU0>0) THEN
         OPEN(UNIT=IO%IU0,FILE=DIR_APP(1:DIR_LEN)//'stdout',STATUS='UNKNOWN')
      ENDIF
!----------------------------------------------------------------------
! try to go to subdir INCAR's for the rest
!----------------------------------------------------------------------
! TODO fix images for HDF5
      IF (DIR_LEN > 0) THEN
         INQUIRE(FILE=DIR_APP(1:DIR_LEN)//INCAR, EXIST=INCAR_SUBDIR)
         IF (INCAR_SUBDIR) THEN
            INCAR=DIR_APP(1:DIR_LEN)//INCAR
            IF (IO%IU0>=0) WRITE(IO%IU0,*) 'using from now: ',INCAR
            INCAR_F = INCAR_FROM_FILE(INCAR)
         ELSE
            INCAR_F%IMAGE = COMM_CHAIN%NODE_ME
         ENDIF
      ENDIF


!
! if no INCAR, try to close vaspin.h5
!
      IF (.NOT. INCAR_FOUND .AND. HDF5_FOUND) THEN
        IH5ERR = VH5_GROUP_CLOSE(IH5ININPUTGROUP_ID)
        IH5ERR = VH5_FILE_CLOSE(IH5INFILEID)
        IH5ERR = VH5_END()
      ENDIF

# 372

      RETURN
      END SUBROUTINE

!***********************************************************************
!
!> make the  directory entry which is used for fileio
!
!***********************************************************************

      SUBROUTINE MAKE_DIR_APP(node)
      INTEGER node

! in principle (1._q,0._q) can chose here any string (1._q,0._q) wants to use
! only DIR_LEN must be adjusted
      WRITE (DIR_APP  , "(I1,I1,'/')") MOD(node/10,10),MOD(node,10)
      IF (IMAGES<=0.OR.KIMAGES>0) THEN
         DIR_LEN=0
      ELSE
         DIR_LEN=3
      ENDIF

      END SUBROUTINE MAKE_DIR_APP

!***********************************************************************
!
!> once unit 6 is open write number of nodes and all other parameters
!> with additional output for @ref openmp
!
!***********************************************************************
      SUBROUTINE WRT_DISTR(IU6)

!! USE moffload_struct_def, ONLY : OFFLOAD_NUM_DEVICES

      USE mopenmp_struct_def, ONLY : omp_nthreads
      INTEGER IU6

      IF (IU6>=0) THEN
# 413

        WRITE(IU6,'(" running ",I4," mpi-ranks, on ",I4," nodes")') &
                        COMM_WORLD%NCPU,COMM_inter_node_world%NCPU

        IF (IMAGES>0 ) &
        WRITE(IU6,'(" each image running on ",I4," cores")') COMM%NCPU
        WRITE(IU6,'(" distrk:  each k-point on ",I4," cores, ",I4," groups")') &
                        COMM_KIN%NCPU,COMM_KINTER%NCPU
        WRITE(IU6,'(" distr:  one band on NCORE=",I4," cores, ",I4," groups")') &
                        COMM_INB%NCPU,COMM_INTER%NCPU

!!   WRITE(IU6,'(" Offloading initialized ... ",I4," GPUs detected")') OFFLOAD_NUM_DEVICES
      ENDIF
# 428

      END SUBROUTINE WRT_DISTR



!***********************************************************************
!> VASP DEMO
!***********************************************************************
      SUBROUTINE VASP_DEMO( COMM_WORLD, IU0 )
      USE tutor, ONLY: vtutor
      TYPE( communic )  :: COMM_WORLD
      INTEGER           :: IU0
! local
      INTEGER, PARAMETER :: NMAX = 4
!$    INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

! immeadiate return if not DEMO

      RETURN

! write a banner to indicate that demo is executed
      IF ( COMM_WORLD%NODE_ME == 1  ) THEN
         WRITE(IU0,*)'-----------------------------------------------------------------------------'
# 453

         IF( COMM_WORLD%NCPU > NMAX ) THEN

         WRITE(IU0,*)'         This                                                                 '
         ENDIF
         WRITE(IU0,*)'                                                                              '
         WRITE(IU0,*)'          V    V   AA    SSSS  PPPPP    DDDDD  EEEEEE  MM MM   OOOO           '
         WRITE(IU0,*)'          V    V  A  A  S    S P    P   D    D E      M  M  M O    O          '
         WRITE(IU0,*)'          V    V A    A  SSS   PPPP     D    D EEEEE  M  M  M O    O          '
         WRITE(IU0,*)'          V    V AAAAAA      S P        D    D E      M  M  M O    O          '
         WRITE(IU0,*)'           V  V  A    A S    S P        D    D E      M     M O    O          '
         WRITE(IU0,*)'            VV   A    A  SSSS  P        DDDDD  EEEEEE M     M  OOOO           '
         WRITE(IU0,*)'                                                                              '
      ENDIF

! if requirements are not fullfiled, stop here
# 471

      IF( COMM_WORLD%NCPU > NMAX ) THEN

         IF ( COMM_WORLD%NODE_ME == 1  ) THEN
            WRITE(IU0,'(A,I2,A,I2,A)') '          is restricted to maximally',NMAX,&
                                              ' MPI Ranks and',NMAX,' OpenMP threads!'
            WRITE(IU0,*)'-----------------------------------------------------------------------------'
         ENDIF
         CALL vtutor%stopCode()
      ENDIF
      IF ( COMM_WORLD%NODE_ME == 1  ) &
         WRITE(IU0,*)'-----------------------------------------------------------------------------'
      END SUBROUTINE


      END MODULE
