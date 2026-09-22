# 1 "greens_orbital.F"
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


# 2 "greens_orbital.F" 2 
!*********************************************************************
!
!
! this module implements handling of the Green's function
! in the orbital basis
! 1 is used to store the matrices
!
!
!*********************************************************************
MODULE greens_orbital
  USE prec
  USE mpimy
  USE GG_base
  USE chi_glb
  USE wave_struct_def
  USE tutor
  USE string, ONLY: str
  USE kpoints_change
  USE mathtools, ONLY : INVERT_REAL_MATRIX
  USE minimax_struct, ONLY : imag_grid_handle
  USE scala 

  IMPLICIT NONE

!
! values discriminating occupied from unoccupied states
! in the reference Green's function  subtracted to accelerate
! Matsubara sum (the value must NOT matter)
!
  REAL(q), PRIVATE :: FERMI_OCC=0.50_q

  CONTAINS

   SUBROUTINE FT_G_TEST(W, GDES_MAT, GU_MAT, GO_MAT, GRIDS, IO )
      TYPE (greens_mat_des), POINTER :: GDES_MAT   ! out: pointer to generated descriptor
      TYPE (wavespin)          :: W
      COMPLEX(q)                     :: GU_MAT(:,:,:,:)
      COMPLEX(q)                     :: GO_MAT(:,:,:,:)
      TYPE( imag_grid_handle ) :: GRIDS
      TYPE( in_struct )        :: IO 
! local
      INTEGER                  :: NKPTS_IRZ
      INTEGER                  :: NTAU_ROOT, NOMEGA_ROOT
      INTEGER                  :: ISP, NK1
      COMPLEX(q), POINTER, CONTIGUOUS :: SIGMAW_MAT(:,:,:,:)=>NULL() ! frequency dependent self-energy

      NKPTS_IRZ = SIZE( GU_MAT, 3 )

      NULLIFY(SIGMAW_MAT)
      CALL ALLOCATE_GREENS_MAT(GDES_MAT, SIGMAW_MAT, NKPTS_IRZ, W%WDES%ISPIN, &
      GRIDS%T%NPOINTS_IN_ROOT_GROUP)
      SIGMAW_MAT=(0._q,0._q)
      
      IF ( IO%IU0>=0 ) OPEN(99,FILE='SIGMA_TAU.dat',STATUS='REPLACE')
      IF ( IO%IU0>=0 ) OPEN(100,FILE='SIGMA_OMEGA.dat',STATUS='REPLACE')
      IF ( IO%IU0>=0 ) OPEN(101,FILE='SIGMA_INV.dat',STATUS='REPLACE')
      DO ISP = 1, W%WDES%ISPIN
         DO NK1 = 1, NKPTS_IRZ
            DO NTAU_ROOT = 1, GRIDS%T%NPOINTS_IN_ROOT_GROUP
               CALL SETUP_INDICES( NTAU_ROOT, GRIDS ) 
! dump in serial mode only
               IF( GRIDS%T%NPOINTS == GRIDS%T%NPOINTS_IN_GROUP )THEN
                  IF( IO%IU0>=0 ) WRITE(99,*)NTAU_ROOT,NK1,ISP,GRIDS%TAU(NTAU_ROOT)
                 CALL DUMP_GREENS_HAM_GDEF( "SIGMAU_TAU", GU_MAT(:,NK1,ISP,NTAU_ROOT),&
                    W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0, 99)
                 CALL DUMP_GREENS_HAM_GDEF( "SIGMAO_TAU", GO_MAT(:,NK1,ISP,NTAU_ROOT),&
                    W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0, 99)
               ENDIF
! perform forward Fourier transform
               CALL FT_SIGMA_MAT( W, GO_MAT(:, NK1, ISP, NTAU_ROOT), GU_MAT(:, NK1, ISP, NTAU_ROOT), & 
                    SIGMAW_MAT(:, NK1, ISP,:), NTAU_ROOT, GRIDS, GDES_MAT )
            ENDDO

! dump  fourier coefficients
            DO NOMEGA_ROOT = 1, GRIDS%F%NPOINTS_IN_ROOT_GROUP
               IF( IO%IU0>=0 )  WRITE(100,*)NOMEGA_ROOT,NK1,ISP,GRIDS%FER_RE(NOMEGA_ROOT)
               CALL DUMP_GREENS_HAM( "SIGMA_OMEGA", SIGMAW_MAT(:,NK1,ISP,NOMEGA_ROOT),&
                 W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0, 100)
            ENDDO

 
! perform inverse transform as a check
            DO NTAU_ROOT=1, GRIDS%T%NPOINTS_IN_ROOT_GROUP
               CALL SETUP_INDICES( NTAU_ROOT, GRIDS )
! transform the correlated part of Greens function to imaginary time
! occ and unocc components SIGMAW_MAT -> GU_MAT and GO_MAT
               CALL INV_FT_SIGMA_MAT(W, W%WDES%NB_TOT, GO_MAT(:, NK1, ISP, NTAU_ROOT),&
                 GU_MAT(:, NK1, ISP, NTAU_ROOT), SIGMAW_MAT(:, NK1, ISP,:), &
                 NTAU_ROOT, GRIDS, GDES_MAT )

               IF( IO%IU0>=0 )  WRITE(101,*)NTAU_ROOT,NK1,ISP,GRIDS%TAU(NTAU_ROOT)
               CALL DUMP_GREENS_HAM_GDEF( "SIGMU_INV", GU_MAT(:,NK1,ISP,NTAU_ROOT),&
                  W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0, 101)
              CALL DUMP_GREENS_HAM_GDEF( "SIGMO_INV", GO_MAT(:,NK1,ISP,NTAU_ROOT),&
                 W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0, 101)
            ENDDO
         ENDDO
      ENDDO
      IF( IO%IU0>=0 ) CLOSE(99)
      IF( IO%IU0>=0 ) CLOSE(100)
      IF( IO%IU0>=0 ) CLOSE(101)

      CONTAINS 

      SUBROUTINE SETUP_INDICES( NTAU_ROOT, GRIDS )
         INTEGER                  :: NTAU_ROOT
         TYPE( imag_grid_handle ) :: GRIDS
!if the local number of tau points in the group coincides with the
!(1._q,0._q) of the root group, the tau loop variable is the same
         GRIDS%T%NPOINTSC = NTAU_ROOT
!other groups may have less tau points locally
         GRIDS%T%LDO_POINT_LOCAL = .TRUE.
         IF ( GRIDS%T%NPOINTSC > GRIDS%T%NPOINTS_IN_GROUP ) THEN
! set to last tau point in group to avoid going over array bounds
            GRIDS%T%NPOINTSC = GRIDS%T%NPOINTS_IN_GROUP
! but actually little work is 1._q since most expensive calls are bypassed
            GRIDS%T%LDO_POINT_LOCAL = .FALSE.
         ENDIF    
! set current imaginary time tau point
         GRIDS%T%POINT_CURRENT = GRIDS%T%POINTS_LOCAL(GRIDS%T%NPOINTSC)
!determine current smallest TAU value of all groups
!this is always the tau value of the first tau group
!this weird construction here is since M_min_d is not available
         GRIDS%T%POINTMIN_CURRENT = -GRIDS%T%POINT_CURRENT
         CALL M_max_d(GRIDS%T%COMM, GRIDS%T%POINTMIN_CURRENT, 1)
         GRIDS%T%POINTMIN_CURRENT = -GRIDS%T%POINTMIN_CURRENT
      END SUBROUTINE SETUP_INDICES

      SUBROUTINE FORWARD_FT( GU_I, GO_I, G_RE, G_IM, GRIDS, N_LOCAL, IO )
         REAL(q)                  :: GU_I
         REAL(q)                  :: GO_I
         COMPLEX(q)               :: G_RE(:)
         COMPLEX(q)               :: G_IM(:)
         TYPE( imag_grid_handle ) :: GRIDS
         INTEGER                  :: N_LOCAL
         TYPE( in_struct )        :: IO
! local
         COMPLEX(q)               :: FT(GRIDS%NOMEGA,GRIDS%NOMEGA)
         INTEGER                  :: I, J 


         DO I = 1, GRIDS%NOMEGA
            FT(I,:)=GRIDS%TO_FER_RE(I,:)+(0._q,1._q)*GRIDS%TO_FER_RE_CONJG(I,:)
         ENDDO


! add forward transform
         DO I = 1, GRIDS%NOMEGA 
            G_RE( I ) = G_RE( I ) + GRIDS%TO_FER_RE( I , N_LOCAL )*GO_I
            G_IM( I ) = G_IM( I ) + GRIDS%TO_FER_RE_CONJG( I , N_LOCAL )*GO_I
         ENDDO

      END SUBROUTINE FORWARD_FT

      SUBROUTINE BACKWARD_FT( GU, GO, G_RE, G_IM, GRIDS, IO )
         COMPLEX(q)               :: GU(:)
         COMPLEX(q)               :: GO(:)
         COMPLEX(q)               :: G_RE(:)
         COMPLEX(q)               :: G_IM(:)
         TYPE( imag_grid_handle ) :: GRIDS
         INTEGER                  :: N_LOCAL
         TYPE( in_struct )        :: IO
! local
         COMPLEX(q)               :: FT(GRIDS%NOMEGA,GRIDS%NOMEGA)
         INTEGER                  :: I, J 


         DO I = 1, GRIDS%NOMEGA
            FT(I,:)=GRIDS%TO_FER_RE(:,I)-(0._q,1._q)*GRIDS%TO_FER_RE_CONJG(:,I)
            FT(I,:)=FT(I,:)*GRIDS%T%NPOINTS
         ENDDO

         IF( IO%IU0>=0 ) THEN
            WRITE(*,*)'Inverse FT MATRIX'
            DO I = 1, GRIDS%NOMEGA
            WRITE(*,'(8F12.6)') REAL(FT(I,:))
            ENDDO
            WRITE(*,*)'imag'
            DO I = 1, GRIDS%NOMEGA
            WRITE(*,'(8F12.6)') IMAG(FT(I,:))
            ENDDO
         ENDIF

! backward transform
         GU=0
         GO=0
         DO J = 1, GRIDS%T%NPOINTS
         DO I = 1, GRIDS%F%NPOINTS
            GU( J ) = GU( J ) + GRIDS%TO_FER_RE(I,J)*GRIDS%T%NPOINTS*REAL(G_RE(I),q)
            GU( J ) = GU( J ) + GRIDS%TO_FER_RE_CONJG(I,J)*GRIDS%T%NPOINTS*REAL(G_IM(I),q)
! imaginary part
            GU( J ) = GU( J ) -(0._q,1._q)*GRIDS%TO_FER_RE_CONJG(I,J)*GRIDS%T%NPOINTS*REAL(G_RE(I),q)
            GU( J ) = GU( J ) -(0._q,1._q)*GRIDS%TO_FER_RE(I,J)*GRIDS%T%NPOINTS*REAL(G_IM(I),q)
           
! occupied and unoccupied green's function are identical
            GO( J ) = GO(J ) + GRIDS%TO_FER_RE(I,J)*GRIDS%T%NPOINTS*REAL(G_RE(I),q)
            GO( J ) = GO(J ) + GRIDS%TO_FER_RE_CONJG(I,J)*GRIDS%T%NPOINTS*REAL(G_IM(I),q)

! imaginary part
            GO( J ) = GO(J ) -(0._q,1._q)*GRIDS%TO_FER_RE_CONJG(I,J)*GRIDS%T%NPOINTS*REAL(G_RE(I),q)
            GO( J ) = GO(J ) -(0._q,1._q)*GRIDS%TO_FER_RE(I,J)*GRIDS%T%NPOINTS*REAL(G_IM(I),q)
         ENDDO
         ENDDO

      END SUBROUTINE BACKWARD_FT

   END SUBROUTINE FT_G_TEST

  SUBROUTINE DUMP_NATURALO_ROT_MATRIX( IUNIT, W, NK1, ISP, CHAM_DISTRI, GDES_MAT, IU)
    USE wave
    USE scala
    IMPLICIT NONE
    INTEGER               :: IUNIT          
    TYPE (wavespin)       :: W
    INTEGER               :: NK1               ! current k-point
    INTEGER               :: ISP               ! current spin
    COMPLEX(q)                  :: CHAM_DISTRI(*)    ! distributed matrix
    TYPE (greens_mat_des) :: GDES_MAT       ! in: descriptor
    INTEGER               :: IU                ! unit to write to (not dump for IU<0)
! local
    INTEGER               :: NB_TOT
    INTEGER               :: NDUMP
    INTEGER               :: NBMAX
    COMPLEX(q),ALLOCATABLE   :: CHAM(:,:)
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW
    LOGICAL                 :: LISOPEN

    NBMAX = 0
    DO N2 = 1, W%WDES%ISPIN
       DO N1 = 1, W%WDES%NKPTS 
          NBMAX = MAX( NBMAX, W%WDES%NB_TOTK( N1, N2 ) )
       ENDDO
    ENDDO

    NB_TOT = W%WDES%NB_TOTK( NK1, ISP ) 
    NDUMP = NB_TOT
    ALLOCATE(CHAM(NB_TOT, NB_TOT))

    CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  GDES_MAT%DESC, 1, NB_TOT)
    CALL M_sum_z(GDES_MAT%COMM_INTRA, CHAM(1,1), SIZE(CHAM))

    IF (IU>=0) THEN
       INQUIRE( UNIT = IUNIT , OPENED = LISOPEN )

       IF ( .NOT. LISOPEN ) THEN
          OPEN( UNIT=IUNIT , FILE='NATURALO_ROT', STATUS='REPLACE')
          WRITE(IUNIT,'(A)')'# This file was created by VASP'
          WRITE(IUNIT,'(4I6)')W%WDES%NKPTS, W%WDES%ISPIN,NBMAX, NBMAX 
          DO N1 = 1, W%WDES%NKPTS
             WRITE(IUNIT,'(I6,3F25.16)')N1,W%WDES%VKPT(:,N1)
          ENDDO 
       ENDIF
       WRITE(IUNIT,'(2I6,3F25.16)')1,NB_TOT, W%WDES%VKPT(:,NK1)
       DO N1 = 1, NB_TOT
          DO N2 = 1, NB_TOT
             WRITE(IUNIT,1)N1,N2,CHAM(N1,N2)
          ENDDO
       ENDDO
    ENDIF

    DEALLOCATE(CHAM)
1   FORMAT(2I6,2F20.12)

  END SUBROUTINE DUMP_NATURALO_ROT_MATRIX

!*********************************************************************
!
! FT SIGMA from time to frequency grid
!  SIGMAW_MAT = (SIGMAO_MAT-SIGMAU_MAT) *  FTCOS
!              -(SIGMAO_MAT+SIGMAU_MAT) *  FTSIN * i
!
!*********************************************************************

   SUBROUTINE FT_SIGMA_MAT( W, SIGMAO_MAT, SIGMAU_MAT, SIGMAW_MAT, & 
        NTAU_ROOT, IMAG_GRIDS, GDES_MAT)
     USE scala
     USE constant
     USE GG_base
     IMPLICIT NONE
     TYPE (wavespin)         :: W
     COMPLEX(q)                    :: SIGMAO_MAT(:)     ! matrix <i| Sigma_o(tau) | a> stored distributed
     COMPLEX(q)                    :: SIGMAU_MAT(:)     ! matrix <i| Sigma_u(tau) | a> stored distributed
     COMPLEX(q)              :: SIGMAW_MAT(:,:)   !in/out self-energy in frequency space
     TYPE(imag_grid_handle)  :: IMAG_GRIDS        !time and frequency grid handle
     INTEGER                 :: NTAU_ROOT         !current tau point
     TYPE (greens_mat_des), POINTER :: GDES_MAT   ! in: descriptor
!local variables
     LOGICAL                 :: LDO_TAU_LOCAL     !determines wether the node holds valid data, that needs to be broadcast
     INTEGER                 :: NTAU_GLOBAL
     INTEGER                 :: NTAU_GLOBAL_RCV
     INTEGER                 :: NOMEGA_GLOBAL
     INTEGER                 :: NT
     COMPLEX(q), ALLOCATABLE       :: SIGMAW(:)
     INTEGER                 :: I

     
     LDO_TAU_LOCAL = IMAG_GRIDS%T%LDO_POINT_LOCAL

     ALLOCATE(SIGMAW(SIZE(SIGMAO_MAT,1)))
     
!this is the global tau point of the node
     NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU_ROOT,&
     IMAG_GRIDS )
!if the node does not hold data set NTAU_GLOBAL to -1
     IF (.NOT. LDO_TAU_LOCAL) NTAU_GLOBAL=-1
     
!
! loop over all time points (distributed among nodes)
!
     DO NT=1,GDES_MAT%COMM_INTER%NCPU
        
        NTAU_GLOBAL_RCV=NTAU_GLOBAL
! first broadcast the frequency point
        CALL M_bcast_i_from(GDES_MAT%COMM_INTER, NTAU_GLOBAL_RCV, 1, NT)
        
! if the time point is valid, bcast the data to all nodes
        IF (NTAU_GLOBAL_RCV>0) THEN
           IF (GDES_MAT%COMM_INTER%NODE_ME==NT) SIGMAW=SIGMAO_MAT-SIGMAU_MAT
           CALL M_bcast_z_from(GDES_MAT%COMM_INTER, SIGMAW(1), SIZE(SIGMAW), NT)

! loop over local frequency points
           DO I=1,SIZE(SIGMAW_MAT,2)
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
! possibly beyond last frequency
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE

              SIGMAW_MAT(:,I)=SIGMAW_MAT(:,I)+SIGMAW(:)* &
                   IMAG_GRIDS%TO_FER_RE(NOMEGA_GLOBAL, NTAU_GLOBAL_RCV)
           ENDDO

           IF (GDES_MAT%COMM_INTER%NODE_ME==NT) SIGMAW=SIGMAO_MAT+SIGMAU_MAT
           CALL M_bcast_z_from(GDES_MAT%COMM_INTER,SIGMAW(1), SIZE(SIGMAW), NT)

! loop over local frequency points
           DO I=1,SIZE(SIGMAW_MAT,2)
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
! possibly beyond last frequency point, then cycle
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE

              SIGMAW_MAT(:,I)=SIGMAW_MAT(:,I)-SIGMAW(:)* &
                   (IMAG_GRIDS%TO_FER_RE_CONJG(NOMEGA_GLOBAL, NTAU_GLOBAL_RCV)*(0.0_q,1.0_q))
           ENDDO
        ENDIF
     ENDDO

     DEALLOCATE(SIGMAW)

     
      
    END SUBROUTINE FT_SIGMA_MAT

!*********************************************************************
! inverse of FT_SIGMA_MAT
! FT SIGMA from time to frequency grid
!
! occupied part:                                   -1
! SIGMAO_MAT = 1/2( SIGMAW_MAT + SIGMAW_MAT*)*FTCOS  +
!                                                  -1
!              i/2( SIGMAW_MAT*- SIGMAW_MAT )*FTSIN
!
! unoccupied part:                                 -1
! SIGMAO_MAT =-1/2( SIGMAW_MAT + SIGMAW_MAT*)*FTCOS  +
!                                                  -1
!              i/2( SIGMAW_MAT*- SIGMAW_MAT )*FTSIN
!
!*********************************************************************


   SUBROUTINE INV_FT_SIGMA_MAT( W, NB_TOT, SIGMAO_MAT, SIGMAU_MAT, SIGMAW_MAT, & 
        NTAU_ROOT, IMAG_GRIDS, GDES_MAT)
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)    :: W
     INTEGER            :: NB_TOT
     COMPLEX(q)               :: SIGMAO_MAT(:)     ! matrix <i| Sigma_o(tau) | a> stored distributed
     COMPLEX(q)               :: SIGMAU_MAT(:)     ! matrix <i| Sigma_u(tau) | a> stored distributed
     COMPLEX(q),INTENT(IN) :: SIGMAW_MAT(:,:)   !in/out self-energy in frequency space
     INTEGER            :: NTAU_ROOT         !current tau point
     TYPE(imag_grid_handle)  :: IMAG_GRIDS        !time and frequency grid handle
     TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
!local variables
     LOGICAL            :: LDO_TAU_LOCAL     !determines wether the node holds valid data, that needs to be broadcast
     INTEGER            :: NTAU_GLOBAL
     INTEGER            :: NOMEGA_GLOBAL_RCV
     INTEGER            :: NOMEGA_GLOBAL
     INTEGER            :: NT
     COMPLEX(q), ALLOCATABLE :: SIGMA1(:)
     COMPLEX(q), ALLOCATABLE :: SIGMA2(:)
     REAL(q),ALLOCATABLE,SAVE :: INVFTCOS(:,:)
     REAL(q),ALLOCATABLE,SAVE :: INVFTSIN(:,:)
     INTEGER            :: I
     COMPLEX(q)  :: FT(IMAG_GRIDS%TIME%N,IMAG_GRIDS%TIME%N )
     COMPLEX(q)  :: INV(IMAG_GRIDS%TIME%N,IMAG_GRIDS%TIME%N )
     COMPLEX(q)  :: MAT(IMAG_GRIDS%TIME%N,IMAG_GRIDS%TIME%N )

     
     LDO_TAU_LOCAL = IMAG_GRIDS%T%LDO_POINT_LOCAL 

     ALLOCATE(SIGMA1(SIZE(SIGMAO_MAT,1)))
     ALLOCATE(SIGMA2(SIZE(SIGMAO_MAT,1)))
! distribution over tau and omega points is governed by NTAUPAR

!initialize arrays
     SIGMAO_MAT=(0._q,0._q)
     SIGMAU_MAT=(0._q,0._q)

!in the first run we allocate inverse transformation matrices
     IF ( NTAU_ROOT == 1 .AND. .NOT.(ALLOCATED(INVFTCOS) )) THEN
        IF ( ALLOCATED(INVFTCOS) ) DEALLOCATE(INVFTCOS)
        IF ( ALLOCATED(INVFTSIN) ) DEALLOCATE(INVFTSIN)
        ALLOCATE(INVFTCOS(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA))
        ALLOCATE(INVFTSIN(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA))
        INVFTCOS=IMAG_GRIDS%TO_FER_RE
        INVFTSIN=IMAG_GRIDS%TO_FER_RE_CONJG
        CALL INVERT_REAL_MATRIX(INVFTCOS,-1)
        CALL INVERT_REAL_MATRIX(INVFTSIN,-1)
     ENDIF 

!the global tau point of each node is
     NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU_ROOT,&
     IMAG_GRIDS )
!if global tau point is not valid -> set it to -1
     IF (.NOT. LDO_TAU_LOCAL) NTAU_GLOBAL=-1

! loop over local frequency points ( in group )
     DO I=1,SIZE(SIGMAW_MAT,2)
!this is the global frequency point of the node (same distribution like tau)
        NOMEGA_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
        IMAG_GRIDS )
!if the node does not hold data set NOMEGA_GLOBAL to -1
        IF (NOMEGA_GLOBAL > IMAG_GRIDS%NOMEGA ) NOMEGA_GLOBAL=-1
!
! loop over all frequency points ( between groups )
!
        DO NT=1,GDES_MAT%COMM_INTER%NCPU
!this decides who gets data
           NOMEGA_GLOBAL_RCV=NOMEGA_GLOBAL
! first broadcast the global frequency point
           CALL M_bcast_i_from(GDES_MAT%COMM_INTER, NOMEGA_GLOBAL_RCV, 1, NT)
! only nodes who have NOMEGA_GLOBAL_RCV > 0 bcast
           IF ( NOMEGA_GLOBAL_RCV > 0 ) THEN
!even part: (cos transform of Sigmaw+Sigmaw*)
! if the frequency point is valid, bcast the data to all groups
              SIGMA1=0
              SIGMA2=0
              IF (GDES_MAT%COMM_INTER%NODE_ME==NT) THEN
                 SIGMA1(:)=SIGMAW_MAT(:,I)
! add the conjugated part (making the matrix Hermitian)
                 CALL PZGEADD( "C", NB_TOT, NB_TOT, &
                      (1.0_q,0.0_q), SIGMA1(1), 1, 1, GDES_MAT%DESC,  &
                      (0.0_q,0.0_q), SIGMA2(1), 1, 1, GDES_MAT%DESC)
                 SIGMA1(:)=SIGMA1(:)+SIGMA2(:)
              ENDIF 
!broadcast to all groups
              CALL M_bcast_z_from(GDES_MAT%COMM_INTER, SIGMA1(1), SIZE(SIGMA1), NT)
               
!and store this piece
              IF (NOMEGA_GLOBAL_RCV>0 .AND. NTAU_GLOBAL>0) THEN
                 SIGMAO_MAT(:)=SIGMAO_MAT(:)+SIGMA1(:)*&
                      INVFTCOS(NTAU_GLOBAL, NOMEGA_GLOBAL_RCV)*0.25_q
                 SIGMAU_MAT(:)=SIGMAU_MAT(:)-SIGMA1(:)*&
                      INVFTCOS(NTAU_GLOBAL, NOMEGA_GLOBAL_RCV)*0.25_q
              ENDIF

! uneven part: (sin transform of Sigmaw-Sigmaw*)
! if the frequency point is valid, bcast the data to all groups
              IF (GDES_MAT%COMM_INTER%NODE_ME==NT) THEN
                 SIGMA1(:)=SIGMAW_MAT(:,I)
! add the conjugated part (making the matrix Hermitian)
                 CALL PZGEADD( "C", NB_TOT, NB_TOT, &
                      (1.0_q,0.0_q), SIGMA1(1), 1, 1, GDES_MAT%DESC,  &
                      (0.0_q,0.0_q), SIGMA2(1), 1, 1, GDES_MAT%DESC)
                 SIGMA1(:)=SIGMA1(:)-SIGMA2(:)
              ENDIF 
!broadcast to all groups
              CALL M_bcast_z_from(GDES_MAT%COMM_INTER, SIGMA1(1), SIZE(SIGMA1), NT)

              IF (NOMEGA_GLOBAL_RCV>0 .AND. NTAU_GLOBAL>0) THEN
                 SIGMAO_MAT(:)=SIGMAO_MAT(:)-SIGMA1(:)* &
                    (INVFTSIN(NTAU_GLOBAL, NOMEGA_GLOBAL_RCV)*(0.0_q,-0.25_q))
                 SIGMAU_MAT(:)=SIGMAU_MAT(:)-SIGMA1(:)* &
                    (INVFTSIN(NTAU_GLOBAL, NOMEGA_GLOBAL_RCV)*(0.0_q,-0.25_q))
              ENDIF
           ENDIF !decide if sending is 1._q
        ENDDO ! groups
     ENDDO !local frequency points

     DEALLOCATE(SIGMA1)
     DEALLOCATE(SIGMA2)

     
      
    END SUBROUTINE INV_FT_SIGMA_MAT

!***********************************************************************
!
! rotate complex Green function or selfenergy according to a given unitary
! transformation matrix U=CHAM_MAT
!    U+ G U
!
!***********************************************************************

  SUBROUTINE ROTATE_GREEN_SIGMA(WDES, CHAM_MAT, G_MAT, IMAG_GRIDS, GDES_MAT)
    IMPLICIT NONE
    TYPE (wavedes)           :: WDES
    COMPLEX(q), POINTER            :: CHAM_MAT(:,:,:,:) ! unitary transformation
    COMPLEX(q), POINTER      :: G_MAT(:,:,:,:)    ! frequency dependent Green function or selfenergy
    TYPE (imag_grid_handle ) :: IMAG_GRIDS        ! imaginary grid descriptor
    TYPE (greens_mat_des), POINTER :: GDES_MAT    ! greens function matrix descriptor
! local
    INTEGER ISP, NK, I, N
    COMPLEX(q), ALLOCATABLE  :: G1(:),CHAM(:)

    ALLOCATE(G1(SIZE(CHAM_MAT,1)), CHAM(SIZE(CHAM_MAT,1)))

    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS
          N= WDES%NB_TOTK(NK, ISP)
! copy CHAM to complex work array
          CHAM=CHAM_MAT(:, NK, ISP, 1)

          DO I=1,SIZE(G_MAT,4)
!this is the global tau point of the node
             IF ( DETERMINE_NTAU_GLOBAL( WDES%COMM_KIN%NODE_ME, I, IMAG_GRIDS)> IMAG_GRIDS%NOMEGA) CYCLE

! G1= U+ G
             CALL PZGEMM( 'C', 'N', N, N, N, (1.0_q,0.0_q),  &
                  CHAM(1), 1, 1, GDES_MAT%DESC, &
                  G_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC, &
                  (0.0_q,0.0_q),  &
                  G1(1), 1, 1, GDES_MAT%DESC )
! G = G1 U
             CALL PZGEMM( 'N', 'N', N, N, N, (1.0_q,0.0_q),  &
                  G1(1), 1, 1, GDES_MAT%DESC, &
                  CHAM(1), 1, 1, GDES_MAT%DESC, &
                  (0.0_q,0.0_q),  &
                  G_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC )
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(G1, CHAM)
  END SUBROUTINE ROTATE_GREEN_SIGMA


!*********************************************************************
!
! this routine determines the change of the density matrix gamma.
! Specifically, the integral
!
!  delta gamma = 1/ pi int_0^infty G_(iw) dw
!
! it is important to subtract the non-interacting part, since
! the frequency integration is otherwise not converging
! this is either the non-interacting DFT Green function
! or the Hartree-Fock Green function
!
!*********************************************************************

! explicit calculation of psi term


   SUBROUTINE GAMMA_FROM_G_W( WMEAN, NKPTS_IRZ, GW, CORR_MAT, & 
        IMAG_GRIDS, GDES_MAT, CHAM_MAT, LUNOCCUPIED )
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)         :: WMEAN
     INTEGER                 :: NKPTS_IRZ
     COMPLEX(q)              :: GW(:,:,:,:)       ! interacting Green function G
     COMPLEX(q)                    :: CORR_MAT(:,:,:,:) ! contribution to density matrix from correlation effects
     TYPE (imag_grid_handle) :: IMAG_GRIDS        ! imaginary grid descriptor
     TYPE (greens_mat_des)   :: GDES_MAT          ! greens function matrix descriptor
     COMPLEX(q)                    :: CHAM_MAT(:,:,:,:) ! Hartree-Fock basis vectors
     LOGICAL, OPTIONAL       :: LUNOCCUPIED       ! obtain holes
! local variables
     INTEGER                 :: I,J                ! frequency points
     INTEGER                 :: NOMEGA_GLOBAL      ! corresponing global frequency point
     INTEGER                 :: NK1                !current q-point
     INTEGER                 :: ISP                !current spin
     COMPLEX(q), ALLOCATABLE       :: GTMP(:)
     COMPLEX(q), ALLOCATABLE       :: GPHI(:)
     COMPLEX(q), ALLOCATABLE       :: GPSI(:)
     LOGICAL                 :: LUNOCCUPIED_

     REAL(q)                 :: THETA( SIZE(WMEAN%CELTOT,1) )
     INTEGER                 :: N


     
!
! Feynman propagator has this form in imaginary time:
!     G(t) = -( G_psi( t ) + sign( t ) * G_phi( t ) )
! where
!                              1  sinh( beta(x-mu)/2( 1-2*t))
!   G_psi(t) = sum_{x in poles}-*---------------------------
!                              2    cosh(beta*(x-mu)/2)
!
!                              1  cosh( beta(x-mu)/2( 1-2*t))
!   G_phi(t) = sum_{x in poles}-*---------------------------
!                              2    cosh(beta*(x-mu)/2)
!
! G_psi(0) and G_phi(0) is set up from G(iw) as follows
!
! G_phi(0) = sum_{n=1}^NOMEGA weight_n * ( G(iw_n) + G(-iw_n) )/2
! G_psi(0) = sum_{n=1}^NOMEGA weight_n * ( G(iw_n) - G(-iw_n) )/2
!

     LUNOCCUPIED_ = .FALSE.
     IF ( LUNOCCUPIED ) THEN
        LUNOCCUPIED_ = LUNOCCUPIED
     ENDIF

     ALLOCATE(GTMP(SIZE(GW,1)))
     GTMP = 0
     ALLOCATE(GPHI(SIZE(GW,1)))
     ALLOCATE(GPSI(SIZE(GW,1)))
     
     CORR_MAT=0

     DO ISP=1,WMEAN%WDES%ISPIN
!test
        DO NK1=1,NKPTS_IRZ
           GPHI = 0
           GPSI = 0
           DO I=1,SIZE(GW,4)  ! loop over local frequency points
!this is the global tau point of the node
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( WMEAN%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE

! transform of even part of density matrix
              GPHI(:) = GPHI(:) + GW(:, NK1, ISP, I)*&
                 IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2

# 648

           ENDDO
           
            
! add the conjugated part (making the matrix Hermitian)
           CALL PZGEADD( 'C', WMEAN%WDES%NB_TOTK(NK1, ISP), WMEAN%WDES%NB_TOTK(NK1, ISP), &
                (1._q,0._q), GPHI(1), 1, 1, GDES_MAT%DESC,  & 
                (0._q,0._q), GTMP(1), 1, 1, GDES_MAT%DESC )
           GPHI( : ) = GPHI( : ) + GTMP( : )
!           CALL DUMP_GREENS_HAM( "G_PHI_0", GPHI(:), WMEAN%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6,IUNIT=101)

           IF ( LFINITE_TEMPERATURE ) THEN
# 674

           THETA=0
           N =  WMEAN%WDES%NB_TOTK(NK1, ISP)
! collective sum between groups is 1._q after this routine
! so prevent double summation
           THETA(1:N)=-0.5_q/IMAG_GRIDS%T%COMM_BETWEEN_GROUPS%NCPU
! determine unitary transformation matrix U
           GTMP(:) =CHAM_MAT(:, NK1, ISP, 1 )
           CALL RIGHT_MULTIPLY_BY_R(N, GTMP, THETA, GDES_MAT%DESC )
           
! calculate U F U+
           CALL PZGEMM( 'N', 'C', N, N, N, (1._q,0._q),  &
                GTMP(1), 1, 1, GDES_MAT%DESC, &
                CHAM_MAT(1,NK1,ISP,1), 1, 1, GDES_MAT%DESC, &
                (0._q,0._q),  &
                GPSI(1), 1, 1, GDES_MAT%DESC )

!           CALL DUMP_GREENS_HAM( "G_PSI_0", GPSI(:), WMEAN%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6,IUNIT=102)

           ENDIF 
!
! Feynman propagator has this form in imaginary time:
!     G(t) = -( G_psi( t ) + sign( t ) * G_phi( t ) )
!
! so the density matrix for electrons (t -> 0-) is
!
! GAM_NEG = G(0-) = G_psi(0) - G_phi( 0 )
!
! and for holes (t->0+) is
!
! GAM_POS = G(0+) = G_psi(0) + G_phi( 0 )
!

           IF ( LUNOCCUPIED_ ) THEN
              CORR_MAT(:, NK1, ISP, 1) =-(GPSI(:) + GPHI( : ))
           ELSE
              CORR_MAT(:, NK1, ISP, 1) =-(GPSI(:) - GPHI( : ))
           ENDIF 

!           CALL DUMP_GREENS_HAM( "GAMMA", CORR_MAT(:,NK1,ISP,1), WMEAN%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6,IUNIT=100)

        ENDDO
     ENDDO

     DEALLOCATE(GTMP)
     DEALLOCATE(GPHI)
     DEALLOCATE(GPSI)
     


     CONTAINS

     SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
! Rational interpolation from NR
        INTEGER, INTENT(IN)    :: N
        REAL(q), INTENT(IN)    :: XA(N)
        REAL(q), INTENT(IN)    :: YA(N)
        REAL(q), INTENT(IN)    :: X
        REAL(q), INTENT(INOUT) :: Y
        REAL(q), INTENT(INOUT) :: DY
! local
        INTEGER, PARAMETER :: NMAX = 32
        REAL(q), PARAMETER :: TINY = 1.E-25_q
!largest expected value NMAX of N, and a small number TINY.
!given arrays XA and YA , each of length N , and given a value of X , this routine returns a
!value of Y and an accuracy estimate DY. The value returned is that of the diagonal rational
!function, evaluated at X , which passes through the N points ( XA_i , YA_i ), i = 1... n .
        INTEGER I,M,NS
        REAL(q) DD,H,HH,T,W,C(NMAX),D(NMAX)

        NS=1
        HH=ABS(X-XA(1))
        DO I=1,N
           H=ABS(X-XA(I))
           IF ( H == 0 )THEN
              Y=YA(I)
              DY=0
              RETURN
           ELSE IF (H < HH) THEN
              NS=I
              HH=H
           ENDIF
           C(I)=YA(I)
           D(I)=YA(I)+TINY
!The TINY part is needed to prevent a rare (0._q,0._q)-over-
!(0._q,0._q) condition.
        ENDDO

        Y=YA(NS)
        NS=NS-1
        DO M=1,N-1
           DO I=1,N-M
              W=C(I+1)-D(I)
              H=XA(I+M)-X
!h will never be (0._q,0._q), since this was tested in the ini-
!tializing loop.
              T=(XA(I)-X)*D(I)/H
              DD=T-C(I+1)
! jF.  this should be never reached if the XA's are all different since then
!      (XA(I)-X)/H=(XA(I)-X)/(XA(I+M)-X) is clearly different from 1. and D(I)
!      which is C(I)+TINY multiplied with this factor different 1 (giving T)
!      should then be sufficiently different from C(I) to avoid DD=0.
!              IF(DD == 0)PAUSE 'FAILURE IN RATINT'
! jF.
!This error condition indicates that the interpolating function has a pole at the re-
!quested value of x.
              DD=W/DD
              D(I)=C(I+1)*DD
              C(I)=T*DD
           ENDDO

           IF (2*NS < N-M)THEN
           DY=C(NS+1)
           ELSE
           DY=D(NS)
           NS=NS-1
           ENDIF
           Y=Y+DY
        ENDDO 
     END SUBROUTINE RATINT

   END SUBROUTINE GAMMA_FROM_G_W


!*********************************************************************
!
! recovers Feynman propagator  at tau=0+ and tau=0-
!
!*********************************************************************
   SUBROUTINE GAM_POS_GAM_NEG_FROM_G_W( W, NKPTS_IRZ, GW, GAM_POS, GAM_NEG,& 
        IMAG_GRIDS, GDES_MAT )
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)         :: W
     INTEGER                 :: NKPTS_IRZ
     COMPLEX(q)              :: GW(:,:,:,:)      ! interacting Green function G
     COMPLEX(q)                    :: GAM_POS(:,:,:,:) ! density matrix for holes
     COMPLEX(q)                    :: GAM_NEG(:,:,:,:) ! density matrix for eletrons
     TYPE (greens_mat_des)   :: GDES_MAT         ! greens function matrix descriptor
     TYPE (imag_grid_handle) :: IMAG_GRIDS       ! imaginary grid descriptor
! local variables
     INTEGER            :: I                  ! frequency points
     INTEGER            :: NOMEGA_GLOBAL      ! corresponing global frequency point
     INTEGER            :: NK1                !current q-point
     INTEGER            :: ISP                !current spin
     INTEGER            :: N
     COMPLEX(q), ALLOCATABLE  :: GTMP(:)
     COMPLEX(q), ALLOCATABLE  :: GPHI(:)
     COMPLEX(q), ALLOCATABLE  :: GPSI(:)
     
!WRITE(*,*)' RE_WEIGHT'
!WRITE(*,'(8E16.9)')IMAG_GRIDS%FER_RE_WEIGHT(:)
!
!WRITE(*,*)' TO_TAU0_CONJG/2'
!WRITE(*,'(8E16.9)')IMAG_GRIDS%TO_TAU0_CONJG(1,:)/2
!
!WRITE(*,*)' TO_TAU0/2'
!WRITE(*,'(8E16.9)')IMAG_GRIDS%TO_TAU0(1,:)/2

     

     ALLOCATE(GTMP(SIZE(GW,1)))
     ALLOCATE(GPHI(SIZE(GW,1)))
     ALLOCATE(GPSI(SIZE(GW,1)))
     
     GAM_POS=0
     GAM_NEG=0
     
     DO ISP=1,W%WDES%ISPIN
!test
        DO NK1=1,NKPTS_IRZ
           GPHI = 0
           GPSI = 0
           N=  W%WDES%NB_TOTK(NK1, ISP)
           DO I=1,SIZE(GW,4)  ! loop over local frequency points
!this is the global tau point of the node
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE

! transform of even part of density matrix
              GPHI(:) = GPHI(:) + GW(:, NK1, ISP, I)*&
                 IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2

! transform odd part of density matrix
! if transform is available
              IF ( ASSOCIATED( IMAG_GRIDS%TO_TAU0 ) ) THEN

                 GPSI(:) = GPSI(:) + GW(:, NK1, ISP, I)*&
                    IMAG_GRIDS%TO_TAU0(1, NOMEGA_GLOBAL)/4
# 873

              ENDIF 
           ENDDO

! add the conjugated part (making the matrix Hermitian)
           CALL PZGEADD( 'C', W%WDES%NB_TOTK(NK1, ISP), W%WDES%NB_TOTK(NK1, ISP), &
                (1._q,0._q), GPHI(1), 1, 1, GDES_MAT%DESC,  & 
                (0._q,0._q), GTMP(1), 1, 1, GDES_MAT%DESC )
           GPHI( : ) = GPHI( : ) + GTMP( : )
!           CALL DUMP_GREENS_HAM( "G_PHI_0", GPHI(:), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6,IUNIT=101)

           IF ( ASSOCIATED( IMAG_GRIDS%TO_TAU0 ) ) THEN
! subtract conjugated part and divide by I
! GTMP <-- GPSI*
              CALL PZGEADD( 'C', W%WDES%NB_TOTK(NK1, ISP), W%WDES%NB_TOTK(NK1, ISP), &
                   (1._q,0._q), GPSI(1), 1, 1, GDES_MAT%DESC,  & 
                   (0._q,0._q), GTMP(1), 1, 1, GDES_MAT%DESC )
! subtract from GPSI: GPSI - GPSI*
              GPSI( : ) = (GPSI( : ) - GTMP( : ))*(0._q,-1._q)
!CALL DUMP_GREENS_HAM( "G_PSI_0", GPSI(:), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6,IUNIT=102)
           ENDIF

!
! Feynman propagator has this form in imaginary time:
!     G(t) = -( G_psi( t ) + sign( t ) * G_phi( t ) )
! where
!
!                              1  sinh( beta(x-mu)/2( 1-2*t))
!   G_psi(t) = sum_{x in poles}-*---------------------------
!                              2   cosh(beta*(x-mu)/2)
!
!                              1   cosh( beta(x-mu)/2( 1-2*t))
!   G_phi(t) = sum_{x in poles}-*----------------------------
!                              2     cosh(beta*(x-mu)/2)
!
!
! so the density matrix for electrons (t -> 0-) is
!
! GAM_NEG = G(0-) = G_psi(0) - G_phi( 0 )
!
! and for holes (t->0+) is
!
! GAM_POS = G(0+) = G_psi(0) + G_phi( 0 )
!

           GAM_POS(:, NK1, ISP, 1) =-(GPSI(:) + GPHI( : ))
           GAM_NEG(:, NK1, ISP, 1) =-(GPSI(:) - GPHI( : ))

!           CALL DUMP_GREENS_HAM( "GAMMA", GAM_POS(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6,IUNIT=100)

        ENDDO
     ENDDO

     DEALLOCATE(GTMP)
     DEALLOCATE(GPHI)

      
! at this point SIGMAW_MAT and CORR_MAT contain the interacting G-G_mean
! and the corresponding density matrix
     CALL M_sum_z(GDES_MAT%COMM_INTER, GAM_POS, SIZE(GAM_POS))
     CALL M_sum_z(GDES_MAT%COMM_INTER, GAM_NEG, SIZE(GAM_NEG))
     
     
   END SUBROUTINE GAM_POS_GAM_NEG_FROM_G_W

!*********************************************************************
!
! same routine as above but this version is designed to
! work e.g. only on the diagonals of the Green's function
!
! because we use 1 above, where the row and column indices
! are merged into (1._q,0._q) common index, differences to previous
! routine are small (except for addition of complex conjugated)
!
!*********************************************************************

   SUBROUTINE GAMMA_FROM_G_W_DIAG( W, NKPTS_IRZ, GW, GAM_NEG, IMAG_GRIDS )
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)         :: W
     INTEGER                 :: NKPTS_IRZ
     COMPLEX(q)              :: GW(:,:,:,:)    !diagonal part of interacting Green function G
     COMPLEX(q)                    :: GAM_NEG(:,:,:) !diagonal part of correlation contr. to density matrix
! for each state, k-point and spin index
     TYPE( imag_grid_handle) :: IMAG_GRIDS     ! imaginary grids

! local variables
     INTEGER            :: I                  ! frequency points
     INTEGER            :: NOMEGA_GLOBAL      ! corresponing global frequency point
     INTEGER            :: NK1                !current q-point
     INTEGER            :: ISP                !current spin
     COMPLEX(q), ALLOCATABLE  :: GPHI(:)
     COMPLEX(q), ALLOCATABLE  :: GPSI(:)

     
     
! initialize negative time gamma matrix (electrons)
     GAM_NEG=0

     ALLOCATE(GPHI(SIZE(GW,1)))
     ALLOCATE(GPSI(SIZE(GW,1)))
     
     DO ISP=1,W%WDES%ISPIN
        DO NK1=1,NKPTS_IRZ
           GPHI = 0 
           GPSI = 0 

           DO I=1,SIZE(GW,4)  ! loop over local frequency points
              
!this is the global tau point of the node
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE

! add to density matrix
!GAM_NEG(:, NK1, ISP)=GAM_NEG(:, NK1, ISP)+GW(:, NK1, ISP,I)*&
!IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2

! transform of even part of density matrix
              GPHI(:) = GPHI(:) + GW(:, NK1, ISP, I)*&
                 IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2

!#undef 
# 1015

           ENDDO
! compute psi term directly from CHAM
! The PSI term is the exchange part
            IF ( LFINITE_TEMPERATURE ) THEN
               IF ( LHOLEGF ) THEN
                  GPSI(:) = +0.5_q/IMAG_GRIDS%T%COMM_BETWEEN_GROUPS%NCPU
               ELSE
                  GPSI(:) = -0.5_q/IMAG_GRIDS%T%COMM_BETWEEN_GROUPS%NCPU
               ENDIF             
            ENDIF


           GPHI( : ) = GPHI( : ) + CONJG(GPHI( : ))
!
! Feynman propagator has this form in imaginary time:
!     G(t) = -( G_psi( t ) + sign( t ) * G_phi( t ) )
! where
!
!                              1  sinh( beta(x-mu)/2( 1-2*t))
!   G_psi(t) = sum_{x in poles}-*---------------------------
!                              2   cosh(beta*(x-mu)/2)
!
!                              1   cosh( beta(x-mu)/2( 1-2*t))
!   G_phi(t) = sum_{x in poles}-*----------------------------
!                              2     cosh(beta*(x-mu)/2)
!
!
! so the density matrix for electrons (t -> 0-) is
!
! GAM_NEG = G(0-) = G_psi(0) - G_phi( 0 )
!
! and for holes (t->0+) is
!
! GAM_POS = G(0+) = G_psi(0) + G_phi( 0 )
!

! electrons
           GAM_NEG(:, NK1, ISP)= GPHI( : ) 
              
! holes
!GAM_POS(:, NK1, ISP)=-GPHI( : )

! in case of metallic compound this term is present as well
           IF ( ASSOCIATED( IMAG_GRIDS%TO_TAU0 ) ) THEN

! electrons
              GAM_NEG(:, NK1, ISP ) = GAM_NEG( :, NK1, ISP) - GPSI( : ) 

! holes
!GAM_POS(:, NK1, ISP ) = GAM_POS( :, NK1, ISP) - GPSI( : )

           ENDIF
        ENDDO
     ENDDO

     DEALLOCATE( GPHI )
     DEALLOCATE( GPSI )
     
   END SUBROUTINE GAMMA_FROM_G_W_DIAG

!*********************************************************************
!
! this routine determines the interacting Green function
! from the full selfenergy for all spin, k-points and
! locally stored frequency points
!
!
!  G(w)= (w - T+ V_H+ V_x + Sigma_c(w))^-1
!
! is calculated e.g. the first order change of the Green function
! in this case, CHAM_MAT is irrelevant
!
! this routine also determines the Galitskii-Migdal Energy:
! E_GM = int dw Tr( G(iw) Sigma_c( iw ))
!
!*********************************************************************

   SUBROUTINE CALCULATE_G_FROM_SIGMA( W, NKPTS_IRZ, EFERMI, SIGMAW_MAT, & 
        IMAG_GRIDS, GDES_MAT, CHAM_MAT, E, IO)
     USE base
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)    :: W!MEAN
     INTEGER            :: NKPTS_IRZ
     REAL(q)            :: EFERMI(:)          ! in: fermi-level
     COMPLEX(q)         :: SIGMAW_MAT(:,:,:,:)! in: matrix <i| Sigma(omega) | a> stored distributed
! out: interacting Green function
     COMPLEX(q)               :: CHAM_MAT(:,:,:,:)  ! in: frequency independent part of self-energy (T +V_x + V_H)
     TYPE( energy )     :: E 
     TYPE( in_struct )  :: IO                 ! output to OUTCAR
     TYPE (imag_grid_handle ) :: IMAG_GRIDS   ! imaginary grids
     TYPE (greens_mat_des) :: GDES_MAT        ! descriptor for distributed matrix
! local variables
     INTEGER                 :: I                  ! frequency points
     INTEGER                 :: N                  ! bands for each k-points
     INTEGER                 :: NOMEGA_GLOBAL      ! corresponing global frequency point
     INTEGER                 :: NK1                !current q-point
     INTEGER                 :: ISP                !current spin
     COMPLEX(q)              :: AUX( SIZE(W%CELTOT,1))
     COMPLEX(q), ALLOCATABLE :: G1(:), G2(:), CHAM(:)
     INTEGER                 :: INFO
     INTEGER                 :: IPIV(W%WDES%NB_TOT)
     INTEGER                 :: LWORK, LIWORK
     INTEGER, ALLOCATABLE    :: IWORK(:)
     COMPLEX(q),ALLOCATABLE  :: WORK(:)
     INTEGER                 :: NB                 ! band index
     COMPLEX(q),ALLOCATABLE  :: SIGMA(:)! auxillary
     COMPLEX(q),ALLOCATABLE  :: GS(:)! auxillary
     COMPLEX(q),ALLOCATABLE  :: GX(:)! auxillary
     COMPLEX(q),ALLOCATABLE  :: GSSG(:)! auxillary
     COMPLEX(q)              :: E_GM
     COMPLEX(q)              :: E_KLGG0,E_KLLOG ! Klein terms
     REAL(q)                 :: E_TR
     REAL(q),ALLOCATABLE     :: EIGGS(:)! auxillary
     COMPLEX(q),ALLOCATABLE  :: CDIAG(:)! auxillary
! functional type
     LOGICAL  :: LLEEUWEN=.TRUE.    ! see PRA73,012511(2006)  Eq. 3
     LOGICAL  :: LWELDEN=.FALSE.    ! see JCP 145, 204106 (2016) Eq. 18
      

     ALLOCATE(G1(SIZE(SIGMAW_MAT,1)), G2(SIZE(SIGMAW_MAT,1)), CHAM(SIZE(SIGMAW_MAT,1)))

! querry work-space and allocate optimal workspace
     ALLOCATE(WORK(1), IWORK(1))
     LWORK=-1
     LIWORK=-1
     CALL PZGETRI( W%WDES%NB_TOT, G1(1), 1, 1, GDES_MAT%DESC, IPIV, WORK, LWORK, &
          IWORK, LIWORK, INFO )
     LWORK=MAX(INT(WORK(1)),1)
     LIWORK=MAX(IWORK(1),1)
     DEALLOCATE(WORK, IWORK)
     ALLOCATE(WORK(LWORK), IWORK(LIWORK))

! for Galitskii-Migdal, two auxillary arrays are required
! (1._q,0._q) to store the self-energy and (1._q,0._q) to store G*Sigma
! this has to be 1._q only at (1._q,0._q) intermediate frequency point
     ALLOCATE(SIGMA(SIZE(SIGMAW_MAT,1)))
     ALLOCATE(GS(SIZE(SIGMAW_MAT,1)))
     ALLOCATE(GX(SIZE(SIGMAW_MAT,1)))
     ALLOCATE(GSSG(SIZE(SIGMAW_MAT,1)))
     SIGMA=0
     GS=0
     GX=0
     GSSG=0
     E%ELOG1G=0
     E%EKLLOG = 0
     E%EKLGG0 = 0
     DO ISP=1,W%WDES%ISPIN
        DO NK1=1,NKPTS_IRZ
           N=  W%WDES%NB_TOTK(NK1, ISP)
           ALLOCATE(EIGGS(N))
           ALLOCATE(CDIAG(N))
           EIGGS=0
           CDIAG=0

! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
! set CHAM to  Hartree-Fock Hamiltonian minus Fermi-level
!#define diag_cham
# 1190

! subtract Fermi-level from diagonal
           CHAM=0
           AUX(1:N)=W%CELTOT(1:N,NK1,ISP)-EFERMI(ISP)
           CALL ADD_TO_DIAGONALE(N, CHAM, AUX, GDES_MAT%DESC )

!              CALL DUMP_GREENS_HAM( "mean HF", CHAM(:), N, GDES_MAT, 6)

! non-interacting grand potential is -Tr Log ( G_0^-1 )
! see Luttinger, Ward PhzRev118,1417 Eq.54
! however, we determine Helmholtz free energy from grand potential
! by adding +\mu * N to \Omega, so there we only evaluate
! -Tr Log ( G_0^-1 ) -> \sum_{i=1..occ} (e^i_x - \mu) at T->0
           E%ELOGG0 = E%ELOGG0 + TR_LOG_G0_TERM ( W, AUX, ISP, NK1, N )

           SIGMA(:)=0
           E_TR=0 
           E_KLGG0=0 
           E_KLLOG=0 
           E_GM=0
           DO I=1,SIZE(SIGMAW_MAT,4)  ! loop over local frequency points

!this is the global tau point of the node
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE


! optionally, use only exchange part of self-energy
              IF ( LGREENHF ) SIGMAW_MAT(:, NK1, ISP, I) = 0 
! store Self-energy (correlation part only)
              SIGMA(:)=SIGMAW_MAT(:, NK1, ISP, I)

!              CALL LUTTINGER_WARD_POTENTIAL_TERMS( NK1, ISP, I, N, E_TR, E_GM, E_KLGG0,E_KLLOG )

! add mean field Hamiltonian
              SIGMAW_MAT(:, NK1, ISP, I)=SIGMAW_MAT(:, NK1, ISP, I)+CHAM

              SIGMAW_MAT(:, NK1, ISP, I)=-SIGMAW_MAT(:, NK1, ISP, I)
!  calculate i w +e_fermi - Sigma_c - (T +V_x + V_H)
              AUX=CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q)
              CALL ADD_TO_DIAGONALE(N, SIGMAW_MAT(:, NK1, ISP, I), AUX, GDES_MAT%DESC)

!
! SIGAMW_MAT contains now ( G_HF ^ -1 - \Sigma_c ) = G^-1
!
!CALL DUMP_GREENS_HAM_GDEF( "SIGMAW_MAT", SIGMA(:), N, GDES_MAT, IO%IU0, 100)
               
! invert selfenergy
! first LU decomposition
              CALL PZGETRF( N, N, & 
                   SIGMAW_MAT(1, NK1, ISP, I), 1, 1, GDES_MAT%DESC, IPIV, INFO )
              IF (INFO/=0) THEN
                 CALL vtutor%bug("PZGETRF in CALCULATE_G_FROM_SIGMA, INFO= " // str(INFO), "greens_orbital.F", 1243)
              ENDIF
! then inversion
              CALL PZGETRI( N, SIGMAW_MAT(1, NK1, ISP, I), 1, 1, GDES_MAT%DESC, IPIV, WORK, LWORK, &
                   IWORK, LIWORK, INFO )
              IF (INFO/=0) THEN
                 CALL vtutor%bug("PZGETRI in CALCULATE_G_FROM_SIGMA, INFO= " // str(INFO), "greens_orbital.F", 1249)
              ENDIF

! evaluate Grand potential terms
! compute Tr { G*S_c}  and Tr ln( 1 - GX * S_c )
! computes also Klein functional terms Tr{Ln( G*G_0^-1)}, Tr{ G*G_0^-1 - 1 }
              CALL LUTTINGER_WARD_POTENTIAL_TERMS( NK1, ISP, I, N, E_TR, E_GM, E_KLGG0,E_KLLOG ) 
           ENDDO
           CALL M_sum_d(GDES_MAT%COMM_INTER, E_TR, 1)
           CALL M_sum_z(GDES_MAT%COMM_INTER, E_GM, 1)
           CALL M_sum_z(GDES_MAT%COMM_INTER, E_KLGG0, 1)
           CALL M_sum_z(GDES_MAT%COMM_INTER, E_KLLOG, 1)

! accumulate Tr Log( 1- G_x S_c )
           E%ELOG1G = E%ELOG1G + E_TR*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

! accumulate Galitskii-Migdal correlation eenergy
           E%ECGWGM = E%ECGWGM + REAL(E_GM,q)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN
           E%EKLGG0 = E%EKLGG0 + REAL(E_KLGG0,q)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN
           E%EKLLOG = E%EKLLOG + REAL(E_KLLOG,q)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN
           DEALLOCATE(EIGGS)
           DEALLOCATE(CDIAG)

        ENDDO
     ENDDO

     DEALLOCATE(WORK, IWORK)
     DEALLOCATE(G1, G2, CHAM)
     DEALLOCATE(GS, GX, SIGMA, GSSG)

      
     CONTAINS

!=======================================================================
! calculate Tr{ ln( -G_0^-1 ) }
!=======================================================================
     FUNCTION TR_LOG_G0_TERM( W, AUX, ISP, NK1, N )
        REAL(q)                :: TR_LOG_G0_TERM
        TYPE (wavespin)        :: W
        INTEGER, INTENT(IN)    :: ISP, NK1, N
        COMPLEX(q), INTENT(IN) :: AUX( : ) ! e_r - \mu
! local
        REAL(q)                :: TR
        INTEGER                :: J
        REAL(q)                :: F, DF, X
        REAL(q),PARAMETER      :: MAXEXP=35._q

        TR= 0
        DO J = 1 , N 
           X = IMAG_GRIDS%BETA*REAL(AUX(J),q)
!
! free propagator is
!
!   G_0(i\omega_l) = \sum_r 1/( i\omega_l  - (e_r - \mu ) )
!
! Luttinger-Ward: PysRev.118.1417.(1960) Eq. 54
!
! -Tr ln( -G_0^-1 ) = -1/\beta \sum_l exp(-i\omega_l 0-) x
!                             \sum_r ln( e_r - \mu - i\omega_l )
!                   =-1/\beta \sum_r ln{ 1+e^[-\beta(e_r -\mu)] }
!                   =1/\beta \sum_r ln{ 1 - f_r }
!
! - ln(1+exp(-x)) =-(  -x/2 + ln(2) + ln(cosh(x/2)) )
!
!
! strongly occupied  ln( 1 + e^(x>30) ) ~ ln( e^(x>30) ) = x
!IF (X < -MAXEXP) THEN
!   TR = TR + AUX(J)
!ELSE IF ( -MAXEXP < X .AND. X < MAXEXP ) THEN
!   TR = TR - LOG( 1 + EXP(-X) )/IMAG_GRIDS%BETA
!END IF
! strongly occupied
           IF (X < -MAXEXP) THEN
              TR = TR + AUX(J)
           ELSE IF ( -MAXEXP < X .AND. X < MAXEXP ) THEN
              TR = TR + AUX(J)/(1+EXP(X))
           END IF
!
! entropy term is actually
!    1/b [ -f ln(f) - (1-f) ln(1-f) ]           =
!    1/b [ -f ln( f/(1-f) ) - ln(1-f) ]         =
!    1/b [ -f ln( e^-beta{e_r-mu} ) - ln(1-f) ] =
!    f {e_r-mu} - 1/b ln(1-f)
!
        ENDDO
        TR_LOG_G0_TERM = TR*W%WDES%RSPIN*W%WDES%WTKPT(NK1)
     END FUNCTION TR_LOG_G0_TERM

     FUNCTION LOGCOSH(X)
        REAL(q) :: X
        REAL(q) :: LOGCOSH
        REAL(q) :: R, X2
       
        R = ABS(X) 
         
! add + LOG(1 + EXP(-2 * ABS(X)))
        X2=EXP(-2*ABS(X) )
        IF( X2 > -1._q ) THEN 
           IF ( ABS( X2 ) > 1.E-6_q ) THEN
              R = R+ LOG( 1._q + X2 ) 
! use mercator series
           ELSE
! log(1+x) = x - x^2/2 + x^3/3 - x^4/4
              R = R + X2 - X2**2/2 + X2**3/3
           ENDIF
        ENDIF

        LOGCOSH = R 
     END FUNCTION LOGCOSH


! accumulates G*S for Galitskii-Migdal energy
! and other terms for the grand potential
! Tr( G * S + Log(1-G*S):w

     SUBROUTINE LUTTINGER_WARD_POTENTIAL_TERMS( NK1, ISP, I, N, E_TR, E_GM,E_KLGG0, E_KLLOG )
        REAL(q)   , INTENT(INOUT) :: E_TR
        COMPLEX(q), INTENT(INOUT) :: E_GM
        COMPLEX(q), INTENT(INOUT) :: E_KLGG0,E_KLLOG
        INTEGER, INTENT( IN )     :: I, NK1, ISP, N
        INTEGER  :: J     ! diagonal index
        COMPLEX(q) :: GS_A(N), GS_B(N)
!#define unitTest
# 1375

   
! the problem at T=0  is that there is no common grid
! for even and odd basis functions in frequency domain
! as compromise we solve the dyson equation on cos grid
! in order to evaluate the w-integral, (1._q,0._q) would require
! to map G and Sigma to tau-domain and integrate there
! similar as we have 1._q it in the ADD_G_SIGMA_G routines
! for now, we do not evaluate the GM and LW terms at T=0
        IF ( .NOT. LFINITE_TEMPERATURE ) RETURN 
! accumulation is 1._q by calculating
! UU = ( G + G* ) x ( S + S* ) =: (GS + G*S*) + (G*S + GS*)
! VV = ( G - G* ) x ( S - S* ) =: (GS + G*S*) - (G*S + GS*)
!
! UU can be integrated with FER_RE grid and FER_RE_WEIGHTS,
! because all terms are proportional to U^2, where U = X/(X^2 + W^2)
! VV can be integrated with FER_IM grid and FER_IM_WEIGHTS,
! because all terms are proportional to V^2, where V = W/(X^2 + W^2)
!
! so the idea is to store GS_A = GS + G*S* and GS_B = G*S + GS*
! an to integrate GS_A + GS_B with FER_RE_WEIGHTS, while

!
! evaluate Galitskii-Migdal term: Tr { G.S_c }
!
! use interacting Green's function for Galtiskii-Migdal term
! see PRA73,012511 (2006) Equ.3
        GX(:) = SIGMAW_MAT(:,NK1,ISP,I)
! original Luttinger-Ward expression in PhyRev.118.1417(1960) Eq. 47
! has a typo. This is obvious by inspecting Eq. 56 in the same paper
! thus G_0 must be replaced by G in Eq. 47.
# 1409


! determine G S_c
        CALL PZGEMM( 'N', 'N', N, N, N, (1._q,0.0_q),  &
             SIGMA(1), 1, 1, GDES_MAT%DESC, &
             GX(1), 1, 1, GDES_MAT%DESC, &
             (0.0_q,0.0_q),  &
             GS(1), 1, 1, GDES_MAT%DESC )
# 1420

! determine S*G*
        CALL PZGEMM( 'C', 'C', N, N, N, (1._q,0.0_q),  &
             GX(1), 1, 1, GDES_MAT%DESC, &
             SIGMA(1), 1, 1, GDES_MAT%DESC, &
             (1.0_q,0.0_q),  &
             GS(1), 1, 1, GDES_MAT%DESC )
! obtain diagonal
        CALL DETERMINE_DIAGONALE(N, GS(:), CDIAG(:), GDES_MAT%DESC)
# 1431

        IF ( LFINITE_TEMPERATURE ) THEN
! accumulate Tr ( G*SIGMA_c)
           DO J = 1 , N
! drop a warining for complex energies, this should certainly never happen
              IF ( AIMAG( CDIAG(J) ) > 1.E-6_q ) THEN
                 WRITE(*,3)AIMAG( CDIAG(J) ) 
              ENDIF
              E_GM = E_GM +CDIAG(J)*IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2
           ENDDO
        ELSE
           GS_A = CDIAG 
! determine G S_c
           CALL PZGEMM( 'C', 'N', N, N, N, (1._q,0.0_q),  &
                SIGMA(1), 1, 1, GDES_MAT%DESC, &
                GX(1), 1, 1, GDES_MAT%DESC, &
                (0.0_q,0.0_q),  &
                GS(1), 1, 1, GDES_MAT%DESC )
! determine S*G*
           CALL PZGEMM( 'C', 'N', N, N, N, (1._q,0.0_q),  &
                GX(1), 1, 1, GDES_MAT%DESC, &
                SIGMA(1), 1, 1, GDES_MAT%DESC, &
                (1.0_q,0.0_q),  &
                GS(1), 1, 1, GDES_MAT%DESC )
! obtain diagonal
           CALL DETERMINE_DIAGONALE(N, GS(:), GS_B(:), GDES_MAT%DESC)

! accumulate Tr ( G*SIGMA_c)
! first term converges, but second (1._q,0._q) does not converge, because
! the term is not given on the proper sine grid, really nasty
! no real way around this, than to integrate in time
           DO J = 1 , N
              E_GM = E_GM &
                 + ( GS_A(J)+GS_B(J) )*IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/4 &
                 + ( GS_A(J)-GS_B(J) )*IMAG_GRIDS%FER_IM_WEIGHT(NOMEGA_GLOBAL)/4
           ENDDO

        ENDIF
3    FORMAT(' WARNING Galiskii-Migdal energy is complex: ',F20.10)
!
! evaluate ln|1 - G_x Sigma_c|^2
!

! Leeuwen and Luttinger agree on, how to evaluate following term
! however, Welden does not
        IF ( LWELDEN ) THEN
! use interacting Green's function for ln( 1-GS)
! see JCP 145, 204106 (2016) Eq. 20
           GX(:) = SIGMAW_MAT(:,NK1,ISP,I)
        ELSE
! i w - ( T + V_ext + V_h + V_x - E_fermi )
! use non-interacting Green's function for ln( 1-GS)
           GX(:)=0
           AUX(1:N)=1/(-W%CELTOT(1:N,NK1,ISP)+EFERMI(ISP)+&
              CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q))
           CALL ADD_TO_DIAGONALE(N,GX(:), AUX, GDES_MAT%DESC)
        ENDIF 

! see Equ. B3 in PRA 73, 012511 (2006)
! = S_c G_x + G_x* S_c* - S_c G_x G_x* S_c*
! determine S_c G_x  -> GS
        CALL PZGEMM( 'N', 'N', N, N, N, (1._q,0.0_q),  &
             SIGMA(1), 1, 1, GDES_MAT%DESC, &
             GX(1), 1, 1, GDES_MAT%DESC, &
             (0.0_q,0.0_q),  &
             GS(1), 1, 1, GDES_MAT%DESC )

# 1528

! form -1 + S_c . G_x
        AUX(1:N)=-1
        CALL ADD_TO_DIAGONALE(N,GS(:), AUX, GDES_MAT%DESC)
! form 1 - S_c . G_x
        GX(:) = -GS(:)
        GS(:) = GX(:)
        CALL PZGEMM( 'N', 'C', N, N, N, (1._q,0.0_q),  &
             GX(1), 1, 1, GDES_MAT%DESC, &
             GS(1), 1, 1, GDES_MAT%DESC, &
             (0.0_q,0.0_q),  &
             GSSG(1), 1, 1, GDES_MAT%DESC )
! diagonalize | 1 - S_c . G_x |^2
        CALL PZHEEVD_DESC( GSSG(:), EIGGS(1), N, GDES_MAT%DESC, GDES_MAT%COMM_INTRA)
! accumulate Tr Log | 1-G S|^2
        DO J = 1 , N
           IF ( EIGGS(J) < 0 ) THEN
              WRITE(*,'(A,4I4,F20.10)')" WARNING, negative argument of Log in &
                 &CALCULATE_G_FROM_SIGMA", ISP, NK1, I, J, EIGGS(J)
              CYCLE
           ENDIF
           E_TR = E_TR + LOG( EIGGS(J) )*&
              IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2
        ENDDO

!
! Klein functional Terms
!
!evaluate Tr{ GX^-1 G - 1 }
!
! quadrature has been created for positive frequencies
! so we need to rewrite as follows
!
! sum_n=-inf^inf ( GX^-1(iw)G(iw) - 1) =
! sum_n=0^inf    ( GX^-1(iw)G(iw) - 1)+(G(-iw)GX^-1(-iw) - 1)
! sum_n=0^inf    ( GX^-1(iw)G(iw) + G(-iw)GX^-1(-iw) - 2 )
!
! set  GX^-1
        GX(:)=0
        AUX(1:N)=(-W%CELTOT(1:N,NK1,ISP)+EFERMI(ISP)+&
           CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q))
        CALL ADD_TO_DIAGONALE(N,GX(:), AUX, GDES_MAT%DESC)
! build GX^-1 . G
        CALL PZGEMM( 'N', 'N', N, N, N, (1._q,0.0_q),  &
             GX(1), 1, 1, GDES_MAT%DESC, &
             SIGMAW_MAT(1,NK1,ISP,I), 1, 1, GDES_MAT%DESC, &
             (0.0_q,0.0_q),  &
             GS(1), 1, 1, GDES_MAT%DESC )
! GSSG <- GX^-1.G
        GSSG(:) = GS(:)  
# 1580

! add G* GX*^-1
        CALL PZGEMM( 'C', 'C', N, N, N, (1._q,0.0_q),  &
             SIGMAW_MAT(1,NK1,ISP,I), 1, 1, GDES_MAT%DESC, &
             GX(1), 1, 1, GDES_MAT%DESC, &
             (1.0_q,0.0_q),  &
             GS(1), 1, 1, GDES_MAT%DESC )
# 1592

        CDIAG = 0 
        CALL DETERMINE_DIAGONALE(N, GS(:), CDIAG(:), GDES_MAT%DESC)
        DO J = 1 , N
           E_KLGG0 = E_KLGG0 + (CDIAG(J)-2)*&
              IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2
        ENDDO
! evaluate -Tr Ln( GX^-1 G ) = -Ln( GX^-1 G * G^+ GX^-1+)
! form GX^-1.G . G*.GX*^-1
        GX(:) = GSSG(:)  
        CALL PZGEMM( 'N', 'C', N, N, N, (1._q,0.0_q),  &
             GSSG(1), 1, 1, GDES_MAT%DESC, &
             GX(1), 1, 1, GDES_MAT%DESC, &
             (0.0_q,0.0_q),  &
             GS(1), 1, 1, GDES_MAT%DESC )
! diagonalize GX^-1.G . G*.GX*^-1
        CALL PZHEEVD_DESC( GS(:), EIGGS(1), N, GDES_MAT%DESC, GDES_MAT%COMM_INTRA)
        DO J = 1 , N
           IF ( EIGGS(J) < 0 ) THEN
              WRITE(*,'(A,4I4,F20.10)')" WARNING, ignoring negative argument of Log in &
                 &Klein functional: CALCULATE_G_FROM_SIGMA", ISP, NK1, I, J, EIGGS(J)
              CYCLE
           ENDIF
           E_KLLOG = E_KLLOG - LOG(EIGGS(J))*&
              IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2
        ENDDO
     END SUBROUTINE LUTTINGER_WARD_POTENTIAL_TERMS

   END SUBROUTINE CALCULATE_G_FROM_SIGMA


!*********************************************************************
!
! this subroutine determines and stores the mean field density matrix
! input/ output  CHAM_MAT
!
!  1/(i 2 pi) int_-infty^infty G_mean(w) e^ (i w eta) dw
!
! poles of states below the Fermi-level above real axis
!
! G_mean = ( w- H_mean) ^-1  ; H_mean = T+ V_H+ V_x
!
! with  H_mean = epsilon_HF  this simply becomes
!  Theta(e_fermi-epsilon_HF)
!
! TODO: old version
! with  H_mean = U epsilon_HF U+  this simply becomes
!  U Theta(e_fermi-epsilon_HF) U+
! where U is the matrix diagonalizing H_mean, and Theta is the
! step function
!
! if the Green function is passed down (GU_MAT and GO_MAT) the
! mean field Green function is also added to the two matrices
!
! Important note: this function must be exactly compatible to
! HF_DENSITY_MATRIX that determines the mean field G from frequency
! integration
!
!*********************************************************************

   SUBROUTINE SUBTRACT_GHF_FROM_G( W, NKPTS_IRZ, EFERMI, GW_MAT, & 
        IMAG_GRIDS, GDES_MAT, CHAM_MAT)
     USE base
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)         :: W
     INTEGER                 :: NKPTS_IRZ
     REAL(q)                 :: EFERMI(:)          ! in: fermi-level
     COMPLEX(q)              :: GW_MAT(:,:,:,:)    ! in/out: interacting Green function
     COMPLEX(q),OPTIONAL           :: CHAM_MAT(:,:,:,:)  ! in: frequency independent part of self-energy (T +V_x + V_H)
     TYPE (imag_grid_handle) :: IMAG_GRIDS         ! imaginary grids
     TYPE (greens_mat_des)   :: GDES_MAT           ! descriptor
! local variables
     INTEGER                 :: I                  ! frequency points
     INTEGER                 :: N                  ! bands for each k-points
     INTEGER                 :: NOMEGA_GLOBAL      ! global frequency point
     INTEGER                 :: ISP, NK1, NB       ! current spin, k-point and band index
     COMPLEX(q)              :: AUX( SIZE(W%CELTOT,1))
     COMPLEX(q), ALLOCATABLE :: G1(:), G2(:)

     ALLOCATE(G1(SIZE(GW_MAT,1)), G2(SIZE(GW_MAT,1)))

     DO ISP=1,W%WDES%ISPIN
        DO NK1=1,NKPTS_IRZ
           N=  W%WDES%NB_TOTK(NK1, ISP)
           DO I=1,SIZE(GW_MAT,4)  ! loop over local frequency points
              
!this is the global frequency point of the node
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE
!
! now subtract the "mean" field  Greens function - (i w - H_mean )^-1
!
              AUX(1:N)=W%CELTOT(1:N,NK1,ISP)-EFERMI(ISP)
! this is actually necessary only for T=0 integration
              IF ( .NOT.LFINITE_TEMPERATURE ) THEN
                 DO NB=1, N
                    IF (W%FERTOT(NB,NK1,ISP)>FERMI_OCC) THEN
! "occupied" shift mean field eigenvalue to below -OMEGAMIN=IMAG_GRIDS%X1
                       AUX(NB)=MIN( REAL(AUX(NB),q),  -IMAG_GRIDS%X1)
                    ELSE
! "un-occupied" shift mean field eigenvalue to above OMEGAMIN
                       AUX(NB)=MAX( REAL(AUX(NB),q),   IMAG_GRIDS%X1)
                    ENDIF
                 ENDDO
              ENDIF
              AUX(1:N)=-1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q)-AUX(1:N))
!
! AUX is now the exchange Green's function in frequency domain
! in the Fock basis and it thus diagonal

! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
# 1718

! TODO: remove G1 and G2
! subtract 1/(i w -epsilon_HF) from the Green's function
              CALL ADD_TO_DIAGONALE(N, GW_MAT(:, NK1, ISP, I), AUX, GDES_MAT%DESC )

           ENDDO
        ENDDO
     ENDDO

     DEALLOCATE(G1, G2 )
   END SUBROUTINE SUBTRACT_GHF_FROM_G


!*********************************************************************
!
! this routine optimizes the EFERMI to conserve the
! number of electrons form the full self-energy
! a merger of several other routines
! including CALCULATE_G_FROM_SIGMA
!           GAMMA_FROM_G_W
! to avoid excessive storage, most things are 1._q on the fly without
! storing intermediate values
! TODO: CHAM could be made COMPLEX(q)
!
!*********************************************************************

   SUBROUTINE ITERATE_EFERMI_G_FROM_SIGMA( W, NKPTS_IRZ, NELECT_IN, NUP_DOWN, EFERMI, EFERMI_SLOPE, SIGMAW_MAT, & 
        IMAG_GRIDS, GDES_MAT, CHAM_MAT, IO)
     USE base
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)    :: W
     INTEGER            :: NKPTS_IRZ
     REAL(q)            :: NELECT_IN, NUP_DOWN! number of electrons, and spin-polarization
     REAL(q)            :: EFERMI(:)          ! in: fermi-level
     REAL(q)            :: EFERMI_SLOPE(2)    ! change of E-fermi/ change of density from ITERATE_EFERMI_FROM_SIGMA_DIAG
     COMPLEX(q)         :: SIGMAW_MAT(:,:,:,:)! in: matrix <i| Sigma_c(omega) | a> stored distributed
! out: interacting Green function
     COMPLEX(q)               :: CHAM_MAT(:,:,:,:)  ! in: frequency independent part of self-energy (T +V_x + V_H)
     TYPE (imag_grid_handle) :: IMAG_GRIDS    ! imaginary grids

     TYPE (greens_mat_des) :: GDES_MAT        ! greens function matrix descriptor
     TYPE( in_struct )  :: IO 
! local variables
     REAL(q)            :: NELECT             ! desired number of electrons
     INTEGER            :: I,J                ! frequency points
     INTEGER            :: N                  ! bands for each k-points
     INTEGER            :: NOMEGA_GLOBAL      ! corresponing global frequency point
     INTEGER            :: ISP, NK, NB        ! current spin, q-point and band index
     COMPLEX(q)         :: AUX( SIZE(W%CELTOT,1))
     COMPLEX(q), ALLOCATABLE :: G1(:), G2(:), CHAM(:), GMAT(:)
     COMPLEX(q), ALLOCATABLE :: GW(:,:,:,:)   ! diagonal part of interacting Green's function at each frequency
     COMPLEX(q), ALLOCATABLE  :: GAMMA(:,:,:)       ! occupancy of each states
     REAL(q) :: CHARGE_MEAN, CHARGE_CORR
     INTEGER            :: INFO
     INTEGER            :: IPIV(W%WDES%NB_TOT)
     INTEGER            :: LWORK, LIWORK
     INTEGER, ALLOCATABLE   :: IWORK(:)
     COMPLEX(q),ALLOCATABLE :: WORK(:)
     REAL(q)            :: RSPIN
     INTEGER            :: ISPIN_MAX, ISPIN, ISP_BASE
     LOGICAL :: LBOUND(2)                     ! determines whether upper/lower bound for EFERMI was found
     REAL(q) :: CHARGE_DIFF(2)                ! error in the number of electrons for upper/lower bound
     REAL(q) :: CHARGE_DIFF_MEAN(2)           ! error in the number of electrons for upper/lower bound for mean field part
     REAL(q) :: EFERMI_BOUND(2)               ! upper and lower bounds for EFERMI
     REAL(q) :: EFERMI_STEP=0.001              ! default step for Fermi-level
     INTEGER :: MAXSTEP=50                    ! maximum number of interations on Fermi-level
     REAL(q) :: CHARGE_THRESHHOLD=1E-9_q
     INTEGER :: IU0
     REAL(q) :: XBETA

     

     ALLOCATE(G1(SIZE(SIGMAW_MAT,1)), G2(SIZE(SIGMAW_MAT,1)), CHAM(SIZE(SIGMAW_MAT,1)), GMAT(SIZE(SIGMAW_MAT,1)))
     ALLOCATE(GW(W%WDES%NB_TOT, NKPTS_IRZ, W%WDES%ISPIN, SIZE(SIGMAW_MAT,4)), & 
              GAMMA(W%WDES%NB_TOT, NKPTS_IRZ,  W%WDES%ISPIN))

     CALL QUERRY_WORK_SPACE

! do we need outer loop since spin polarization is fixed
     IF (NUP_DOWN >= 0 .AND. W%WDES%ISPIN>1 ) THEN
        ISPIN_MAX=1
        ISPIN    =1    ! inner loop maximum is 1
        RSPIN    =2    ! spin multiplicity
     ELSE
        ISPIN_MAX=0
        ISPIN=W%WDES%ISPIN
        RSPIN=W%WDES%RSPIN
     ENDIF
!
! outer loop, required if spin up and spin down are treated seperately
!
out: DO ISP_BASE=0,ISPIN_MAX

! set desired number of electrons
     IF ( ISPIN_MAX==1 ) THEN
        IF (ISP_BASE==0) NELECT=NELECT_IN+NUP_DOWN
        IF (ISP_BASE==1) NELECT=NELECT_IN-NUP_DOWN
     ELSE
        NELECT=NELECT_IN
     ENDIF

     LBOUND=.FALSE. ! no bounds yet found

! loop over the maximum number of steps
     DO J=1,MAXSTEP

! determine diagonal part of correlated Green's function
! note that frequency points are blocked here
     CALL DETERMINE_DIAG_G_FROM_SIGMA(W, EFERMI,G1,G2,CHAM_MAT,GMAT, GDES_MAT, GW)
! G_mean is not subtracted at T>0
! get diagonal of corresponding density matrix
     CALL GAMMA_FROM_G_W_DIAG( W, NKPTS_IRZ, GW, GAMMA, IMAG_GRIDS )
! sum over all distributed time points
     CALL M_sum_z(GDES_MAT%COMM_INTER, GAMMA, SIZE(GAMMA))

! sum charge from all diagonal elements, and add mean field occupancies
     CHARGE_CORR=0
     CHARGE_MEAN=0
     DO ISP=1,ISPIN
     DO NK=1,NKPTS_IRZ
! update FERMI weights in mean-field function
        IF ( LFINITE_TEMPERATURE ) THEN
           DO I = 1 , W%WDES%NB_TOTK( NK, ISP+ISP_BASE ) 
              XBETA = (W%CELTOT( I, NK, ISP+ISP_BASE )-EFERMI(ISP+ISP_BASE)) *IMAG_GRIDS%T%BETA
              IF ( XBETA < -40.0_q ) THEN
                 W%FERTOT(I,NK,ISP+ISP_BASE) = 1.0_q
              ELSE IF ( XBETA > 40.0_q ) THEN
                 W%FERTOT(I,NK,ISP+ISP_BASE) = 0.0_q
              ELSE
                 W%FERTOT(I,NK,ISP+ISP_BASE) = 1/( 1+EXP( XBETA ) )
              ENDIF
             
              GAMMA(I,NK,ISP)=GAMMA(I,NK,ISP) - W%FERTOT(I,NK, ISP+ISP_BASE)
           END DO
        ENDIF
! sum all bands
        CHARGE_CORR=CHARGE_CORR+SUM(GAMMA(1:W%WDES%NB_TOTK(NK,ISP+ISP_BASE),NK, ISP+ISP_BASE))*KPOINTS_ORIG%WTKPT(NK)*RSPIN
        DO NB=1,W%WDES%NB_TOTK(NK, ISP+ISP_BASE)
           IF ( LFINITE_TEMPERATURE ) THEN
              GAMMA(NB,NK, ISP)=GAMMA(NB,NK, ISP)+W%FERTOT(NB,NK,ISP+ISP_BASE)
              CHARGE_MEAN=CHARGE_MEAN+KPOINTS_ORIG%WTKPT(NK)*RSPIN*W%FERTOT(NB,NK,ISP+ISP_BASE)
           ELSE
              IF (W%FERTOT(NB,NK,ISP+ISP_BASE)>FERMI_OCC) THEN
                 GAMMA(NB,NK, ISP)=GAMMA(NB,NK, ISP)+1
                 CHARGE_MEAN=CHARGE_MEAN+KPOINTS_ORIG%WTKPT(NK)*RSPIN
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     ENDDO

!WRITE(200,'(I4, 2L3, 8F14.6)')J,LBOUND, CHARGE_CORR, CHARGE_MEAN,EFERMI,EFERMI_BOUND
! subtract the desired electron number
     CHARGE_CORR=CHARGE_CORR+CHARGE_MEAN-NELECT
     CHARGE_MEAN=CHARGE_MEAN-NELECT
!WRITE(201,'(I4, 2L3, 8F14.6)')J,LBOUND, CHARGE_CORR, CHARGE_MEAN,EFERMI,EFERMI_BOUND

     IF (J==1) THEN
! use slope determined by ITERATE_EFERMI_FROM_SIGMA_DIAG
        EFERMI_STEP=ABS(EFERMI_SLOPE(1+ISP_BASE)*CHARGE_CORR)
!CHARGE_THRESHHOLD=MAX(CHARGE_THRESHHOLD, ABS(CHARGE_CORR)/100000)
     ENDIF

     IF (IO%NWRITE>3 ) THEN
        IF (IO%IU0 >=0) &
        WRITE(*,'(" E-Fermi and density",6F14.7)') EFERMI(1+ISP_BASE), CHARGE_MEAN, CHARGE_CORR
     ENDIF

     IF (ABS(CHARGE_CORR)<CHARGE_THRESHHOLD) EXIT

     IF (CHARGE_CORR>0) THEN
! found a new upper bound
        LBOUND(1)=.TRUE.
        CHARGE_DIFF(1)=CHARGE_CORR
        CHARGE_DIFF_MEAN(1)=CHARGE_MEAN
        EFERMI_BOUND(1)=EFERMI(1+ISP_BASE)
     ELSE
! found a new lower bound
        LBOUND(2)=.TRUE.
        CHARGE_DIFF(2)=CHARGE_CORR
        CHARGE_DIFF_MEAN(2)=CHARGE_MEAN
        EFERMI_BOUND(2)=EFERMI(1+ISP_BASE)
     ENDIF
     IF (LBOUND(1) .AND. LBOUND(2)) THEN
        IF (CHARGE_DIFF_MEAN(1)==CHARGE_DIFF_MEAN(2) .AND.  ABS(EFERMI_BOUND(1)-EFERMI_BOUND(2))<IMAG_GRIDS%X1*4) THEN
! if mean field electron number is identical and estimates are close
! try to interpolate Fermi-level to obtain correct number
! TODO: we might need some fall back to intervall bisectioning when stagnetion sets in (Brent's algorithm)
           EFERMI(1+ISP_BASE)=EFERMI_BOUND(1)-(EFERMI_BOUND(1)-EFERMI_BOUND(2))/(CHARGE_DIFF(1)-CHARGE_DIFF(2))*CHARGE_DIFF(1)
        ELSE
! use intervall bisectioning
           EFERMI(1+ISP_BASE)=(EFERMI_BOUND(1)+EFERMI_BOUND(2))/2
        ENDIF
     ELSE IF (LBOUND(1)) THEN
! only upper bound found
        EFERMI(1+ISP_BASE)=EFERMI(1+ISP_BASE)-EFERMI_STEP
     ELSE
! only lower  bound found
        EFERMI(1+ISP_BASE)=EFERMI(1+ISP_BASE)+EFERMI_STEP
     ENDIF
! finally fix fermi-level for second spin channel (if outer loop is NOT used)
     IF (ISPIN_MAX==0.AND.SIZE(EFERMI)>1) EFERMI(2)=EFERMI(1)

     ENDDO
   ENDDO out

     DEALLOCATE(GW, GAMMA)
     DEALLOCATE(WORK, IWORK)
     DEALLOCATE(G1, G2, CHAM, GMAT)

     

   CONTAINS

!=====================================================================
! querry 1 work-space and allocate optimal workspace for PZGETRI
! used below
!=====================================================================

  SUBROUTINE QUERRY_WORK_SPACE
     ALLOCATE(WORK(1), IWORK(1))
     LWORK=-1
     LIWORK=-1
     CALL PZGETRI( W%WDES%NB_TOT, G1(1), 1, 1, GDES_MAT%DESC, IPIV, WORK, LWORK, &
          IWORK, LIWORK, INFO )
     LWORK=MAX(INT(WORK(1)),1)
     LIWORK=MAX(IWORK(1),1)
     DEALLOCATE(WORK, IWORK)
     ALLOCATE(WORK(LWORK), IWORK(LIWORK))
   END SUBROUTINE QUERRY_WORK_SPACE

!=====================================================================
! determine diagonal part of correlated Green's function
!=====================================================================

   SUBROUTINE DETERMINE_DIAG_G_FROM_SIGMA(W, EFERMI, G1, G2, CHAM_MAT,&
      GMAT, GDES_MAT, GW  )
     TYPE (wavespin)    :: W
     REAL(q)            :: EFERMI(:)          ! in: fermi-level
     COMPLEX(q)         :: G1(:)
     COMPLEX(q)         :: G2(:)
     COMPLEX(q)               :: CHAM_MAT(:,:,:,:)  ! in: frequency independent part of self-energy (T +V_x + V_H)
     COMPLEX(q)         :: GMAT(:)
     TYPE (greens_mat_des) :: GDES_MAT        ! greens function matrix descriptor
     COMPLEX(q)         :: GW(:,:,:,:)
     INTEGER      :: NK, NB, ISP, N, I, NOMEGA_GLOBAL
     COMPLEX(q)   :: AUX( SIZE(W%CELTOT,1))

     GW=0
     DO ISP=1,ISPIN
        DO NK=1,NKPTS_IRZ
           N=  W%WDES%NB_TOTK(NK, ISP+ISP_BASE)
! this is H_HF - MU
           AUX(1:N)=W%CELTOT(1:N,NK,ISP+ISP_BASE)-EFERMI(ISP+ISP_BASE)
! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
# 1988

! subtract Fermi-level from diagonal Hartree-Fock Hamiltonian
           CHAM=0
           CALL ADD_TO_DIAGONALE(N, CHAM, AUX, GDES_MAT%DESC )

! now CHAM contains H_HF - MU

           DO I=1,SIZE(SIGMAW_MAT,4)  ! loop over local frequency points
!this is the global omega point of the node
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE

              IF ( LGREENHF ) THEN
! add mean field HF and change sign -Sigma_c - H_HF + MU
                 GMAT=-CHAM
              ELSE
! add mean field HF and change sign -Sigma_c - H_HF + MU
                 GMAT=-(SIGMAW_MAT(:, NK, ISP+ISP_BASE, I)+CHAM)
              ENDIF
!  calculate i w +e_fermi - Sigma_c - (T +V_x + V_H)
              AUX=CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q)
! add explicit frequency dependence
! this make G^-1 complete
              CALL ADD_TO_DIAGONALE(N, GMAT, AUX, GDES_MAT%DESC)
              
! CALL DUMP_GREENS_HAM( "G", G(:), N, GDES_MAT, 6)
              
!
! invert full Green's function
!
! first LU decomposition
              CALL PZGETRF( N, N, & 
                   GMAT, 1, 1, GDES_MAT%DESC, IPIV, INFO )
              IF (INFO/=0) THEN
                 CALL vtutor%bug("PZGETRF in CALCULATE_G_FROM_SIGMA, INFO= " // str(INFO), "greens_orbital.F" , 2023 )
              ENDIF
! then inversion
              CALL PZGETRI( N, GMAT, 1, 1, GDES_MAT%DESC, IPIV, WORK, LWORK, &
                   IWORK, LIWORK, INFO )
              IF (INFO/=0) THEN
                 CALL vtutor%bug("PZGETRI in CALCULATE_G_FROM_SIGMA, INFO= " // str(INFO), "greens_orbital.F", 2029 )
              ENDIF

              CALL DETERMINE_DIAGONALE( W%WDES%NB_TOTK(NK,ISP+ISP_BASE), GMAT, GW(:, NK, ISP+ISP_BASE, I), GDES_MAT%DESC)
              CALL M_sum_z(GDES_MAT%COMM_INTRA, GW(:, NK, ISP+ISP_BASE, I), W%WDES%NB_TOTK(NK,ISP+ISP_BASE) )
! now subtract a reference mean field Greens function - (i w - H_mean )^-1
              AUX(1:N)=W%CELTOT(1:N,NK,ISP+ISP_BASE)-EFERMI(ISP+ISP_BASE)
! to obtain convergent Matsubara sum (does not really matter much
! what is subtracted, as long as (1._q,0._q) does not get too close to OMEGAMIN)
! this is actually necessary only for T=0 integration
              IF ( .NOT.LFINITE_TEMPERATURE ) THEN
                 DO NB=1, N
                    IF (W%FERTOT(NB,NK,ISP+ISP_BASE)>FERMI_OCC) THEN
! occupied shift mean field eigenvalue to below -OMEGAMIN
                       AUX(NB)=MIN( REAL(AUX(NB),q),-IMAG_GRIDS%X1)
                    ELSE
! un-occupied shift mean field eigenvalue to above OMEGAMIN
                       AUX(NB)=MAX( REAL(AUX(NB),q), IMAG_GRIDS%X1)
                    ENDIF
                 ENDDO
                 GW(1:N, NK, ISP+ISP_BASE, I)=GW(1:N, NK, ISP+ISP_BASE, I)- &
                      1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q)-AUX(1:N))
              ENDIF
           ENDDO
        ENDDO
     ENDDO
   END SUBROUTINE DETERMINE_DIAG_G_FROM_SIGMA
   
   END SUBROUTINE ITERATE_EFERMI_G_FROM_SIGMA


!*********************************************************************
!
! this routine optimizes EFERMI from the diagonal entries in
! the self-energy to conserve NELECT
!
!*********************************************************************

   SUBROUTINE ITERATE_EFERMI_FROM_SIGMA_DIAG(W, NKPTS_IRZ, NELECT_IN, NUP_DOWN, EFERMI, EFERMI_SLOPE, & 
        IMAG_GRIDS, G_COS, G_SIN, GDES_MAT, IO)
     USE constant
     USE wave
     IMPLICIT NONE
     TYPE (wavespin) :: W                  ! orbital
     INTEGER         :: NKPTS_IRZ          ! number of k-points in IRZ
     REAL(q)         :: NELECT_IN, NUP_DOWN! number of electrons, and spin-polarization
     REAL(q)         :: EFERMI(:)          ! in: fermi-level
     REAL(q)         :: EFERMI_SLOPE(2)    ! initial error in charge density, on exit: change of E-fermi/ change of density
     TYPE (imag_grid_handle) :: IMAG_GRIDS ! imaginary grids
     COMPLEX(q)      :: G_COS(:,:,:,:)     !diagonal of cos transformed Sigma in orbital basis
     COMPLEX(q)      :: G_SIN(:,:,:,:)     !diagonal of sin transformed Sigma in orbital basis
     TYPE (greens_mat_des) :: GDES_MAT     !descriptor
     TYPE( in_struct )  :: IO 
! local variables

     REAL(q)            :: NELECT             ! desired number of electrons
     INTEGER            :: I,J                ! frequency points
     INTEGER            :: N                  ! bands for each k-points
     INTEGER            :: ISP, NK, NB        ! current spin, q-point and band index
     REAL(q)            :: CHARGE_MEAN, CHARGE_CORR
     COMPLEX(q), ALLOCATABLE  :: GAMMA(:,:,:)       ! occupancy of each states
     COMPLEX(q), ALLOCATABLE :: GW(:,:,:,:)   ! diagonal part of interacting Green's function at each frequency
     REAL(q)            :: RSPIN
     INTEGER            :: ISPIN_MAX, ISPIN, ISP_BASE
     LOGICAL            :: LBOUND(2)          ! determines whether upper/lower bound for EFERMI was found
     REAL(q) :: CHARGE_DIFF(2)                ! error in the number of electrons for upper/lower bound
     REAL(q) :: CHARGE_DIFF_MEAN(2)           ! error in the number of electrons for upper/lower bound for mean field part
     REAL(q) :: EFERMI_BOUND(2)               ! upper and lower bounds for EFERMI
     REAL(q) :: EFERMI_IN(2)                  ! initial Fermi-leve
     REAL(q) :: EFERMI_STEP=0.05              ! default step for Fermi-level
     INTEGER :: MAXSTEP=50                    ! maximum number of interations on Fermi-level
     REAL(q) :: CHARGE_THRESHHOLD=1E-7_q
     INTEGER :: IU0
     REAL(q) :: XBETA

     

     ALLOCATE(GW(W%WDES%NB_TOT, NKPTS_IRZ, W%WDES%ISPIN, SIZE(G_COS,1)), & 
          GAMMA(W%WDES%NB_TOT, NKPTS_IRZ,  W%WDES%ISPIN))

! do we need outer loop since spin polarization is fixed
     IF (NUP_DOWN >= 0 .AND. W%WDES%ISPIN>1 ) THEN
        ISPIN_MAX=1
        ISPIN    =1    ! inner loop maximum is 1
        RSPIN    =2    ! spin multiplicity
     ELSE
        ISPIN_MAX=0
        ISPIN=W%WDES%ISPIN
        RSPIN=W%WDES%RSPIN
     ENDIF

     EFERMI_IN(1)=EFERMI(1)
     EFERMI_IN(2)=EFERMI(SIZE(EFERMI))
!
! outer loop, required if spin up and spin down are treated seperately
!
out: DO ISP_BASE=0,ISPIN_MAX

! set desired number of electrons
     IF ( ISPIN_MAX==1 ) THEN
        IF (ISP_BASE==0) NELECT=NELECT_IN+NUP_DOWN
        IF (ISP_BASE==1) NELECT=NELECT_IN-NUP_DOWN
     ELSE
        NELECT=NELECT_IN
     ENDIF

     LBOUND=.FALSE. ! no bounds yet found
! loop over the maximum number of steps
     DO J=1,MAXSTEP

! determine diagonal part of correlated Green's function
     CALL DETERMINE_GW_FROM_SIGMA_DIAG(J)
! G_mean is not subtracted at T>0
! get correlation contribution to diagonal of density matrix
     CALL GAMMA_FROM_G_W_DIAG_NONBLOCKED( W, NKPTS_IRZ, GW, GAMMA, IMAG_GRIDS )
! subtract meanfield density matrix

! sum charge from all diagonal elements, and add mean field occupancies
     CHARGE_CORR=0
     CHARGE_MEAN=0
     DO ISP=1,ISPIN
     DO NK=1,NKPTS_IRZ
! first update FERMI weights of mean-field part
! because they determine mean-field charge
        IF ( LFINITE_TEMPERATURE ) THEN
           DO I = 1 , W%WDES%NB_TOTK( NK, ISP+ISP_BASE ) 
              XBETA = (W%CELTOT( I, NK, ISP+ISP_BASE )-EFERMI(ISP+ISP_BASE)) *IMAG_GRIDS%T%BETA
              IF ( XBETA < -40.0_q ) THEN
                 W%FERTOT(I,NK,ISP+ISP_BASE) = 1.0_q
              ELSE IF ( XBETA > 40.0_q ) THEN
                 W%FERTOT(I,NK,ISP+ISP_BASE) = 0.0_q
              ELSE
                 W%FERTOT(I,NK,ISP+ISP_BASE) = 1/( 1+EXP( XBETA ) )
              ENDIF
              GAMMA(I,NK,ISP)=GAMMA(I,NK,ISP) - W%FERTOT(I,NK, ISP+ISP_BASE)
           END DO
        ENDIF

! obtain correlated charge contribution as sum of all bands
! including k-point integration and spin multiplicity
        CHARGE_CORR=CHARGE_CORR+SUM(GAMMA(1:W%WDES%NB_TOTK(NK,ISP+ISP_BASE),NK, ISP+ISP_BASE))*KPOINTS_ORIG%WTKPT(NK)*RSPIN

! now do the same integration of the mean-field part
        DO NB=1,W%WDES%NB_TOTK(NK, ISP+ISP_BASE)
           IF ( LFINITE_TEMPERATURE ) THEN
! at finite temperature we take the correct fermi weights
! (determined above) into account
              GAMMA(NB,NK, ISP)=GAMMA(NB,NK, ISP)+W%FERTOT(NB,NK,ISP+ISP_BASE)
              CHARGE_MEAN=CHARGE_MEAN+KPOINTS_ORIG%WTKPT(NK)*RSPIN*W%FERTOT(NB,NK,ISP+ISP_BASE)
           ELSE
! at (0._q,0._q) temperature weights larger as FERMI_OCC ( that is typically 0.5)
! are counted as occupied, so integer occupancies are forced
! for the mean-field part. This causes an error for metals at T=0
              IF (W%FERTOT(NB,NK,ISP+ISP_BASE)>FERMI_OCC) THEN
                 GAMMA(NB,NK, ISP)=GAMMA(NB,NK, ISP)+1
                 CHARGE_MEAN=CHARGE_MEAN+KPOINTS_ORIG%WTKPT(NK)*RSPIN
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     ENDDO

!WRITE(100,'(I4, 2L3, 8F14.6)')J,LBOUND, CHARGE_CORR, CHARGE_MEAN,EFERMI,EFERMI_BOUND
! subtract the desired electron number
     CHARGE_CORR=CHARGE_CORR+CHARGE_MEAN-NELECT
     CHARGE_MEAN=CHARGE_MEAN-NELECT
!WRITE(101,'(I4, 2L3, 8F14.6)')J,LBOUND, CHARGE_CORR, CHARGE_MEAN, EFERMI,EFERMI_BOUND


     IF (IO%NWRITE>3 ) THEN
        IF (IO%IU0 >=0) &
        WRITE(*,'(" E-Fermi and density DIAG ",6F14.7)') EFERMI(1+ISP_BASE), CHARGE_MEAN, CHARGE_CORR
     ENDIF

! set break condition to charge violation upon entry /1000
! (since the routine is very fast, we can invest some compute time)
     IF (J==1) THEN
        EFERMI_SLOPE(1+ISP_BASE)=CHARGE_CORR
        CHARGE_THRESHHOLD=MAX(1E-8_q, ABS(CHARGE_CORR)/1000)
     ENDIF

!     EFERMI(1+ISP_BASE)=EFERMI(1+ISP_BASE)+0.0001
!     CYCLE

     IF (ABS(CHARGE_CORR)<CHARGE_THRESHHOLD) EXIT

     IF (CHARGE_CORR>0) THEN
! found a new upper bound
        LBOUND(1)=.TRUE.
        CHARGE_DIFF(1)=CHARGE_CORR
        CHARGE_DIFF_MEAN(1)=CHARGE_MEAN
        EFERMI_BOUND(1)=EFERMI(1+ISP_BASE)
     ELSE
! found a new lower bound
        LBOUND(2)=.TRUE.
        CHARGE_DIFF(2)=CHARGE_CORR
        CHARGE_DIFF_MEAN(2)=CHARGE_MEAN
        EFERMI_BOUND(2)=EFERMI(1+ISP_BASE)
     ENDIF
     IF (LBOUND(1) .AND. LBOUND(2)) THEN
! estimate from linear interpolation (not always robust)
! EFERMI(1+ISP_BASE)=EFERMI_BOUND(1)-(EFERMI_BOUND(1)-EFERMI_BOUND(2))/(CHARGE_DIFF(1)-CHARGE_DIFF(2))*CHARGE_DIFF(1)
! use intervall bisectioning
        EFERMI(1+ISP_BASE)=(EFERMI_BOUND(1)+EFERMI_BOUND(2))/2
     ELSE IF (LBOUND(1)) THEN
! only upper bound found
        EFERMI(1+ISP_BASE)=EFERMI(1+ISP_BASE)-EFERMI_STEP
     ELSE
! only lower  bound found
        EFERMI(1+ISP_BASE)=EFERMI(1+ISP_BASE)+EFERMI_STEP
     ENDIF
! finally fix fermi-level for second spin channel (if outer loop is NOT used)
     IF (ISPIN_MAX==0.AND.SIZE(EFERMI)>1) EFERMI(2)=EFERMI(1)

     ENDDO
   ENDDO out

! determine slope: required change of fermi-energy/ initial charge error
   DO ISP=1,W%WDES%ISPIN
      EFERMI_SLOPE(ISP)=(EFERMI_IN(ISP)-EFERMI(ISP))/EFERMI_SLOPE(ISP)
   ENDDO

     IF (IO%NWRITE>3 ) THEN
        IF (IO%IU0 >=0) &
      WRITE(*,'(" E-Fermi slope ",6F14.7)') EFERMI_SLOPE
   ENDIF


   DEALLOCATE(GW, GAMMA)
   
   CONTAINS
!=====================================================================
! determine diagonal part of correlated Green's function from diagonal
! part of self energy
!=====================================================================

   SUBROUTINE DETERMINE_GW_FROM_SIGMA_DIAG(J)
     INTEGER :: J
     REAL(q) :: E
     COMPLEX(q) :: SIGMA,EPS_MU
     INTEGER :: ISP, NK, I 
     INTEGER :: N

     GW=0
!loop over all orbitals
     DO ISP=1,ISPIN
        DO NK=1,NKPTS_IRZ
           N= W%WDES%NB_TOTK(NK,ISP+ISP_BASE)
           
! frequency points are not blocked here
           DO I=1,SIZE(G_COS,1)  ! loop over all frequency points
              
! loop over all bands
            DO NB=1, N
               IF ( LGREENHF ) THEN
! Hartree-Fock greens function is determined
                  SIGMA = 0
               ELSE
! Sigma_c = self-energy
                  SIGMA = G_COS(I,NB,NK,ISP+ISP_BASE)+G_SIN(I,NB,NK,ISP+ISP_BASE)*(0._q,1._q)
               ENDIF
! H_HF - MU = epsilon - mu
                
               EPS_MU = W%CELTOT(NB,NK,ISP+ISP_BASE)-EFERMI(ISP+ISP_BASE)
! diagonal of Green function  G(iw) = 1/(i w - ( H_HF - mu ) + Sigma_c )
               GW(NB, NK, ISP+ISP_BASE, I) = 1.0_q/((0._q,1._q)*IMAG_GRIDS%FER_RE(I)+SIGMA-EPS_MU)
           
! subtract mean field Green function at T=0
               IF ( .NOT. LFINITE_TEMPERATURE ) THEN
                  E = W%CELTOT(NB,NK,ISP+ISP_BASE)-EFERMI(ISP+ISP_BASE)
!, but shift poles if it is too close for numerical integration
! this is actually necessary only for T=0 integration
                  IF ( .NOT.LFINITE_TEMPERATURE ) THEN
                     IF (W%FERTOT(NB,NK,ISP+ISP_BASE)>FERMI_OCC) THEN
! occupied shift mean field eigenvalue to below -OMEGAMIN
                        E=MIN(E,-IMAG_GRIDS%X1)
                     ELSE
                        E=MAX(E, IMAG_GRIDS%X1)
                     ENDIF
                  ENDIF
                  GW(NB, NK, ISP+ISP_BASE, I)=GW(NB, NK, ISP+ISP_BASE, I)-&
                  1.0_q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-E)
               ENDIF 
            ENDDO
         ENDDO
      ENDDO  !k-point
      ENDDO  !spin

    END SUBROUTINE DETERMINE_GW_FROM_SIGMA_DIAG


    SUBROUTINE GAMMA_FROM_G_W_DIAG_NONBLOCKED( W, NKPTS_IRZ, GW, GAM_NEG, IMAG_GRIDS )
      USE scala
      USE constant
      IMPLICIT NONE
      TYPE (wavespin)         :: W
      INTEGER                 :: NKPTS_IRZ
      COMPLEX(q)              :: GW(:,:,:,:)        !diagonal part of interacting Green function G
      COMPLEX(q)                    :: GAM_NEG(:,:,:)    !diagonal part of correlation contr. to density matrix
! for each state, k-point and spin index
      TYPE( imag_grid_handle) :: IMAG_GRIDS         ! imaginary grids
 
! local variables
      INTEGER            :: I                  ! frequency points
      INTEGER            :: NOMEGA_GLOBAL      ! corresponing global frequency point
      INTEGER            :: NK1                !current q-point
      INTEGER            :: ISP                !current spin
      COMPLEX(q), ALLOCATABLE  :: GPHI(:)
      COMPLEX(q), ALLOCATABLE  :: GPSI(:)
 
! initialize negative time gamma matrix (electrons)
      GAM_NEG=0
 
      ALLOCATE(GPHI(SIZE(GW,1)))
      ALLOCATE(GPSI(SIZE(GW,1)))
      
      DO ISP=1,W%WDES%ISPIN
         DO NK1=1,NKPTS_IRZ
            GPHI = 0 
            GPSI = 0 
 
            DO I=1,SIZE(GW,4)  ! loop over local frequency points
! this loop is non-blocked
               NOMEGA_GLOBAL = I 
! transform of even part of density matrix
               GPHI(:) = GPHI(:) + GW(:, NK1, ISP, I)*&
                  IMAG_GRIDS%FER_RE_WEIGHT(NOMEGA_GLOBAL)/2
 
!#undef 
# 2377

            ENDDO
            GPSI = 0 
! compute psi term directly from CHAM
! The PSI term is the exchange part
            IF ( LFINITE_TEMPERATURE ) THEN
               IF ( LHOLEGF ) THEN
                  GPSI(1:W%WDES%NB_TOTK(NK1, ISP)) = +0.5_q
               ELSE
                  GPSI(1:W%WDES%NB_TOTK(NK1, ISP)) = -0.5_q
               ENDIF             
            ENDIF


 
            GPHI( : ) = GPHI( : ) + CONJG(GPHI( : ))

!
! Feynman propagator has this form in imaginary time:
!     G(t) = -( G_psi( t ) + sign( t ) * G_phi( t ) )
! where
!
!                              1  sinh( beta(x-mu)/2( 1-2*t))
!   G_psi(t) = sum_{x in poles}-*---------------------------
!                              2   cosh(beta*(x-mu)/2)
!
!                              1   cosh( beta(x-mu)/2( 1-2*t))
!   G_phi(t) = sum_{x in poles}-*----------------------------
!                              2     cosh(beta*(x-mu)/2)
!
!
! so the density matrix for electrons (t -> 0-) is
!
! GAM_NEG = G(0-) = G_psi(0) - G_phi( 0 )
!
! and for holes (t->0+) is
!
! GAM_POS = G(0+) = G_psi(0) + G_phi( 0 )
!
 
! electrons
            GAM_NEG(:, NK1, ISP)= GPHI( : ) 
               
! holes
!GAM_POS(:, NK1, ISP)=-GPHI( : )
 
! in case of metallic compound this term is present as well
            IF ( ASSOCIATED( IMAG_GRIDS%TO_TAU0 ) ) THEN
 
! electrons
               GAM_NEG(:, NK1, ISP ) = GAM_NEG( :, NK1, ISP) - GPSI(:)
! holes
!GAM_POS(:, NK1, ISP ) = GAM_POS( :, NK1, ISP) - GPSI( : )
 
            ENDIF
!IF (IO%IU0>=0 ) THEN
!WRITE(99,'(3I4,4F20.10)')(NK1,ISP,I, GPHI(I),GPSI(I), I=1,SIZE( GPSI ))
!ENDIF
         ENDDO
      ENDDO
 
      DEALLOCATE( GPHI )
      DEALLOCATE( GPSI )
    END SUBROUTINE GAMMA_FROM_G_W_DIAG_NONBLOCKED


    END SUBROUTINE ITERATE_EFERMI_FROM_SIGMA_DIAG

!*********************************************************************
!
! this subroutine determines and stores the mean field density matrix
! input/ output  CHAM_MAT
!
!  1/(i 2 pi) int_-infty^infty G_mean(w) e^ (i w eta) dw
!
! poles of states below the Fermi-level above real axis
!
! G_mean = ( w- H_mean) ^-1  ; H_mean = T+ V_H+ V_x
!
! with  H_mean = epsilon_HF  this simply becomes
!  Theta(e_fermi-epsilon_HF)
!
! TODO: old version
! with  H_mean = U epsilon_HF U+  this simply becomes
!  U Theta(e_fermi-epsilon_HF) U+
! where U is the matrix diagonalizing H_mean, and Theta is the
! step function
!
! if the Green function is passed down (GU_MAT and GO_MAT) the
! mean field Green function is also added to the two matrices
!
! Important note: this function must be exactly compatible to
! SUBTRACT_GHF_FROM_G taking off the mean field G before
! the frequency to time transformation
!
!*********************************************************************

   SUBROUTINE HF_DENSITY_MATRIX( W, NK1, ISP, CHAM_MAT, GDES_MAT, IMAG_GRIDS,&
     EFERMI, GU_MAT, GO_MAT, LUNOCCUPIED)
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)         :: W
     INTEGER                 :: NK1, ISP      ! k-point and spin index
     COMPLEX(q)                    :: CHAM_MAT(:)   ! in: frequency independent part of self-energy (T +V_x + V_H)
! out: mean field density matrix
     TYPE (greens_mat_des)   :: GDES_MAT      ! greens function matrix descriptor
     TYPE (imag_grid_handle) :: IMAG_GRIDS    ! imaginary grid descriptor
     REAL(q)                 :: EFERMI        ! Fermi-level
     COMPLEX(q), POINTER, OPTIONAL :: GU_MAT(:,:,:,:), GO_MAT(:,:,:,:) ! Green function in imaginary time
     LOGICAL, OPTIONAL  :: LUNOCCUPIED        ! obtain holes
! local variables
     INTEGER                 :: N, NB, I, NTAU_GLOBAL
     REAL(q)                 :: R( SIZE(W%CELTOT,1)), THETA(  SIZE(W%CELTOT,1))
     REAL(q)                 :: EMIN, EMAX
     INTEGER                 :: IERROR
     LOGICAL            :: LUNOCCUPIED_
      

! decide if electron or hole Green's function is determined
     LUNOCCUPIED_ = .FALSE.
     IF ( LUNOCCUPIED ) THEN
        LUNOCCUPIED_ = LUNOCCUPIED
     ENDIF

! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
# 2509


     N=  W%WDES%NB_TOTK(NK1, ISP)
     R= W%CELTOT(:, NK1, ISP)

! calculate G1= U Theta(e_fermi- epsilon_i)
     THETA=0
     DO NB=1,N
! at finite temperature states can be partially filled
        IF ( LFINITE_TEMPERATURE ) THEN
           IF ( LUNOCCUPIED_ ) THEN
              THETA(NB)=-(1-W%FERTOT(NB,NK1,ISP))
           ELSE
              THETA(NB)=-W%FERTOT(NB,NK1,ISP)
           ENDIF 
        ELSE
! larger then 0.5 treat as occupied, otherwise empty
           IF (W%FERTOT(NB,NK1,ISP)>FERMI_OCC) THEN
              IF ( LUNOCCUPIED_ ) THEN
                 THETA(NB)=0
              ELSE
                 THETA(NB)=1
              ENDIF 
           ELSE
              IF ( LUNOCCUPIED_ ) THEN
                 THETA(NB)=1
              ELSE
                 THETA(NB)=0
              ENDIF 
           ENDIF
        ENDIF
     ENDDO
# 2553

     CHAM_MAT=0
     CALL ADD_TO_DIAGONALE_REAL(N, CHAM_MAT, THETA,  GDES_MAT%DESC )

     
! add mean field Green function exp tau H_mean to
! the correlation part of the Green function
! this is simply  U exp (tau (e-efermi)) U+
     IF (PRESENT(GO_MAT) .AND. .NOT. LFINITE_TEMPERATURE) THEN
        IF (ASSOCIATED(GO_MAT)) THEN
           DO I=1,SIZE(GO_MAT,4)
              NTAU_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS )
! occupied Green function
! calculate G1= U exp( tau (e-EFERMI))
              THETA=0
              DO NB=1, N
! for finite temperatures this is the negative time GF
                 IF ( LFINITE_TEMPERATURE ) THEN
                    THETA(NB) = GF_WEIGHT( .TRUE.,  &
                       IMAG_GRIDS%TAU(NTAU_GLOBAL), &
                       IMAG_GRIDS%T%BETA, &
                       EFERMI, R(NB), W%FERTOT( NB, NK1, ISP ) )
                 ELSE
                    IF (W%FERTOT(NB,NK1,ISP)>FERMI_OCC) THEN
                       IF (R(NB)-EFERMI<0.0_q) THEN
                          THETA(NB)=EXP(IMAG_GRIDS%TAU(NTAU_GLOBAL)*(R(NB)-EFERMI))
                       ELSE
                          THETA(NB)=EXP(-IMAG_GRIDS%TAU(NTAU_GLOBAL)*IMAG_GRIDS%X1)
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
# 2599

              CALL ADD_TO_DIAGONALE_REAL(N, GO_MAT(:, NK1, ISP, I), THETA,  GDES_MAT%DESC )


! unoccupied Green function
! calculate G1= U exp( tau (EFERMI-e))
              THETA=0
              DO NB=1,N
! for finite temperatures this is the positive time GF
                 IF ( LFINITE_TEMPERATURE ) THEN
                    THETA(NB) = GF_WEIGHT( .FALSE.,  &
                       IMAG_GRIDS%TAU(NTAU_GLOBAL), &
                       IMAG_GRIDS%T%BETA, &
                       EFERMI, R(NB), W%FERTOT( NB, NK1, ISP ) )
                 ELSE
                    IF (W%FERTOT(NB,NK1,ISP)<=FERMI_OCC) THEN
                       IF (R(NB)-EFERMI>=0.0_q) THEN
                          THETA(NB)=EXP(IMAG_GRIDS%TAU(NTAU_GLOBAL)*(EFERMI-R(NB)))
                       ELSE
                          THETA(NB)=EXP(-IMAG_GRIDS%TAU(NTAU_GLOBAL)*IMAG_GRIDS%X1)
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
# 2636

              CALL ADD_TO_DIAGONALE_REAL(N, GU_MAT(:, NK1, ISP, I), THETA,  GDES_MAT%DESC )


           ENDDO
        ENDIF
     ENDIF

! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
# 2649


!     CALL DUMP_GREENS_HAM( "G1", G1(:), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6)
!     CALL DUMP_GREENS_HAM( "G2", G2(:), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6)
!     CALL DUMP_GREENS_HAM( "CHAM", CHAM_MAT(:), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, 6)
      
        
   END SUBROUTINE HF_DENSITY_MATRIX

!*********************************************************************
!
! determine HF Hamiltonian from unitary rotation matrix
! and supplied eigenvalues
!
!*********************************************************************


   SUBROUTINE HF_MATRIX( W, NK1, ISP, CHAM_MAT, GDES_MAT )
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)    :: W
     INTEGER            :: NK1, ISP           ! k-point and spin index
     COMPLEX(q)               :: CHAM_MAT(:)        ! in: frequency independent part of self-energy (T +V_x + V_H)
! out: Hartree-Fock density matrix
     TYPE (greens_mat_des) :: GDES_MAT        ! descriptor
! local variables
     INTEGER            :: N
     REAL(q)            :: AUX( SIZE(W%CELTOT,1))
# 2680

     INTEGER :: IERROR

     N=  W%WDES%NB_TOTK(NK1, ISP)
! if diag_cham is set, the self-energy is always stored in the NO basis
! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
!   the second option is much more efficient
# 2702

! TODO remove G1 and G2
     AUX(1:N)=W%CELTOT(1:N,NK1,ISP)
     CHAM_MAT=0
     CALL ADD_TO_DIAGONALE_REAL(N, CHAM_MAT, AUX,  GDES_MAT%DESC )


   END SUBROUTINE HF_MATRIX

!*********************************************************************
!
! determine rotation matrix from density matrix using
! perturbation theory
! this routine can be used to diagonalize any matrix
! in principle
!
!*********************************************************************

   SUBROUTINE ITERATIVE_DIAG( N, NK1, ISP, GAMMA, R, GDES_MAT, ITER)
     USE scala
     USE constant
     IMPLICIT NONE
     INTEGER            :: N                  ! matrix size
     INTEGER            :: NK1, ISP           ! k-point and spin index
     COMPLEX(q)               :: GAMMA(:)           ! in: frequency independent part of self-energy (T +V_x + V_H)
! out: Hartree-Fock density matrix
     COMPLEX(q)               :: R(:)               ! approximate eigenvalues of GAMMA
     TYPE (greens_mat_des) :: GDES_MAT        ! descriptor
     INTEGER :: ITER                          ! number of iterations performed for rotation matrix
! local variables
     INTEGER :: I
     COMPLEX(q)               :: DIAGONAL( N )
     COMPLEX(q), ALLOCATABLE  :: ROT_CURRENT(:), TMP(:), ROT(:)
     INTEGER :: IERROR

     ALLOCATE(ROT_CURRENT(SIZE(GAMMA,1)), TMP(SIZE(GAMMA,1)), ROT(SIZE(GAMMA,1)))

! initialize rotation matrix to unit matrix 1
     ROT=0
     CALL ADD_TO_DIAGONALE_REAL(N, ROT, (/ (1.0_q, I=1,N) /), GDES_MAT%DESC )

! iterate rotation matrix
     DO I=1,ITER
 
        ROT_CURRENT=GAMMA
! get diagonal elements of density matrix
        CALL DETERMINE_DIAGONALE_GDEF(N, ROT_CURRENT(:), DIAGONAL, GDES_MAT%DESC)
        CALL M_sum_z(GDES_MAT%COMM_INTRA, DIAGONAL, N)

! perturbation theory to obtain an approximately unitary rotation matrix
        CALL DIVIDE_BY_EIGENVALUEDIFF_SIGNED( N, ROT_CURRENT(:), DIAGONAL,&
           GDES_MAT%DESC, LFINITE_TEMPERATURE)
! make it exactly unitary by Gramm-Schmidt orthogonalization
        CALL MAKE_UNITARY_DESC( ROT_CURRENT(:), N, GDES_MAT%DESC)

        IF (I<ITER) THEN
! except for last iteration, rotate density matrix by multiplication
! TMP= U+ GAMMA
           CALL PZGEMM( 'C', 'N', N, N, N, (1._q,0._q),  &
                ROT_CURRENT(1), 1, 1, GDES_MAT%DESC, &
                GAMMA(1), 1, 1, GDES_MAT%DESC, &
                (0._q,0._q),  &
                TMP(1), 1, 1, GDES_MAT%DESC )
! GAMMA = TMP U
           CALL PZGEMM( 'N', 'N', N, N, N, (1._q,0._q),  &
                TMP(1), 1, 1, GDES_MAT%DESC, &
                ROT_CURRENT(1), 1, 1, GDES_MAT%DESC, &
                (0._q,0._q),  &
                GAMMA(1), 1, 1, GDES_MAT%DESC )
        ELSE
           CALL DETERMINE_DIAGONALE_GDEF(N, GAMMA(:), R, GDES_MAT%DESC)
           CALL M_sum_z(GDES_MAT%COMM_INTRA, R, N )
        ENDIF
! accumulate rotation matrix
        CALL PZGEMM( 'N', 'N', N, N, N, (1._q,0._q),  &
             ROT(1), 1, 1, GDES_MAT%DESC, &
             ROT_CURRENT(1), 1, 1, GDES_MAT%DESC, &
             (0._q,0._q),  &
             TMP(1), 1, 1, GDES_MAT%DESC )
        ROT=TMP
     ENDDO
     
! return final rotation matrix
     GAMMA=ROT

     DEALLOCATE(ROT_CURRENT, TMP, ROT)

   END SUBROUTINE

!*********************************************************************
!
! this subroutine calculates the HF singles contribution
! i.e. the change of the mean field energy when a unitary
! transformation of the orbitals is performed
! this is essentially the Trace of the occupied part of
! the HF matrix before and after the unitary rotation
!
!  Tr [F U+ H U] - Tr [F H]
!
! CHAM_MAT stores the HF mean field matrix
! upon return it is rotated U+ H U
!
!*********************************************************************


   SUBROUTINE HF_SINGLES( W, NK1, ISP, CUNI, CHAM_MAT, GDES_MAT, EBANDSTR )
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)    :: W
     INTEGER            :: NK1, ISP           ! k-point and spin index
     COMPLEX(q)               :: CHAM_MAT(:)        ! in: Hartree-Fock mean field Hamiltonian
! out: rotated Hamiltonian
     COMPLEX(q)               :: CUNI    (:)        ! in: unitary transformation
     TYPE (greens_mat_des) :: GDES_MAT        ! descriptor for matrices
     REAL(q)            :: EBANDSTR
! local variables
     INTEGER :: NBANDSK, N
     COMPLEX(q), ALLOCATABLE  :: CTMP(:)
     COMPLEX(q)               :: HFEIG(W%WDES%NB_TOT) ! HF eigenvalues
     INTEGER :: IERROR

     
     ALLOCATE(CTMP(SIZE(CHAM_MAT,1)))

     NBANDSK=  W%WDES%NB_TOTK(NK1, ISP)

! TMP = U+ H
     CALL PZGEMM( 'C', 'N', NBANDSK, NBANDSK, NBANDSK, (1._q,0._q), &
          CUNI(1), 1, 1, GDES_MAT%DESC, &
          CHAM_MAT(1), 1, 1, GDES_MAT%DESC, &
          (0._q,0._q),  &
          CTMP(1), 1, 1, GDES_MAT%DESC )
! H' = TMP U
     CALL PZGEMM( 'N', 'N', NBANDSK, NBANDSK, NBANDSK, (1._q,0._q), &
          CTMP(1), 1, 1, GDES_MAT%DESC, &
          CUNI(1), 1, 1, GDES_MAT%DESC, &
          (0._q,0._q),  &
          CHAM_MAT(1), 1, 1, GDES_MAT%DESC )

     CALL DETERMINE_DIAGONALE_GDEF( W%WDES%NB_TOTK(NK1,ISP), CHAM_MAT(:), HFEIG(:), GDES_MAT%DESC)
     CALL M_sum_z(GDES_MAT%COMM_INTRA, HFEIG(1), SIZE(HFEIG,1))

! subtract
     DO N=1,W%WDES%NB_TOTK(NK1,ISP)
        EBANDSTR=EBANDSTR+W%WDES%RSPIN* REAL( HFEIG(N) ,KIND=q)*W%FERTOT(N,NK1,ISP)*W%WDES%WTKPT(NK1)
     ENDDO
     
     DEALLOCATE( CTMP)
        
     
   END SUBROUTINE HF_SINGLES

!*********************************************************************
!
! this routine performs an integration/summation test
! and checks wether the number of electrons is conserved
! especially if following identity is satisfied:
! Tr ( G_pos(t=0) ) =  N_elect
! Tr ( G_neg(t=0) ) =  N_holes
!
!*********************************************************************

   SUBROUTINE PERFORM_INTEGRATION_TEST( W, SIGMAW_MAT, IMAG_GRIDS, GDES_MAT,&
      NKPTS_IRZ, EFERMI, E, IO )  
      USE scala 
      USE tutor, ONLY: vtutor
      TYPE (wavespin)         :: W                   ! wave function

      COMPLEX(q)              :: SIGMAW_MAT(:,:,:,:) ! frequency dependent self-energy
! in orbital basis

      TYPE (imag_grid_handle) :: IMAG_GRIDS          ! imaginary grids handle

      TYPE (greens_mat_des), POINTER :: GDES_MAT     ! descriptor for
! functions in orbital basis
      INTEGER                 :: NKPTS_IRZ
      REAL(q)                 :: EFERMI(:)
      TYPE (energy)           :: E
      TYPE( in_struct )       :: IO 
! local variables
      COMPLEX(q), POINTER, CONTIGUOUS :: GAM_NEG(:,:,:,:)=>NULL() !density matrix (tau<0)
      COMPLEX(q), POINTER, CONTIGUOUS :: GAM_POS(:,:,:,:)=>NULL() !density matrix (tau>0)
      COMPLEX(q), POINTER, CONTIGUOUS :: GF_MAT(:,:,:,:)=>NULL()!green function in orbital basis
      INTEGER    :: ISP, NK1, N, I
      INTEGER    :: NOMEGA_GLOBAL
      COMPLEX(q)       :: TR_GNEG( W%WDES%ISPIN )
      COMPLEX(q)       :: TR_GPOS( W%WDES%ISPIN )
      COMPLEX(q)       :: TR_GALL( W%WDES%ISPIN )
      REAL(q)    :: MU_I( W%WDES%ISPIN )  ! chemical potential
      INTEGER, PARAMETER :: MAX_ITER = 1000 
      REAL(q), PARAMETER :: CHEM_SHIFT = 30._q

! density matrix for states propagating negatively in time
      CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAM_NEG, NKPTS_IRZ, W%WDES%ISPIN, 1)
      GAM_NEG=(0._q,0._q)

! density matrix for states propagating negatively in time
      CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAM_POS, NKPTS_IRZ, W%WDES%ISPIN, 1)
      GAM_POS=(0._q,0._q)

! Green's function matrix in imaginary frequency domain
      CALL ALLOCATE_GREENS_MAT(GDES_MAT, GF_MAT, NKPTS_IRZ, W%WDES%ISPIN,&
       IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP)
      GF_MAT=0
 
! scan for different fermi occupancies
      IF ( IO%IU0>=0 ) WRITE(IO%IU0,10) EFERMI 
      IF ( IO%IU6>=0 ) OPEN( UNIT=100, FILE='dos.dat', STATUS='REPLACE' )

      DO I = 1, MAX_ITER      
! N_el( mu_i ) is calculated for chemical potential mu_i
! with mu_i in [ EFERMI - CHEM_SHIFT, EFERMI + CHEM_SHIFT ]
         MU_I(:) = EFERMI(:)-CHEM_SHIFT + ((I-1)*(2*CHEM_SHIFT))/(MAX_ITER-1)
! test
! set self-energy to Green's function in imaginary frequency
! TEST
!         CALL SET_TO_DIAGONAL_GF( MU_I )

! determine green's function at mu_i, that is sum_{nk} 1/(iw - eps_nk + mu_i
         CALL SET_NON_INTERACTING_PROPAGATOR( SIGMAW_MAT, MU_I, IMAG_GRIDS )

!         CALL CALCULATE_G_FROM_SIGMA( W, NKPTS_IRZ, MU_I, SIGMAW_MAT, &
!             IMAG_GRIDS, GDES_MAT, GF_MAT, E, IO )
   
! calculate density matrix als upper and lower limit t->0
         CALL GAM_POS_GAM_NEG_FROM_G_W( W, NKPTS_IRZ, SIGMAW_MAT, GAM_POS, GAM_NEG, &
              IMAG_GRIDS, GDES_MAT )
      
         CALL OBTAIN_NUMBER_OF_STATES( GAM_NEG, GAM_POS, TR_GNEG, TR_GPOS, TR_GALL )

         IF ( IO%IU0>=0 ) & 
         WRITE(IO%IU0,11) ( MU_I(ISP), TR_GNEG(ISP), TR_GPOS(ISP), TR_GALL(ISP), ISP=1,W%WDES%ISPIN)

         IF ( IO%IU6>=0 ) & 
         WRITE(100,11) ( MU_I(ISP), TR_GNEG(ISP), TR_GPOS(ISP), TR_GALL(ISP), ISP=1,W%WDES%ISPIN)
 
      ENDDO  !mu step
      IF ( IO%IU6>=0 ) CLOSE( 100 )

      CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAM_NEG)
      CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAM_POS)
      CALL DEALLOCATE_GREENS_MAT(GDES_MAT, GF_MAT)

      CALL vtutor%stopCode()

10    FORMAT( "E-fermi: ", 2E14.6,/, " mu_i(up), N_el( mu_i ), N_h( mu_i ), N_tot( mu_i ), mu_i(dn) ..." )
11    FORMAT( 14E14.6 )


      CONTAINS 

      SUBROUTINE SET_TO_DIAGONAL_GF( MU_I ) 
         INTEGER :: ISP, NK1, N, I, NOMEGA_GLOBAL 
         REAL(q) :: MU_I(:)
         REAL(q) :: MU
         COMPLEX(q) :: AUX( SIZE(W%CELTOT,1))

         DO ISP=1,W%WDES%ISPIN
            MU = MU_I(ISP)
            DO NK1=1,NKPTS_IRZ
               N=  W%WDES%NB_TOTK(NK1, ISP)
               DO I=1,SIZE(SIGMAW_MAT,4)  ! loop over local frequency points
!this is the global tau point of the node
                  NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
                  IMAG_GRIDS)
                  SIGMAW_MAT(:, NK1, ISP, I)=0
!  calculate i w +mu_x -(T +V_x + V_H)
                  AUX(1:N) = 0
                  AUX(1:N) = CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q)
                  AUX(1:N) = AUX(1:N) - W%CELTOT(1:N,NK1,ISP)+MU_I(ISP)
                  AUX(1:N) = 1/AUX(1:N) 
                  CALL ADD_TO_DIAGONALE(N, SIGMAW_MAT(:, NK1, ISP, I), AUX, GDES_MAT%DESC)
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE SET_TO_DIAGONAL_GF

      SUBROUTINE SET_NON_INTERACTING_PROPAGATOR( GF_MAT, MU_I, IMAG_GRIDS,&
         LCONSTANT ) 
         COMPLEX(q)              :: GF_MAT(:,:,:,:)
         REAL(q)                 :: MU_I(:)
         TYPE (imag_grid_handle) :: IMAG_GRIDS          ! imaginary grids handle
         LOGICAL, OPTIONAL       :: LCONSTANT
! local
         INTEGER :: ISP, NK1, N, I, NOMEGA_GLOBAL 
         REAL(q) :: MU
         COMPLEX(q) :: AUX( SIZE(W%CELTOT,1))
         LOGICAL    :: LCONSTANT_
         
         LCONSTANT_=.FALSE.
         IF ( PRESENT( LCONSTANT) ) THEN
            LCONSTANT_ = LCONSTANT 
         ENDIF

         GF_MAT = 0
         DO ISP=1,W%WDES%ISPIN
            MU = MU_I(ISP)
            DO NK1=1,NKPTS_IRZ
               N=  W%WDES%NB_TOTK(NK1, ISP)
               DO I=1,SIZE(GF_MAT,4)  ! loop over local frequency points
!this is the global tau point of the node
                  NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
                  IMAG_GRIDS)
                  IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE

!  calculate i w +mu_x -(T +V_x + V_H)
                  AUX(1:N) = 0
                  IF ( LCONSTANT_ ) THEN
                     AUX(1:N) = 0.5_q
                  ELSE
                     AUX(1:N) = CMPLX(0._q,IMAG_GRIDS%FER_RE(NOMEGA_GLOBAL),q)
                     AUX(1:N) = AUX(1:N) - W%CELTOT(1:N,NK1,ISP)+MU_I(ISP)
                     AUX(1:N) = 1/AUX(1:N) 
                  ENDIF
                  CALL ADD_TO_DIAGONALE(N, GF_MAT(:, NK1, ISP, I), AUX, GDES_MAT%DESC)
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE SET_NON_INTERACTING_PROPAGATOR

      SUBROUTINE OBTAIN_NUMBER_OF_STATES( GAM_NEG, GAM_POS, TR_GNEG, TR_GPOS, TR_GALL )
         COMPLEX(q)       :: GAM_POS(:,:,:,:) ! density matrix for holes
         COMPLEX(q)       :: GAM_NEG(:,:,:,:) ! density matrix for eletrons
         COMPLEX(q)       :: TR_GNEG(:)
         COMPLEX(q)       :: TR_GPOS(:)
         COMPLEX(q)       :: TR_GALL(:)

         COMPLEX(q)       :: R( SIZE(W%CELTOT,1))

         DO ISP=1,W%WDES%ISPIN
            TR_GNEG(ISP) = 0
            TR_GPOS(ISP) = 0
            TR_GALL(ISP) = 0
            DO NK1=1,NKPTS_IRZ
! number of electrons from FERMI function
               TR_GALL(ISP) = TR_GALL(ISP) + SUM(W%FERTOT(1:W%WDES%NB_TOTK(NK1,ISP),NK1,ISP))*&
                  KPOINTS_ORIG%WTKPT(NK1)*W%WDES%RSPIN

! obtain trace of G_NEG
               CALL DETERMINE_DIAGONALE_GDEF(W%WDES%NB_TOTK(NK1,ISP), &
                  GAM_NEG(:,NK1,ISP,1), R, GDES_MAT%DESC)
               CALL M_sum_z(GDES_MAT%COMM_INTRA, R, W%WDES%NB_TOTK(NK1, ISP) )
               TR_GNEG(ISP) = TR_GNEG(ISP) + SUM(R(1:W%WDES%NB_TOTK(NK1,ISP)))*&
                  KPOINTS_ORIG%WTKPT(NK1)*W%WDES%RSPIN

! obtain trace of G_POS
               CALL DETERMINE_DIAGONALE_GDEF(W%WDES%NB_TOTK(NK1,ISP), &
                  GAM_POS(:,NK1,ISP,1), R, GDES_MAT%DESC)
               CALL M_sum_z(GDES_MAT%COMM_INTRA, R, W%WDES%NB_TOTK(NK1, ISP) )
               TR_GPOS(ISP) = TR_GPOS(ISP) + SUM(R(1:W%WDES%NB_TOTK(NK1,ISP)))*&
                  KPOINTS_ORIG%WTKPT(NK1)*W%WDES%RSPIN
            ENDDO
         ENDDO
      
     END SUBROUTINE OBTAIN_NUMBER_OF_STATES


   END SUBROUTINE PERFORM_INTEGRATION_TEST


!*********************************************************************
!
! merge diagonal elements of self-energy from all nodes
! and write it to file
!
!*********************************************************************

   SUBROUTINE WRITE_SELFENERGY_W( W, NKPTS_IRZ, SIGMAW_MAT, IMAG_GRIDS, GDES_MAT, IU0) 
     USE scala
     USE constant
     IMPLICIT NONE
     TYPE (wavespin)         :: W
     INTEGER                 :: NKPTS_IRZ
     COMPLEX(q),POINTER      :: SIGMAW_MAT(:,:,:,:)  !in self-energy in frequency space
     TYPE (imag_grid_handle) :: IMAG_GRIDS       ! imaginary grid descriptor
     TYPE (greens_mat_des)   :: GDES_MAT          ! in: descriptor
     INTEGER                 :: IU0                  ! unit for dump
! local variables
     INTEGER                 :: NOMEGA_GLOBAL
     INTEGER                 :: ISP, NK1, N, I
     COMPLEX(q),ALLOCATABLE  :: R(:,:)
!COMPLEX(q)         :: R(W%WDES%NB_TOT, NOMEGA)
     
     ALLOCATE( R( W%WDES%NB_TOT, IMAG_GRIDS%NOMEGA ) )
     R=(0._q,0._q)

     IF (IU0>=0) OPEN( 200, FILE="SIGMA.DAT", STATUS="REPLACE" )
   
     DO ISP=1,W%WDES%ISPIN
        DO NK1=1,NKPTS_IRZ
           IF (IU0>=0) WRITE(200,'(A,I3,A,I5,3F14.6)')'Spin',ISP,'k-point:',NK1, W%WDES%VKPT(:,NK1)
           R=0
           DO I=1,SIZE(SIGMAW_MAT,4)
              NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
              IMAG_GRIDS)
              IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE
              CALL DETERMINE_DIAGONALE(W%WDES%NB_TOTK(NK1,ISP), SIGMAW_MAT(:,NK1,ISP,I), R(:,NOMEGA_GLOBAL), GDES_MAT%DESC)
           ENDDO
! sum over all nodes (between tau and in tau)
           CALL M_sum_z(GDES_MAT%COMM, R, SIZE(R))
           IF (IU0>=0) THEN
              DO N=1,MIN(W%WDES%NB_TOTK(NK1,ISP),32)
                 WRITE(200,'(A,I6)')'Band-no.:',N
                 DO I=1,IMAG_GRIDS%NOMEGA
                    WRITE(200,'(3F20.7)') IMAG_GRIDS%FER_RE(I), R(N, I)
                 ENDDO
                 WRITE(200,*)
              ENDDO
           END IF
        ENDDO
     ENDDO

     IF (IU0>=0) CLOSE(200)

     DEALLOCATE( R )
     
   END SUBROUTINE WRITE_SELFENERGY_W

!***********************************************************************
!
! determine diagonal of self-energy in HF basis
! write self-energy and QP energies from Pade fit
!
!***********************************************************************

  SUBROUTINE PADE_FIT_GW0( W, WMEAN, IMAG_GRIDS, GDES_MAT, SIGMAW_MAT, NKPTS_IRZ, NBANDSGW, IO )
    USE pade_fit, ONLY: QP_PADE_FIT_FROM_TAU
    USE scala, ONLY: DETERMINE_DIAGONALE
    TYPE(wavespin)                 :: W          !wave function
    TYPE(wavespin)                 :: WMEAN      !identical to W, except for CELTOT and FERTOT
    TYPE (imag_grid_handle)        :: IMAG_GRIDS ! imaginary grids handle
    TYPE (greens_mat_des), POINTER :: GDES_MAT   ! out: pointer to generated descriptor
    COMPLEX(q)                     :: SIGMAW_MAT(:,:,:,:)
    INTEGER                        :: NKPTS_IRZ
    INTEGER                        :: NBANDSGW
    TYPE(in_struct)                :: IO
!Local
    INTEGER :: ISP, NK, I, NOMEGA_GLOBAL
    COMPLEX(q), POINTER :: SIGMAO_TAU(:,:,:,:)     ! self-energy in time (occupied)
    COMPLEX(q), POINTER :: SIGMAU_TAU(:,:,:,:)     ! self-energy in time (unoccupied)
    COMPLEX(q), POINTER :: SIGMAW_DIA(:,:,:,:)  ! diagonal part of SIGMAW
    COMPLEX(q) :: R(W%WDES%NB_TOT)

    ALLOCATE(SIGMAW_DIA(NOMEGA,W%WDES%NB_TOT,NKPTS_IRZ,W%WDES%ISPIN))
    ALLOCATE(SIGMAO_TAU(NOMEGA,W%WDES%NB_TOT,NKPTS_IRZ,W%WDES%ISPIN))
    ALLOCATE(SIGMAU_TAU(NOMEGA,W%WDES%NB_TOT,NKPTS_IRZ,W%WDES%ISPIN))
    SIGMAW_DIA=(0._q,0._q)
    SIGMAO_TAU=(0._q,0._q)
    SIGMAU_TAU=(0._q,0._q)
    IF (IO%IU6>=0) WRITE(IO%IU6,'(1X,A)', ADVANCE="NO")"QP shifts evaluated in HF basis"

    DO ISP=1,W%WDES%ISPIN
       DO NK=1,NKPTS_IRZ
          DO I=1,SIZE(SIGMAW_MAT,4)  ! loop over local frequency points
             NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, I,&
             IMAG_GRIDS)
             IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE
             R=(0._q,0._q)
             CALL DETERMINE_DIAGONALE(W%WDES%NB_TOT, SIGMAW_MAT(:,NK,ISP,I), R, GDES_MAT%DESC)
! copy all diagonal entries at once
             SIGMAW_DIA( NOMEGA_GLOBAL, :, NK, ISP ) = R
          ENDDO

       ENDDO
    ENDDO
!
! collective sum of diagonal entries, every core has the complete frequency dependent array
    CALL M_sum_z(GDES_MAT%COMM, SIGMAW_DIA, SIZE(SIGMAW_DIA))
!
! inverse transformation to time is needed here for pade fit below
!
    CALL INV_FT_G_OR_SIGMA(W, SIGMAW_DIA, SIGMAO_TAU, SIGMAU_TAU, IMAG_GRIDS, NKPTS_IRZ )
    CALL QP_PADE_FIT_FROM_TAU(W,WMEAN,IMAG_GRIDS,SIGMAU_TAU,SIGMAO_TAU,NKPTS_IRZ,NBANDSGW,IO)

!CALL WRITE_SELF_ENERGY(WMEAN, SIGMA_COS, SIGMA_SIN, NKPTS_IRZ,IMAG_GRIDS, NBANDSGW, IO%IU6,NELM,'HF basis')
!    CALL WRITE_SELF_ENERGY(WMEAN, SIGMA_COS, SIGMA_SIN, NKPTS_IRZ,IMAG_GRIDS, NBANDSGW, IO%IU6)
!    CALL QP_PADE_FIT(W,WMEAN,IMAG_GRIDS%FER_RE,SIGMA_COS,SIGMA_SIN,NKPTS_IRZ,NBANDSGW,IO%IU6,IO%IU0)
   
    DEALLOCATE( SIGMAW_DIA )
    DEALLOCATE( SIGMAO_TAU , SIGMAU_TAU)

  END SUBROUTINE PADE_FIT_GW0

END MODULE greens_orbital
