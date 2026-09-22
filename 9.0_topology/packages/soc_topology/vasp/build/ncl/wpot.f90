# 1 "wpot.F"
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


# 2 "wpot.F" 2 
!*********************************************************************
!
!> This module reads and writes the screened potential
!>
!> Constructing the potential at q points that are not in
!> the IRZ is managed through a handle wpothandle
!
!*********************************************************************

MODULE wpot
  USE chi_base
  USE wave_high
  USE kpoints_change
  USE minimax_struct, ONLY: loop_des
  USE mpimy
  IMPLICIT NONE

  TYPE wpothandle
     LOGICAL, POINTER :: LSET(:) =>NULL()
     INTEGER :: IU0, IU6
     TYPE (wavedes), POINTER ::       WGW =>NULL()
     TYPE (skpoints_trans),POINTER :: KPOINTS_TRANS =>NULL()
     COMPLEX(q), POINTER :: C(:,:) => NULL()
     TYPE (responsefunction), POINTER :: WPOT(:) => NULL()
     INTEGER :: NKREDLFX, NKREDLFY, NKREDLFZ

! frequency distribution descriptor
     TYPE( loop_des ) :: O
     INTEGER :: ID_NOMEGAPAR  = -1        !< which OMEGAPAR group is treated currently
     INTEGER :: NOMEGA_LOCAL  = -1        !< which OMEGA point is treated locally
     INTEGER :: NOMEGA_GLOBAL = -1
     
     INTEGER :: NQ_IN_IRZ = 0        !< number of q-points in irreducible BZ
     INTEGER, POINTER :: NQ_DIST(:,:) =>NULL()!< distribution of q-points
     REAL(q), POINTER :: OMEGA(:)     =>NULL()!< frequency points
     REAL(q), POINTER :: OMEGA_TAU_WEIGHT(:)  =>NULL()!< fourier transformation from frequency to time
     INTEGER :: NOMEGA    = -1       !< number of frequency points
     INTEGER :: NQ_LOCAL = 0         !< number of local frequency points
     REAL(q) :: FSG0=0               !< HF singularity correction

     INTEGER :: NQ_IN_FULL = 0       !< current q-point in full BZ
     LOGICAL :: LFULL = .FALSE.      !< is irz. BZ -> full BZ  necessary

  END TYPE wpothandle

  CONTAINS

!*********************************************************************
!
!> initialize the potential handle
!
!*********************************************************************

  SUBROUTINE INIT_WPOT_HANDLE(WH, WGW, NQ , IU6, IU0, NKREDLFX, NKREDLFY, NKREDLFZ)
    USE kpoints_change
    IMPLICIT NONE
    TYPE(wpothandle), POINTER :: WH
    TYPE (wavedes), TARGET :: WGW
    INTEGER :: NQ
    INTEGER :: IU6, IU0

    INTEGER :: NKREDLFX, NKREDLFY, NKREDLFZ
! local
    INTEGER :: NK

    ALLOCATE(WH)
    ALLOCATE(WH%LSET(NQ), WH%C(WGW%NGDIM, NQ), WH%KPOINTS_TRANS)
    NULLIFY(WH%WPOT)
    WH%WGW=> WGW
    WH%IU6=IU6
    WH%IU0=IU0
    WH%C=0
    WH%LSET=.FALSE.
    WH%NKREDLFX=NKREDLFX
    WH%NKREDLFY=NKREDLFY
    WH%NKREDLFZ=NKREDLFZ
    CALL GENERATE_KPOINTS_TRANS(WGW%GRID, KPOINTS_ORIG%NKPTS, WGW, KPOINTS_FULL_ORIG, WH%KPOINTS_TRANS)

! allocate the potential handle
    ALLOCATE(WH%WPOT(KPOINTS_FULL_ORIG%NKPTS))
    DO NK=1,KPOINTS_FULL_ORIG%NKPTS
       NULLIFY(WH%WPOT(NK)%RESPONSEFUN)
    ENDDO
    
  END SUBROUTINE INIT_WPOT_HANDLE

!*********************************************************************
!
!> initializes the potential handle,
!> compatible if POTENTIAL is distributed amound NOMEGAPAR groups in
!> a column major fashion
!
!*********************************************************************

  SUBROUTINE INIT_WPOT_HANDLE_DISTRIBUTED(WH, IMAG_GRIDS, NUMBER_OF_NQ , FSG0, WGW, NQ,&
    IU6, IU0, NKREDLFX, NKREDLFY, NKREDLFZ)
    USE kpoints_change
    USE minimax_struct,  ONLY: imag_grid_handle
    USE mathtools, ONLY: INVERT_REAL_MATRIX
    IMPLICIT NONE
    TYPE(wpothandle), POINTER :: WH
    TYPE( imag_grid_handle )  :: IMAG_GRIDS
    INTEGER, INTENT(IN) :: NUMBER_OF_NQ
    REAL(q)                   :: FSG0
    TYPE (wavedes), TARGET    :: WGW
    INTEGER :: NQ
    INTEGER :: IU6, IU0

    INTEGER :: NKREDLFX, NKREDLFY, NKREDLFZ
! local
    INTEGER :: NK
    INTEGER :: I,J
    REAL(q) :: FT(IMAG_GRIDS%T%NPOINTS,IMAG_GRIDS%B%NPOINTS)

    ALLOCATE(WH)
    ALLOCATE(WH%LSET(NQ), WH%C(WGW%NGDIM, NQ), WH%KPOINTS_TRANS)
    NULLIFY(WH%WPOT)
    WH%WGW=> WGW
    WH%IU6=IU6
    WH%IU0=IU0
    WH%C=0
    WH%LSET=.FALSE.
    WH%NKREDLFX=NKREDLFX
    WH%NKREDLFY=NKREDLFY
    WH%NKREDLFZ=NKREDLFZ
    CALL GENERATE_KPOINTS_TRANS(WGW%GRID, KPOINTS_ORIG%NKPTS, WGW, KPOINTS_FULL_ORIG, WH%KPOINTS_TRANS)

! allocate the potential handle, only two q-points necessary for distributed WPOT
! WPOT(1), WPOT(2) is used for acutual computations in dmft.F
    ALLOCATE(WH%WPOT(2))
    NULLIFY(WH%WPOT(1)%RESPONSEFUN)
    NULLIFY(WH%WPOT(2)%RESPONSEFUN)
    
! also add info about frequency distribution
    WH%O%NPOINTS               = IMAG_GRIDS%B%NPOINTS
    WH%O%NPOINTSC              = IMAG_GRIDS%B%NPOINTSC
    WH%O%NPOINTS_IN_GROUP      = IMAG_GRIDS%B%NPOINTS_IN_GROUP
    WH%O%NPOINTS_IN_ROOT_GROUP = IMAG_GRIDS%B%NPOINTS_IN_ROOT_GROUP
    WH%O%POINT_CURRENT         = IMAG_GRIDS%B%POINT_CURRENT
    WH%O%POINTMIN_CURRENT      = IMAG_GRIDS%B%POINTMIN_CURRENT
    WH%O%LDO_POINT_LOCAL       = IMAG_GRIDS%B%LDO_POINT_LOCAL
! copy all commuicators
    WH%O%COMM_BETWEEN_GROUPS   = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS
    WH%O%COMM_IN_GROUP         = IMAG_GRIDS%B%COMM_IN_GROUP
    WH%O%COMM                  = IMAG_GRIDS%B%COMM

! copy frequency points
    NULLIFY( WH%OMEGA )
    ALLOCATE( WH%OMEGA( IMAG_GRIDS%NOMEGA ) )
    WH%OMEGA = IMAG_GRIDS%BOS_RE

! construct integration weights that transform
! from frequency to time and integrate time representation
! to obtain the omega=0 value of WPOT
    FT = IMAG_GRIDS%TO_BOS_RE
    CALL INVERT_REAL_MATRIX(FT,-1)
    NULLIFY( WH%OMEGA_TAU_WEIGHT )
    ALLOCATE( WH%OMEGA_TAU_WEIGHT( IMAG_GRIDS%NOMEGA ) )
    DO J = 1, IMAG_GRIDS%B%NPOINTS
       WH%OMEGA_TAU_WEIGHT(J) = 0 
       DO I= 1, IMAG_GRIDS%T%NPOINTS
            WH%OMEGA_TAU_WEIGHT(J) =  WH%OMEGA_TAU_WEIGHT(J) + 2*IMAG_GRIDS%TAU_WEIGHT(I) * FT( I, J )  
       ENDDO
    ENDDO

    NULLIFY( WH%O%DISTRIBUTION )
    NULLIFY( WH%O%GROUP_OF_NODE )

    ALLOCATE( WH%O%DISTRIBUTION( IMAG_GRIDS%B%COMM_BETWEEN_GROUPS%NCPU , 2 ) )
    ALLOCATE( WH%O%GROUP_OF_NODE( WH%O%COMM%NCPU ) )
# 175

    WH%O%DISTRIBUTION = IMAG_GRIDS%B%DISTRIBUTION
    WH%O%GROUP_OF_NODE = IMAG_GRIDS%B%GROUP_OF_NODE 

    WH%ID_NOMEGAPAR  = -1 ! which OMEGAPAR group
    WH%NOMEGA_LOCAL  = -1 ! which local OMEGA point
    WH%NOMEGA_GLOBAL = -1 ! which global OMEGA point
    WH%NOMEGA = IMAG_GRIDS%NOMEGA 
 
! number of NQ in full BZ
    WH%NQ_IN_IRZ = NUMBER_OF_NQ
! set map of full BZ to IRZ
    NULLIFY( WH%NQ_DIST )
    ALLOCATE( WH%NQ_DIST( WGW%NB_PAR, KPOINTS_FULL_ORIG%NKPTS ) )
    WH%NQ_DIST(:,:) = 0
    WH%FSG0 = FSG0 
  END SUBROUTINE INIT_WPOT_HANDLE_DISTRIBUTED

!*********************************************************************
!
!> destroy the potential handle
!
!*********************************************************************

  SUBROUTINE DESTROY_WPOT_HANDLE(WH)
    IMPLICIT NONE
    TYPE(wpothandle), POINTER :: WH
! local
    INTEGER :: NK

    IF (.NOT. ASSOCIATED(WH)) RETURN

    CALL DEALLOCATE_KPOINTS_TRANS(WH%KPOINTS_TRANS)
    IF (ASSOCIATED(WH%WPOT)) THEN
! deallocate the potential handle
       DO NK=1,KPOINTS_FULL_ORIG%NKPTS
          IF (ASSOCIATED(WH%WPOT(NK)%RESPONSEFUN)) THEN
             CALL DEALLOCATE_RESPONSEFUN(WH%WPOT(NK) )
          ENDIF
       ENDDO
       DEALLOCATE(WH%WPOT)
       NULLIFY(WH%WPOT)
    ENDIF

    DEALLOCATE(WH%LSET, WH%C, WH%KPOINTS_TRANS)
    NULLIFY( WH%LSET, WH%C, WH%KPOINTS_TRANS)
    DEALLOCATE(WH)
   

  END SUBROUTINE DESTROY_WPOT_HANDLE

!*********************************************************************
!
!> destroy the distributed potential handle
!
!*********************************************************************

  SUBROUTINE DESTROY_DISTRIBUTED_WPOTH(WH)
    IMPLICIT NONE
    TYPE(wpothandle), POINTER :: WH
! local
    INTEGER :: NK

    IF (.NOT. ASSOCIATED(WH)) RETURN

    CALL DEALLOCATE_KPOINTS_TRANS(WH%KPOINTS_TRANS)
    IF (ASSOCIATED(WH%WPOT)) THEN
! deallocate the potential handle
       DO NK=1,2
          IF (ASSOCIATED(WH%WPOT(NK)%RESPONSEFUN)) THEN
             CALL DEALLOCATE_RESPONSEFUN(WH%WPOT(NK) )
          ENDIF
       ENDDO
       DEALLOCATE(WH%WPOT)
       NULLIFY(WH%WPOT)
    ENDIF

    DEALLOCATE(WH%LSET, WH%C, WH%KPOINTS_TRANS)

! get rid of all imaginary time grid info
    IF ( ASSOCIATED( WH%OMEGA) ) THEN
       DEALLOCATE( WH%OMEGA ) 
    ENDIF
    NULLIFY( WH%OMEGA ) 

    IF ( ASSOCIATED( WH%OMEGA_TAU_WEIGHT) ) THEN
       DEALLOCATE( WH%OMEGA_TAU_WEIGHT ) 
    ENDIF
    NULLIFY( WH%OMEGA_TAU_WEIGHT ) 

    IF ( ASSOCIATED( WH%O%DISTRIBUTION) ) THEN
       DEALLOCATE( WH%O%DISTRIBUTION ) 
    ENDIF
    NULLIFY( WH%O%DISTRIBUTION ) 

    IF ( ASSOCIATED( WH%O%GROUP_OF_NODE) ) THEN
       DEALLOCATE( WH%O%GROUP_OF_NODE) 
    ENDIF
    NULLIFY( WH%O%GROUP_OF_NODE ) 

    IF ( ASSOCIATED( WH%NQ_DIST) ) THEN
       DEALLOCATE( WH%NQ_DIST) 
    ENDIF
    NULLIFY( WH%NQ_DIST ) 

    WH%FSG0 = 0
    WH%ID_NOMEGAPAR  = -1 ! which OMEGAPAR group
    WH%NOMEGA_LOCAL  = -1 ! which local OMEGA point
    WH%NOMEGA_GLOBAL = -1 ! which global OMEGA point
    WH%NOMEGA = -1

    WH%O%NPOINTS               = -1
    WH%O%NPOINTSC              = -1
    WH%O%NPOINTS_IN_GROUP      = -1
    WH%O%NPOINTS_IN_ROOT_GROUP = -1
    WH%O%POINT_CURRENT         = -1
    WH%O%POINTMIN_CURRENT      = -1

    DEALLOCATE(WH)

  END SUBROUTINE DESTROY_DISTRIBUTED_WPOTH

!*********************************************************************
!
!> determine screened two electron potential
!> at an arbitrary q point supplied by the calling routine
!>
!> some comments:
!>@param WPOT is of the type of a response function array
!> a the gamma point the entries
!> are interpreted as sin and cosine transforms
!> and the potential is stored as a real valued function
!> in RESPONSER
!
!*********************************************************************

  SUBROUTINE GET_WPOT( WH, NQ_IN_FULL, POTFAK, LFULL )
    USE sym_prec
    USE lattice
    USE tutor, ONLY: vtutor, isError, newArgument, Wmissing
    IMPLICIT NONE
    TYPE(wpothandle), POINTER :: WH
    INTEGER :: NQ_IN_FULL
    REAL(q) :: POTFAK(WH%WGW%NGDIM)
! local
    INTEGER :: NQ, NQ_IN_FULL_ORIG
    COMPLEX(q) :: C(WH%WGW%NGDIM)
    LOGICAL :: LFULL
    INTEGER :: IERR, I, NP, ierror
    INTEGER, SAVE :: ICOUNTER=0
    TYPE (wavedes1) WGWQ

    
    
    CALL SETWDES(WH%WGW, WGWQ, NQ_IN_FULL)

    IF (.NOT. WH%LSET(NQ_IN_FULL)) THEN
       IF (.NOT. ASSOCIATED(WH%WPOT)) THEN
          CALL vtutor%bug("internal error in GET_WPOT: the WH%WPOT handle is not allocated", "wpot.F", 333)
       ENDIF

! map this k-point to corresponding k-point in the KPOINTS_FULL_ORIG
       NQ_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(KPOINTS_FULL%VKPT(:,NQ_IN_FULL),KPOINTS_FULL_ORIG)
! corresponding k-point in IRZ
       NQ=KPOINTS_FULL_ORIG%NEQUIV(NQ_IN_FULL_ORIG)
       IF (ABS(KPOINTS_FULL_ORIG%VKPT(1,NQ_IN_FULL_ORIG)-KPOINTS_FULL%VKPT(1,NQ_IN_FULL))>TINY .OR. & 
           ABS(KPOINTS_FULL_ORIG%VKPT(2,NQ_IN_FULL_ORIG)-KPOINTS_FULL%VKPT(2,NQ_IN_FULL))>TINY .OR. &
           ABS(KPOINTS_FULL_ORIG%VKPT(3,NQ_IN_FULL_ORIG)-KPOINTS_FULL%VKPT(3,NQ_IN_FULL))>TINY) THEN
          CALL vtutor%bug("internal error in GET_WPOT shift: " // str(KPOINTS_FULL_ORIG%VKPT(:,&
             NQ_IN_FULL_ORIG) - KPOINTS_FULL%VKPT(:,NQ_IN_FULL)), "wpot.F", 344)
       ENDIF

       IF (WH%WGW%NGVECTOR(NQ) /=WH%WGW%NGVECTOR(NQ_IN_FULL)) THEN
          CALL vtutor%bug("internal error in GET_WPOT G-vector: " // str(WH%WGW%NGVECTOR(NQ)) // " " &
             // str(WH%WGW%NGVECTOR(NQ_IN_FULL)), "wpot.F", 349)
       ENDIF

       LFULL=.FALSE.
       NP=WGWQ%NGVECTOR

! allocate response function array at the k-point in the IRZ and read the potential
       IF (.NOT. ASSOCIATED(WH%WPOT(NQ)%RESPONSEFUN)) THEN

! allocate response function array for gamma-point only WH%WPOT%LREALSTORE=.TRUE.
! if the file W????.tmp is read than LFULL is set to .FALSE. and WH%WPOT(NQ) is
! deallocated before leaving the routine
# 364

          CALL ALLOCATE_RESPONSEFUN(WH%WPOT(NQ), WH%WGW%NGDIM, WH%WGW%LGAMMA, WH%WGW%LGAMMA, 1)


          CALL READ_WPOT(WH%WPOT(NQ), WGWQ, NQ, LFULL, WH%IU0, IERR)
# 377

! divide by number of grid points to get proper integration weight
! ok, this is really "stupid" and not very clean
          CALL SCALE_RESPONSE_RESPONSE(WH%WPOT(NQ), 1.0_q/GRIDHF%NPLWV)

! scale convergence corrections by NKREDLFX
! this is only correct if reduction in x y and z is identical
! obviously this will require some carefull reconsideration
! probably I need the full epsilon with a full recalculation ...
          IF ( WH%NKREDLFX /= WH%NKREDLFY .OR. WH%NKREDLFX /= WH%NKREDLFZ) THEN
             CALL vtutor%bug("internal error in GET_WPOT G-vector: NKREDLF must be identical in all " &
                // "directions", "wpot.F", 388)
          ENDIF
          IF (  ABS(SUM(WGWQ%VKPT*WGWQ%VKPT))<=G2ZERO ) THEN
! scale the "head" (1,1) element of the response function array
! in principle the first version is the correct version
! empirically the first version is more accurate for LiF and the second (1._q,0._q) for C
             CALL SCALE_RESPONSE_0(WH%WPOT(NQ), 1.0_q/WH%NKREDLFX/WH%NKREDLFX) 
!            CALL SCALE_RESPONSE_0(WH%WPOT(NQ), 1.0_q/WH%NKREDLFX/WH%NKREDLFY/WH%NKREDLFZ)
          ENDIF

          IF (IERR/=0) THEN
             CALL vtutor%write(isError, Wmissing, newArgument(ival = [NQ]))
          ENDIF
       ELSE
! potential at k-point in IRZ already associated, hence full potential is
! available
          LFULL=.TRUE.
       ENDIF
!==========================================================================
! create potential at corresponding point in full BZ
!==========================================================================
       IF (LFULL) THEN
          IF (.NOT. ASSOCIATED(WH%WPOT(NQ_IN_FULL)%RESPONSEFUN)) THEN
# 414

             CALL ALLOCATE_RESPONSEFUN(WH%WPOT(NQ_IN_FULL), WH%WGW%NGDIM, WH%WGW%LGAMMA, WH%WGW%LGAMMA, 1)


             IF (WH%WPOT(NQ_IN_FULL)%LLEAD) &
             CALL ROTATE_WPOT(NP,WH%WPOT(NQ_IN_FULL)%RESPONSEFUN(1,1,1), WH%WPOT(NQ)%RESPONSEFUN(1,1,1), SIZE( WH%WPOT(NQ)%RESPONSEFUN,1 ), & 
                     WH%KPOINTS_TRANS%CPHASE(1,NQ_IN_FULL), WH%KPOINTS_TRANS%NINDPW(1,NQ_IN_FULL),  &
                     WH%KPOINTS_TRANS%LINV(NQ_IN_FULL), WH%KPOINTS_TRANS%LSHIFT(NQ_IN_FULL))
# 426


          ENDIF
! here (1._q,0._q) can dump the potential to a file to compare it to e.g. calculation
! without symmetry
!          WRITE(77,'(I10,3F8.3)') NQ_IN_FULL, KPOINTS_FULL%VKPT(:,NQ_IN_FULL)
!          DO I=1,NP
!             WRITE(77,'(32F8.3)') WH%WPOT(NQ_IN_FULL)%RESPONSEFUN(1:NP,I,1)*GRIDHF%NPLWV
!          ENDDO
       ENDIF
!==========================================================================
! this part is only required if the full potential has not been read in
! it generates the diagonal part of the potential
! it can be executed anyhow, if (1._q,0._q) desires to set POTFAK correctly
!==========================================================================
       IF (.NOT. LFULL .OR. .TRUE.) THEN
          IF (WH%WPOT(NQ)%LREALSTORE) THEN
! the real value in WH%C stores the potential factor for the cosine transforms
! whereas the imag value in WH%C stores that for the sine transform
             DO I=1,NP
                C(I)=CMPLX(WH%WPOT(NQ)%RESPONSER  ((I-1)*2+1,(I-1)*2+1,1), WH%WPOT(NQ)%RESPONSER  ((I-1)*2+2,(I-1)*2+2,1),q)
             ENDDO
! at this point, the real value in WH%C stores the potential factor for the cosine transforms
! whereas the imag value in WH%C stores that for the sin transform
! but only the real part is passed back in POTFAK
! so better average over both components (except for G=0)
             DO I=2,NP
                C(I)=(WH%WPOT(NQ)%RESPONSER  ((I-1)*2+1,(I-1)*2+1,1)+WH%WPOT(NQ)%RESPONSER  ((I-1)*2+2,(I-1)*2+2,1))/2
             ENDDO
          ELSE
             DO I=1,NP
                C(I)=WH%WPOT(NQ)%RESPONSEFUN(I,I,1)
             ENDDO
          ENDIF
          
! re-index to new k-point
          DO I=1,WGWQ%NGVECTOR
             WH%C(WH%KPOINTS_TRANS%NINDPW(I,NQ_IN_FULL),NQ_IN_FULL)=C(I)
          ENDDO

! deallocate response function at the k-point in the IRZ
          IF (.NOT. LFULL) CALL DEALLOCATE_RESPONSEFUN( WH%WPOT(NQ))
       ENDIF
    ENDIF

    WH%LSET(NQ_IN_FULL)=.TRUE.
    POTFAK(1:WGWQ%NGVECTOR)=WH%C(1:WGWQ%NGVECTOR, NQ_IN_FULL)

    IF (ASSOCIATED(WH%WPOT(NQ_IN_FULL)%RESPONSEFUN)) THEN
       LFULL=.TRUE.
    ELSE
       LFULL=.FALSE.
    ENDIF

    
  END SUBROUTINE GET_WPOT

!*********************************************************************
!> determines screened two electron potential, where the frequency
!> points are distributed among @param NOMEGAPAR groups
!>
!> at an arbitrary q point for a specific frequency group and local
!> frequency points supplied by the calling routine
!> some comments:
!> WPOT_DIST is of the type of a response function array
!> at the gamma point the entries are interpreted as sin and cosine
!> transforms and the potential is stored as a real valued function
!> in RESPONSER
!
!*********************************************************************
  SUBROUTINE GATHER_DISTRIBUTED_WPOT( WH, WPOT_DIST, NQ_IN_FULL, POTFAK )
    USE sym_prec
    USE lattice
    IMPLICIT NONE
    TYPE(wpothandle), POINTER :: WH
    TYPE (responsefunction) WPOT_DIST(:)   !< screened potential in frequency domain
    INTEGER :: NQ_IN_FULL
    REAL(q) :: POTFAK(WH%WGW%NGDIM)
! local
    INTEGER :: NQ, NQ_IN_FULL_ORIG
    COMPLEX(q) :: C(WH%WGW%NGDIM)
    INTEGER :: I, NP, ierror
    INTEGER, SAVE :: ICOUNTER = 0 
    TYPE (wavedes1) WGWQ

    

! scale convergence corrections by NKREDLFX
! this is only correct if reduction in x y and z is identical
! obviously this will require some carefull reconsideration
! probably I need the full epsilon with a full recalculation ...
    IF ( WH%NKREDLFX /= WH%NKREDLFY .OR. WH%NKREDLFX /= WH%NKREDLFZ) THEN
       CALL vtutor%bug("internal error in GET_WPOT G-vector: NKREDLF must be identical in all " &
          // "directions", "wpot.F", 519)
    ENDIF

! on entry NQ_IN_FULL >=0, because all 1 ranks have to communicate
! however, WH%NQ_IN_FULL == -1 indicates if the rank does not need to receive data

! set current q-point (in full BZ)
    CALL SETWDES(WH%WGW, WGWQ, NQ_IN_FULL)

! map this k-point to corresponding k-point in the KPOINTS_FULL_ORIG
    NQ_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(KPOINTS_FULL%VKPT(:,NQ_IN_FULL),KPOINTS_FULL_ORIG)
! corresponding k-point in IRZ
    NQ=KPOINTS_FULL_ORIG%NEQUIV(NQ_IN_FULL_ORIG)
    IF (ABS(KPOINTS_FULL_ORIG%VKPT(1,NQ_IN_FULL_ORIG)-KPOINTS_FULL%VKPT(1,NQ_IN_FULL))>TINY .OR. & 
        ABS(KPOINTS_FULL_ORIG%VKPT(2,NQ_IN_FULL_ORIG)-KPOINTS_FULL%VKPT(2,NQ_IN_FULL))>TINY .OR. &
        ABS(KPOINTS_FULL_ORIG%VKPT(3,NQ_IN_FULL_ORIG)-KPOINTS_FULL%VKPT(3,NQ_IN_FULL))>TINY) THEN
       CALL vtutor%bug("internal error in GET_WPOT shift: " // str(KPOINTS_FULL_ORIG%VKPT(:,&
          NQ_IN_FULL_ORIG) - KPOINTS_FULL%VKPT(:,NQ_IN_FULL)), "wpot.F", 536)
    ENDIF

! internal consistency test, should not be visible to the user
    IF (WH%WGW%NGVECTOR(NQ) /=WH%WGW%NGVECTOR(NQ_IN_FULL)) THEN
       CALL vtutor%bug("internal error in GET_WPOT G-vector: " // str(WH%WGW%NGVECTOR(NQ)) // " " &
          // str(WH%WGW%NGVECTOR(NQ_IN_FULL)), "wpot.F", 542)
    ENDIF

! number of plane waves to take into account
    NP=WGWQ%NGVECTOR

! allocate storage for screened potential
    IF (.NOT. ASSOCIATED(WH%WPOT(1)%RESPONSEFUN)) THEN
! allocate response function array for gamma-point only WH%WPOT%LREALSTORE=.TRUE.
       CALL ALLOCATE_RESPONSEFUN(WH%WPOT(1), WH%WGW%NGDIM, WH%WGW%LGAMMA, WH%WGW%LGAMMA, 1)
    ENDIF
!
! scatters current q-point of screened potential to the 1 rank holding that q-point
! set correct k-point
    CALL COLLECT_WPOT( WH, WGWQ, WPOT_DIST, NQ, WH%IU0 ) 

! this decides whether NQ_IN_FULL is inside IRZ or not
    WH%LFULL=.FALSE.

! in case rank participates only in collective communication
! clean response function for these ranks and return
    IF (  WH%NQ_IN_FULL < 0 ) THEN
       WH%WPOT(1)%HEAD=0
       IF (WH%WPOT(1)%LREALSTORE) THEN
          WH%WPOT(1)%WINGR=0
          WH%WPOT(1)%CWINGR=0
          WH%WPOT(1)%RESPONSER=0
       ELSE
          WH%WPOT(1)%WING=0
          WH%WPOT(1)%CWING=0
          WH%WPOT(1)%RESPONSEFUN=0
       ENDIF
       NQ_IN_FULL = -1
       
       RETURN
    ENDIF
! divide by number of grid points to get proper integration weight
! ok, this is really "stupid" and not very clean
    CALL SCALE_RESPONSE_RESPONSE(WH%WPOT(1), 1.0_q/GRIDHF%NPLWV)

! scale convergence corrections by NKREDLFX
! this is only correct if reduction in x y and z is identical
! obviously this will require some carefull reconsideration
! probably I need the full epsilon with a full recalculation ...
    IF ( WH%NKREDLFX /= WH%NKREDLFY .OR. WH%NKREDLFX /= WH%NKREDLFZ) THEN
       CALL vtutor%bug("internal error in GET_WPOT G-vector: NKREDLF must be identical in all " &
          // "directions", "wpot.F", 588)
    ENDIF
! scale the "head" (1,1) element of the response function array
! in principle the first version is the correct version
! empirically the first version is more accurate for LiF and the second (1._q,0._q) for C
    IF (  ABS(SUM(WGWQ%VKPT*WGWQ%VKPT))<=G2ZERO ) THEN
       CALL SCALE_RESPONSE_0(WH%WPOT(1), 1.0_q/WH%NKREDLFX/WH%NKREDLFX) 
    ENDIF

! rotate potential only if point is outside irreducible BZ
    IF ( NQ_IN_FULL > WH%NQ_IN_IRZ ) WH%LFULL = .TRUE.
!==========================================================================
! create potential at corresponding point in full BZ
!==========================================================================
    IF ( WH%LFULL ) THEN
       IF (.NOT. ASSOCIATED(WH%WPOT(2)%RESPONSEFUN)) THEN
# 607

          CALL ALLOCATE_RESPONSEFUN(WH%WPOT(2), WH%WGW%NGDIM, WH%WGW%LGAMMA, WH%WGW%LGAMMA, 1)

       ENDIF
!
! store q-point outside of irreducible BZ into WPOT(2) temporarily
!
       IF (WH%WPOT(2)%LLEAD) THEN
          CALL ROTATE_WPOT(NP,WH%WPOT(2)%RESPONSEFUN(1,1,1), WH%WPOT(1)%RESPONSEFUN(1,1,1), &
                  SIZE( WH%WPOT(1)%RESPONSEFUN,1 ), WH%KPOINTS_TRANS%CPHASE(1,NQ_IN_FULL), &
                  WH%KPOINTS_TRANS%NINDPW(1,NQ_IN_FULL), WH%KPOINTS_TRANS%LINV(NQ_IN_FULL), &
                  WH%KPOINTS_TRANS%LSHIFT(NQ_IN_FULL))
! TODO: skip this and use WPOT(2) directly
! point to WPOT(1), since this is the (1._q,0._q) used later on
          WH%WPOT(1)%RESPONSEFUN = WH%WPOT(2)%RESPONSEFUN
       ENDIF
       CALL DEALLOCATE_RESPONSEFUN( WH%WPOT(2) )
# 628

    ENDIF

!==========================================================================
! this part is only required if the full potential has not been read in
! it generates the diagonal part of the potential
! it can be executed anyhow, if (1._q,0._q) desires to set POTFAK correctly
!==========================================================================
    IF (.NOT. WH%LFULL .OR. .TRUE.) THEN
       IF (WH%WPOT(1)%LREALSTORE) THEN
! the real value in WH%C stores the potential factor for the cosine transforms
! whereas the imag value in WH%C stores that for the sine transform
          DO I=1,NP
             C(I)=CMPLX(WH%WPOT(1)%RESPONSER((I-1)*2+1,(I-1)*2+1,1), WH%WPOT(1)%RESPONSER((I-1)*2+2,(I-1)*2+2,1),q)
          ENDDO
! at this point, the real value in WH%C stores the potential factor for the cosine transforms
! whereas the imag value in WH%C stores that for the sin transform
! but only the real part is passed back in POTFAK
! so better average over both components (except for G=0)
          DO I=2,NP
             C(I)=(WH%WPOT(1)%RESPONSER((I-1)*2+1,(I-1)*2+1,1)+WH%WPOT(1)%RESPONSER((I-1)*2+2,(I-1)*2+2,1))/2
          ENDDO
       ELSE
          DO I=1,NP
             C(I)=WH%WPOT(1)%RESPONSEFUN(I,I,1)
          ENDDO
       ENDIF
       
! re-index to new k-point
       DO I=1,WGWQ%NGVECTOR
          WH%C(WH%KPOINTS_TRANS%NINDPW(I,NQ_IN_FULL),NQ_IN_FULL)=C(I)
       ENDDO
    ENDIF

! set diagonal potential here to POTFAK
    WH%LSET(NQ_IN_FULL)=.TRUE.
    POTFAK(1:WGWQ%NGVECTOR)=WH%C(1:WGWQ%NGVECTOR, NQ_IN_FULL)

    IF (ASSOCIATED(WH%WPOT(1)%RESPONSEFUN)) THEN
       WH%LFULL=.TRUE.
    ELSE
       WH%LFULL=.FALSE.
    ENDIF
    
  END SUBROUTINE GATHER_DISTRIBUTED_WPOT

!*********************************************************************
!
!> remove the G=0 component from the stored WPOT and return it
!> to the calling routine
!
!*********************************************************************

  FUNCTION GET_WPOT_GZERO(WH)
    IMPLICIT NONE
    TYPE(wpothandle), POINTER :: WH
    REAL(q) :: GET_WPOT_GZERO
! local
    INTEGER :: NK

    IF (.NOT. ASSOCIATED(WH)) RETURN

! deallocate the potential handle
    DO NK=1,KPOINTS_FULL_ORIG%NKPTS
       IF (  ABS(SUM(KPOINTS_FULL_ORIG%VKPT(:,NK)*KPOINTS_FULL_ORIG%VKPT(:,NK)))<=G2ZERO ) THEN
          IF (ASSOCIATED(WH%WPOT)) THEN
             GET_WPOT_GZERO=WH%WPOT(NK)%RESPONSEFUN(1,1,1)
             WH%WPOT(NK)%RESPONSEFUN(1,1,1)=0
          ELSE IF (ASSOCIATED(WH%C)) THEN
             GET_WPOT_GZERO=WH%C(1, 1)
             WH%C(1, 1)=0
          ENDIF
       ENDIF
    ENDDO


! ok this scaling is really odd, but thats used in SET_GFAC_WAVEFUN
    GET_WPOT_GZERO=GET_WPOT_GZERO*GRIDHF%NPLWV &
           *KPOINTS_FULL_ORIG%WTKPT(1)*WH%NKREDLFX*WH%NKREDLFY*WH%NKREDLFZ

  END FUNCTION GET_WPOT_GZERO

!*********************************************************************
!
!> write the screened potential to a file
!>
!> a slight complication is caused by the possibly subtracted
!> bare Hartree Fock potential for q/=0 the bare kernel is stored
!> in DATAKE
!> the bare singularity correction at q=0 is passed by the calling
!> routine (FSG0)
!
!*********************************************************************

  SUBROUTINE WRITE_WPOT(WPOT, WGWQ, FSG0, NQ, LFOCK_SUBTRACT, FILEBASE_ )

    IMPLICIT NONE
    TYPE (responsefunction) :: WPOT  !< screened two electron potential-HF
    TYPE (wavedes1)         :: WGWQ  !< basis set descriptor
    REAL (q)                :: FSG0  !< HF singularity correction
    INTEGER              NQ          !< q-point for which local field correction are required
    LOGICAL :: LFOCK_SUBTRACT
    CHARACTER(*), OPTIONAL :: FILEBASE_
! local
    INTEGER              IERR        ! error status
    CHARACTER (4) :: APP
    INTEGER       :: IU=72, NOMEGA, I, NP
    COMPLEX(q)    :: W( WPOT%NP2, WPOT%NP2)
    CHARACTER(32) :: FILEBASE

    IERR=0
    IF (WPOT%NOMEGA_LOW==1) THEN
! WRITE(*,*) 'node ',WGWQ%COMM_INTER%NODE_ME,' writing ',NQ
       NOMEGA=1
       WRITE (APP  , "(4I1)") MOD(NQ/1000,10),MOD(NQ/100,10),MOD(NQ/10,10), MOD(NQ,10)
       FILEBASE="W"
       IF(PRESENT(FILEBASE_)) FILEBASE=TRIM(FILEBASE_)
       OPEN( UNIT=IU, FILE=TRIM(FILEBASE)//APP//".tmp", IOSTAT=IERR, FORM='UNFORMATTED')

       NP=WGWQ%NGVECTOR
       IF (WGWQ%LGAMMA) NP=NP*2
       IF (NP>WPOT%NP2) THEN
          CALL vtutor%bug("internal error in WRITE_WPOT: NP>WPOT%NP2 " // str(NP) // " " // &
             str(WPOT%NP2), "wpot.F", 751)
       ENDIF

       IF (IERR==0) WRITE(IU) NP, 0
       IF (IERR==0) WRITE(IU) WPOT%HEAD(:,:,NOMEGA)
       IF (WPOT%LREALSTORE) THEN
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%WINGR(1:NP,:,NOMEGA)
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%CWINGR(1:NP,:,NOMEGA)
       ELSE
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%WING(1:NP,:,NOMEGA)
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%CWING(1:NP,:,NOMEGA)
       ENDIF

       IF (WPOT%LREALSTORE) THEN
          W=WPOT%RESPONSER  (:,:,NOMEGA)
       ELSE
          W=WPOT%RESPONSEFUN(:,:,NOMEGA)
       ENDIF


       IF (LFOCK_SUBTRACT) THEN
          IF (WPOT%LGAMMA) THEN
             DO I=2,NP
                W(I,I)=W(I,I)+WGWQ%DATAKE(I,1)
             ENDDO
             W(1,1)=W(1,1)+FSG0
          ELSE
             DO I=1,NP
                W(I,I)=W(I,I)+WGWQ%DATAKE(I,1)
             ENDDO
          ENDIF
       ENDIF
       IF (IERR==0) THEN
          IF (WPOT%LREAL) THEN
             WRITE(IU, IOSTAT=IERR) (REAL(W(I,I),q),I=1, NP)
          ELSE
             WRITE(IU, IOSTAT=IERR) (W(I,I),I=1, NP)
          ENDIF
       ENDIF

       CLOSE(IU)
    ENDIF
    CALL M_sum_i(WGWQ%COMM_INTER, IERR, 1)

  END SUBROUTINE WRITE_WPOT

!*********************************************************************
!
!> write the screened potential to a file
!> this version writes the entire matrix @param WPOT(G,G')
!
!*********************************************************************

  SUBROUTINE WRITE_WPOT_FULL(WPOT, WGWQ, FSG0, NQ, LFOCK_SUBTRACT ) 

    IMPLICIT NONE
    TYPE (responsefunction) :: WPOT  !< screened two electron potential-HF
    TYPE (wavedes1)         :: WGWQ  !< basis set descriptor
    REAL (q)                :: FSG0  !< HF singularity correction
    INTEGER              NQ          !< q-point for which local field correction are required
    LOGICAL :: LFOCK_SUBTRACT
! local
    INTEGER              IERR        ! error status
    CHARACTER (4) :: APP
    INTEGER       :: IU=72, NOMEGA, I, NP
    COMPLEX(q)    :: W( WPOT%NP2, WPOT%NP2)

    IERR=0
! write full potential at every k-point
    IF (WPOT%NOMEGA_LOW==1) THEN
       NOMEGA=1
       WRITE (APP  , "(4I1)") MOD(NQ/1000,10),MOD(NQ/100,10),MOD(NQ/10,10), MOD(NQ,10)
       OPEN( UNIT=IU, FILE="WFULL"//APP//".tmp", IOSTAT=IERR, FORM='UNFORMATTED')

       NP=WGWQ%NGVECTOR
       IF (WGWQ%LGAMMA) NP=NP*2
       IF (NP>WPOT%NP2) THEN
          CALL vtutor%bug("internal error in WRITE_WPOT_FULL: NP>WPOT%NP2 " // str(NP) // " " // &
             str(WPOT%NP2), "wpot.F", 829)
       ENDIF

       IF (IERR==0) WRITE(IU) NP, NP
       IF (IERR==0) WRITE(IU) WPOT%HEAD(:,:,NOMEGA)
       IF (WPOT%LREALSTORE) THEN
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%WINGR(1:NP,:,NOMEGA)
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%CWINGR(1:NP,:,NOMEGA)
       ELSE
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%WING(1:NP,:,NOMEGA)
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%CWING(1:NP,:,NOMEGA)

# 854

       ENDIF

       IF (WPOT%LREALSTORE) THEN
          W=WPOT%RESPONSER  (:,:,NOMEGA)
       ELSE
          W=WPOT%RESPONSEFUN(:,:,NOMEGA)
       ENDIF

       IF (LFOCK_SUBTRACT) THEN
          IF (WPOT%LGAMMA) THEN
             DO I=2,NP
                W(I,I)=W(I,I)+WGWQ%DATAKE(I,1)
             ENDDO
             W(1,1)=W(1,1)+FSG0
          ELSE
             DO I=1,NP
                W(I,I)=W(I,I)+WGWQ%DATAKE(I,1)
             ENDDO
          ENDIF
       ENDIF
       IF (WPOT%LREAL) THEN
! WPOT%LREAL indicates that we are using the gamma point only version
! in this case the response function is necessarily real valued
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) REAL(W(1:NP,1:NP),q)
       ELSE
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) W(1:NP,1:NP)
# 886

       ENDIF
       CLOSE(IU)
    ENDIF
    CALL M_sum_i(WGWQ%COMM_INTER, IERR, 1)

  END SUBROUTINE WRITE_WPOT_FULL

!*********************************************************************
!
!> read screened potential from the file
!>
!> this version always reads from WFULLXXXX.tmp but is capable
!> to read both WXXXX.tmp
!> so if you want to read WXXXX.tmp copy it to WFULLXXXX.tmp
!
!*********************************************************************

  SUBROUTINE READ_WPOT(WPOT, WGWQ, NQ, LFULL, IU0, IERR)

    IMPLICIT NONE
    TYPE (responsefunction) :: WPOT  !< screened two electron potential
    TYPE (wavedes1)         :: WGWQ  !< basis set descriptor
    INTEGER              NQ         !< q-point for which local field correction are required
    LOGICAL              LFULL      !< all elements read
    INTEGER              IU0        !< stderr
    INTEGER              IERR       !< error status
! local
    CHARACTER (4) :: APP
    INTEGER       :: N1, N2, NP
    INTEGER :: IU=72, NOMEGA, I, NALLOC
    LOGICAL :: LREAD
    INTEGER :: ierror

    IERR=0
    NOMEGA=1
    WRITE (APP  , "(4I1)") MOD(NQ/1000,10),MOD(NQ/100,10),MOD(NQ/10,10), MOD(NQ,10)

    OPEN( UNIT=IU, FILE="WFULL"//APP//".tmp", STATUS="OLD", IOSTAT=IERR, FORM='UNFORMATTED')
    IF (IERR/=0) THEN
       IF (IU0>=0) WRITE(IU0,*) 'reading now W'//APP//".tmp"
       OPEN( UNIT=IU, FILE="W"//APP//".tmp", STATUS="OLD", IOSTAT=IERR, FORM='UNFORMATTED')
    ELSE
       IF (IU0>=0) WRITE(IU0,*) 'reading now WFULL'//APP//".tmp"
    ENDIF

! stop if neither WFULLXXXX.tmp nor WXXXX.tmp files are read
    IF (IERR/=0) CALL vtutor%error("Could not read WFULLXXXX.tmp/WXXXX.tmp files. &
                                    Make sure that these files are present in the directory.")

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2
    IF (NP>WPOT%NP2) THEN
       CALL vtutor%bug("internal error in READ_WPOT: NP>WPOT%NP2 " // str(NP) // " " // str(WPOT%NP2), "wpot.F", 939)
    ENDIF

    IF (IERR==0) READ(IU) N1, N2
! there are N1 data on the actual file
! either these match WPOT%NP2 (old version)
! or they match NP
! in both cases we can read the file without problems, since
! all arrays are allocted as WPOT%NP2 and NP is smaller than WPOT%NP2

! aT & gK: 3.04.24 this condition should be sufficinet to catch ENCUTGW mismatch
    IF (N1 /= NP) THEN
       CALL vtutor%error("Number of G-vectors in WFULLXXXX.tmp/WXXXX.tmp (" &
           // str(N1) //") does not match " &
           // "the basis set of the response function ("// str(NP) // "). \n" &
           // "Make sure that ENCUTGW is the same in both calculations.")
       IERR=1
    ENDIF

! WPOT%NP2 = NGDIM in WDES and N1 is the actual number of Gs in WFULL, which corresponds to NGVECTOR
! the condition N1 /= WPOT%NP2 is not useful
!IF ( N1/= WPOT%NP2 .AND. N1/= NP) THEN
!   IERR=1
!   IF (IU0>=0) WRITE(IU0,'(" ",A,2I10)') "N1 differs from file",N1,WPOT%NP2
!   IF (N1<=WPOT%NP2) THEN
!      IF (IU0>=0) WRITE(IU0,'(" ",A,2I10)') "trying to continue anyhow"
!      IERR=0
!   ENDIF
!ENDIF

    IF (IERR==0) READ(IU) WPOT%HEAD(:,:,NOMEGA)
    IF (WPOT%LREALSTORE) THEN
       IF (IERR==0) READ(IU, IOSTAT=IERR) WPOT%WINGR(1:N1,:,NOMEGA)
       IF (IERR==0) READ(IU, IOSTAT=IERR) WPOT%CWINGR(1:N1,:,NOMEGA)
    ELSE
        IF (IERR==0) READ(IU, IOSTAT=IERR) WPOT%WING (1:N1,:,NOMEGA)
        IF (IERR==0) READ(IU, IOSTAT=IERR) WPOT%CWING(1:N1,:,NOMEGA)
    ENDIF

    IF (IERR==0) THEN
       IF (N2==0) THEN
          LFULL=.FALSE.
       ELSE
          LFULL=.TRUE.
       ENDIF
    ENDIF

    IF (IERR==0 .AND. WPOT%LLEAD) THEN
       IF (N2==0) THEN
! read only diagonal components
          WPOT%RESPONSEFUN=0
          IF (WPOT%LREALSTORE) THEN
             READ(IU, IOSTAT=IERR) (WPOT%RESPONSER  (I,I,NOMEGA),I=1,N1)
          ELSE
             READ(IU, IOSTAT=IERR) (WPOT%RESPONSEFUN(I,I,NOMEGA),I=1,N1)
          ENDIF
          IF (IERR/=0 .AND. IU0>=0) THEN
             WRITE(IU0,'(" ",A)') "error upon READ of WXXXX.tmp file"
          ENDIF
       ELSE
! read all components
          IF (WPOT%LREALSTORE) THEN
             WPOT%RESPONSER(:,:,NOMEGA)=0
             READ(IU, IOSTAT=IERR)  WPOT%RESPONSER  (1:N1,1:N1,NOMEGA)
          ELSE
             WPOT%RESPONSEFUN(:,:,NOMEGA)=0
             READ(IU, IOSTAT=IERR)  WPOT%RESPONSEFUN(1:N1,1:N1,NOMEGA)
          ENDIF
          IF (IERR/=0 .AND. IU0>=0) THEN
             WRITE(IU0,'(" ",A)') "error upon READ of WXXXX.tmp file"
          ENDIF
       ENDIF
    ENDIF

# 1020

    CLOSE(IU)

  END SUBROUTINE READ_WPOT

!*********************************************************************
!> collects response function among frequency groups into @param WH\%WPOT(1)
!> for the currently passed q-point NQ in the first BZ
!*********************************************************************
  SUBROUTINE COLLECT_WPOT(WH, WGWQ, WPOT_DIST, NQ, IU0)
    IMPLICIT NONE
    TYPE( wpothandle )      :: WH
    TYPE (wavedes1)         :: WGWQ  !< basis set descriptor
    TYPE (responsefunction) :: WPOT_DIST(:)  !< screened two electron potential
    INTEGER              NQ         !< q-point for which local field correction are required
    INTEGER              IU0        !< stderr
! local
    INTEGER              :: NP, NP1, NP2, I
    INTEGER              :: IRANK_SEND ! rank who scatters data to the other ranks
    INTEGER              :: IRANK_IN_GROUP ! local rank in omega group
    INTEGER              :: IPOS_GLOBAL
    INTEGER              :: IRANK_RECV
    INTEGER, ALLOCATABLE :: NQ_IBZ_DIST(:)
    INTEGER              :: IERROR

    

    CALL SET_RESPONSE_KPOINT(WH%WPOT(1), WH%WGW%VKPT(:,NQ), NQ)

! size consistency checks
    NP = WGWQ%NGVECTOR
    IF ( WGWQ%LGAMMA ) NP=NP*2
    IF (NP>WH%WPOT(1)%NP2) THEN
       CALL vtutor%bug("internal error in COLLECT_WPOT: NP>WPOT%NP2 " // str(NP) // " " // &
          str(WH%WPOT(1)%NP2), "wpot.F", 1054)
    ENDIF

    NP1 = WPOT_DIST(NQ)%NP1   ! blocking size of row
    NP2 = WPOT_DIST(NQ)%NP2   ! full column size of response function
! depends on # of CPUS in (1._q,0._q) OMEGA group


! determine which rank handles which irreducible q-point
    ALLOCATE( NQ_IBZ_DIST( WH%O%COMM%NCPU )  ) 
    NQ_IBZ_DIST = 0 
    NQ_IBZ_DIST( WH%O%COMM%NODE_ME ) = NQ 
    CALL M_sum_i( WH%O%COMM, NQ_IBZ_DIST, WH%O%COMM%NCPU )

! clean screened potential
    WH%WPOT(1)%HEAD=0 
    IF (WH%WPOT(1)%LREALSTORE) THEN
       WH%WPOT(1)%WINGR=0
       WH%WPOT(1)%CWINGR=0
       WH%WPOT(1)%RESPONSER=0
    ELSE
       WH%WPOT(1)%WING=0
       WH%WPOT(1)%CWING=0
       WH%WPOT(1)%RESPONSEFUN=0
    ENDIF
          
! WPOT_DIST for current q-point is gathered into WH%WPOT( 1 )
! every local 1 rank in the NOMEGAPAR group holding the
! given WH%NOMEGA_GLOBAL frequency point is sending data
! of the appropriate q-point to the global 1 ranks.
! the appropriate Q-point for each rank is stored in NQ_IBZ_DIST( NODE_ME )
! Important:
! All global 1 ranks are participating in this resdistribution
! however, ranks with current NQ=1 point and WH%NQ_IN_FULL=0
! do not receive data.

! loop over ranks holding current frequency point and which is sending
    send: DO IRANK_IN_GROUP = 1, WH%O%DISTRIBUTION( WH%ID_NOMEGAPAR, 1 )

! find group and ID of rank who is scattering data
       IRANK_SEND = 0 
! global index in response function (adress in receiving buffer)
       IPOS_GLOBAL = 0 
! determine sending node
       IF ( WH%ID_NOMEGAPAR == WH%O%COMM_BETWEEN_GROUPS%NODE_ME .AND. &
          IRANK_IN_GROUP == WH%O%COMM_IN_GROUP%NODE_ME ) THEN
          IRANK_SEND = WH%O%COMM%NODE_ME
          IPOS_GLOBAL = (WH%O%COMM_IN_GROUP%NODE_ME-1)*NP2 + 1 
       ENDIF
       CALL M_sum_i( WH%O%COMM, IRANK_SEND, 1 ) 
       CALL M_sum_i( WH%O%COMM, IPOS_GLOBAL, 1 ) 

       recv: DO IRANK_RECV = 1 , WH%O%COMM%NCPU
! send stripe hold by rank IRANK_SEND to all other ranks
! send data directly to every node
! advantage over using MPI_SCATTER: no need to set up SENDBUF
          IF ( WH%O%COMM%NODE_ME == IRANK_SEND ) THEN
!
! do not send if data is on the same node
!
             IF ( IRANK_RECV == IRANK_SEND ) THEN
!WRITE(*,'(A,5I4)')'on node',WH%O%COMM%NODE_ME, IRANK_SEND, IRANK_RECV, WH%NQ_IN_FULL, NQ_IBZ_DIST( IRANK_RECV )
!
! gamma only version works with double precision arrays
!
                IF ( WH%WPOT(1)%LREALSTORE ) THEN
                   WH%WPOT(1)%RESPONSER(1:2*NP1,IPOS_GLOBAL:IPOS_GLOBAL+2*NP2-1,1) = &
                      WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%RESPONSER( 1: 2*NP1, 1:2*NP2, WH%NOMEGA_LOCAL)
                   WH%WPOT(1)%CWINGR(1:2*NP1,1:3,1) = &
                      WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%CWINGR( 1:2*NP1,1:3,WH%NOMEGA_LOCAL)
                   WH%WPOT(1)%WINGR(1:2*NP1,1:3,1) = &
                      WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%WINGR( 1:2*NP1,1:3,WH%NOMEGA_LOCAL)
!
! full-k version works with double complex arrays
!
                ELSE
                   WH%WPOT(1)%RESPONSEFUN(1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1,1) = &
                      WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%RESPONSEFUN( 1: NP1, 1:NP2, WH%NOMEGA_LOCAL)
                   WH%WPOT(1)%CWING(1:NP1,1:3,1) = &
                      WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%CWING( 1:NP1,1:3,WH%NOMEGA_LOCAL)
                   WH%WPOT(1)%WING(1:NP1,1:3,1) = &
                      WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%WING( 1:NP1,1:3,WH%NOMEGA_LOCAL)
                ENDIF
                WH%WPOT(1)%HEAD(1:3,1:3,1) = &
                   WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%HEAD( 1:3, 1:3, WH%NOMEGA_LOCAL)
!
! otherwise send data to rank IRANK_RECV
!
             ELSE
!WRITE(*,'(A,5I4)')'sending',WH%O%COMM%NODE_ME, IRANK_SEND, IRANK_RECV, WH%NQ_IN_FULL, NQ_IBZ_DIST( IRANK_RECV )
!
! gamma only version works with double precision arrays
!
                IF ( WH%WPOT(1)%LREALSTORE ) THEN
                   CALL M_send_d( WH%O%COMM, IRANK_RECV, WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%RESPONSER( 1: 2*NP1, 1:2*NP2, WH%NOMEGA_LOCAL), 4*NP1*NP2 )
                   CALL M_send_d( WH%O%COMM, IRANK_RECV, WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%CWINGR( 1:2*NP1,1:3,WH%NOMEGA_LOCAL), 6*NP1 )
                   CALL M_send_d( WH%O%COMM, IRANK_RECV, WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%WINGR( 1:2*NP1,1:3,WH%NOMEGA_LOCAL), 6*NP1 )
!
! full-k version works with double complex arrays
!
                ELSE
                   CALL M_send_z( WH%O%COMM, IRANK_RECV, WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%RESPONSEFUN( 1: NP1, 1:NP2, WH%NOMEGA_LOCAL), NP1*NP2 )
                   CALL M_send_z( WH%O%COMM, IRANK_RECV, WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%CWING( 1:NP1,1:3,WH%NOMEGA_LOCAL), 3*NP1 )
                   CALL M_send_z( WH%O%COMM, IRANK_RECV, WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%WING( 1:NP1,1:3,WH%NOMEGA_LOCAL), 3*NP1 )
                ENDIF
                CALL M_send_z( WH%O%COMM, IRANK_RECV, WPOT_DIST( NQ_IBZ_DIST( IRANK_RECV ) )%HEAD( 1:3, 1:3, WH%NOMEGA_LOCAL), 9 )
             ENDIF
          ENDIF

          IF ( WH%O%COMM%NODE_ME == IRANK_RECV .AND. IRANK_RECV /= IRANK_SEND ) THEN
!WRITE(*,'(A,5I4)')'receiving',WH%O%COMM%NODE_ME, IRANK_SEND, IRANK_RECV, WH%NQ_IN_FULL, NQ_IBZ_DIST( IRANK_RECV )
!
! gamma only version works with double precision arrays
!
             IF ( WH%WPOT(1)%LREALSTORE ) THEN
                CALL M_recv_d( WH%O%COMM, IRANK_SEND, WH%WPOT(1)%RESPONSER(1,IPOS_GLOBAL,1), 4*NP1*NP2 )
                CALL M_recv_d( WH%O%COMM, IRANK_SEND, WH%WPOT(1)%CWINGR(1,1,1), 6*NP1 )
                CALL M_recv_d( WH%O%COMM, IRANK_SEND, WH%WPOT(1)%WINGR(1,1,1), 6*NP1 )
!
! full-k version works with double complex arrays
!
             ELSE
                CALL M_recv_z( WH%O%COMM, IRANK_SEND, WH%WPOT(1)%RESPONSEFUN(1,IPOS_GLOBAL,1), NP1*NP2 )
                CALL M_recv_z( WH%O%COMM, IRANK_SEND, WH%WPOT(1)%CWING(1,1,1), 3*NP1 )
                CALL M_recv_z( WH%O%COMM, IRANK_SEND, WH%WPOT(1)%WING(1,1,1), 3*NP1 )
             ENDIF
             CALL M_recv_z( WH%O%COMM, IRANK_SEND, WH%WPOT(1)%HEAD(1,1,1), 9 )
          ENDIF

          CALL M_barrier( WH%O%COMM )

       ENDDO recv
    ENDDO send
# 1198

   
! add V_GG'(q) to the diagonal of the matrix, that is V_x
    IF (WH%WPOT(1)%LREALSTORE) THEN
       DO I=2,NP
          WH%WPOT(1)%RESPONSER(I,I,1) = &
             WH%WPOT(1)%RESPONSER(I,I,1) + WH%WGW%DATAKE(I,1,NP)
       ENDDO
       WH%WPOT(1)%RESPONSER(1,1,1) = &
          WH%WPOT(1)%RESPONSER(1,1,1) + WH%FSG0
    ELSE
! special treatment for gamma point
       IF ( WH%WPOT(1)%LGAMMA) THEN
          WH%WPOT(1)%RESPONSEFUN(1,1,1) = &
             WH%WPOT(1)%RESPONSEFUN(1,1,1) + WH%FSG0
          DO I=2,NP
             WH%WPOT(1)%RESPONSEFUN(I,I,1) = &
                WH%WPOT(1)%RESPONSEFUN(I,I,1) + WH%WGW%DATAKE(I,1,NQ)
          ENDDO
! for non-gamma points head is obtained from V_GG
       ELSE
          DO I=1,NP
             WH%WPOT(1)%RESPONSEFUN(I,I,1) = &
                WH%WPOT(1)%RESPONSEFUN(I,I,1) + WH%WGW%DATAKE(I,1,NQ)
          ENDDO
       ENDIF
    ENDIF

    
  END SUBROUTINE COLLECT_WPOT

!*********************************************************************
!> collects response function among frequency groups into @param WH\%WPOT(1)
!> for the currently passed q-point NQ in the first BZ
!*********************************************************************
  SUBROUTINE ACCUMULATE_WPOT(WH, WGWQ, WPOT_DIST, NQ, COMM_KIN)
    IMPLICIT NONE
    TYPE( wpothandle )      :: WH
    TYPE (wavedes1)         :: WGWQ          !< basis set descriptor
    TYPE (responsefunction) :: WPOT_DIST     !< screened two electron potential
    INTEGER                 :: NQ            !< q-point for which local field correction are required
    TYPE( communic )        :: COMM_KIN      !< communicator
! local
    INTEGER                 :: NP, NP1, NP2
    INTEGER                 :: IERROR

    INTEGER ID_NOMEGA
    INTEGER NOMEGA_LOCAL         ! which OMEGAPAR group is treated currently
    INTEGER NOMEGA_GLOBAL
    INTEGER NCPU , NPOINTS_IN_ROOT_GROUP , NPOINTS

    INTEGER              :: IRANK_SEND ! rank who scatters data to the other ranks
    INTEGER              :: IRANK_IN_GROUP ! local rank in omega group
    INTEGER              :: IPOS_GLOBAL
    INTEGER              :: IRANK_RECV, I
    INTEGER, ALLOCATABLE :: NQ_IBZ_DIST(:)
    COMPLEX(Q), ALLOCATABLE :: RECV_CBUFF1(:,:)
    COMPLEX(Q), ALLOCATABLE :: RECV_CBUFF2(:,:)
    COMPLEX(Q), ALLOCATABLE :: RECV_CBUFF3(:,:)
    COMPLEX(Q), ALLOCATABLE :: RECV_CBUFF4(:,:)
    REAL(q), ALLOCATABLE :: RECV_RBUFF1(:,:)
    REAL(q), ALLOCATABLE :: RECV_RBUFF2(:,:)
    REAL(q), ALLOCATABLE :: RECV_RBUFF3(:,:)
    REAL(q) :: WEIGHT = 1._q
# 1268

    

!clear the response array
    CALL CLEAR_RESPONSE(WH%WPOT(1))

    CALL SET_RESPONSE_KPOINT(WH%WPOT(1), WH%WGW%VKPT(:,NQ), NQ)

! size consistency checks
    NP = WGWQ%NGVECTOR
    NP1 = WPOT_DIST%NP1   ! blocking size of row
    NP2 = WPOT_DIST%NP2   ! full column size of response function
    IF ( WGWQ%LGAMMA ) THEN
       NP=NP*2
       NP1 = 2*NP1
    ENDIF

! depends on # of CPUS in (1._q,0._q) OMEGA group

! determine which rank handles which irreducible q-point
    ALLOCATE( NQ_IBZ_DIST( COMM_KIN%NCPU )  ) 
    NQ_IBZ_DIST = 0 
!NQ_IBZ_DIST( COMM_KIN%NODE_ME ) = NQ
! only (1._q,0._q) q-point is handles at a time
    NQ_IBZ_DIST( COMM_KIN%NODE_ME ) = 1
    CALL M_sum_i( COMM_KIN, NQ_IBZ_DIST, COMM_KIN%NCPU )

    NPOINTS_IN_ROOT_GROUP = WH%O%NPOINTS_IN_ROOT_GROUP
    NCPU = WH%O%COMM_BETWEEN_GROUPS%NCPU
    NPOINTS = WH%O%NPOINTS
# 1302


!> some auxillary mpi buffers where received data is stored
    IF ( COMM_KIN%NODE_ME == 1  ) THEN 
       IF ( WH%WGW%LGAMMA ) THEN
          ALLOCATE(RECV_RBUFF1( NP1,NP2 ) )
          ALLOCATE(RECV_RBUFF2( NP1,3 ) )
          ALLOCATE(RECV_RBUFF3( NP1,3 ) )
       ELSE
          ALLOCATE(RECV_CBUFF1( NP1, NP2 ) )
          ALLOCATE(RECV_CBUFF2( NP1, 3 ) )
          ALLOCATE(RECV_CBUFF3( NP1, 3 ) )
       ENDIF
       ALLOCATE(RECV_CBUFF4( 3,3 ) )
    ENDIF


    freq_in_group: DO NOMEGA_LOCAL = 1, NPOINTS_IN_ROOT_GROUP
    freq_groups: DO ID_NOMEGA = 1, NCPU
! loop over frequency groups
       NOMEGA_GLOBAL = ID_NOMEGA + ( NOMEGA_LOCAL - 1 ) * NCPU
       IF ( NOMEGA_GLOBAL < 0 .OR. NOMEGA_GLOBAL > NPOINTS ) CYCLE
! set current frequency point update in wpot handle
       WH%ID_NOMEGAPAR = ID_NOMEGA
       WH%NOMEGA_LOCAL = NOMEGA_LOCAL
       WH%NOMEGA_GLOBAL = NOMEGA_GLOBAL 

!WRITE(100+WGWQ%COMM%NODE_ME,'(L2, I3,A,10I4)') WH%WPOT(1)%LREALSTORE, NOMEGA_GLOBAL, "1-"//TRIM(ADJUSTL(APP)), WH%ID_NOMEGAPAR , WH%O%DISTRIBUTION( WH%ID_NOMEGAPAR, : )
! gather specific frequency point
       IF  ( WH%NOMEGA_GLOBAL /= NOMEGA_DUMP .AND. NOMEGA_DUMP > 0  ) CYCLE

! in case (1._q,0._q) specific point is gathered
       IF  ( NOMEGA_DUMP > 0  ) WEIGHT = 1._q 
! in case integration over omega is 1._q
       IF  ( NOMEGA_DUMP == 0  ) WEIGHT = WH%OMEGA_TAU_WEIGHT( WH%NOMEGA_GLOBAL ) 

!WRITE(100+WGWQ%COMM%NODE_ME,FM) ( REAL(WPOT_DIST%RESPONSEFUN(I,1:NP2,WH%NOMEGA_LOCAL) ), I=1,NP1 )

! WPOT_DIST for current q-point is gathered into WH%WPOT( 1 )
! every local 1 rank in the NOMEGAPAR group holding the
! given WH%NOMEGA_GLOBAL frequency point is sending data
! of the appropriate q-point to the global 1 ranks.
! the appropriate Q-point for each rank is stored in NQ_IBZ_DIST( NODE_ME )
! Important:
! All global 1 ranks are participating in this resdistribution
! however, ranks with current NQ=1 point and WH%NQ_IN_FULL=0
! do not receive data.

! loop over ranks holding current frequency point and which is sending
    send: DO IRANK_IN_GROUP = 1, WH%O%DISTRIBUTION( WH%ID_NOMEGAPAR, 1 )
! find group and ID of rank who is scattering data
       IRANK_SEND = 0 
! global index in response function (adress in receiving buffer)
       IPOS_GLOBAL = 0 
! determine sending node
       IF ( WH%ID_NOMEGAPAR == WH%O%COMM_BETWEEN_GROUPS%NODE_ME .AND. &
          IRANK_IN_GROUP == WH%O%COMM_IN_GROUP%NODE_ME ) THEN
          IRANK_SEND = COMM_KIN%NODE_ME
          IPOS_GLOBAL = (WH%O%COMM_IN_GROUP%NODE_ME-1)*NP2 + 1 
       ENDIF
       CALL M_sum_i( COMM_KIN, IRANK_SEND, 1 ) 
       CALL M_sum_i( COMM_KIN, IPOS_GLOBAL, 1 ) 

!WRITE(*,*) WGWQ%LGAMMA, WH%WPOT(1)%LREALSTORE, COMM_KIN%NODE_ME, WH%NOMEGA_LOCAL, WH%NOMEGA_GLOBAL, IRANK_SEND, IPOS_GLOBAL,IRANK_IN_GROUP
       IRANK_RECV = 1
! send stripe hold by rank IRANK_SEND to all other ranks
! send data directly to every node
! advantage over using MPI_SCATTER: no need to set up SENDBUF
       IF ( COMM_KIN%NODE_ME == IRANK_SEND ) THEN
!
! do not send if data is on the same node
!
          IF ( IRANK_RECV == IRANK_SEND ) THEN
!WRITE(*,'(A,5I5)')'on node',COMM_KIN%NODE_ME, IRANK_SEND, IRANK_RECV, IPOS_GLOBAL
!
! gamma only version works with double precision arrays
!
             IF ( WH%WGW%LGAMMA ) THEN
                WH%WPOT(1)%RESPONSER(1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1,1) = &
                   WH%WPOT(1)%RESPONSER(1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1,1) + &
                   WPOT_DIST%RESPONSER( 1: NP1, 1:NP2, WH%NOMEGA_LOCAL)* WEIGHT
                WH%WPOT(1)%CWINGR(1:NP1,1:3,1) = &
                   WH%WPOT(1)%CWINGR(1:NP1,1:3,1) + &
                   WPOT_DIST%CWINGR( 1:NP1,1:3,WH%NOMEGA_LOCAL)* WEIGHT
                WH%WPOT(1)%WINGR(1:NP1,1:3,1) = &
                   WH%WPOT(1)%WINGR(1:NP1,1:3,1) + &
                   WPOT_DIST%WINGR( 1:NP1,1:3,WH%NOMEGA_LOCAL)* WEIGHT
!
! full-k version works with double complex arrays
!
             ELSE
                WH%WPOT(1)%RESPONSEFUN(1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1,1) = &
                   WH%WPOT(1)%RESPONSEFUN(1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1,1) + &
                   WPOT_DIST%RESPONSEFUN( 1: NP1, 1:NP2, WH%NOMEGA_LOCAL)* WEIGHT
                WH%WPOT(1)%CWING(1:NP1,1:3,1) = &
                   WH%WPOT(1)%CWING(1:NP1,1:3,1) + &
                   WPOT_DIST%CWING( 1:NP1,1:3,WH%NOMEGA_LOCAL)* WEIGHT
                WH%WPOT(1)%WING(1:NP1,1:3,1) = &
                   WH%WPOT(1)%WING(1:NP1,1:3,1) + &
                   WPOT_DIST%WING( 1:NP1,1:3,WH%NOMEGA_LOCAL)* WEIGHT
             ENDIF
             WH%WPOT(1)%HEAD(1:3,1:3,1) = &
                WH%WPOT(1)%HEAD(1:3,1:3,1) + &
                WPOT_DIST%HEAD( 1:3, 1:3, WH%NOMEGA_LOCAL)* WEIGHT
!
! otherwise send data to rank IRANK_RECV
!
          ELSE

!
! gamma only version works with double precision arrays
!
             IF ( WH%WGW%LGAMMA ) THEN
!WRITE(*,'(A,10I6)')'sendingR',COMM_KIN%NODE_ME, IRANK_SEND, IRANK_RECV, IPOS_GLOBAL, SIZE(WPOT_DIST%RESPONSER,1),SIZE(WPOT_DIST%RESPONSER,2), 2*NP1,2*NP2
                CALL M_send_d( COMM_KIN, IRANK_RECV, WPOT_DIST%RESPONSER( 1, 1, WH%NOMEGA_LOCAL), NP1*NP2 )
                CALL M_send_d( COMM_KIN, IRANK_RECV, WPOT_DIST%CWINGR( 1,1,WH%NOMEGA_LOCAL), 3*NP1 )
                CALL M_send_d( COMM_KIN, IRANK_RECV, WPOT_DIST%WINGR( 1,1,WH%NOMEGA_LOCAL), 3*NP1 )
!
! full-k version works with double complex arrays
!
             ELSE
!WRITE(*,'(A,10I6)')'sending',COMM_KIN%NODE_ME, IRANK_SEND, IRANK_RECV, IPOS_GLOBAL,  SIZE(WPOT_DIST%RESPONSEFUN,1),SIZE(WPOT_DIST%RESPONSEFUN,2), NP1,NP2
                CALL M_send_z( COMM_KIN, IRANK_RECV, WPOT_DIST%RESPONSEFUN( 1, 1, WH%NOMEGA_LOCAL), NP1*NP2 )
                CALL M_send_z( COMM_KIN, IRANK_RECV, WPOT_DIST%CWING( 1,1,WH%NOMEGA_LOCAL), 3*NP1 )
                CALL M_send_z( COMM_KIN, IRANK_RECV, WPOT_DIST%WING( 1,1,WH%NOMEGA_LOCAL), 3*NP1 )
             ENDIF
             CALL M_send_z( COMM_KIN, IRANK_RECV, WPOT_DIST%HEAD( 1, 1, WH%NOMEGA_LOCAL), 9 )
          ENDIF
       ENDIF
       IF ( COMM_KIN%NODE_ME == IRANK_RECV .AND. IRANK_RECV /= IRANK_SEND ) THEN
!
! gamma only version works with double precision arrays
!
          IF ( WH%WGW%LGAMMA ) THEN
!WRITE(*,*) "HERE",COMM_KIN%NODE_ME
             CALL M_recv_d( COMM_KIN, IRANK_SEND, RECV_RBUFF1(1,1), NP1*NP2 )
             CALL M_recv_d( COMM_KIN, IRANK_SEND, RECV_RBUFF2(1,1), NP1*3 )
             CALL M_recv_d( COMM_KIN, IRANK_SEND, RECV_RBUFF3(1,1), NP1*3 )
! accumulate
! accumulate
             WH%WPOT(1)%RESPONSER( 1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1, 1 ) =  &
                WH%WPOT(1)%RESPONSER( 1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1, 1 ) + RECV_RBUFF1(1:NP1,1:NP2) * WEIGHT
             WH%WPOT(1)%CWING( 1:NP1,1:3, 1 ) =  &
                WH%WPOT(1)%CWING( 1:NP1,1:3, 1 ) + RECV_RBUFF2(1:NP1,1:3) * WEIGHT
             WH%WPOT(1)%WING( 1:NP1,1:3, 1 ) =  &
                WH%WPOT(1)%WING( 1:NP1,1:3, 1 ) + RECV_RBUFF3(1:NP1,1:3) * WEIGHT
!
! full-k version works with double complex arrays
!
          ELSE
!WRITE(*,'(A,8I6)')'receiving',COMM_KIN%NODE_ME, IRANK_SEND, IRANK_RECV, IPOS_GLOBAL, NP1,NP2
!WRITE(*,'(A,I4,6I6)')'allocated',COMM_KIN%NODE_ME, SIZE(RECV_CBUFF,1), SIZE( RECV_CBUFF,2), NP1, NP2
             CALL M_recv_z( COMM_KIN, IRANK_SEND, RECV_CBUFF1(1,1), NP1*NP2 )
             CALL M_recv_z( COMM_KIN, IRANK_SEND, RECV_CBUFF2(1,1), NP1*3 )
             CALL M_recv_z( COMM_KIN, IRANK_SEND, RECV_CBUFF3(1,1), NP1*3 )

! accumulate
             WH%WPOT(1)%RESPONSEFUN( 1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1, 1 ) =  &
                WH%WPOT(1)%RESPONSEFUN( 1:NP1,IPOS_GLOBAL:IPOS_GLOBAL+NP2-1, 1 ) + RECV_CBUFF1(1:NP1,1:NP2) * WEIGHT
             WH%WPOT(1)%CWING( 1:NP1,1:3, 1 ) =  &
                WH%WPOT(1)%CWING( 1:NP1,1:3, 1 ) + RECV_CBUFF2(1:NP1,1:3) * WEIGHT
             WH%WPOT(1)%WING( 1:NP1,1:3, 1 ) =  &
                WH%WPOT(1)%WING( 1:NP1,1:3, 1 ) + RECV_CBUFF3(1:NP1,1:3) * WEIGHT
          ENDIF

! receive head buffer
          CALL M_recv_z( COMM_KIN, IRANK_SEND, RECV_CBUFF4, 9 )
! accumulate  head
          WH%WPOT(1)%HEAD( 1:3,1:3, 1 ) =  &
             WH%WPOT(1)%HEAD( 1:3,1:3, 1 ) + RECV_CBUFF4(1:3,1:3) * WEIGHT
       ENDIF
       CALL M_barrier( COMM_KIN )
    ENDDO send
# 1482

    ENDDO freq_groups
    ENDDO freq_in_group
    
!> get rid of auxillary arrays
    IF ( ALLOCATED( RECV_RBUFF1 ) ) DEALLOCATE(RECV_RBUFF1 )
    IF ( ALLOCATED( RECV_RBUFF2 ) ) DEALLOCATE(RECV_RBUFF2 )
    IF ( ALLOCATED( RECV_RBUFF3 ) ) DEALLOCATE(RECV_RBUFF3 )
    IF ( ALLOCATED( RECV_CBUFF1 ) ) DEALLOCATE(RECV_CBUFF1 )
    IF ( ALLOCATED( RECV_CBUFF2 ) ) DEALLOCATE(RECV_CBUFF2 )
    IF ( ALLOCATED( RECV_CBUFF3 ) ) DEALLOCATE(RECV_CBUFF3 )
    IF ( ALLOCATED( RECV_CBUFF4 ) ) DEALLOCATE(RECV_CBUFF4 )
    

!> decide who is writing
    WH%WPOT(1)%NOMEGA_LOW = -1
    IF( COMM_KIN%NODE_ME == 1 ) THEN
       WH%WPOT(1)%NOMEGA_LOW = 1
# 1505

    ENDIF 
!write potential
!CALL WRITE_WPOT_FULL(WH%WPOT(1), WGWQ, WH%FSG0, NQ, LFOCK_SUBTRACT )
    CALL WRITE_CONJG_WPOT_FULL(WH%WPOT(1), WGWQ, WH%FSG0, NQ, LFOCK_SUBTRACT )


    CONTAINS 

!*********************************************************************
!
!> write complex conjugated of the screened potential to a file
!> this version writes the entire matrix @param WPOT(G,G')
!
!*********************************************************************

  SUBROUTINE WRITE_CONJG_WPOT_FULL(WPOT, WGWQ, FSG0, NQ, LFOCK_SUBTRACT ) 

    IMPLICIT NONE
    TYPE (responsefunction) :: WPOT  !< screened two electron potential-HF
    TYPE (wavedes1)         :: WGWQ  !< basis set descriptor
    REAL (q)                :: FSG0  !< HF singularity correction
    INTEGER              NQ          !< q-point for which local field correction are required
    LOGICAL :: LFOCK_SUBTRACT
! local
    INTEGER              IERR        ! error status
    CHARACTER (4) :: APP
    INTEGER       :: IU=72, NOMEGA, I, NP
    COMPLEX(q)    :: W( WPOT%NP2, WPOT%NP2)

    IERR=0
! write full potential at every k-point
    IF (WPOT%NOMEGA_LOW==1) THEN
       NOMEGA=1
       WRITE (APP  , "(4I1)") MOD(NQ/1000,10),MOD(NQ/100,10),MOD(NQ/10,10), MOD(NQ,10)
       OPEN( UNIT=IU, FILE="WFULL"//APP//".tmp", IOSTAT=IERR, FORM='UNFORMATTED')

       NP=WGWQ%NGVECTOR
       IF (WGWQ%LGAMMA) NP=NP*2
       IF (NP>WPOT%NP2) THEN
          CALL vtutor%bug("internal error in WRITE_CONJG_WPOT_FULL: NP>WPOT%NP2 " // str(NP) // " " // &
             str(WPOT%NP2), "wpot.F", 1546)
       ENDIF

       IF (IERR==0) WRITE(IU) NP, NP
       IF (IERR==0) WRITE(IU) CONJG( WPOT%HEAD(:,:,NOMEGA) )
       IF (WPOT%LREALSTORE) THEN
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%WINGR(1:NP,:,NOMEGA)
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) WPOT%CWINGR(1:NP,:,NOMEGA)
       ELSE
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) CONJG( WPOT%WING(1:NP,:,NOMEGA) )
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) CONJG( WPOT%CWING(1:NP,:,NOMEGA) )

# 1572

       ENDIF

       IF (WPOT%LREALSTORE) THEN
          W=WPOT%RESPONSER  (:,:,NOMEGA)
       ELSE
          W=CONJG( WPOT%RESPONSEFUN(:,:,NOMEGA) )
       ENDIF

       IF (LFOCK_SUBTRACT) THEN
          IF (WPOT%LGAMMA) THEN
             DO I=2,NP
                W(I,I)=W(I,I)+WGWQ%DATAKE(I,1)
             ENDDO
             W(1,1)=W(1,1)+FSG0
          ELSE
             DO I=1,NP
                W(I,I)=W(I,I)+WGWQ%DATAKE(I,1)
             ENDDO
          ENDIF
       ENDIF
       IF (WPOT%LREAL) THEN
! WPOT%LREAL indicates that we are using the gamma point only version
! in this case the response function is necessarily real valued
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) REAL(W(1:NP,1:NP),q)
       ELSE
          IF (IERR==0) WRITE(IU, IOSTAT=IERR) W(1:NP,1:NP)
# 1604

       ENDIF
       CLOSE(IU)
    ENDIF
    CALL M_sum_i(WGWQ%COMM_INTER, IERR, 1)

  END SUBROUTINE WRITE_CONJG_WPOT_FULL

  END SUBROUTINE ACCUMULATE_WPOT

END MODULE wpot
