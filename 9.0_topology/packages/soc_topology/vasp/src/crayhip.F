#include "symbol.inc"
#ifdef CRAYHIP
#ifdef OFFLOADING
!***********************************************************************
!
!***********************************************************************
#ifdef MPI
    MODULE mdevice_info

      USE main_mpi
      USE iso_c_binding

      PRIVATE
      PUBLIC :: GET_NUM_DEVICES_GLOBAL, GET_NUM_RANKS_ON_DEVICE, GET_FREE_MEM_DEVICE
      PUBLIC :: CHECK_MPI

      INTERFACE
         SUBROUTINE hip_device_get_uuid(id, device_num) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr),    VALUE :: id
            INTEGER(c_int), VALUE :: device_num
         END SUBROUTINE hip_device_get_uuid


         SUBROUTINE hip_mem_get_info(free, total, err) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr), VALUE :: free, total, err
         END SUBROUTINE hip_mem_get_info
      END INTERFACE

    CONTAINS

!****************** FUNCTION GET_NUM_DEVICES_GLOBAL ********************
!
!***********************************************************************

      INTEGER FUNCTION GET_NUM_DEVICES_GLOBAL(MY_DEVICE_NUM)

      INTEGER :: MY_DEVICE_NUM

      CHARACTER(KIND=c_char), TARGET :: uuidchar(16)
      INTEGER :: uuidint(16)
      INTEGER, ALLOCATABLE :: uuidints(:,:)
      INTEGER, ALLOCATABLE :: uuidunique(:,:)
      INTEGER :: totalGPUs, I, J
      LOGICAL :: isUnique

      CALL hip_device_get_uuid(c_loc(uuidchar), MY_DEVICE_NUM)
      uuidint = ICHAR(uuidchar(:))

      ALLOCATE(uuidints(16, COMM_intra_node_world%NCPU))

      uuidints = 0; uuidints(:, COMM_intra_node_world%NODE_ME) = uuidint(:)
      CALL M_sum_i(COMM_intra_node_world, uuidints(1,1), 16*COMM_intra_node_world%NCPU)

      ALLOCATE(uuidunique(16, COMM_intra_node_world%NCPU))

      uuidunique = 0 ; totalGPUs = 0
      DO I=1, COMM_intra_node_world%NCPU
         isUnique = .TRUE.
         DO J=1,totalGPUs
            IF ( ALL( uuidints(:,I) == uuidunique(:,J) ) ) THEN
               isUnique = .FALSE.
               EXIT
            ENDIF
         ENDDO
         IF ( isUnique ) THEN
            totalGPUs = totalGPUs + 1
            uuidunique(:,totalGPUs) = uuidints(:,I)
         ENDIF
      ENDDO

      DEALLOCATE(uuidints, uuidunique)

      CALL M_sum_i(COMM_inter_node_world, totalGPUs, 1)

      GET_NUM_DEVICES_GLOBAL = totalGPUs

      END FUNCTION GET_NUM_DEVICES_GLOBAL


!****************** FUNCTION GET_NUM_RANKS_ON_DEVICE *******************
!
!***********************************************************************

      INTEGER FUNCTION GET_NUM_RANKS_ON_DEVICE(MY_DEVICE_NUM)

      INTEGER :: MY_DEVICE_NUM

      CHARACTER(KIND=c_char), TARGET :: uuidchar(16)
      INTEGER :: uuidint(16)
      INTEGER, ALLOCATABLE :: uuidints(:,:)
      INTEGER :: numranks, I

      CALL hip_device_get_uuid(c_loc(uuidchar), MY_DEVICE_NUM)
      uuidint = ICHAR(uuidchar(:))

      ALLOCATE(uuidints(16, COMM_intra_node_world%NCPU))

      uuidints = 0; uuidints(:, COMM_intra_node_world%NODE_ME) = uuidint(:)
      CALL M_sum_i(COMM_intra_node_world, uuidints(1,1), 16*COMM_intra_node_world%NCPU)

      numranks = 0

      DO I=1, COMM_intra_node_world%NCPU
         IF ( ALL( uuidints(:,I) == uuidint(:) ) ) numranks = numranks + 1
      ENDDO

      DEALLOCATE(uuidints)

      GET_NUM_RANKS_ON_DEVICE = numranks

      END FUNCTION GET_NUM_RANKS_ON_DEVICE


!****************** FUNCTION GET_FREE_MEM_DEVICE ***********************
!
!***********************************************************************

      FUNCTION GET_FREE_MEM_DEVICE(MFREE) RESULT(IERR)

      REAL(q)         :: MFREE
      INTEGER, TARGET :: IERR
      ! local variables
      INTEGER(c_size_t), TARGET :: FREE, TOTAL

      IERR = -1; MFREE = 0

      CALL hip_mem_get_info(c_loc(FREE), c_loc(TOTAL), c_loc(IERR))

      MFREE = REAL(FREE, KIND=q)

      END FUNCTION GET_FREE_MEM_DEVICE


!****************** SUBROUTINE CHECK_MPI *******************************
!
!***********************************************************************

      SUBROUTINE CHECK_MPI()
      RETURN
      END SUBROUTINE CHECK_MPI

    END MODULE mdevice_info
#endif // MPI

!***********************************************************************
!
!***********************************************************************

    MODULE mrocfft

      USE prec
      USE moffload_struct_def
      USE mgrid_struct_def, ONLY : grid_3d

GPACC USE openacc, ONLY : ACC_ASYNC_SYNC, ACC_ASYNC_NOVAL, acc_get_cuda_stream

      USE tutor
      USE string

      USE iso_c_binding

      IMPLICIT NONE

      PUBLIC :: FFT_OFFLOAD_MAKEPLANS, FFT_OFFLOAD_MAKEPLAN, FFT_OFFLOAD_DESTROYPLAN, FFT_OFFLOAD_GETPLAN
      PUBLIC :: fft_offload_plan

      PUBLIC :: FFTBAS_PLAN_OFFLOAD, FFTBRC_PLAN_OFFLOAD
      PUBLIC :: FFTBAS_PLAN_MU_OFFLOAD, FFTBRC_PLAN_MU_OFFLOAD
      PUBLIC :: FFTBAS_PLAN_MPI_OFFLOAD, FFTBRC_PLAN_MPI_OFFLOAD

      TYPE, BIND(C):: plan_type
         TYPE(c_ptr) :: plan_f, plan_i
         TYPE(c_ptr) :: desc_f, desc_i
         TYPE(c_ptr) :: info_f, info_i
         TYPE(c_ptr) :: work
         INTEGER :: Asize, Bsize
      END TYPE plan_type

      TYPE fft_offload_plan
         TYPE(grid_3d), POINTER :: GRID
         INTEGER :: stream
         INTEGER :: nbatch
         INTEGER :: nstride
         CHARACTER(LEN=:), ALLOCATABLE :: tag
         ! plans for serial FFTs (3D)
         TYPE(plan_type) :: FFT_3D_Z2Z
         TYPE(plan_type) :: FFT_3D_D2Z
         ! plans for MPI distributed FFTs (1D)
         TYPE(plan_type) :: FFT_RC_Z2Z
         TYPE(plan_type) :: FFT_IN_Z2Z
         TYPE(plan_type) :: FFT_RL_Z2Z
         TYPE(plan_type) :: FFT_RL_D2Z
      END TYPE fft_offload_plan

      TYPE(fft_offload_plan), ALLOCATABLE :: plans(:)

      INTEGER :: createdPlans = 0

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE planAssign
      END INTERFACE

      INTERFACE OPERATOR(==)
         MODULE PROCEDURE planCompare
      END INTERFACE

      INTERFACE
#ifdef OMP_OFFLOAD
         SUBROUTINE hip_rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist, &
                                      plan_f,plan_i,desc_f,desc_i,info_f,info_i,work) BIND(C, NAME="hip_rocfft_create_omp")
            USE iso_c_binding
            INTEGER(c_int), VALUE :: ttype,tdim,tnum,idist,odist
            TYPE(c_ptr), VALUE :: tsizes,istrides,ostrides
            TYPE(c_ptr) :: plan_f,plan_i,desc_f,desc_i,info_f,info_i
            TYPE(c_ptr) :: work
         END SUBROUTINE hip_rocfft_create
#endif
#ifdef ACC_OFFLOAD
         SUBROUTINE hip_rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist, &
                                      plan_f,plan_i,desc_f,desc_i,info_f,info_i,work,stream) BIND(C, NAME="hip_rocfft_create_acc")
            USE iso_c_binding
            INTEGER(c_int), VALUE :: ttype,tdim,tnum,idist,odist
            TYPE(c_ptr), VALUE :: tsizes,istrides,ostrides
            TYPE(c_ptr) :: plan_f,plan_i,desc_f,desc_i,info_f,info_i
            TYPE(c_ptr) :: work
            TYPE(c_ptr) :: stream
         END SUBROUTINE hip_rocfft_create
#endif

         SUBROUTINE hip_rocfft_destroy(plan_f,plan_i,desc_f,desc_i,info_f,info_i,work) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr) :: plan_f,plan_i,desc_f,desc_i,info_f,info_i
            TYPE(c_ptr) :: work
         END SUBROUTINE hip_rocfft_destroy


         SUBROUTINE hip_rocfft_exec(plan,A,B,info) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr) :: plan,info
            TYPE(c_ptr), VALUE :: A,B
         END SUBROUTINE hip_rocfft_exec
      END INTERFACE

    CONTAINS

      SUBROUTINE rocfft_exec_fwd(plan,A,B)
      TYPE(plan_type) :: plan
      COMPLEX(q), TARGET :: A(*)
      COMPLEX(q), TARGET :: B(*)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B)
GPACC !$ACC HOST_DATA USE_DEVICE(A,B)
      CALL hip_rocfft_exec(plan%plan_f,c_loc(A),c_loc(B),plan%info_f)
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
      END SUBROUTINE rocfft_exec_fwd


      SUBROUTINE rocfft_exec_bck(plan,A,B)
      TYPE(plan_type) :: plan
      COMPLEX(q), TARGET :: A(*)
      COMPLEX(q), TARGET :: B(*)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B)
GPACC !$ACC HOST_DATA USE_DEVICE(A,B)
      CALL hip_rocfft_exec(plan%plan_i,c_loc(A),c_loc(B),plan%info_i)
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
      END SUBROUTINE rocfft_exec_bck


      SUBROUTINE rocfft_exec_d2z(plan,A,B)
      TYPE(plan_type) :: plan
      REAL(q)   , TARGET :: A(*)
      COMPLEX(q), TARGET :: B(*)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B)
GPACC !$ACC HOST_DATA USE_DEVICE(A,B)
      CALL hip_rocfft_exec(plan%plan_f,c_loc(A),c_loc(B),plan%info_f)
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
      END SUBROUTINE rocfft_exec_d2z


      SUBROUTINE rocfft_exec_z2d(plan,A,B)
      TYPE(plan_type) :: plan
      COMPLEX(q), TARGET :: A(*)
      REAL(q)   , TARGET :: B(*)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B)
GPACC !$ACC HOST_DATA USE_DEVICE(A,B)
      CALL hip_rocfft_exec(plan%plan_i,c_loc(A),c_loc(B),plan%info_i)
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
      END SUBROUTINE rocfft_exec_z2d


      SUBROUTINE rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist,istream,plan)
      INTEGER :: ttype,tdim,tnum,idist,odist,tsizes(3),istrides(3),ostrides(3),istream
      TYPE(plan_type) :: plan
#ifdef OMP_OFFLOAD
      CALL hip_rocfft_create(ttype, tdim, tnum, c_loc(tsizes), c_loc(istrides), idist, c_loc(ostrides), odist, &
                             plan%plan_f,plan%plan_i,plan%desc_f,plan%desc_i,plan%info_f,plan%info_i,plan%work)
#endif
#ifdef ACC_OFFLOAD
      CALL hip_rocfft_create(ttype, tdim, tnum, c_loc(tsizes), c_loc(istrides), idist, c_loc(ostrides), odist, &
                             plan%plan_f,plan%plan_i,plan%desc_f,plan%desc_i,plan%info_f,plan%info_i,plan%work, &
                             acc_get_cuda_stream(istream))
#endif
      END SUBROUTINE rocfft_create


      SUBROUTINE rocfft_destroy(plan)
      TYPE(plan_type) :: plan
      CALL hip_rocfft_destroy(plan%plan_f,plan%plan_i,plan%desc_f,plan%desc_i,plan%info_f,plan%info_i,plan%work)
      END SUBROUTINE rocfft_destroy


!****************** SUBROUTINE FFT_OFFLOAD_MAKEPLANS *******************
!
!***********************************************************************

      SUBROUTINE FFT_OFFLOAD_MAKEPLANS(GRID, NSTREAMS)
      TYPE(grid_3d), INTENT(IN), TARGET :: GRID
      INTEGER, INTENT(IN), OPTIONAL :: NSTREAMS
      ! local variables
GPACC INTEGER ::N, I

      PROFILING_START('fft_offload_makeplans')

      IF (.NOT.OFFLOAD_ALLOW) THEN
         PROFILING_STOP('fft_offload_makeplans')
         RETURN
      ENDIF
#ifdef ACC_OFFLOAD
      ! Always create a plan for ACC_ASYNC_SYNC
      CALL FFT_OFFLOAD_MAKEPLAN(GRID,INT(ACC_ASYNC_SYNC))

      ! And for ACC_ASYNC_ASYNC
      CALL FFT_OFFLOAD_MAKEPLAN(GRID,ACC_ASYNC_ASYNC)

      ! Then for additional streams if requested
      N=1; IF (PRESENT(NSTREAMS)) N=NSTREAMS
      DO I=1,N
         CALL FFT_OFFLOAD_MAKEPLAN(GRID,I)
      ENDDO

      ! And last but not least for ACC_ASYNC_NOVAL
      CALL FFT_OFFLOAD_MAKEPLAN(GRID,INT(ACC_ASYNC_NOVAL))
#endif
#ifdef OMP_OFFLOAD
      CALL FFT_OFFLOAD_MAKEPLAN(GRID,-1)
#endif
      PROFILING_STOP('fft_offload_makeplans')

      END SUBROUTINE FFT_OFFLOAD_MAKEPLANS


!****************** SUBROUTINE FFT_OFFLOAD_MAKEPLAN ********************
!
!***********************************************************************

      SUBROUTINE FFT_OFFLOAD_MAKEPLAN(GRID, ISTREAM, NBATCH, NSTRIDE, TAG)
      TYPE(grid_3d), INTENT(IN), TARGET      :: GRID
      INTEGER, INTENT(IN)                    :: ISTREAM
      INTEGER, INTENT(IN), OPTIONAL          :: NBATCH
      INTEGER, INTENT(IN), OPTIONAL          :: NSTRIDE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: TAG

      ! local variables
      TYPE(fft_offload_plan) :: NEW

      PROFILING_START('fft_offload_makeplan')

      IF (.NOT.OFFLOAD_ALLOW) THEN
         PROFILING_STOP('fft_offload_makeplan')
         RETURN
      ENDIF

      NEW%GRID => GRID

GPOMP NEW%STREAM = -1
GPACC NEW%STREAM = ISTREAM

      NEW%NBATCH = 1 ; IF (PRESENT(NBATCH)) NEW%NBATCH = MAX(NBATCH,1)
      NEW%NSTRIDE=-1 ; IF (PRESENT(NSTRIDE)) NEW%NSTRIDE = NSTRIDE

      ! when we already have a plan for this GRID+ISTREAM+NBATCH we may exit here
      IF (planFind(NEW)/=0) THEN
         PROFILING_STOP('fft_offload_makeplan')
         RETURN
      ENDIF

      NEW%TAG = "" ; IF (PRESENT(TAG)) NEW%TAG = TAG

      CALL MAKEPLAN(NEW)

      CALL planAdd(NEW)

      PROFILING_STOP('fft_offload_makeplan')

      END SUBROUTINE FFT_OFFLOAD_MAKEPLAN


      SUBROUTINE MAKEPLAN(PLAN)
      IMPLICIT NONE
      TYPE(fft_offload_plan), TARGET :: PLAN
      ! local variables
      TYPE(plan_type), POINTER :: ptype
      INTEGER :: tsizes(3),istrides(3),ostrides(3),idist,odist,tdim,tnum,ttype,istream
      INTEGER :: NG(3), BATCHES(3), LD

GPOMP IF (PLAN%NBATCH>1) &
GPOMP    CALL vtutor%bug("MAKEPLAN: batched FFTs not implemented yet", __FILE__, __LINE__)

      istream = PLAN%STREAM

      NG = PLAN%GRID%NGPTAR

      ! serial FFT (3D)
      IF (PLAN%GRID%RL%NFAST==1 .AND. PLAN%GRID%RL_FFT%NFAST==1) THEN

         ! single FFT
         IF (PLAN%NBATCH==1) THEN
            ! real-2-complex
            IF (PLAN%GRID%REAL2CPLX) THEN
               ! FFT_3D_D2Z: 3d real-2-complex/complex-2-real
               ptype => PLAN%FFT_3D_D2Z
               ttype = 1; tdim = 3; tnum = 1;
               tsizes(1) = NG(1); tsizes(2) = NG(2); tsizes(3) = NG(3)
               idist = 0;
               istrides(1) = 1;
               istrides(2) = NG(1)+2;
               istrides(3) = (NG(1)+2)*NG(2)
               odist = 0;
               ostrides(1) = 1;
               ostrides(2) = (NG(1)+2)/2;
               ostrides(3) = (NG(1)+2)/2*NG(2)
               ptype%Asize = ((NG(1)+2)*NG(2)*NG(3))
               ptype%Bsize = ((NG(1)+2)/2*NG(2)*NG(3))
            ! complex-2-complex
            ELSE
               ! FFT_3D_Z2Z: 3d complex-2-complex
               ptype => PLAN%FFT_3D_Z2Z
               ttype = 0; tdim = 3; tnum = 1;
               tsizes(1) = NG(1); tsizes(2) = NG(2); tsizes(3) = NG(3)
               idist = 0;
               istrides(1) = 1;
               istrides(2) = NG(1);
               istrides(3) = NG(1)*NG(2)
               odist = 0;
               ostrides(1) = 1;
               ostrides(2) = NG(1);
               ostrides(3) = NG(1)*NG(2)
               ptype%Asize = NG(1)*NG(2)*NG(3)
               ptype%Bsize = NG(1)*NG(2)*NG(3)
            ENDIF
         ! batched FFT
         ELSE
!!            CALL vtutor%bug("MAKEPLAN: batched serial FFTs not implemented yet", __FILE__, __LINE__)
            LD = PLAN%GRID%MPLWV; IF (PLAN%NSTRIDE > 0) LD = PLAN%NSTRIDE

            IF (PLAN%GRID%REAL2CPLX) THEN
               ptype => PLAN%FFT_3D_D2Z
               ttype = 1; tdim = 3; tnum = PLAN%NBATCH
               tsizes(1) = NG(1); tsizes(2) = NG(2); tsizes(3) = NG(3)
               idist = 2*LD;
               istrides(1) = 1;
               istrides(2) = NG(1);
               istrides(3) = NG(1)*NG(2)
               odist = 2*LD;
               ostrides(1) = 1;
               ostrides(2) = (NG(1)+2)/2;
               ostrides(3) = (NG(1)+2)/2*NG(2)
               ptype%Asize = NG(1)*NG(2)*NG(3)
               ptype%Bsize = (NG(1)+2)/2*NG(2)*NG(3)
            ELSE
               ptype => PLAN%FFT_3D_Z2Z
               ttype = 0; tdim = 3; tnum = PLAN%NBATCH
               tsizes(1) = NG(1); tsizes(2) = NG(2); tsizes(3) = NG(3)
               idist = LD;
               istrides(1) = 1;
               istrides(2) = NG(1);
               istrides(3) = NG(1)*NG(2)
               odist = LD;
               ostrides(1) = 1;
               ostrides(2) = NG(1);
               ostrides(3) = NG(1)*NG(2)
               ptype%Asize = NG(1)*NG(2)*NG(3)
               ptype%Bsize = NG(1)*NG(2)*NG(3)
            ENDIF
         ENDIF

         CALL rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist,istream,ptype)

      ! parallel FFT
      ELSE
         IF (PLAN%NBATCH/=1) THEN
            CALL vtutor%bug("MAKEPLAN: batched parallel ffts not supported yet",__FILE__,__LINE__)
         ENDIF

         BATCHES(1)=PLAN%GRID%RC%NCOL
         BATCHES(2)=PLAN%GRID%IN%NCOL
         BATCHES(3)=PLAN%GRID%RL_FFT%NCOL

         ! FFT_RC_Z2Z: 1d complex-2-complex FFTs on data with reciprocal space layout
         IF (BATCHES(1)>0) THEN
            ptype => PLAN%FFT_RC_Z2Z
            ttype = 0; tdim = 1; tnum = BATCHES(1);
            tsizes(1) = NG(1);
            idist = NG(1);
            istrides(1) = 1;
            odist = NG(1);
            ostrides(1) = 1;
            ptype%Asize = idist*tnum
            ptype%Bsize = odist*tnum
            CALL rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist,istream,ptype)
         ENDIF
         ! FFT_IN_Z2Z: 1d complex-2-complex FFTs on data with intermediate layout
         IF (BATCHES(2)>0) THEN
             ptype => PLAN%FFT_IN_Z2Z
             ttype = 0; tdim = 1; tnum = BATCHES(2);
             tsizes(1) = NG(2);
             idist = 1;
             istrides(1) = batches(2);
             odist = 1;
             ostrides(1) = batches(2);
             ptype%Asize = NG(2)*tnum
             ptype%Bsize = NG(2)*tnum
             CALL rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist,istream,ptype)
         ENDIF
         ! 1d FFTs on data with real space layout
         IF (BATCHES(3)>0) THEN
            ! FFT_RL_D2Z: real-2-complex/complex-2-real
            IF (PLAN%GRID%REAL2CPLX) THEN
               ptype => PLAN%FFT_RL_D2Z
               ttype = 1; tdim = 1; tnum = BATCHES(3);
               tsizes(1) = NG(3);
               idist = NG(3)+2;
               istrides(1) = 1;
               odist = (NG(3)+2)/2;
               ostrides(1) = 1;
               ptype%Asize = (NG(3)+2)*tnum
               ptype%Bsize = (NG(3)+2)/2*tnum
               CALL rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist,istream,ptype)
            ! FFT_RL_Z2Z: complex-2-complex
            ELSE
               ptype => PLAN%FFT_RL_Z2Z
               ttype = 0; tdim = 1; tnum = BATCHES(3);
               tsizes(1) = NG(3);
               idist = NG(3);
               istrides(1) = 1;
               odist = NG(3);
               ostrides(1) = 1;
               ptype%Asize = NG(3)*tnum
               ptype%Bsize = NG(3)*tnum
               CALL rocfft_create(ttype,tdim,tnum,tsizes,istrides,idist,ostrides,odist,istream,ptype)
            ENDIF
         ENDIF
      ENDIF
      END SUBROUTINE MAKEPLAN


!****************** SUBROUTINE FFT_OFFLOAD_DESTROYPLAN *****************
!
!***********************************************************************

      SUBROUTINE FFT_OFFLOAD_DESTROYPLAN(GRID,ISTREAM,NBATCH,TAG)
      TYPE(grid_3d), INTENT(IN), TARGET      :: GRID
      INTEGER, INTENT(IN), OPTIONAL          :: ISTREAM
      INTEGER, INTENT(IN), OPTIONAL          :: NBATCH
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: TAG
      ! local variables
      TYPE(fft_offload_plan), ALLOCATABLE :: tempPlans(:)
      INTEGER :: i,j
      LOGICAL :: keepPlan

      IF (.NOT.OFFLOAD_ALLOW) RETURN

      j=0

      ALLOCATE(tempPlans(createdPlans))

      DO i=1,createdPlans
         keepPlan = .FALSE.

         IF (.NOT.ASSOCIATED(plans(i)%GRID,GRID)) keepPLan = .TRUE.

         IF (PRESENT(ISTREAM)) THEN
            IF (plans(i)%stream/=ISTREAM) keepPlan = .TRUE.
         ENDIF

         IF (PRESENT(NBATCH)) THEN
            IF (plans(i)%NBATCH/=MAX(NBATCH,1)) keepPlan = .TRUE.
         ENDIF

         IF (PRESENT(TAG)) keepPlan = (plans(i)%TAG/=TAG)

         IF (keepPlan) THEN
            j=j+1; tempPlans(j)=plans(i)
         ELSE
            CALL DESTROYPLAN(plans(i))
         ENDIF
      ENDDO

      DEALLOCATE(plans)

      IF (j/=0) THEN
         ALLOCATE(plans(j))
         DO i=1,j
            plans(i)=tempPlans(i)
         ENDDO
      ENDIF

      DEALLOCATE(tempPlans)

      createdPlans=j

      END SUBROUTINE FFT_OFFLOAD_DESTROYPLAN


      SUBROUTINE DESTROYPLAN(PLAN)
      IMPLICIT NONE
      TYPE(fft_offload_plan) :: PLAN

      IF (PLAN%GRID%RL%NFAST==1) THEN
         ! serial fft
         IF (PLAN%GRID%REAL2CPLX) THEN
            CALL rocfft_destroy(PLAN%FFT_3D_D2Z)
         ELSE
            CALL rocfft_destroy(PLAN%FFT_3D_Z2Z)
         ENDIF
      ELSE
         IF (PLAN%GRID%RC%NCOL>0) CALL rocfft_destroy(PLAN%FFT_RC_Z2Z)
         IF (PLAN%GRID%IN%NCOL>0) CALL rocfft_destroy(PLAN%FFT_IN_Z2Z)
         IF (PLAN%GRID%RL_FFT%NCOL>0) THEN
            IF (PLAN%GRID%REAL2CPLX) THEN
               CALL rocfft_destroy(PLAN%FFT_RL_D2Z)
            ELSE
               CALL rocfft_destroy(PLAN%FFT_RL_Z2Z)
            ENDIF
         ENDIF
      ENDIF
      END SUBROUTINE DESTROYPLAN


!****************** SUBROUTINE FFT_OFFLOAD_GETPLAN *********************
!
!***********************************************************************

      FUNCTION FFT_OFFLOAD_GETPLAN(GRID,ISTREAM,NBATCH,RES) RESULT(ISTAT)
      TYPE(grid_3d), INTENT(IN), TARGET :: GRID
      INTEGER, INTENT(IN) :: ISTREAM
      INTEGER, INTENT(IN) :: NBATCH
      TYPE(fft_offload_plan), INTENT(OUT) :: RES
      INTEGER :: ISTAT
      ! local variables
      TYPE(fft_offload_plan) :: q
      INTEGER :: i

      q%GRID   =>GRID
      q%stream = ISTREAM
      q%nbatch = NBATCH

      ISTAT=0

      ! search matching plan
      i=planFind(q)

      IF (i/=0) THEN
         ! return matching plan
         RES=plans(i)
      ELSE
         ! no match found: return in error
         ISTAT=-1
      ENDIF
      END FUNCTION FFT_OFFLOAD_GETPLAN


!****************** SUBROUTINE PLANADD *********************************
!
!***********************************************************************

      SUBROUTINE planAdd(a)
      TYPE(fft_offload_plan), INTENT(IN) :: a
      ! local variables
      INTEGER :: i

      ! Early exit if plan already exists
      IF (planFind(a)/=0) RETURN

      ! Find storage space for a new plan
      i=planAlloc()

      ! Add plan
      plans(i)=a
      END SUBROUTINE planAdd


!****************** FUNCTION PLANALLOC *********************************
!
!***********************************************************************

      FUNCTION planAlloc()
      INTEGER :: planAlloc
      ! local variables
      TYPE(fft_offload_plan), ALLOCATABLE :: tempPlans(:)
      INTEGER :: i,oldSize

      ! on first call the plans array is not allocated
      IF (.NOT.ALLOCATED(plans)) THEN
         ALLOCATE(plans(1))
         NULLIFY(plans(1)%grid)
         createdPlans=0

      ! when the plans array is full
      ELSEIF(createdPlans==SIZE(plans)) THEN
         oldSize=SIZE(plans)
         ALLOCATE(tempPlans(oldsize))
         DO i=1,oldSize
            tempPlans(i)=plans(i)
         ENDDO
         DEALLOCATE(plans)
         ALLOCATE(plans(oldsize+1))
         DO i=1,oldSize
            plans(i)=tempPlans(i)
         ENDDO
         DEALLOCATE(tempPlans)
         NULLIFY(plans(oldsize+1)%grid)

      ! error .. this may not happen
      ELSEIF(createdPlans>SIZE(plans)) THEN
         STOP
      ENDIF

      ! advance the number of created plans by one
      createdPlans=createdPlans+1

      ! and pass the index back
      planAlloc=createdPlans
      END FUNCTION planAlloc


!****************** FUNCTION PLANFIND **********************************
!
!***********************************************************************

      FUNCTION planFind(a)
      TYPE(fft_offload_plan), INTENT(IN) :: a
      INTEGER :: planFind
      ! local variables
      INTEGER :: i

      planFind=0
      DO i=1,createdPlans
         IF (plans(i)==a) THEN
            planFind=i; RETURN
         ENDIF
      ENDDO
      END FUNCTION planFind


!****************** SUBROUTINE PLANASSIGN ******************************
!
!***********************************************************************

      SUBROUTINE planAssign(a,b)
      TYPE(fft_offload_plan), INTENT(OUT) :: a
      TYPE(fft_offload_plan), INTENT(IN)  :: b

      a%GRID          =>b%GRID
      a%stream        = b%stream
      a%nbatch        = b%nbatch
      a%nstride       = b%nstride
      a%tag           = b%tag

      a%FFT_3D_Z2Z    = b%FFT_3D_Z2Z
      a%FFT_3D_D2Z    = b%FFT_3D_D2Z

      a%FFT_RC_Z2Z    = b%FFT_RC_Z2Z
      a%FFT_IN_Z2Z    = b%FFT_IN_Z2Z
      a%FFT_RL_Z2Z    = b%FFT_RL_Z2Z
      a%FFT_RL_D2Z    = b%FFT_RL_D2Z
      END SUBROUTINE planAssign


!****************** FUNCTION PLANCOMPARE *******************************
!
!***********************************************************************

      FUNCTION planCompare(a,b) RESULT(res)
      TYPE(fft_offload_plan), INTENT(IN) :: a,b
      LOGICAL :: res

      res=.FALSE.

      IF (.NOT.ASSOCIATED(a%GRID,b%GRID)) RETURN
      IF (a%stream/=b%stream) RETURN
      IF (a%nbatch/=b%nbatch) RETURN

      res=.TRUE.
      END FUNCTION planCompare


!************************* SUBROUTINE FFTBAS_PLAN_OFFLOAD **************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_OFFLOAD(C,GRID,ISIGN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      PROFILING_START('fftbas_plan_offload')

GPOMP IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
GPACC IERR = FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBAS_PLAN_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      IF (ISIGN<=0) THEN
         CALL rocfft_exec_fwd(PLAN%FFT_3D_Z2Z,C,C)
      ELSE
         CALL rocfft_exec_bck(PLAN%FFT_3D_Z2Z,C,C)
      ENDIF

      PROFILING_STOP('fftbas_plan_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_OFFLOAD


!************************* SUBROUTINE FFTBRC_PLAN_OFFLOAD **************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_OFFLOAD(C,GRID,ISIGN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q), TARGET :: C(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      REAL(q), POINTER :: A(:)

      PROFILING_START('fftbrc_plan_offload')

      CALL c_f_pointer(c_loc(C),A,[2*GRID%MPLWV])

GPOMP IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
GPACC IERR = FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBRC_PLAN_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      IF (ISIGN<=0) THEN
         CALL rocfft_exec_d2z(PLAN%FFT_3D_D2Z,A,C)
      ELSE
         CALL rocfft_exec_z2d(PLAN%FFT_3D_D2Z,C,A)
      ENDIF

      PROFILING_STOP('fftbrc_plan_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_OFFLOAD


!****************** SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD ******************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD(N,C,LDC,GRID,ISIGN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

!!      CALL vtutor%bug("FFTBRC_PLAN_MU_OFFLOAD: not implemented yet", __FILE__, __LINE__)

      PROFILING_START('fftbas_plan_mu_offload')

GPOMP IERR=FFT_OFFLOAD_GETPLAN(GRID,-1,N,PLAN)
GPACC IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),N,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBAS_PLAN_MU_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      IF (ISIGN<=0) THEN
         CALL rocfft_exec_fwd(PLAN%FFT_3D_Z2Z,C,C)
      ELSE
         CALL rocfft_exec_bck(PLAN%FFT_3D_Z2Z,C,C)
      ENDIF

      PROFILING_STOP('fftbas_plan_mu_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_MU_OFFLOAD


!****************** SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD ******************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD(N,C,LDC,GRID,ISIGN)
      USE tutor, ONLY : vtutor

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q), TARGET :: C(*)
      INTEGER :: N,LDC,ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      REAL(q), POINTER :: A(:)

!!      CALL vtutor%bug("FFTBRC_PLAN_OFFLOAD: not implemented yet", __FILE__, __LINE__)

      PROFILING_START('fftbrc_plan_mu_offload')

GPOMP IERR=FFT_OFFLOAD_GETPLAN(GRID,-1,N,PLAN)
GPACC IERR=FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),N,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBRC_PLAN_MU_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      CALL c_f_pointer(c_loc(C),A,[2*LDC*N])

      IF (ISIGN<=0) THEN
         CALL rocfft_exec_d2z(PLAN%FFT_3D_D2Z,A,C)
      ELSE
         CALL rocfft_exec_z2d(PLAN%FFT_3D_D2Z,C,A)
      ENDIF

      PROFILING_STOP('fftbrc_plan_mu_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_MU_OFFLOAD


!****************** SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD *****************
!
!***********************************************************************

      SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD(A,GRID,ISIGN)
      USE smart_allocate
      USE tutor, ONLY : vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      REAL(q), TARGET :: A(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      COMPLEX(q), POINTER :: C(:)
      COMPLEX(q), POINTER, SAVE :: RCVBUF(:),SNDBUF(:)

      PROFILING_START('fftbas_plan_mpi_offload')

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,GRID%MPLWV)
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,GRID%MPLWV)

GPOMP IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
GPACC IERR = FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBAS_PLAN_MPI_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      CALL c_f_pointer(c_loc(A),C,[GRID%MPLWV])

!=======================================================================
! transformation (q->r)
!=======================================================================
      IF (ISIGN == 1) THEN
         IF (GRID%RC%NCOL>0) CALL rocfft_exec_bck(PLAN%FFT_RC_Z2Z,C,C)

         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%IN%NCOL>0) CALL rocfft_exec_bck(PLAN%FFT_IN_Z2Z,C,C)

         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%RL_FFT%NCOL>0) CALL rocfft_exec_bck(PLAN%FFT_RL_Z2Z,C,C)

!=======================================================================
! transformation (r->q)
!=======================================================================
      ELSE
         IF (GRID%RL_FFT%NCOL>0)  CALL rocfft_exec_fwd(PLAN%FFT_RL_Z2Z,C,C)

         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%IN%NCOL>0) CALL rocfft_exec_fwd(PLAN%FFT_IN_Z2Z,C,C)

         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%RC%NCOL>0) CALL rocfft_exec_fwd(PLAN%FFT_RC_Z2Z,C,C)
      ENDIF

      PROFILING_STOP('fftbas_plan_mpi_offload')

      RETURN
      END SUBROUTINE FFTBAS_PLAN_MPI_OFFLOAD


!****************** SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD *****************
!
!***********************************************************************

      SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD(A,GRID,ISIGN)
      USE smart_allocate
      USE tutor, ONLY : vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      REAL(q), TARGET :: A(*)
      INTEGER :: ISIGN

      ! local variables
      TYPE(fft_offload_plan) :: PLAN
      INTEGER :: IERR

      COMPLEX(q), POINTER :: C(:)
      COMPLEX(q), POINTER, SAVE :: RCVBUF(:),SNDBUF(:)

      PROFILING_START('fftbrc_plan_mpi_offload')

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,GRID%MPLWV)
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,GRID%MPLWV)

GPOMP IERR = FFT_OFFLOAD_GETPLAN(GRID,-1,1,PLAN)
GPACC IERR = FFT_OFFLOAD_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,PLAN)
      IF (IERR/=0) &
         CALL vtutor%bug("FFTBRC_PLAN_MPI_OFFLOAD: no FFT plan found", __FILE__, __LINE__)

      CALL c_f_pointer(c_loc(A),C,[GRID%MPLWV])

!=======================================================================
! transformation (q->r)
!=======================================================================
      IF (ISIGN == 1) THEN
         IF (GRID%RC%NCOL>0) CALL rocfft_exec_bck(PLAN%FFT_RC_Z2Z,C,C)

         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%IN%NCOL>0) CALL rocfft_exec_bck(PLAN%FFT_IN_Z2Z,C,C)

         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%RL_FFT%NCOL>0) CALL rocfft_exec_z2d(PLAN%FFT_RL_D2Z,C,A)

!=======================================================================
! transformation (r->q)
!=======================================================================
      ELSE
         IF (GRID%RL_FFT%NCOL>0)  CALL rocfft_exec_d2z(PLAN%FFT_RL_D2Z,A,C)

         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

         IF (GRID%IN%NCOL>0) CALL rocfft_exec_fwd(PLAN%FFT_IN_Z2Z,C,C)

         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

         IF (GRID%RC%NCOL>0) CALL rocfft_exec_fwd(PLAN%FFT_RC_Z2Z,C,C)
      ENDIF

      PROFILING_STOP('fftbrc_plan_mpi_offload')

      RETURN
      END SUBROUTINE FFTBRC_PLAN_MPI_OFFLOAD

    END MODULE mrocfft


!***********************************************************************
!
!***********************************************************************

    MODULE mrocblas
      USE prec
      USE moffload_struct_def

GPACC USE openacc, ONLY : ACC_ASYNC_SYNC, ACC_ASYNC_NOVAL, acc_get_cuda_stream

      USE tutor, ONLY : vtutor

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(c_ptr) :: ROCBLAS_HANDLE

      INTERFACE

         SUBROUTINE HIP_INIT_ROCBLAS(PTR) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
         END SUBROUTINE HIP_INIT_ROCBLAS


         SUBROUTINE HIP_FINISH_ROCBLAS(PTR) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
         END SUBROUTINE HIP_FINISH_ROCBLAS


GPACC    SUBROUTINE HIP_ROCBLAS_SET_STREAM(PTR, STREAM) BIND(C)
GPACC       USE iso_c_binding
GPACC       TYPE(c_ptr)                      :: PTR, STREAM
GPACC    END SUBROUTINE HIP_ROCBLAS_SET_STREAM


         SUBROUTINE HIP_DCOPY(PTR, N, DX, INCX, DY, INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            TYPE(c_ptr), VALUE               :: DX, DY
         END SUBROUTINE HIP_DCOPY


         SUBROUTINE HIP_DCOPY_STRIDED_BATCHED(PTR, N, DX, INCX, SX, DY, &
                INCY, SY, BC) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY, SX, SY, BC
            TYPE(c_ptr), VALUE               :: DX, DY
         END SUBROUTINE HIP_DCOPY_STRIDED_BATCHED


         SUBROUTINE HIP_ZCOPY_STRIDED_BATCHED(PTR, N, ZX, INCX, SX, ZY, &
                INCY, SY, BC) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY, SX, SY, BC
            TYPE(c_ptr), VALUE               :: ZX, ZY
         END SUBROUTINE HIP_ZCOPY_STRIDED_BATCHED


         SUBROUTINE HIP_ZCOPY(PTR, N, ZX, INCX, ZY, INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            TYPE(c_ptr), VALUE               :: ZX, ZY
         END SUBROUTINE HIP_ZCOPY


         SUBROUTINE HIP_SCOPY(PTR, N, DX, INCX, DY, INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            TYPE(c_ptr), VALUE               :: DX, DY
         END SUBROUTINE HIP_SCOPY


         SUBROUTINE HIP_CCOPY(PTR, N, DX, INCX, DY, INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            TYPE(c_ptr), VALUE               :: DX, DY
         END SUBROUTINE HIP_CCOPY


         SUBROUTINE HIP_DAXPY(PTR,N,DA,DX,INCX,DY,INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            REAL(c_double), VALUE            :: DA
            TYPE(c_ptr), VALUE               :: DX, DY
         END SUBROUTINE HIP_DAXPY


         SUBROUTINE HIP_ZAXPY(PTR,N,ZA,ZX,INCX,ZY,INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int),VALUE             :: N, INCX, INCY
            COMPLEX(c_double_complex), VALUE :: ZA
            TYPE(c_ptr),VALUE                :: ZX, ZY
         END SUBROUTINE HIP_ZAXPY


         SUBROUTINE HIP_DSCAL(PTR,N,DA,DX,INCX) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX
            REAL(c_double), VALUE            :: DA
            TYPE(c_ptr), VALUE               :: DX
         END SUBROUTINE HIP_DSCAL


         SUBROUTINE HIP_ZSCAL(PTR,N,ZA,ZX,INCX) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX
            COMPLEX(c_double_complex), VALUE :: ZA
            TYPE(c_ptr), VALUE               :: ZX
         END SUBROUTINE HIP_ZSCAL


         SUBROUTINE HIP_ZDSCAL(PTR,N,DA,ZX,INCX) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX
            REAL(c_double), VALUE            :: DA
            TYPE(c_ptr), VALUE               :: ZX
         END SUBROUTINE HIP_ZDSCAL


         SUBROUTINE HIP_DDOT(PTR,N,DX,INCX,DY,INCY,RESULT) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            TYPE(c_ptr), VALUE               :: DX,DY,RESULT
         END SUBROUTINE HIP_DDOT


         SUBROUTINE HIP_ZDOTC(PTR,N,ZX,INCX,ZY,INCY,RESULT) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            TYPE(c_ptr), VALUE               :: ZX,ZY,RESULT
         END SUBROUTINE HIP_ZDOTC


         SUBROUTINE HIP_CDOTC(PTR,N,ZX,INCX,ZY,INCY,RESULT) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: N, INCX, INCY
            TYPE(c_ptr), VALUE               :: ZX,ZY,RESULT
         END SUBROUTINE HIP_CDOTC


         SUBROUTINE HIP_DGEMV(PTR,MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: MODE,M,N,LDA,INCX,INCY
            REAL(c_double), VALUE            :: ALPHA,BETA
            TYPE(c_ptr), VALUE               :: A,X,Y
         END SUBROUTINE HIP_DGEMV


         SUBROUTINE HIP_ZGEMV(PTR,MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: MODE,M,N,LDA,INCX,INCY
            COMPLEX(c_double_complex), VALUE :: ALPHA,BETA
            TYPE(c_ptr), VALUE               :: A,X,Y
         END SUBROUTINE HIP_ZGEMV


         SUBROUTINE HIP_DGEMM(PTR,TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: TRANSA,TRANSB,M,N,K,LDA,LDB,LDC
            REAL(c_double), VALUE            :: ALPHA,BETA
            TYPE(c_ptr), VALUE               :: A,B,C
         END SUBROUTINE HIP_DGEMM


         SUBROUTINE HIP_ZGEMM(PTR,TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: TRANSA,TRANSB,M,N,K,LDA,LDB,LDC
            COMPLEX(c_double_complex), VALUE :: ALPHA,BETA
            TYPE(c_ptr), VALUE               :: A,B,C
         END SUBROUTINE HIP_ZGEMM


         SUBROUTINE HIP_SGEMM(PTR,TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: TRANSA,TRANSB,M,N,K,LDA,LDB,LDC
            REAL(c_float), VALUE             :: ALPHA,BETA
            TYPE(c_ptr), VALUE               :: A,B,C
         END SUBROUTINE HIP_SGEMM


         SUBROUTINE HIP_CGEMM(PTR,TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: TRANSA,TRANSB,M,N,K,LDA,LDB,LDC
            COMPLEX(c_float_complex), VALUE  :: ALPHA,BETA
            TYPE(c_ptr), VALUE               :: A,B,C
         END SUBROUTINE HIP_CGEMM


         SUBROUTINE HIP_DTRMM(PTR,SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE (c_ptr)                     :: PTR
            INTEGER(c_int), VALUE            :: SIDE,UPLO,TRANSA,DIAG,M,N,LDA,LDB
            REAL(c_double), VALUE            :: ALPHA
            TYPE(c_ptr), VALUE               :: A,B
         END SUBROUTINE HIP_DTRMM


         SUBROUTINE HIP_ZTRMM(PTR,SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE (c_ptr)                     :: PTR
            INTEGER(c_int), VALUE            :: SIDE,UPLO,TRANSA,DIAG,M,N,LDA,LDB
            COMPLEX(c_double_complex), VALUE :: ALPHA
            TYPE(c_ptr), VALUE               :: A,B
         END SUBROUTINE HIP_ZTRMM

      END INTERFACE

    CONTAINS

      INTEGER FUNCTION COLNUM(MODE,N1,N2)
      CHARACTER :: MODE
      INTEGER   :: N1,N2
      SELECT CASE (MODE)
         CASE ('N', 'n')
             COLNUM = N1
         CASE ('L', 'l')
             COLNUM = N1
         CASE DEFAULT
             COLNUM = N2
      END SELECT
      END FUNCTION COLNUM


      INTEGER FUNCTION CHAR_TO_OP(MODE)
      CHARACTER :: MODE
      SELECT CASE (MODE)
          CASE ('N', 'n')
              CHAR_TO_OP = 111 !OCBLAS_OPERATION_NONE
          CASE ('T', 't')
              CHAR_TO_OP = 112 !ROCBLAS_OEPRATION_TRANSPOSE
          CASE ('C', 'c')
              CHAR_TO_OP = 113 !ROCBLAS_OPERATION_CONJUGATE_TRANSPOSE
          CASE DEFAULT
              WRITE(*,*) "INVALID CHAR_TO_OP", MODE
              CALL EXIT(1)
      END SELECT
      END FUNCTION CHAR_TO_OP


      INTEGER FUNCTION CHAR_TO_FILL(MODE)
      CHARACTER :: MODE
      SELECT CASE (MODE)
          CASE ('U','u')
              CHAR_TO_FILL = 121
          CASE ('L','l')
              CHAR_TO_FILL = 122
          !CASE ('F','f')
          !  CHAR_TO_FILL = 123
          CASE DEFAULT
              WRITE(*,*) "INVALID CHAR_TO_FILL", MODE
              CALL EXIT(1)
      END SELECT
      END FUNCTION CHAR_TO_FILL


      INTEGER FUNCTION CHAR_TO_SIDE(MODE)
      CHARACTER :: MODE
      SELECT CASE (MODE)
          CASE ('L','l')
              CHAR_TO_SIDE = 141
          CASE ('R','r')
              CHAR_TO_SIDE = 142
          CASE DEFAULT
              WRITE(*,*) "INVALID CHAR_TO_SIDE", MODE
              CALL EXIT(1)
      END SELECT
      END FUNCTION CHAR_TO_SIDE


      INTEGER FUNCTION CHAR_TO_DIAG(MODE)
      CHARACTER :: MODE
      SELECT CASE (MODE)
          CASE ('N','n')
              CHAR_TO_DIAG = 131
          CASE ('U','u')
              CHAR_TO_DIAG = 132
          CASE DEFAULT
              WRITE(*,*) "INVALID CHAR_TO_DIAG", MODE
              CALL EXIT(1)
      END SELECT
      END FUNCTION CHAR_TO_DIAG


      INTEGER FUNCTION CHAR_TO_JOBZ(MODE)
      CHARACTER :: MODE
      SELECT CASE (MODE)
          CASE ('N','m')
              CHAR_TO_JOBZ = 213
          CASE DEFAULT
              CHAR_TO_JOBZ = 211
      END SELECT
      END FUNCTION CHAR_TO_JOBZ


      SUBROUTINE ROCBLAS_INIT()
      IF (.NOT.OFFLOAD_ALLOW) RETURN
      CALL HIP_INIT_ROCBLAS(ROCBLAS_HANDLE)
      END SUBROUTINE ROCBLAS_INIT


      SUBROUTINE ROCBLAS_FINISH()
      IF (.NOT.OFFLOAD_ALLOW) RETURN
      CALL HIP_FINISH_ROCBLAS(ROCBLAS_HANDLE)
      END SUBROUTINE ROCBLAS_FINISH


      SUBROUTINE DCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
      IMPLICIT NONE
      REAL(q), TARGET :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(DX,DY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(DX,DY)

      CALL HIP_DCOPY(ROCBLAS_HANDLE,N,c_loc(DX),INCX,c_loc(DY),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE DCOPY_OFFLOAD


      SUBROUTINE SCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
      IMPLICIT NONE
      REAL(qs), TARGET :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(DX,DY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(DX,DY)

      CALL HIP_SCOPY(ROCBLAS_HANDLE,N,c_loc(DX),INCX,c_loc(DY),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE SCOPY_OFFLOAD


      SUBROUTINE CCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
      IMPLICIT NONE
      COMPLEX(qs), TARGET :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX,ZY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(ZX,ZY)

      CALL HIP_CCOPY(ROCBLAS_HANDLE,N,c_loc(ZX),INCX,c_loc(ZY),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE CCOPY_OFFLOAD


      SUBROUTINE ZCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX,ZY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(ZX,ZY)

      CALL HIP_ZCOPY(ROCBLAS_HANDLE,N,c_loc(ZX),INCX,c_loc(ZY),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZCOPY_OFFLOAD


      SUBROUTINE DAXPY_OFFLOAD(N,DA,DX,INCX,DY,INCY)
      IMPLICIT NONE
      REAL(q), TARGET :: DX(*),DY(*)
      REAL(q) :: DA
      INTEGER :: N,INCX,INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(DX,DY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(DX,DY)

      CALL HIP_DAXPY(ROCBLAS_HANDLE,N,DA,c_loc(DX),INCX,c_loc(DY),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE DAXPY_OFFLOAD


      SUBROUTINE ZAXPY_OFFLOAD(N,ZA,ZX,INCX,ZY,INCY)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: ZX(*),ZY(*)
      COMPLEX(q) :: ZA
      INTEGER :: N,INCX,INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX,ZY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(ZX,ZY)

      CALL HIP_ZAXPY(ROCBLAS_HANDLE,N,ZA,c_loc(ZX),INCX,c_loc(ZY),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZAXPY_OFFLOAD


      SUBROUTINE DSCAL_OFFLOAD(N,DA,DX,INCX)
      IMPLICIT NONE
      REAL(q), TARGET :: DX(*)
      REAL(q) :: DA
      INTEGER :: N,INCX

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(DX)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(DX)

      CALL HIP_DSCAL(ROCBLAS_HANDLE,N,DA,c_loc(DX),INCX)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE DSCAL_OFFLOAD


      SUBROUTINE ZSCAL_OFFLOAD(N,ZA,ZX,INCX)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: ZX(*)
      COMPLEX(q) :: ZA
      INTEGER :: N,INCX

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(ZX)

      CALL HIP_ZSCAL(ROCBLAS_HANDLE,N,ZA,c_loc(ZX),INCX)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZSCAL_OFFLOAD


      SUBROUTINE ZDSCAL_OFFLOAD(N,DA,ZX,INCX)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: ZX(*)
      REAL(q) :: DA
      INTEGER :: N,INCX

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(ZX)

      CALL HIP_ZDSCAL(ROCBLAS_HANDLE,N,DA,c_loc(ZX),INCX)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZDSCAL_OFFLOAD


      FUNCTION DDOT_OFFLOAD(N,DX,INCX,DY,INCY)
      IMPLICIT NONE
      REAL(q), TARGET :: DDOT_OFFLOAD
      REAL(q), TARGET :: DX(*),DY(*)
      INTEGER :: N, INCX, INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(DX,DY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(DX,DY)

      CALL HIP_DDOT(ROCBLAS_HANDLE,N,c_loc(DX),INCX,c_loc(DY),INCY,c_loc(DDOT_OFFLOAD))

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END FUNCTION DDOT_OFFLOAD


      FUNCTION ZDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: ZDOTC_OFFLOAD
      COMPLEX(q), TARGET :: ZX(*),ZY(*)
      INTEGER :: N, INCX, INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX,ZY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(ZX,ZY)

      CALL HIP_ZDOTC(ROCBLAS_HANDLE,N,c_loc(ZX),INCX,c_loc(ZY),INCY,c_loc(ZDOTC_OFFLOAD))

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END FUNCTION ZDOTC_OFFLOAD


      FUNCTION CDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
      IMPLICIT NONE
      COMPLEX(qs), TARGET :: CDOTC_OFFLOAD
      COMPLEX(qs), TARGET :: ZX(*),ZY(*)
      INTEGER :: N, INCX, INCY

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX,ZY)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(ZX,ZY)

      CALL HIP_CDOTC(ROCBLAS_HANDLE,N,c_loc(ZX),INCX,c_loc(ZY),INCY,c_loc(CDOTC_OFFLOAD))

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END FUNCTION CDOTC_OFFLOAD


      SUBROUTINE DGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*),X(*),Y(*)
      REAL(q) :: ALPHA,BETA
      INTEGER :: M,N,LDA,INCX,INCY
      CHARACTER(1) :: MODE

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(X,Y,A)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(X,Y,A)

      CALL HIP_DGEMV(ROCBLAS_HANDLE,CHAR_TO_OP(MODE),M,N,ALPHA,c_loc(A),LDA,c_loc(X),INCX,BETA,c_loc(Y),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE DGEMV_OFFLOAD


      SUBROUTINE ZGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*),X(*),Y(*)
      COMPLEX(q) :: ALPHA,BETA
      INTEGER :: M,N,LDA,INCX,INCY
      CHARACTER(1) :: MODE

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,X,Y)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,X,Y)

      CALL HIP_ZGEMV(ROCBLAS_HANDLE,CHAR_TO_OP(MODE),M,N,ALPHA,c_loc(A),LDA,c_loc(X),INCX,BETA,c_loc(Y),INCY)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZGEMV_OFFLOAD


      SUBROUTINE DGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*),B(LDB,*),C(LDC,*)
      REAL(q) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B,C)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,B,C)

      CALL HIP_DGEMM(ROCBLAS_HANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,c_loc(A),LDA,c_loc(B),LDB,BETA,c_loc(C),LDC)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE DGEMM_OFFLOAD


      SUBROUTINE ZGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*),B(LDB,*),C(LDC,*)
      COMPLEX(q) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B,C)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,B,C)

      CALL HIP_ZGEMM(ROCBLAS_HANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,c_loc(A),LDA,c_loc(B),LDB,BETA,c_loc(C),LDC)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZGEMM_OFFLOAD


      SUBROUTINE SGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
      REAL(qs), TARGET :: A(LDA,*),B(LDB,*),C(LDC,*)
      REAL(qs) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B,C)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,B,C)

      CALL HIP_SGEMM(ROCBLAS_HANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,c_loc(A),LDA,c_loc(B),LDB,BETA,c_loc(C),LDC)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE SGEMM_OFFLOAD


      SUBROUTINE CGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
      COMPLEX(qs), TARGET :: A(LDA,*),B(LDB,*),C(LDC,*)
      COMPLEX(qs) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B,C)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,B,C)

      CALL HIP_CGEMM(ROCBLAS_HANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,c_loc(A),LDA,c_loc(B),LDB,BETA,c_loc(C),LDC)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE CGEMM_OFFLOAD


      SUBROUTINE DTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*),B(LDB,*)
      REAL(q) :: ALPHA
      INTEGER :: M,N,LDA,LDB
      CHARACTER(1) :: SIDE,UPLO,TRANSA,DIAG

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,B)

      CALL HIP_DTRMM(ROCBLAS_HANDLE,CHAR_TO_SIDE(SIDE),CHAR_TO_FILL(UPLO),CHAR_TO_OP(TRANSA),CHAR_TO_DIAG(DIAG),M,N,ALPHA,c_loc(A),LDA,c_loc(B),LDB)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE DTRMM_OFFLOAD


      SUBROUTINE ZTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*),B(LDB,*)
      COMPLEX(q) :: ALPHA
      INTEGER :: M,N,LDA,LDB
      CHARACTER(1) :: SIDE,UPLO,TRANSA,DIAG

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B)

GPACC CALL HIP_ROCBLAS_SET_STREAM(ROCBLAS_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,B)

      CALL HIP_ZTRMM(ROCBLAS_HANDLE,CHAR_TO_SIDE(SIDE),CHAR_TO_FILL(UPLO),CHAR_TO_OP(TRANSA),CHAR_TO_DIAG(DIAG),M,N,ALPHA,c_loc(A),LDA,c_loc(B),LDB)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZTRMM_OFFLOAD


      SUBROUTINE DLACPY_OFFLOAD(UPLO,M,N,A,LDA,B,LDB)
      IMPLICIT NONE
      REAL(q) :: A(LDA,*),B(LDB,*)
      INTEGER :: LDA,LDB,M,N
      CHARACTER(1) :: UPLO

      CALL vtutor%bug("DLACPY_OFFLOAD: not implemented yet",__FILE__,__LINE__)

      IF (UPLO/='F'.AND.UPLO/='f') THEN
         CALL vtutor%bug("DLACPY_OFFLOAD: only UPLO=F (full matrix) is supported",__FILE__,__LINE__)
      ENDIF

      END SUBROUTINE DLACPY_OFFLOAD

    END MODULE mrocblas


!***********************************************************************
!
!***********************************************************************

    MODULE mrocsolver
      USE prec
      USE moffload_struct_def

      USE mrocblas, ONLY : CHAR_TO_JOBZ, CHAR_TO_FILL, CHAR_TO_DIAG, CHAR_TO_OP, &
                           HIP_DCOPY_STRIDED_BATCHED, HIP_ZCOPY_STRIDED_BATCHED

GPACC USE openacc, ONLY : ACC_ASYNC_SYNC, ACC_ASYNC_NOVAL, acc_get_cuda_stream

      USE iso_c_binding

      IMPLICIT NONE

      TYPE(c_ptr) :: ROCSOLVER_HANDLE

      INTERFACE

         SUBROUTINE HIP_INIT_ROCSOLVER(PTR) BIND(C, NAME="hip_init_rocblas")
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
         END SUBROUTINE HIP_INIT_ROCSOLVER


         SUBROUTINE HIP_FINISH_ROCSOLVER(PTR) BIND(C, NAME="hip_finish_rocblas")
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
         END SUBROUTINE HIP_FINISH_ROCSOLVER


GPACC    SUBROUTINE HIP_ROCSOLVER_SET_STREAM(PTR, STREAM) BIND(C, NAME="hip_rocblas_set_stream")
GPACC       USE iso_c_binding
GPACC       TYPE(c_ptr)                      :: PTR, STREAM
GPACC    END SUBROUTINE HIP_ROCSOLVER_SET_STREAM


         SUBROUTINE HIP_DPOTRF(PTR,UPLO,N,A,LDA,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: UPLO,N,LDA
            TYPE(c_ptr), VALUE               :: A,INFO
         END SUBROUTINE HIP_DPOTRF


         SUBROUTINE HIP_ZPOTRF(PTR,UPLO,N,A,LDA,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: UPLO,N,LDA
            TYPE(c_ptr), VALUE               :: A,INFO
         END SUBROUTINE HIP_ZPOTRF


         SUBROUTINE HIP_DGETRF(PTR,M,N,A,LDA,IPIV,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: M,N,LDA
            TYPE(c_ptr), VALUE               :: A,IPIV,INFO
         END SUBROUTINE HIP_DGETRF


         SUBROUTINE HIP_ZGETRF(PTR,M,N,A,LDA,IPIV,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: M,N,LDA
            TYPE(c_ptr), VALUE               :: A,IPIV,INFO
         END SUBROUTINE HIP_ZGETRF


         SUBROUTINE HIP_DGETRS(PTR,TRANS,N,NRHS,A,LDA,IPIV,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: TRANS,N,NRHS,LDA,LDB
            TYPE(c_ptr), VALUE               :: A,B,IPIV
         END SUBROUTINE HIP_DGETRS


         SUBROUTINE HIP_ZGETRS(PTR,TRANS,N,NRHS,A,LDA,IPIV,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: TRANS,N,NRHS,LDA,LDB
            TYPE(c_ptr), VALUE               :: A,B,IPIV
         END SUBROUTINE HIP_ZGETRS


         SUBROUTINE HIP_DTRTRI(PTR,FILL,DIAG,N,A,LDA,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: FILL,DIAG,N,LDA,LDB
            TYPE(c_ptr), VALUE               :: A,B
         END SUBROUTINE HIP_DTRTRI


         SUBROUTINE HIP_ZTRTRI(PTR,FILL,DIAG,N,A,LDA,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: FILL,DIAG,N,LDA,LDB
            TYPE(c_ptr), VALUE               :: A,B
         END SUBROUTINE HIP_ZTRTRI


         SUBROUTINE HIP_DTRSM(PTR,SIDE,FILL,TRANS,DIAG,M,N,SCAL,A,LDA,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: SIDE,FILL,TRANS,DIAG,M,N,LDA,LDB
            REAL(c_double), VALUE            :: SCAL
            TYPE(c_ptr), VALUE               :: A,B
         END SUBROUTINE HIP_DTRSM


         SUBROUTINE HIP_ZTRSM(PTR,SIDE,FILL,TRANS,DIAG,M,N,SCAL,A,LDA,B,LDB) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: SIDE,FILL,TRANS,DIAG,M,N,LDA,LDB
            COMPLEX(c_double_complex), VALUE :: SCAL
            TYPE(c_ptr), VALUE               :: A,B
         END SUBROUTINE HIP_ZTRSM


         SUBROUTINE HIP_DSYGV(PTR,ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: ITYPE,JOBZ,UPLO,N,LDA,LDB
            TYPE(c_ptr), VALUE               :: A,B,W,WORK,INFO
         END SUBROUTINE HIP_DSYGV


         SUBROUTINE HIP_ZHEGV(PTR,ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: ITYPE,JOBZ,UPLO,N,LDA,LDB
            TYPE(c_ptr), VALUE               :: A,B,W,WORK,INFO
         END SUBROUTINE HIP_ZHEGV


         SUBROUTINE HIP_ZHEGVJ(PTR,ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,RES,N_SW,W,WORK,INFO_) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: ITYPE,JOBZ,UPLO,N,LDA,LDB
            TYPE(c_ptr), VALUE               :: A,B,W,WORK,INFO_,RES,N_SW
         END SUBROUTINE HIP_ZHEGVJ


         SUBROUTINE HIP_DSYEV(PTR,JOBZ,UPLO,N,A,LDA,W,WORK,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: JOBZ,UPLO,N,LDA
            TYPE(c_ptr), VALUE               :: A,W,WORK,INFO
         END SUBROUTINE HIP_DSYEV


         SUBROUTINE HIP_ZHEEV(PTR,JOBZ,UPLO,N,A,LDA,W,WORK,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: JOBZ,UPLO,N,LDA
            TYPE(c_ptr), VALUE               :: A,W,WORK,INFO
         END SUBROUTINE HIP_ZHEEV


         SUBROUTINE HIP_SSYEVX(PTR,JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,NEV,W,Z,LDZ,IFAIL,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: JOBZ,RANGE,UPLO,N,LDA,IL,IU,LDZ
            REAL(c_float), VALUE             :: VL,VU,ABSTOL
            TYPE(c_ptr), VALUE               :: A,NEV,W,Z,IFAIL,INFO
         END SUBROUTINE HIP_SSYEVX


         SUBROUTINE HIP_DSYEVX(PTR,JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,NEV,W,Z,LDZ,IFAIL,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: JOBZ,RANGE,UPLO,N,LDA,IL,IU,LDZ
            REAL(c_double), VALUE            :: VL,VU,ABSTOL
            TYPE(c_ptr), VALUE               :: A,NEV,W,Z,IFAIL,INFO
         END SUBROUTINE HIP_DSYEVX


         SUBROUTINE HIP_ZHEEVX(PTR,JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,NEV,W,Z,LDZ,IFAIL,INFO) BIND(C)
            USE iso_c_binding
            TYPE(c_ptr)                      :: PTR
            INTEGER(c_int), VALUE            :: JOBZ,RANGE,UPLO,N,LDA,IL,IU,LDZ
            REAL(c_double), VALUE            :: VL,VU,ABSTOL
            TYPE(c_ptr), VALUE               :: A,NEV,W,Z,IFAIL,INFO
         END SUBROUTINE HIP_ZHEEVX

      END INTERFACE

    CONTAINS

      INTEGER FUNCTION CHAR_TO_ERANGE(MODE)
      CHARACTER :: MODE
      SELECT CASE (MODE)
          CASE ('A','a')
              CHAR_TO_ERANGE = 231
          CASE ('V','v')
              CHAR_TO_ERANGE = 232
          CASE ('I','i')
              CHAR_TO_ERANGE = 233
          CASE DEFAULT
              WRITE(*,*) "INVALID CHAR_TO_ERANGE", MODE
              CALL EXIT(1)
      END SELECT
      END FUNCTION CHAR_TO_ERANGE


      SUBROUTINE ROCSOLVER_INIT()
      IF (.NOT.OFFLOAD_ALLOW) RETURN
      CALL HIP_INIT_ROCSOLVER(ROCSOLVER_HANDLE)
      END SUBROUTINE ROCSOLVER_INIT


      SUBROUTINE ROCSOLVER_FINISH()
      IF (.NOT.OFFLOAD_ALLOW) RETURN
      CALL HIP_FINISH_ROCSOLVER(ROCSOLVER_HANDLE)
      END SUBROUTINE ROCSOLVER_FINISH


      SUBROUTINE DPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*)
      INTEGER, TARGET :: INFO
      INTEGER :: N,LDA
      CHARACTER(1) :: UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(TO:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,INFO)

      CALL HIP_DPOTRF(ROCSOLVER_HANDLE,CHAR_TO_FILL(UPLO),N,c_loc(A),LDA,c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE DPOTRF_OFFLOAD


      SUBROUTINE ZPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*)
      INTEGER, TARGET :: INFO
      INTEGER :: N,LDA
      CHARACTER(1) :: UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(TO:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,INFO)

      CALL HIP_ZPOTRF(ROCSOLVER_HANDLE,CHAR_TO_FILL(UPLO),N,c_loc(A),LDA,c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE ZPOTRF_OFFLOAD


      SUBROUTINE DGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*)
      INTEGER, TARGET :: IPIV(*),INFO
      INTEGER :: M,N,LDA

GPOMP !$OMP TARGET ENTER DATA MAP(TO:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,IPIV,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,IPIV,INFO)

      CALL HIP_DGETRF(ROCSOLVER_HANDLE,M,N,c_loc(A),LDA,c_loc(IPIV),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE DGETRF_OFFLOAD


      SUBROUTINE ZGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*)
      INTEGER, TARGET :: IPIV(*),INFO
      INTEGER :: M,N,LDA

GPOMP !$OMP TARGET ENTER DATA MAP(TO:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,IPIV,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,IPIV,INFO)

      CALL HIP_ZGETRF(ROCSOLVER_HANDLE,M,N,c_loc(A),LDA,c_loc(IPIV),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE ZGETRF_OFFLOAD


      SUBROUTINE DGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*),B(LDB,*)
      INTEGER, TARGET :: IPIV(*)
      INTEGER :: N,NRHS,LDA,LDB,INFO
      CHARACTER(1) :: TRANS

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,IPIV,B)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,IPIV,B)

      CALL HIP_DGETRS(ROCSOLVER_HANDLE,CHAR_TO_OP(TRANS),N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE DGETRS_OFFLOAD


      SUBROUTINE ZGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*),B(LDB,*)
      INTEGER, TARGET :: IPIV(*)
      INTEGER :: N,NRHS,LDA,LDB,INFO
      CHARACTER(1) :: TRANS

GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B,IPIV)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC HOST_DATA USE_DEVICE(A,B,IPIV)

      CALL HIP_ZGETRS(ROCSOLVER_HANDLE,CHAR_TO_OP(TRANS),N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)

GPACC !$ACC END HOST_DATA

GPOMP !$OMP END TARGET DATA
      END SUBROUTINE ZGETRS_OFFLOAD


      SUBROUTINE DTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO,DIAG
      ! local variables
      REAL(q), ALLOCATABLE, TARGET :: WORK(:)

      ALLOCATE(WORK(N*N))

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WORK)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,WORK)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(WORK) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,WORK)

      CALL HIP_DTRTRI(ROCSOLVER_HANDLE,CHAR_TO_FILL(UPLO),CHAR_TO_DIAG(DIAG),N,c_loc(A),LDA,c_loc(WORK),N)
      CALL HIP_DCOPY_STRIDED_BATCHED(ROCSOLVER_HANDLE,N,c_loc(WORK),1,N,c_loc(A),1,LDA,N)

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC WAIT(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:WORK)

      DEALLOCATE(WORK)
      END SUBROUTINE DTRTRI_OFFLOAD


      SUBROUTINE ZTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO,DIAG
      ! local variables
      COMPLEX(q), ALLOCATABLE, TARGET :: WORK(:)

      ALLOCATE(WORK(N*N))

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WORK)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,WORK)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(WORK) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,WORK)

      CALL HIP_ZTRTRI(ROCSOLVER_HANDLE,CHAR_TO_FILL(UPLO),CHAR_TO_DIAG(DIAG),N,c_loc(A),LDA,c_loc(WORK),N)
      CALL HIP_ZCOPY_STRIDED_BATCHED(ROCSOLVER_HANDLE,N,c_loc(WORK),1,N,c_loc(A),1,LDA,N)

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC WAIT(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:WORK)

      DEALLOCATE(WORK)
      END SUBROUTINE ZTRTRI_OFFLOAD


      SUBROUTINE DSYGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*),B(LDB,*),W(*),WORK(*)
      INTEGER, TARGET :: INFO
      INTEGER :: ITYPE,N,LDA,LDB,LWORK
      CHARACTER(1) :: JOBZ,UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B,W,WORK,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,B,W,WORK,INFO)

      CALL HIP_DSYGV(ROCSOLVER_HANDLE,ITYPE,CHAR_TO_JOBZ(JOBZ),CHAR_TO_FILL(UPLO), &
                     N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),c_loc(WORK),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE DSYGV_OFFLOAD


      SUBROUTINE ZHEGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*),B(LDB,*)
      REAL(q), TARGET :: W(*),RWORK(*)
      INTEGER, TARGET :: INFO
      COMPLEX(q) :: WORK(*)
      INTEGER :: ITYPE,N,LDA,LDB,LWORK
      CHARACTER(1) :: JOBZ,UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,B,W,RWORK,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,B,W,RWORK,INFO)

      CALL HIP_ZHEGV(ROCSOLVER_HANDLE,ITYPE,CHAR_TO_JOBZ(JOBZ),CHAR_TO_FILL(UPLO), &
                     N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),c_loc(RWORK),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE ZHEGV_OFFLOAD


      SUBROUTINE DSYEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*),W(*),WORK(*)
      INTEGER, TARGET :: INFO
      INTEGER :: N,LDA,LWORK
      CHARACTER(1) :: JOBZ,UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,W,WORK,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,W,WORK,INFO)

       CALL HIP_DSYEV(ROCSOLVER_HANDLE,CHAR_TO_JOBZ(JOBZ),CHAR_TO_FILL(UPLO), &
                      N,c_loc(A),LDA,c_loc(W),c_loc(WORK),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE DSYEV_OFFLOAD


      SUBROUTINE ZHEEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*)
      REAL(q), TARGET :: W(*),RWORK(*)
      INTEGER, TARGET :: INFO
      COMPLEX(q) :: WORK(*)
      INTEGER :: N,LDA,LWORK
      CHARACTER(1) :: JOBZ,UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,W,RWORK,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,W,RWORK,INFO)

      CALL HIP_ZHEEV(ROCSOLVER_HANDLE,CHAR_TO_JOBZ(JOBZ),CHAR_TO_FILL(UPLO), &
                     N,c_loc(A),LDA,c_loc(W),c_loc(RWORK),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE ZHEEV_OFFLOAD


      SUBROUTINE SSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
      IMPLICIT NONE
      REAL(qs), TARGET :: A(LDA,*),W(*),WORK(*),Z(LDZ,*)
      INTEGER, TARGET :: IFAIL(*),INFO
      REAL(qs) :: ABSTOL,VL,VU
      INTEGER :: IWORK(*),N,LDA,IL,IU,M,LDZ,LWORK
      CHARACTER(1) :: JOBZ,RANGE,UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,W,WORK,Z,IFAIL,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,W,WORK,Z,IFAIL,INFO)

      CALL HIP_SSYEVX(ROCSOLVER_HANDLE,CHAR_TO_JOBZ(JOBZ),CHAR_TO_ERANGE(RANGE),CHAR_TO_FILL(UPLO), &
                     N,c_loc(A),LDA,VL,VU,IL,IU,ABSTOL,c_loc(W),c_loc(WORK),c_loc(Z),LDZ,c_loc(IFAIL),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE SSYEVX_OFFLOAD


      SUBROUTINE DSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
      IMPLICIT NONE
      REAL(q), TARGET :: A(LDA,*),W(*),WORK(*),Z(LDZ,*)
      INTEGER, TARGET :: IFAIL(*),INFO
      REAL(q) :: VL,VU,ABSTOL
      INTEGER :: IWORK(*),N,LDA,IL,IU,M,LDZ,LWORK
      CHARACTER(1) :: JOBZ,RANGE,UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,W,WORK,Z,IFAIL,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,W,WORK,Z,IFAIL,INFO)

      CALL HIP_DSYEVX(ROCSOLVER_HANDLE,CHAR_TO_JOBZ(JOBZ),CHAR_TO_ERANGE(RANGE),CHAR_TO_FILL(UPLO), &
                     N,c_loc(A),LDA,VL,VU,IL,IU,ABSTOL,c_loc(W),c_loc(WORK),c_loc(Z),LDZ,c_loc(IFAIL),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE DSYEVX_OFFLOAD


      SUBROUTINE ZHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
      IMPLICIT NONE
      COMPLEX(q), TARGET :: A(LDA,*),Z(LDZ,*),WORK(*)
      REAL(q), TARGET :: W(*)
      INTEGER, TARGET :: IFAIL(*),INFO
      REAL(q) :: VL,VU,ABSTOL,RWORK(*)
      INTEGER :: IWORK(*),N,LDA,IL,IU,M,LDZ,LWORK
      CHARACTER(1) :: JOBZ,RANGE,UPLO

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:INFO)
GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(A,W,WORK,Z,IFAIL,INFO)

GPACC CALL HIP_ROCSOLVER_SET_STREAM(ROCSOLVER_HANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
GPACC !$ACC ENTER DATA CREATE(INFO) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(A,W,WORK,Z,IFAIL,INFO)

      CALL HIP_ZHEEVX(ROCSOLVER_HANDLE,CHAR_TO_JOBZ(JOBZ),CHAR_TO_ERANGE(RANGE),CHAR_TO_FILL(UPLO), &
                     N,c_loc(A),LDA,VL,VU,IL,IU,ABSTOL,c_loc(W),c_loc(WORK),c_loc(Z),LDZ,c_loc(IFAIL),c_loc(INFO))

GPACC !$ACC END HOST_DATA
GPACC !$ACC EXIT DATA COPYOUT(INFO) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP END TARGET DATA
GPOMP !$OMP TARGET UPDATE FROM(INFO)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:INFO)
      END SUBROUTINE ZHEEVX_OFFLOAD

    END MODULE mrocsolver
#endif // OFFLOADING
#endif // CRAYHIP
