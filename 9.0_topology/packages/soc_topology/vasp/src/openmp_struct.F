#include "symbol.inc"
     MODULE mopenmp_struct_def
      ! overall threading strategy

      !> @ref openmp :
      !> number of threads to be used in \"general\":
      !> by default this is set to the maximum number of available threads,
      !> mopenmp::omp_nthreads=omp_get_max_threads()
      INTEGER :: omp_nthreads = 1

      ! specific routines

      !> @ref openmp :
      !> number of threads used in ::m_alltoall_d_omp
      !> (experimental and unused for now).
      INTEGER :: omp_nthreads_alltoall = 1

      !> @ref openmp :
      !> number of threads used in the domain decomposition of the
      !> real space projection operators.
      !> Default: mopenmp::omp_nthreads_nonlr_rspace=mopenmp::omp_nthreads
      INTEGER :: omp_nthreads_nonlr_rspace = 1

      !> @ref openmp :
      !> number of threads to be used in connection with OpenACC:
      !> in most instances this involves a distribution over orbitals
      INTEGER :: omp_nthreads_acc = 1

      !> @ref openmp :
      !> if set to true the domain of the real space projection operators
      !> is decomposed over (x,y)-planes, in a round-robin fashion, if set
      !> to false the projectors are distributed over (z,y)-columns.
      LOGICAL :: omp_nonlr_planewise = .TRUE.

      !> @ref openmp :
      !> if set to true ::fftmakeplan or ::fftmakeplan_mpi will call
      !> dfftw_init_threads and set mopenmp::omp_dfftw_init_threads=.FALSE.
      LOGICAL :: omp_dfftw_init_threads = .TRUE.

#ifdef OMP_OFFLOAD
      !> The device number to which this MPI-rank connects itself
      !> by means of OMP_SET_DEFAULT_DEVICE in INIT_OFFLOAD.
      INTEGER :: OMP_MY_DEVICE_NUM
#endif
     END MODULE mopenmp_struct_def
