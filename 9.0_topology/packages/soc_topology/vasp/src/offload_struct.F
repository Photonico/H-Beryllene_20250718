#include "symbol.inc"
#ifdef OFFLOADING
    MODULE moffload_struct_def

GPACC USE mopenacc_struct_def
GPOMP USE mopenmp_struct_def

      !> Flags whether device-sided execution is allowed in principle
      LOGICAL, PUBLIC :: OFFLOAD_ALLOW=.FALSE.
      !> Flags whether device-sided execution is momentarily activated
      LOGICAL, PUBLIC :: OFFLOAD_ON=.FALSE.

      !> Total number of devices detected by the OpenACC runtime
      INTEGER, PUBLIC :: OFFLOAD_NUM_DEVICES
      !> Maximum number of MPI-ranks per device
      INTEGER, PUBLIC :: OFFLOAD_MAX_RANKS_PER_DEVICE

    END MODULE moffload_struct_def
#endif // OFFLOADING
