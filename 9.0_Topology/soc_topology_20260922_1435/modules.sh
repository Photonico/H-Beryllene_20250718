set +u
source /etc/profile.d/modules.sh
module load pbspro oneapi-2024.2/compiler/2024.2.1 oneapi-2024.2/mkl/latest oneapi-2024.2/mpi/latest hdf/5/1.14.1-2_intel2021
set -u
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
