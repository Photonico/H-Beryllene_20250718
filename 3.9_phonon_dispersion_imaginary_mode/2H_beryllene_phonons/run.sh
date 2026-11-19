#!/bin/bash --login

#SBATCH --job-name=2H_Beryllene_SC-12x12
#SBATCH --account=pawsey0328
#SBATCH --partition=work
#SBATCH --ntasks=768
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=230G      #Indicate the amount of memory per node when asking for share resources
#SBATCH --time=23:59:00
#SBATCH --output=output.txt

#Modules
module load hdf5/1.14.3-api-v112
module load netlib-scalapack/2.2.0
module load fftw/3.3.10

# OpenMP settings
export OMP_NUM_THREADS=1   #To define the number of threads, in this case will be 32
#export OMP_PLACES=cores     #To bind threads to cores
#export OMP_PROC_BIND=close  #To bind (fix) threads (allocating them as close as possible). This option works together with the "places" indicated above, then: allocates threads in closest cores.
#Settings for MPI jobs
export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_VERBOSE=1

srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c $OMP_NUM_THREADS -m block:block:block /home/oconquest/vasp6.6.0/vasp_std
E=`tail -1 OSZICAR`
echo $E >> ENERGY.dat
