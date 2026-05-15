#!/bin/csh
#PBS -N Bands_2H_Beryllene
#PBS -q cmt
#PBS -j oe
#PBS -l select=1:ncpus=84:mem=250GB
#PBS -l walltime=24:00:00
#PBS -m a
#PBS -M oliver.conquest@sydney.edu.au

cd "$PBS_O_WORKDIR"

module purge
module load pbspro
module load oneapi-2024.2/compiler-rt32/latest
module load oneapi-2024.2/compiler/2024.2.1 
module load oneapi-2024.2/mkl/latest
module load oneapi-2024.2/mpi/latest    
module load hdf/5/1.14.1-2_intel2021

mpirun -n 84 /cmt2/ocon2505/VASP/vasp.6.5.1/bin/vasp_std > vasp.out