#!/bin/csh
#PBS -N geo_b_third_th2
#PBS -q cmt
#PBS -j oe
#PBS -l select=1:ncpus=24:mpiprocs=24:mem=400GB
#PBS -l walltime=48:00:00
#PBS -m a
#PBS -M luke.niu@sydney.edu.au

cd "$PBS_O_WORKDIR"

module load pbspro
module load oneapi-2024.2/compiler-rt32/latest
module load oneapi-2024.2/mkl/latest
module load oneapi-2024.2/mpi/latest    
module load hdf/5/1.14.1-2_intel2021

set VASP=/cmt2/ocon2505/VASP/vasp.6.5.0/bin/vasp_std
set BIN=/cmt2/ocon2505/VASP/vasp.6.5.0/bin/vasp_std

mpirun -np 24 $VASP > vasp_cluster.out