#!/bin/csh
#PBS -N ph0_b_p3m1
#PBS -q cmt
#PBS -j oe
#PBS -l select=1:ncpus=24:mpiprocs=24:mem=100GB
#PBS -l walltime=24:00:00
#PBS -J 1-2
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

set i = `printf "%03d" $PBS_ARRAY_INDEX`
cd "$PBS_O_WORKDIR/disp-$i"
mpirun -np 24 $VASP > vasp_physics_cluster.out
