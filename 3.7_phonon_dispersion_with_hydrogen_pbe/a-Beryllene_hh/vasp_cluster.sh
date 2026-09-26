#!/bin/csh
#PBS -N ph_Beryllene_hh
#PBS -q cmt
#PBS -j oe
#PBS -l select=1:ncpus=24:mpiprocs=24:mem=100GB
#PBS -l walltime=4:00:00
#PBS -J 1-6
#PBS -m a
#PBS -M luke.niu@sydney.edu.au

set i = `printf "%03d" $PBS_ARRAY_INDEX`
cd "$PBS_O_WORKDIR/disp-$i"

module load pbspro
module load oneapi-2024.2/compiler-rt32/latest
module load oneapi-2024.2/mkl/latest
module load oneapi-2024.2/mpi/latest    
module load hdf/5/1.14.1-2_intel2021

set VASP=/cmt2/ocon2505/VASP/vasp.6.5.0/bin/vasp_std
set BIN=/cmt2/ocon2505/VASP/vasp.6.5.0/bin/vasp_std

mpirun -np 24 $VASP > vasp_physics_cluster.out