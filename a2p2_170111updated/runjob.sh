#!/bin/sh
#PBS -N JOB
#PBS -q full
#PBS -l nodes=1:ppn=10:full

VASP_EXE='/usr/local/vasp/vasp-5.3.5-vtst_i13_mkl'
NPROCS=`cat $PBS_NODEFILE|wc -l`
cd $PBS_O_WORKDIR

mpirun -machinefile $PBS_NODEFILE -np $NPROCS $VASP_EXE > log 2>&1
