#!/bin/bash
#
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=24:mpiprocs=24:ompthreads=1:ngpus=0:mem=12gb
#PBS -N sn2_string

NODES=24

module load lang/intel-parallel-studio-xe/2019.u3
export I_MPI_COMPATIBILITY=3

cd $PBS_O_WORKDIR
source /home/vn18144/soft/amber18/amber.sh 

mkdir results
./in.sh $NODES
mpiexec sander.MPI -ng $NODES -groupfile string.groupfile

