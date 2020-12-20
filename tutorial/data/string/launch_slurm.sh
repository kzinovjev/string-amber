#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --job-name=sn2_string

NODES=24

module load languages/intel/2017.01
export I_MPI_COMPATIBILITY=3
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

cd $SLURM_SUBMIT_DIR
source /mnt/storage/home/vn18144/amber18/amber.sh 

mkdir results
./in.sh $NODES
srun sander.MPI -ng $NODES -groupfile string.groupfile

