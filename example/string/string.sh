#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --job-name=string_test

# Number of task does not necessarily have to be 28, any number around 30 is OK

# If needed, load the proper modules
# module load languages/intel/2017.01

export I_MPI_COMPATIBILITY=3 # Needed for new versions of MPI

# Make sure $AMBERHOME is set or source amber.sh directly
source $AMBERHOME/amber.sh 

STRING_TEST_DIR= # Path to the directory with example files
cd $STRING_TEST_DIR/string 
rm -rf results
mkdir results
./in.sh 28
cp STOP_STRING results/ # Can be done while job is running
srun sander.MPI -ng 28 -groupfile string.groupfile

