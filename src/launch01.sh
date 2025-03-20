#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH --account=NN9561K

# The name of the script is myjob
#SBATCH -J a7p5p.sh

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 03-23:59:00

# Number of nodes
#SBATCH -N 8
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=128
# Number of MPI processes.
#SBATCH -n 1024
####SBATCH -n 2048

#SBATCH -e sb_error.e
#SBATCH -o sb_output.o

# Run the executable named myexe 
# and write the output into my_output_file
#
#module unload PrgEnv-intel
#module load PrgEnv-gnu
#
#module unload PrgEnv-gnu
#module swap PrgEnv-cray PrgEnv-intel
#export CRAYTYPE_LINK_TYPE=dynamic
#export CRAY_ROOFTFS=DSL
#
ulimit -s unlimited
#
module restore system   # Restore loaded modules to the default
#module load intel/2019b
#module load OpenMPI/3.1.4-GCC-8.3.0
module load OpenMPI/4.0.3-GCC-9.3.0
#module load OpenMPI/4.0.5-GCC-10.2.0 # seems to work better on scalability
#
#srun -n 1 make veryclean > makeveryclean.log 2>&1
#srun -n 1 make libraries > makelib.log 2>&1m -n 1 make clean > makeclean.log 2>&1
mpirun -n 1 make clean > makeclean.log 2>&1
mpirun -n 1 make > make.log 2>&1
mpirun -n 1024 ./cans > out.log 2>&1
