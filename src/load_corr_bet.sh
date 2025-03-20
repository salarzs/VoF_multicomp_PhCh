export CXX=CC
export CC=cc
export FC=mpif90
#
module restore system   # Restore loaded modules to the default
#module load intel/2019b
#module load OpenMPI/3.1.4-GCC-8.3.0 # doesn't work
module load OpenMPI/4.0.3-GCC-9.3.0 # works but above 4096 not good
#module load OpenMPI/4.0.5-GCC-10.2.0 # seems to work better on scalability
#
#make veryclean
make clean
#make libraries
make

