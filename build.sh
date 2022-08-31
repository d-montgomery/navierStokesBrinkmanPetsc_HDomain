#!/bin/bash
# Load compilers for mio then make executable file
module load PrgEnv/devtoolset-9
module load MPI/openmpi/4.1.1/gcc
module load Lib/petsc/3.13-gnu-ompi-dbg

make CLEAN
make run_main
rm *.mod *.o
