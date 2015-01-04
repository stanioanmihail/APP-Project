#!/bin/bash
module load libraries/openmpi-1.6-gcc-4.6.3
module load compilers/gnu-4.6.3 
module load compilers/solarisstudio-12.3

mpirun -np 10 /opt/tools/compilers/solarisstudio12.3/lib/analyzer/lib/../../../bin/collect -o mpi-test.1.er -p on -S on -A on /export/home/acs/stud/i/ioan.stan/APP/APP-Project/MPI/mpi_canny /export/home/acs/stud/i/ioan.stan/APP/APP-Project/images/alps.bmp 

