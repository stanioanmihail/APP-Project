#!/bin/bash
image=../images/alps.bmp
module load libraries/openmpi-1.6-gcc-4.6.3
module load compilers/gnu-4.6.3 

mpirun -n 10  ./mpi_canny "../images/$image"  

