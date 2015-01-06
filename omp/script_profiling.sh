#!/bin/bash
module load libraries/openmpi-1.6-gcc-4.6.3
module load compilers/gnu-4.6.3 
module load compilers/solarisstudio-12.3
export OMP_NUM_THREADS=6
/opt/tools/compilers/solarisstudio12.3/lib/analyzer/lib/../../../bin/collect -o /export/home/acs/stud/i/ioan.stan/Public/APP/APP-Project/omp/prof/omp-test.1.er -p on -S on -A on /export/home/acs/stud/i/ioan.stan/Public/APP/APP-Project/omp/omp_canny /export/home/acs/stud/i/ioan.stan/Public/APP/APP-Project/plus_size_img_gen/alps_x16.bmp 

