#!/bin/bash
module load libraries/openmpi-1.6-gcc-4.6.3
module load compilers/gnu-4.6.3 
module load compilers/solarisstudio-12.3

/opt/tools/compilers/solarisstudio12.3/lib/analyzer/lib/../../../bin/collect -o test.1.er -d /export/home/acs/stud/i/ioan.stan/Public/APP/APP-Project/PThreads/profiling -p on -S on -A on /export/home/acs/stud/i/ioan.stan/Public/APP/APP-Project/PThreads/pthreads_canny /export/home/acs/stud/i/ioan.stan/Public/APP/APP-Project/plus_size_img_gen/alps_x16.bmp 6 
