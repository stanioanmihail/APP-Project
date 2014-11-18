#!/bin/bash

if [ $# != 1 ]
then
	echo "Usage : ./qscript1.sh option_number"
elif [ $1 == 1 ]
then
	qsub -q ibm-opteron.q -cwd ./script1.sh omp_canny_edge earth.pgm 1.80 0.40 0.80
elif [ $1 == 2 ]
then
	qsub -q ibm-opteron.q -cwd ./script1.sh debug_omp_canny_edge earth.pgm 1.80 0.40 0.80
elif [ $1 == 3 ]
then
	qsub -q ibm-opteron.q -cwd ./sec_script.sh debug_omp_canny_edge earth.pgm 1.80 0.40 0.80
else
	echo "The f**k is this?"
fi
