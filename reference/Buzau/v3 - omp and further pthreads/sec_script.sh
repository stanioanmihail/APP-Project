#!/bin/bash

module load compilers/solarisstudio-12.3

if [ $# == 5 ]
then
	export OMP_NUM_THREADS=4
	collect -o five_test.er ./$1 $2 $3 $4 $5
elif [ $# == 6 ]
then
	export OMP_NUM_THREADS=4
	collect -o six_test.er ./$1 $2 $3 $4 $5 $6
else
	echo "Usage: ./sec_script.sh exec_name pic_name sigma tlow thigh [dirname]"
fi
