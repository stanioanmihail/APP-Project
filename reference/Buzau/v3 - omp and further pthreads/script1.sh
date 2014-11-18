#!/bin/bash

if [ $# == 5 ]
then
	export OMP_NUM_THREADS=4
	time ./$1 $2 $3 $4 $5
elif [ $# == 6 ]
then
	export OMP_NUM_THREADS=4
	time ./$1 $2 $3 $4 $5 $6
else
	echo "Usage: ./script1.sh exec_name pic_name sigma tlow thigh [dirname]"
fi
