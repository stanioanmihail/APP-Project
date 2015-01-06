#!/bin/bash
module load libraries/openmpi-1.6-gcc-4.6.3
module load compilers/gnu-4.6.3

files=(color50000_3 color100000_3 color150000_3 color500000_3 color1000000_3 color3000000_3 color5000000_3)
input=.txt
output=.rez

for i in ${files[*]}
do
	fisier_out=$i$output
	rm -rf $fisier_out
done

for i in ${files[*]}
do
	fisier_in=$i$input
	fisier_out=$i$output
	temp_file=testing
	touch $temp_file
	echo "****$fisier_out****"
	for j in 1 2 4 8
	do
		#(time mpirun -n $j ./mpi_main $fisier_in 3 0.001) >> $fisier_out 2>&1
		(time mpirun -n $j ./mpi_main $fisier_in 3 0.001) &> $temp_file
		echo -ne "$j\t" >> $fisier_out
		cat $temp_file | head -2 | tail -1| cut -f 2 | tr -d "ms" | tr -d "0" >> $fisier_out
	done 
done
