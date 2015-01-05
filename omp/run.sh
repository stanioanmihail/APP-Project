#!/bin/bash
module load libraries/openmpi-1.6-gcc-4.6.3
module load compilers/gnu-4.6.3

files=(../plus_size_img_gen/alps_x16 ../plus_size_img_gen/alps_x18 ../plus_size_img_gen/alps_x20)
input=.bmp
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
		export OMP_NUM_THREADS=$j
		(time ./omp_canny $fisier_in) &> $temp_file
		echo -ne "$j\t" >> $fisier_out
		cat $temp_file | head -2 | tail -1| cut -f 2 >> $fisier_out
	done 
done
