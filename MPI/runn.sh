#!/bin/bash
module load libraries/openmpi-1.6-gcc-4.6.3
module load compilers/gnu-4.6.3

#files=(../plus_size_img_gen/alps_x16 ../plus_size_img_gen/alps_x18 ../plus_size_img_gen/alps_x20)
#input=.bmp
image_location=../plus_size_img_gen/
#image_location=../images/
output=.mpirez

#for i in ${files[*]}
#do
#	fisier_out=$i$output
#	rm -rf $fisier_out
#done

rm -rf *.$output
#make clean && make

#for i in ${files[*]}
#do
#	fisier_in=$i$input
#	fisier_out=$i$output
#	temp_file=testing
#	touch $temp_file
#	echo "****$fisier_out****"
#	for j in 1 2 4 8
#	do
#		#(time mpirun -n $j ./mpi_main $fisier_in 3 0.001) >> $fisier_out 2>&1
#		(time mpirun -n $j ./mpi_canny $fisier_in) &> $temp_file
#		echo -ne "$j\t" >> $fisier_out
#		cat $temp_file | head -2 | tail -1| cut -f 2 | tr -d "ms" | tr -d "0" >> $fisier_out
#	done 
#done
for bmp in $(ls "$image_location" | grep .bmp); do
    echo $bmp
    file_in="$image_location$bmp"
    file_out="$bmp$output"
    file_temp=testing
    for cores in $(seq 2 2 10); do
        echo "File $file_in processed by $cores cores"
        (time mpirun -n $cores ./mpi_canny $file_in $cores) &> $file_temp
        echo -ne "$cores\t" >> $file_out
		cat $file_temp | head -2 | tail -1| cut -f 2 | tr -d "ms" | tr -d "0" >> $file_out
        cat $file_out
    done
done
