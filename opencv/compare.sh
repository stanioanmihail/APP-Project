#!/bin/bash

#set -x
THOLD=40

rm -f out.bmp pthreads_out.bmp mpi_out.bmp
if ! test -d refs; then
	mkdir refs
fi

make clean &> /dev/null && make &> /dev/null

for image in $(ls ../images); do
	./DisplayImage "../images/$image" "./refs/opencv_$image" $THOLD &> /dev/null
done

echo "***************Serial Test*********************"
cd ../serial
make clean &> /dev/null && make &> /dev/null
cd ../opencv

for image in $(ls ../images); do
	echo "$image"
	../serial/canny "../images/$image" &> /dev/null
	./CompareImage "out.bmp" "./refs/opencv_$image"
done

cd ../serial
make clean &> /dev/null
cd ../opencv

echo "***************Pthreads Test*********************"
cd ../PThreads
make clean &> /dev/null && make &> /dev/null
cd ../opencv

for image in $(ls ../images); do
	echo "$image"
	../PThreads/pthreads_canny "../images/$image" 4 &> /dev/null
	./CompareImage "pthreads_out.bmp" "./refs/opencv_$image"
done

cd ../PThreads
make clean &> /dev/null
cd ../opencv

echo "***************MPI Test*********************"
cd ../MPI
make clean &> /dev/null && make &> /dev/null
cd ../opencv

for image in $(ls ../images); do
	echo "$image"
	mpirun -n 4 ../MPI/mpi_canny "../images/$image" 4 &> /dev/null
	./CompareImage "mpi_out.bmp" "./refs/opencv_$image"
done

cd ../MPI
make clean &> /dev/null
cd ../opencv

echo "***************OMP Test*********************"
cd ../omp
make clean &> /dev/null && make &> /dev/null
cd ../opencv

for image in $(ls ../images); do
	echo "$image"
	../omp/omp_canny "../images/$image" &> /dev/null
	./CompareImage "out.bmp" "./refs/opencv_$image"
done

cd ../omp
make clean &> /dev/null
cd ../opencv


rm -rf refs
rm -rf *.bmp
exit 0
