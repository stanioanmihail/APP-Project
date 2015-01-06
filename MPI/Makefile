GCC = mpicc
CCFLAGS = -Wall -lm
SOURCE = mpi_canny.c
TARGET = mpi_canny

build: mpi_canny.o file_io.o
	$(GCC) --std=c99 -g $^ -o $(TARGET) $(CCFLAGS)

mpi_canny.o: mpi_canny.c
	$(GCC) --std=c99 -g -c $(SOURCE) $(CCFLAGS)

file_io.o: file_io.c
	$(GCC) --std=c99 -g -c $^ $(CCFLAGS)

run:
	mpirun -n 4 ./$(TARGET) $1
fep:
	qsub -q ibm-nehalem.q -pe openmpi*1 10 -cwd ./runn.sh	
profiling:
	qsub -q ibm-nehalem.q -pe openmpi*1 10 -cwd ./script_profiling.sh	
clean:
	rm -rf *~ *.o *.bmp run.sh.* $(TARGET)
clean-extra:
	rm -rf *.mpirez testing runn.sh.*
