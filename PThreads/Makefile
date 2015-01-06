GCC = gcc
CCFLAGS = -Wall -lm -lpthread
SOURCE = pthreads_canny.c
TARGET = pthreads_canny

build: pthreads_canny.o file_io.o
	gcc --std=c99 -g $^ -o $(TARGET) $(CCFLAGS)

pthreads_canny.o: pthreads_canny.c
	gcc --std=c99 -g -c $(SOURCE) $(CCFLAGS)

file_io.o: file_io.c
	gcc --std=c99 -g -c $^ $(CCFLAGS)

profiling:
	qsub -q ibm-nehalem -cwd ./profiling.sh

run:
	./$(TARGET) $1 $2
clean:
	rm -rf *~ *.o *.bmp $(TARGET)
