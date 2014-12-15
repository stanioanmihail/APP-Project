GCC = gcc
CCFLAGS = -Wall -lm
SOURCE = omp_canny.c
TARGET = omp_canny

build: omp_canny.o
	$(GCC) --std=c99 -g $^ -o $(TARGET) $(CCFLAGS) -fopenmp

omp_canny.o: omp_canny.c
	$(GCC) --std=c99 -g -c $(SOURCE) $(CCFLAGS) -fopenmp

clean:
	rm -rf *~ *.o out.bmp $(TARGET)
