GCC = gcc
CCFLAGS = -Wall -lm
SOURCE = canny.c
TARGET = canny

build: canny.o file_io.o
	$(GCC) --std=c99 -g $^ -o $(TARGET) $(CCFLAGS)

canny.o: canny.c
	$(GCC) --std=c99 -g -c $(SOURCE) $(CCFLAGS)

file_io.o: file_io.c
	$(GCC) --std=c99 -g -c $^ $(CCFLAGS)

clean:
	rm -rf *~ *.o out.bmp $(TARGET)
