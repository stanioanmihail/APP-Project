GCC = gcc
CCFLAGS = -Wall -lm
SOURCE = canny.c
TARGET = canny

build: canny.o
	$(GCC) --std=c99 -g $^ -o $(TARGET) $(CCFLAGS)

canny.o: canny.c
	$(GCC) --std=c99 -g -c $(SOURCE) $(CCFLAGS)

clean:
	rm -rf *~ *.o out.bmp $(TARGET)
