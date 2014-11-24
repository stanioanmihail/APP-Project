CC= gcc
FLAGS= -Wall
SRC= LoadBitmap.c
EXE= LoadBitmap

build:
	$(CC) -o $(EXE) $(SRC) $(FLAGS)

exec:
	./$(EXE) ./images/road.bmp ./output/road.bmp
	./$(EXE) ./images/car.bmp ./output/car.bmp
	./$(EXE) ./images/alps.bmp ./output/alps.bmp

test:
	diff ./images/road.bmp ./output/road.bmp
	diff ./images/car.bmp ./images/car.bmp
	diff ./images/alps.bmp ./output/alps.bmp

clean:
	rm -rf $(EXE) ./output/*

