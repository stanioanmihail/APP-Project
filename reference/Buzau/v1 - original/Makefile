all:	build_normal build_debug

build_normal:
	gcc -o canny_edge canny_edge.c hysteresis.c pgm_io.c -Wall -lm

build_debug:
	gcc -o canny_edge_debug canny_edge.c hysteresis.c pgm_io.c -g -Wall -lm

clean:
	rm -f *~ canny_edge canny_edge_debug
