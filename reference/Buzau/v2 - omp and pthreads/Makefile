all:	build_omp build_omp_debug

build_omp:
	gcc -o omp_canny_edge canny_edge.c hysteresis.c pgm_io.c -fopenmp -Wall -pthread -lm

build_omp_debug:
	gcc -o debug_omp_canny_edge canny_edge.c hysteresis.c pgm_io.c -fopenmp -g -Wall -pthread -lm

clean:
	rm -f *~ omp_canny_edge debug_omp_canny_edge
