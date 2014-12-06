#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include "mpi.h"

#include "file_io.h"
 
#define MAX_BRIGHTNESS 255
 
// C99 doesn't define M_PI (GNU-C99 does)
#define M_PI 3.14159265358979323846264338327

#define _DEBUG 0

// pixel_t *canny_edge_detection(const pixel_t *in,
//                               const bitmap_info_header_t *bmp_ih,
//                               const int tmin, const int tmax,
//                               const float sigma)
// {
//     const int nx = bmp_ih->width;
//     const int ny = bmp_ih->height;
 
//     pixel_t *G = calloc(nx * ny * sizeof(pixel_t), 1);
//     pixel_t *after_Gx = calloc(nx * ny * sizeof(pixel_t), 1);
//     pixel_t *after_Gy = calloc(nx * ny * sizeof(pixel_t), 1);
//     pixel_t *nms = calloc(nx * ny * sizeof(pixel_t), 1);
//     pixel_t *out = malloc(bmp_ih->bmp_bytesz * sizeof(pixel_t));
 
//     if (G == NULL || after_Gx == NULL || after_Gy == NULL ||
//         nms == NULL || out == NULL) {
//         fprintf(stderr, "canny_edge_detection:"
//                 " Failed memory allocation(s).\n");
//         exit(1);
//     }
 
//     gaussian_filter(in, out, nx, ny, sigma);

//     return out;
// }


int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv); //Initialize MPI medium

	int tasks, rank, resultlen;
	// int i,j;
	char proc_name[MPI_MAX_PROCESSOR_NAME];
	// MPI_Status stat;

	MPI_Comm_size(MPI_COMM_WORLD, & tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(proc_name, &resultlen);

	if (_DEBUG) {
		printf("This is thread[%d] from processor %s \n", rank, proc_name);
	}

	bitmap_info_header_t ih; 
	pixel_t* in_bitmap_data;
	// pixel_t* out_bitmap_data;
	int width, height; 
	int start, stop;
	if (rank == 0) { //Master process
		in_bitmap_data = load_bmp(argv[1], &ih);
    		if (in_bitmap_data == NULL) {
        		fprintf(stderr, "main: BMP image not loaded.\n");
        		return 1;
        	}
        printf("Master[%d] read .bmp file\n", rank);
        width = ih.width;
        height = ih.height;
	}

	// Sending width and height
	MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);

	printf("Thread[%d] has width=%d and height=%d\n", rank,width, height);

	int chunck = height / tasks; // Size of chunck for each process
	int reminder = height % tasks; // Reminder if height does not divide with no reminder

	start = rank * chunck; // Starting position for current chunck
	stop  = (rank + 1) * chunck - 1; // Ending position for current chunck
	if (rank == tasks - 1) {
		stop += reminder;
		chunck += reminder;
	}

	
	printf("Thread[%d] has start = %d and stop = %d\n", rank, start, stop);
	
	MPI_Bcast(&ih, sizeof(bitmap_info_header_t), MPI_CHAR, 0, MPI_COMM_WORLD);
	if (rank != 0) {
		in_bitmap_data = (pixel_t*)malloc(sizeof(pixel_t) * width * height);
		assert(in_bitmap_data != NULL);
	}
    MPI_Bcast(in_bitmap_data, sizeof(pixel_t) * width * height, MPI_CHAR, 0, MPI_COMM_WORLD);	
	
	char file[255];
	sprintf(file, "mpi_out_%d.bmp", rank);
    // save_bmp(file, &ih, in_bitmap_data); 

    pixel_t* local_bitmap_data = (pixel_t*)calloc(height * width, sizeof(pixel_t));
    assert(local_bitmap_data != NULL);
    memcpy(local_bitmap_data, &in_bitmap_data[start*width], sizeof(pixel_t) * chunck * width);


	// save_bmp(file, &ih, local_bitmap_data);

	
	free(in_bitmap_data);
	free(local_bitmap_data);
	MPI_Finalize();
	return 0;
}
