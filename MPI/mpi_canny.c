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


// if normalize is true, map pixels to range 0..MAX_BRIGHTNESS
void convolution(const pixel_t *in, pixel_t *out, const float *kernel,
                 const int nx, const int ny, 
                 const int rank, 
                 const int kn,
                 const bool normalize)
{
	assert(rank >= 0);
    assert(kn % 2 == 1);
    assert(nx > kn && ny > kn);
    const int khalf = kn / 2;
    float min = FLT_MAX, max = -FLT_MAX;

    int  start, stop, tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);
    
    printf("Thread[%d] has width=%d and height=%d\n", rank, nx, ny);

    if (rank == 0) {
    	start = khalf;
    	stop  = ny - 4 ;
    } else if (rank == tasks - 1) {
    	start = 4;
    	stop  = ny - khalf;
    } else {
    	start = 4;
    	stop  = ny - 4;
    }
 
    if (normalize)
        for (int m = khalf; m < nx - khalf; m++)
            for (int n = start; n < stop; n++) {
                float pixel = 0.0;
                size_t c = 0;
                for (int j = -khalf; j <= khalf; j++)
                    for (int i = -khalf; i <= khalf; i++) {
                        pixel += in[(n - j) * nx + m - i] * kernel[c];
                        c++;
                    }
                if (pixel < min)
                    min = pixel;
                if (pixel > max)
                    max = pixel;
                }

    // printf("Thread[%d] has start=%d and stop=%d \n", rank, start, stop);

    for (int m = khalf; m < nx - khalf; m++)
        for (int n = start; n < stop; n++) {

            float pixel = 0.0;
            size_t c = 0;
            for (int j = -khalf; j <= khalf; j++)
                for (int i = -khalf; i <= khalf; i++) {
                    pixel += in[(n - j) * nx + m - i] * kernel[c];
                    c++;
                }
 
            if (normalize)
                pixel = MAX_BRIGHTNESS * (pixel - min) / (max - min);
            out[n * nx + m] = (pixel_t)pixel;
        }
}

/*
 * gaussianFilter:
 * http://www.songho.ca/dsp/cannyedge/cannyedge.html
 * determine size of kernel (odd #)
 * 0.0 <= sigma < 0.5 : 3
 * 0.5 <= sigma < 1.0 : 5
 * 1.0 <= sigma < 1.5 : 7
 * 1.5 <= sigma < 2.0 : 9
 * 2.0 <= sigma < 2.5 : 11
 * 2.5 <= sigma < 3.0 : 13 ...
 * kernelSize = 2 * int(2*sigma) + 3;
 */
void gaussian_filter(const pixel_t *in, pixel_t *out,
                     const int nx, const int ny, 
                     const int rank,
                     const float sigma)
{
    const int n = 2 * (int)(2 * sigma) + 3;
    const float mean = (float)floor(n / 2.0);
    float kernel[n * n]; // variable length array
 
    fprintf(stderr, "gaussian_filter: kernel size %d, sigma=%g\n",
            n, sigma);

    size_t c = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            kernel[c] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
                                    pow((j - mean) / sigma, 2.0)))
                        / (2 * M_PI * sigma * sigma);
            c++;
        }
 
    convolution(in, out, kernel, nx, ny, rank, n, true);
}



pixel_t* canny_edge_detection(const pixel_t* in,
                              const bitmap_info_header_t* bmp_ih,
                              int width, int height,
                              int rank,
                              const float sigma)
{
    const int nx = width;
    const int ny = height;
 
    // pixel_t* G = (pixel_t*)calloc(nx * ny, sizeof(pixel_t));
    // assert(G != NULL);
    pixel_t* after_Gx = calloc((bmp_ih->bmp_bytesz/3), sizeof(pixel_t));
    assert(after_Gx != NULL);
    pixel_t* out = (pixel_t*)calloc((bmp_ih->bmp_bytesz/3), sizeof(pixel_t));
    assert(out != NULL);
    pixel_t* local_out = (pixel_t*)calloc(width * height, sizeof(pixel_t));
    if (rank == 0) {
        printf("allocation = %d \n", width * height * sizeof(pixel_t));
        printf("sizeof(local_out) = %d \n", sizeof(local_out));
    }
    assert(local_out != NULL);

    gaussian_filter(in, local_out, nx, ny, rank, sigma);

    memcpy(out, local_out, width * height * sizeof(pixel_t));

    const float Gx[] = {-1, 0, 1,
                        -2, 0, 2,
                        -1, 0, 1};
 
    convolution(out, after_Gx, Gx, nx, ny, rank, 3, false);    

    return out;
}


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
	pixel_t* out_bitmap_data;
	int width, height; 
	int start, stop;
	if (rank == 0) { //Master process
		in_bitmap_data = load_bmp(argv[1], &ih);
    		if (in_bitmap_data == NULL) {
        		fprintf(stderr, "main: BMP image not loaded.\n");
        		return 1;
        	}
        // printf("Master[%d] read .bmp file\n", rank);
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

    // printf("Thread[%d] has chunck_height=%d \n", rank, chunck);

	
	printf("Thread[%d] has start = %d and stop = %d\n", rank, start, stop);
	
	MPI_Bcast(&ih, sizeof(bitmap_info_header_t), MPI_CHAR, 0, MPI_COMM_WORLD);
	if (rank != 0) {
		in_bitmap_data = (pixel_t*)malloc(sizeof(pixel_t) * width * height);
		assert(in_bitmap_data != NULL);
	}
    MPI_Bcast(in_bitmap_data, sizeof(pixel_t) * width * height, MPI_CHAR, 0, MPI_COMM_WORLD);	
	
	char file[255];
	sprintf(file, "mpi_out_blur_test%d.bmp", rank);
    // save_bmp(file, &ih, in_bitmap_data); 

	pixel_t* local_bitmap_data = NULL;
	if (rank == 0) {
    	local_bitmap_data = (pixel_t*)calloc((chunck + 4) * width, sizeof(pixel_t));
    } else if (rank != 0 && rank != tasks - 1) {
    	local_bitmap_data = (pixel_t*)calloc((chunck + 8) * width, sizeof(pixel_t));
    } else {
    	local_bitmap_data = (pixel_t*)calloc((chunck + 4) * width, sizeof(pixel_t));
    }
    assert(local_bitmap_data != NULL);

    if (rank == 0) {
    	memcpy(local_bitmap_data, &in_bitmap_data[start*width], sizeof(pixel_t) * (chunck + 4) * width); // Extra 2 rows at the end
    } else if (rank != 0 && rank != tasks - 1) {
    	memcpy(local_bitmap_data, &in_bitmap_data[(start-4)*width], 
    		   sizeof(pixel_t) * (chunck + 8) * width); // Extra 4 rows(2 up, 2 at the end)
    } else {
    	memcpy(local_bitmap_data, &in_bitmap_data[(start-4)*width], 
    		   sizeof(pixel_t) * (chunck+4) * width); // Extra 2 rows up
    }

	// save_bmp(file, &ih, local_bitmap_data);
	if (rank == 0) {
    	out_bitmap_data = canny_edge_detection(local_bitmap_data, &ih, width, chunck + 4, rank, 1.0f);
    } else if(rank != 0 && rank != tasks - 1) {
    	out_bitmap_data = canny_edge_detection(local_bitmap_data, &ih, width, chunck + 8, rank, 1.0f);
    } else {
    	out_bitmap_data = canny_edge_detection(local_bitmap_data, &ih, width, chunck + 4, rank, 1.0f);
    }

    assert(out_bitmap_data != NULL);
    save_bmp(file, &ih, out_bitmap_data);
	
	free(in_bitmap_data);
	free(local_bitmap_data);
	free(out_bitmap_data);
	MPI_Finalize();
	return 0;
}
