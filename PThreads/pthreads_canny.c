#define _XOPEN_SOURCE 600
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <pthread.h>

#include "file_io.h"

#define MAX_BRIGHTNESS 255

// C99 doesn't define M_PI (GNU-C99 does)
// #define M_PI 3.14159265358979323846264338327


#define _DEBUG 0

typedef struct _threadArgs {
	int tid; // Thread id
	pixel_t* in_bitmap_data; // Input image
	pixel_t* out; // Output image
	int width;	// Image width
	int height;	// Image height
	bitmap_info_header_t* ih;
} threadArgs;


int num_threads;
pthread_mutex_t mux;

// if normalize is true, map pixels to range 0..MAX_BRIGHTNESS
void convolution(const pixel_t *in, pixel_t *out, const float *kernel,
		const int nx, const int ny,
		const int tid,
		const int kn,
		const bool normalize)
{
	assert(tid >= 0);
	assert(kn % 2 == 1);
	assert(nx > kn && ny > kn);
	const int khalf = kn / 2;
	float min = FLT_MAX, max = -FLT_MAX;

	int  start, stop;

	if (tid == 0) {
		if (kn == 7) { // gaussian_filter
			start = khalf;
			stop  = ny - 3;
		} else { // kn == 3, after_Gx && after_Gy
			start = khalf;
			stop  = ny - 1;
		}
	} else if (tid == num_threads - 1) {
		if (kn == 7) { // gaussian_filter
			start = 3;
			stop  = ny - khalf;
		} else { // kn == 3, after_Gx && after_Gy
			start = 1;
			stop = ny - khalf;
		}
	} else {
		if (kn == 7) { // gaussian_filter
			start = 3;
			stop  = ny - 3;
		} else { // kn == 3, after_Gx && after_Gy
			start = 1;
			stop  = ny - 1;
		}
	}

	if (normalize)
		for (int m = khalf; m < nx - khalf; m++)
			for (int n = start; n <= stop; n++) {
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

	for (int m = khalf; m < nx - khalf; m++)
		for (int n = start; n <= stop; n++) {

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
		const int tid,
		const float sigma)
{
	const int n = 2 * (int)(2 * sigma) + 3;
	const float mean = (float)floor(n / 2.0);
	float kernel[n * n]; // variable length array

	if (_DEBUG) {
		fprintf(stderr, "gaussian_filter: kernel size %d, sigma=%g\n",
				n, sigma);
	}

	size_t c = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			kernel[c] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
						pow((j - mean) / sigma, 2.0)))
				/ (2 * M_PI * sigma * sigma);
			c++;
		}

	convolution(in, out, kernel, nx, ny, tid, n, true);
}



pixel_t* canny_edge_detection(const pixel_t* in,
		const bitmap_info_header_t* bmp_ih,
		int width, int height,
		int tid,
		const int tmin, const int tmax,
		const float sigma)
{
	const int nx = width;
	const int ny = height;

	pixel_t* after_Gx = calloc((width * (height+2)), sizeof(pixel_t));
	assert(after_Gx != NULL);

	pixel_t* after_Gy = calloc((width * (height+2)), sizeof(pixel_t));
	assert(after_Gy != NULL);

	pixel_t* local_out = (pixel_t*)calloc(width * (height+2), sizeof(pixel_t));
	assert(local_out != NULL);

	pixel_t *nms = calloc(nx * (ny+2) * sizeof(pixel_t), 1);
	assert(nms != NULL);

	pixel_t *G = calloc(nx * (ny+2) * sizeof(pixel_t), 1);
	assert(G != NULL);

	gaussian_filter(in, local_out, nx, ny, tid, sigma);


	const float Gx[] = {-1, 0, 1,
		-2, 0, 2,
		-1, 0, 1};

	convolution(local_out, after_Gx, Gx, nx, ny, tid, 3, false);

	const float Gy[] = { 1, 2, 1,
		0, 0, 0,
		-1,-2,-1};

	convolution(local_out, after_Gy, Gy, nx, ny, tid, 3, false);

	for (int i = 1; i < nx - 1; i++)
		for (int j = 1; j < ny - 1; j++) {
			const int c = i + nx * j;
			G[c] = (pixel_t)hypot(after_Gx[c], after_Gy[c]);
		}

	// Non-maximum suppression, straightforward implementation.
	for (int i = 1; i < nx - 1; i++)
		for (int j = 1; j < ny - 1; j++) {
			const int c = i + nx * j;
			const int nn = c - nx;
			const int ss = c + nx;
			const int ww = c + 1;
			const int ee = c - 1;
			const int nw = nn + 1;
			const int ne = nn - 1;
			const int sw = ss + 1;
			const int se = ss - 1;

			const float dir = (float)(fmod(atan2(after_Gy[c],
							after_Gx[c]) + M_PI,
						M_PI) / M_PI) * 8;

			if (((dir <= 1 || dir > 7) && G[c] > G[ee] &&
						G[c] > G[ww]) || // 0 deg
					((dir > 1 && dir <= 3) && G[c] > G[nw] &&
					 G[c] > G[se]) || // 45 deg
					((dir > 3 && dir <= 5) && G[c] > G[nn] &&
					 G[c] > G[ss]) || // 90 deg
					((dir > 5 && dir <= 7) && G[c] > G[ne] &&
					 G[c] > G[sw]))   // 135 deg
				nms[c] = G[c];
			else
				nms[c] = 0;
		}

	// Reuse array
	// used as a stack. nx*ny/2 elements should be enough.
	int *edges = (int*) after_Gy;
	memset(local_out, 0, sizeof(pixel_t) * nx * ny);
	memset(edges, 0, sizeof(pixel_t) * nx * ny);

	// Tracing edges with hysteresis . Non-recursive implementation.
	size_t c = 1;
	for (int j = 1; j < ny - 1; j++)
		for (int i = 1; i < nx - 1; i++) {
			if (nms[c] >= tmax && local_out[c] == 0) { // trace edges
				local_out[c] = MAX_BRIGHTNESS;
				int nedges = 1;
				edges[0] = c;

				do {
					nedges--;
					const int t = edges[nedges];

					int nbs[8]; // neighbours
					nbs[0] = t - nx;     // nn
					nbs[1] = t + nx;     // ss
					nbs[2] = t + 1;      // ww
					nbs[3] = t - 1;      // ee
					nbs[4] = nbs[0] + 1; // nw
					nbs[5] = nbs[0] - 1; // ne
					nbs[6] = nbs[1] + 1; // sw
					nbs[7] = nbs[1] - 1; // se

					for (int k = 0; k < 8; k++)
						if (nms[nbs[k]] >= tmin && local_out[nbs[k]] == 0) {
							local_out[nbs[k]] = MAX_BRIGHTNESS;
							edges[nedges] = nbs[k];
							nedges++;
						}
				} while (nedges > 0);
			}
			c++;
		}

	free(after_Gx);
	free(after_Gy);
	free(G);
	free(nms);

	return local_out;
}


void* run(void* arg) {

	threadArgs* currentArgs = (threadArgs*)arg;
	if (_DEBUG) {
		printf("Thread[%d] says Hello to you\n", currentArgs->tid);
	}

	int start, stop, chunck, reminder;
	pixel_t* local_bitmap_data = NULL;
	pixel_t* out_bitmap_data = NULL;
	pixel_t* rc = NULL;

	if (_DEBUG) {
		printf("Thread[%d] says that total threads is %d \n", currentArgs->tid, num_threads);
	}

	chunck = currentArgs->height / num_threads; // Size of chunck for each thread
	reminder = currentArgs->height % num_threads; // Reminder if height does not divide with no reminder

	start = currentArgs->tid * chunck; // Starting position for current chunck
	stop  = (currentArgs->tid + 1) * chunck - 1; // Ending position for current chunck

	if (currentArgs->tid == num_threads - 1) {
		stop += reminder;
		chunck += reminder;
	}

	if (currentArgs->tid == 0) {
		local_bitmap_data = (pixel_t*)calloc((chunck + 4 + 2) * 
									currentArgs->width, sizeof(pixel_t));
	} else if (currentArgs->tid != 0 && currentArgs->tid != num_threads - 1) {
		local_bitmap_data = (pixel_t*)calloc((chunck + 8 + 2) * 
									currentArgs->width, sizeof(pixel_t));
	} else {
		local_bitmap_data = (pixel_t*)calloc((chunck + 4 + 2) * 
									currentArgs->width, sizeof(pixel_t));
	}
	assert(local_bitmap_data != NULL); // Def check

	if (_DEBUG) {
		printf("Thread[%d] done allocating memory \n", currentArgs->tid);
	}

	if (currentArgs->tid == 0) {
		rc = memcpy(local_bitmap_data, &currentArgs->in_bitmap_data[start*currentArgs->width],
				sizeof(pixel_t) * (chunck + 4) * currentArgs->width); // Extra 2 rows at the end
	} else if (currentArgs->tid != 0 && currentArgs->tid != num_threads - 1) {
		rc = memcpy(local_bitmap_data, &currentArgs->in_bitmap_data[(start-4)*currentArgs->width],
				sizeof(pixel_t) * (chunck + 8) * currentArgs->width); // Extra 4 rows(2 up, 2 at the end)
	} else{
		rc = memcpy(local_bitmap_data, &currentArgs->in_bitmap_data[(start-4)*currentArgs->width],
				sizeof(pixel_t) * (chunck + 4) * currentArgs->width); // Extra 2 rows up
	}
	assert(rc != NULL && rc == local_bitmap_data); // Sanity check!!!!

	if (_DEBUG) {
		printf("Thread[%d] done copying \n", currentArgs->tid);
	}

	if (currentArgs->tid == 0) {
		out_bitmap_data = canny_edge_detection(local_bitmap_data, currentArgs->ih,
			currentArgs->width, chunck + 4, currentArgs->tid, 40, 50, 1.0f);
	} else if (currentArgs->tid != 0 && currentArgs->tid != num_threads - 1) {
		out_bitmap_data = canny_edge_detection(local_bitmap_data, currentArgs->ih,
			currentArgs->width, chunck + 8, currentArgs->tid, 40, 50, 1.0f);
	} else {
		out_bitmap_data = canny_edge_detection(local_bitmap_data, currentArgs->ih,
			currentArgs->width, chunck + 4, currentArgs->tid, 40, 50, 1.0f);
	}
	assert(out_bitmap_data != NULL);

	if (_DEBUG) {
		printf("Thread[%d] done algorithm \n", currentArgs->tid);
	}

	pthread_mutex_lock(&mux);
	if (currentArgs->tid == 0) {
		memcpy(&currentArgs->out[0], &out_bitmap_data[0], chunck * currentArgs->width * sizeof(pixel_t));
	} else if (currentArgs->tid != 0 && currentArgs->tid != num_threads - 1) {
		memcpy(&currentArgs->out[chunck * currentArgs->tid * currentArgs->width],
				&out_bitmap_data[4*currentArgs->width], chunck * currentArgs->width * sizeof(pixel_t));
	} else {
		memcpy(&currentArgs->out[(chunck-reminder) * currentArgs->tid * currentArgs->width], 
				&out_bitmap_data[4*currentArgs->width], chunck * currentArgs->width * sizeof(pixel_t));
	}
	pthread_mutex_unlock(&mux);

	if (_DEBUG) {
		printf("Thread[%d] done making final \n", currentArgs->tid);
	}

	free(local_bitmap_data);
	free(out_bitmap_data);

	pthread_exit(NULL);
}


int main(int argc, char* argv[]) {

	int width, height, rc;
	pixel_t* in_bitmap_data = NULL;
	pixel_t* out = NULL;
	bitmap_info_header_t ih;

	if (argc != 3) {
		perror("Wrong arguments: call with ./<exec> <filename> <num_threads");
		exit(-1);
	}

	num_threads = atoi(argv[2]);
	pthread_t threads[num_threads];
	threadArgs arguments[num_threads];

	in_bitmap_data = load_bmp(argv[1], &ih);
	assert(in_bitmap_data != NULL);

	width = ih.width;
	height = ih.height;

	out = (pixel_t*)malloc(width * height * sizeof(pixel_t));
	assert(out != NULL);

	rc  = pthread_mutex_init(&mux, NULL);
	if (rc != 0) {
		fprintf(stderr,"Could not initialize mutex object. Exiting.\n");
		exit(-1);
	}

	for (int i = 0; i < num_threads; i++) {
		arguments[i].tid = i; // [In]
		arguments[i].in_bitmap_data = in_bitmap_data; // [In]
		arguments[i].ih = &ih; // [In]
		arguments[i].width = width; // [In]
		arguments[i].height = height; // [In]	
		arguments[i].out = out; // [Out]
	}

	for (int i = 0; i < num_threads; i++) {
		rc = pthread_create(&threads[i], NULL, run, (void*)(&arguments[i]));
		assert(rc == 0);
	}

	for (int i = 0; i < num_threads; i++) {
		rc = pthread_join(threads[i], NULL);
		assert(rc == 0);
	}

	rc = pthread_mutex_destroy(&mux);
	if (rc != 0) {
		fprintf(stderr,"Could not destroy mutex object. Exiting. \n");
		exit(-1);
	}

	char file[255];
	sprintf(file, "pthreads_out.bmp");
	save_bmp(file, &ih, out);

	free(in_bitmap_data);
	free(out);

	return 0;

}
