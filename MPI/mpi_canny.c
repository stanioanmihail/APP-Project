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

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv); //Initialize MPI medium

	int tasks, rank, resultlen, i, j;
	char proc_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Status stat;

	MPI_Comm_size(MPI_COMM_WORLD, & tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(proc_name, &resultlen);

	//printf("This is thread[%d] from processor %s \n", rank, proc_name);

	bitmap_info_header_t ih; 
	pixel_t* in_bitmap_data;
	pixel_t* out_bitmap_data;
	char* send_buffer; // Use for sending data
	char* recv_buffer; // Use for receiving data
	int width, height; 
	int start, stop;
	if (rank == 0) { //Master process
		in_bitmap_data = load_bmp(argv[1], &ih);
    		if (in_bitmap_data == NULL) {
        		fprintf(stderr, "main: BMP image not loaded.\n");
        		return 1;
		}
		send_buffer = (char*)malloc(sizeof(char) * (ih.bmp_bytesz * sizeof(pixel_t)));
		assert(send_buffer != NULL);

		// Copy from in_bitmap_data to send_buffer
		if (_DEBUG) {
			printf("[master]%d \n", ih.bmp_bytesz);
			printf("[master]%d %d %d\n", ih.width, ih.height, ih.width * ih.height);
			printf("[master] %d\n", ih.bmp_bytesz * sizeof(pixel_t));
		}
		
		memcpy(send_buffer, in_bitmap_data, ih.bmp_bytesz * sizeof(pixel_t));
		
		

		int send[2];

		send[0] = width  = ih.width;
		send[1] = height = ih.height;

		if (_DEBUG) {
			printf("[Master] sending width = %d and height = %d \n", width, height );
		}

		for (i = 1; i < tasks; i++) {
			MPI_Send(&send[0], 2, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

    	} else {
		int recv[2];

		MPI_Recv(&recv, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
		width = recv[0];
		height = recv[1];
		printf("Thread[%d] received width = %d and height = %d\n", rank, width, height);
		in_bitmap_data = (pixel_t *) malloc(sizeof(pixel_t) * width * height * 3);

		
	}

	int chunck = height / tasks; // Size of chunck for each process
	int reminder = height % tasks; // Reminder if height does not divide with no reminder

	start = rank * chunck; // Starting position for current chunck
	stop  = (rank + 1) * chunck - 1; // Ending position for current chunck
	if (rank == tasks - 1) {
		stop += reminder;
	}

	
	printf("Thread[%d] has start = %d and stop = %d\n", rank, start, stop);
	MPI_Bcast(&ih, sizeof(bitmap_info_header_t), MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(in_bitmap_data, sizeof(pixel_t) * width * height * 3, MPI_CHAR, 0, MPI_COMM_WORLD);	
	char file[255];
	sprintf(file, "mpi_out_%d.bmp", rank);
    	save_bmp(file, &ih, in_bitmap_data); 

	
	MPI_Finalize();
	return 0;
}
