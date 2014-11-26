#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include "mpi.h"
 
#define MAX_BRIGHTNESS 255
 
// C99 doesn't define M_PI (GNU-C99 does)
#define M_PI 3.14159265358979323846264338327

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv); //Initialize MPI medium

	int tasks, rank, resultlen;
	char proc_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Comm_size(MPI_COMM_WORLD, & tasks);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Get_processor_name(proc_name, &resultlen);

	//printf("This is thread[%d] from processor %s \n", rank, proc_name);

	MPI_Finalize();

	return 0;
}
