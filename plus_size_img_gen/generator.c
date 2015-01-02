#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "LoadBitmap.h"

int main(int argc, char *argv[]){

	if(argc != 3) return -1;

	char *file;
	char *aux;
	char out_file[100];

	int n = atoi(argv[2]);
	BitmapFileHeader bfileHeader;
	BitmapInfoHeader binfoHeader;
	Pixel** image;

	image = LoadBitmapFile(argv[1], &binfoHeader, &bfileHeader);
	file = aux = argv[1];

	aux = strrchr(aux, '/');
	if(aux != NULL){
		file = aux + 1;
	}

	file = strtok(file, ".");
	sprintf(out_file, "%s_x%d.bmp", file, n); 
	
	SaveBitmapFileNTimes(out_file, &binfoHeader, &bfileHeader, image, n);
	PrintBitmapInfo(&binfoHeader, &bfileHeader);

	

	return 0;
}

