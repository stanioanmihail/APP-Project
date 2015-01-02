#include <stdio.h>
#include <stdlib.h>
#include "LoadBitmap.h"

/*
 * PrintBitmapInfo prints to stdout bitmap file headers
*/
void PrintBitmapInfo(BitmapInfoHeader *bitmapInfoHeader, BitmapFileHeader *bitmapFileHeader){

	printf("=== Bitmap File Header === \n");

	printf("\t\tbFileType = %x \n\
		bFileSize = %d \n\
		bFileReserved1 = %d \n\
		bFileReserved2 = %d \n\
		bOffset = %d \n", bitmapFileHeader->bfileType, bitmapFileHeader->bfileSize, bitmapFileHeader->bfileReserved1,
		bitmapFileHeader->bfileReserved2, bitmapFileHeader->bOffset);

	
	printf("=== Bitmap Info Header === \n");

	printf("\t\tbinfoSize = %d \n\
		binfoWidth = %d \n\
		binfoHeight = %d \n\
		binfoPlanes = %d \n\
		binfoCount = %d \n\
		binfoCompression = %d \n\
		binfoSizeImage = %d \n\
		binfoXPixelsPMeter = %d \n\
		binfoYPixelsPMeter = %d \n\
		binfoNrClrUsed = %d \n\
		binfoNrClrImportant = %d \n", bitmapInfoHeader->binfoSize, bitmapInfoHeader->binfoWidth, bitmapInfoHeader->binfoHeight,
		bitmapInfoHeader->binfoPlanes, bitmapInfoHeader->binfoCount, bitmapInfoHeader->binfoCompression, 
		bitmapInfoHeader->binfoSizeImage, bitmapInfoHeader->binfoXPixelsPMeter, bitmapInfoHeader->binfoYPixelsPMeter,
		bitmapInfoHeader->binfoNrClrUsed, bitmapInfoHeader->binfoNrClrImportant);

}
/* 
* LoadBitmapFile will load a bmp file and put 
* all the data in specific structures for future
* manipulation 
*/
Pixel** LoadBitmapFile(char *filename, BitmapInfoHeader *bitmapInfoHeader, BitmapFileHeader *bfileHeader){

	FILE *bfilePtr;
	Pixel **pixelMatrix;
	DWORD pIndex_x, pIndex_y;
	int ret;

	bfilePtr = fopen(filename, "rb");
	if(bfilePtr == NULL){
		fprintf(stderr, "OpenBitmapFileError\n");
		return NULL;
	}

	// Read file Header 
	ret = fread(bfileHeader, sizeof(BitmapFileHeader), 1, bfilePtr);
	if(ret < 0) return NULL;
	if(bfileHeader->bfileType != 0x4D42){
		fprintf(stderr, "UnknownFormatFileError\n");
		fclose(bfilePtr);
		return NULL;
	}

	// Read image description header
	ret = fread(bitmapInfoHeader, sizeof(BitmapInfoHeader), 1, bfilePtr);
	if(ret < 0) return NULL;
	PrintBitmapInfo(bitmapInfoHeader, bfileHeader);
	fseek(bfilePtr, bfileHeader->bOffset, SEEK_SET);


	pixelMatrix = (Pixel**) malloc(sizeof(Pixel*) * bitmapInfoHeader->binfoHeight);
	if(pixelMatrix == NULL){	
		fprintf(stderr, "UnableToAllocPixelsMatrixError\n");
		free(pixelMatrix);
		fclose(bfilePtr);
		return NULL;

	}
	for(pIndex_x = 0; pIndex_x < bitmapInfoHeader->binfoHeight; pIndex_x++){
		pixelMatrix[pIndex_x] = (Pixel*) malloc(sizeof(Pixel) * bitmapInfoHeader->binfoWidth);
	}

	// The endiannes that describes a RGB pixel is little-endian 
	fseek(bfilePtr, bfileHeader->bOffset, SEEK_SET);
	for(pIndex_x = 0; pIndex_x < bitmapInfoHeader->binfoHeight; pIndex_x++){
		for(pIndex_y =0; pIndex_y < bitmapInfoHeader->binfoWidth; pIndex_y++){
			ret = fread(&pixelMatrix[pIndex_x][pIndex_y].b, sizeof(BYTE), 1, bfilePtr);
			if(ret < 0) return NULL;
			ret = fread(&pixelMatrix[pIndex_x][pIndex_y].g, sizeof(BYTE), 1, bfilePtr);
			if(ret < 0) return NULL;
			ret = fread(&pixelMatrix[pIndex_x][pIndex_y].r, sizeof(BYTE), 1, bfilePtr);
			if(ret < 0) return NULL;
		}
	}

	
	
	fclose(bfilePtr);
	return pixelMatrix;

}

/*
 * SaveBitmapFile will write a bmp file as it is described in the
 * specific structures
*/
void SaveBitmapFile(char *filename, BitmapInfoHeader *bitmapInfoHeader, BitmapFileHeader *bfileHeader, Pixel **pixelMatrix){
	FILE *bfilePtr;
	DWORD pIndex_x, pIndex_y;

	bfilePtr = fopen(filename, "wb");
	if(bfilePtr == NULL){
		fprintf(stderr, "OpenBitmapFileError\n");
		return;
	}

	fwrite(bfileHeader, sizeof(BitmapFileHeader), 1, bfilePtr);
	fwrite(bitmapInfoHeader, sizeof(BitmapInfoHeader), 1, bfilePtr);

	// The endiannes that describes a RGB pixel is little-endian 
	for(pIndex_x = 0; pIndex_x < bitmapInfoHeader->binfoHeight; pIndex_x++){
		for(pIndex_y =0; pIndex_y < bitmapInfoHeader->binfoWidth; pIndex_y++){
			fwrite(&pixelMatrix[pIndex_x][pIndex_y].b, sizeof(BYTE), 1, bfilePtr);
			fwrite(&pixelMatrix[pIndex_x][pIndex_y].g, sizeof(BYTE), 1, bfilePtr);
			fwrite(&pixelMatrix[pIndex_x][pIndex_y].r, sizeof(BYTE), 1, bfilePtr);
		}
	}
	fclose(bfilePtr);
	return;
}

void SaveBitmapFileNTimes(char *filename, BitmapInfoHeader *bitmapInfoHeader, BitmapFileHeader *bfileHeader, Pixel **pixelMatrix, int n){
	FILE *bfilePtr;
	DWORD pIndex_x, pIndex_y;
	int copy_index_x, copy_index_y;

	bfilePtr = fopen(filename, "wb");
	if(bfilePtr == NULL){
		fprintf(stderr, "OpenBitmapFileError\n");
		return;
	}
	
	int orig_img_width, orig_img_height;

	orig_img_width = bitmapInfoHeader->binfoWidth;
	orig_img_height = bitmapInfoHeader->binfoHeight;
	bitmapInfoHeader->binfoWidth *= n;
	bitmapInfoHeader->binfoHeight *= n;
	bitmapInfoHeader->binfoSizeImage *= (n * n);
	bfileHeader->bfileSize = bitmapInfoHeader->binfoSizeImage + bfileHeader->bOffset;
	bitmapInfoHeader->binfoXPixelsPMeter = bitmapInfoHeader->binfoYPixelsPMeter = 0;
	printf("width x heigh = %d x %d = %d => %d \n", bitmapInfoHeader->binfoWidth, bitmapInfoHeader->binfoHeight, bitmapInfoHeader->binfoSizeImage, bfileHeader->bfileSize);

	fwrite(bfileHeader, sizeof(BitmapFileHeader), 1, bfilePtr);
	fwrite(bitmapInfoHeader, sizeof(BitmapInfoHeader), 1, bfilePtr);

	// The endiannes that describes a RGB pixel is little-endian 
	for(copy_index_x = 0; copy_index_x < n; copy_index_x++){
		for(pIndex_x = 0; pIndex_x < orig_img_height; pIndex_x++){
			for(copy_index_y = 0; copy_index_y < n; copy_index_y++){
				for(pIndex_y =0; pIndex_y < orig_img_width; pIndex_y++){
					fwrite(&pixelMatrix[pIndex_x][pIndex_y].b, sizeof(BYTE), 1, bfilePtr);
					fwrite(&pixelMatrix[pIndex_x][pIndex_y].g, sizeof(BYTE), 1, bfilePtr);
					fwrite(&pixelMatrix[pIndex_x][pIndex_y].r, sizeof(BYTE), 1, bfilePtr);
				}
			}
		}
	}
	fclose(bfilePtr);
	return;
}



