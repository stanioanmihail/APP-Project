#ifndef __LOAD_BITMAP
#define __LOAD_BITMAP

#include <stdio.h>
#include <stdlib.h>

typedef unsigned char BYTE;
typedef unsigned short WORD;
typedef unsigned int DWORD;
typedef int LONG;

typedef struct __attribute__ ((__packed__)) __BitmapFileHeader{

	WORD bfileType; //file type
	DWORD bfileSize; // size in bytes of file
	WORD bfileReserved1; //reserved
	WORD bfileReserved2; //reserved
	DWORD bOffset; //offset in bytes from the header bitmap bits

}BitmapFileHeader;

typedef struct __attribute__ ((__packed__)) __BitmapInfoHeader{
	
	DWORD binfoSize; // size in bytes of struct
	DWORD binfoWidth; // width in pixels
	DWORD binfoHeight; // height in pixels
	WORD binfoPlanes; // = 1
	WORD binfoCount; // number of bits per pixel
	DWORD binfoCompression; 
	DWORD binfoSizeImage; //size of image in bytes
	DWORD binfoXPixelsPMeter;
	DWORD binfoYPixelsPMeter;
	DWORD binfoNrClrUsed; //number of colors used
	DWORD binfoNrClrImportant; //num of important colors
	
}BitmapInfoHeader;

typedef struct __attribute__ ((__packed__)) __PixelFormat{

	BYTE r;
	BYTE g;
	BYTE b;

}Pixel;

#endif //__LOAD_BITMAP
