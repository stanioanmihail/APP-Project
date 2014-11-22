#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include "LoadBitmap.h"

#define MAX_BRIGHTNESS 255

float min3(float a, float b, float c){

	float temp = (a < b) ? a : b;
	float min = (c < temp) ? c : temp;

	return min;
}

Pixel** convolution(const Pixel **in, const float *kernel, const int height, const int width,
const int kernel_size, const bool normalize){

	assert(kernel_size % 2 == 1);
	assert(height > kernel_size && width > kernel_size);

	const int khalf = kernel_size / 2;

	float min_float = FLT_MAX;
	float max_float = -FLT_MAX;

	int x, y, c, i, j;
	float pixel_R, pixel_G, pixel_B;
	float min_pixel_component;

	Pixel** out; 

	if(normalize){
		
		for(x = 0; x < height; x++){
			for(y = 0; y < width; y++){

				pixel_R = 0.0;
				pixel_G = 0.0;
				pixel_B = 0.0;

				c = 0;

				for(i = -khalf; i <= khalf; i++){
					for(j = -khalf; j <= khalf; j++){
						
						if(x + i < 0 || x + i >= width || y + j < 0 || y + j >= height){
							c++;
							continue;
						}
						pixel_R += in[x+i][y+j].r * kernel[c];
						pixel_G += in[x+i][y+j].r * kernel[c];
						pixel_B += in[x+i][y+j].r * kernel[c];
						c++;
					}
				}
				min_pixel_component = min3(pixel_R, pixel_G, pixel_B);
				if(min_pixel_component < min_float){
					min_float = min_pixel_component;
				}
				if(min_pixel_component > max_float){
					max_float = min_pixel_component;
				}
			}
		}
	}

	out = (Pixel**) malloc(sizeof(Pixel*) * height);
	for(x = 0; x < height; x++){
		out[x] = (Pixel*) malloc(sizeof(Pixel) * width);
	}

	for(x = 0; x < height; x++){
		for(y = 0; y < width; y++){

			pixel_R = 0.0;
			pixel_G = 0.0;
			pixel_B = 0.0;

			c = 0;

			for(i = -khalf; i <= khalf; i++){
				for(j = -khalf; j <= khalf; j++){
						
					if(x + i < 0 || x + i >= width || y + j < 0 || y + j >= height){
						c++;
						continue;
					}
					pixel_R += 1.0 * in[x+i][y+j].r * kernel[c];
					pixel_G += 1.0 * in[x+i][y+j].r * kernel[c];
					pixel_B += 1.0 * in[x+i][y+j].r * kernel[c];
					c++;
				}
			}
			if(normalize){

				pixel_R = MAX_BRIGHTNESS * (pixel_R - min_float) / (max_float - min_float);
				pixel_G = MAX_BRIGHTNESS * (pixel_G - min_float) / (max_float - min_float);
				pixel_B = MAX_BRIGHTNESS * (pixel_B - min_float) / (max_float - min_float);
			}	
			out[x][y].r = pixel_R;
			out[x][y].g = pixel_G;
			out[x][y].b = pixel_B;
		}
	}

	return out;
}

int main(int argc, char *argv[]){

	BitmapInfoHeader bitmapInfoHeader;
	BitmapFileHeader bitmapFileHeader;
	Pixel **PixelMatrix;

	if(argc != 3){
		fprintf(stdout, "./<exec> <input bmp file> <output bmp file>\n");
		fprintf(stderr, "NumberOfArgumentsError\n");
		return 0;
	}
	PixelMatrix = LoadBitmapFile(argv[1], &bitmapInfoHeader, &bitmapFileHeader);
	SaveBitmapFile(argv[2], &bitmapInfoHeader, &bitmapFileHeader, PixelMatrix);
	return 0;
}
