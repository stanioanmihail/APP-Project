#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
 
#define MAX_BRIGHTNESS 255
 
// C99 doesn't define M_PI (GNU-C99 does)
#define M_PI 3.14159265358979323846264338327

#define DEBUG 0 

/*
 * Loading part taken from
 * http://www.vbforums.com/showthread.php?t=261522
 * BMP info:
 * http://en.wikipedia.org/wiki/BMP_file_format
 *
 * Note: the magic number has been removed from the bmpfile_header_t
 * structure since it causes alignment problems
 *     bmpfile_magic_t should be written/read first
 * followed by the
 *     bmpfile_header_t
 * [this avoids compiler-specific alignment pragmas etc.]
 */
 
typedef struct {
    uint8_t magic[2];
} bmpfile_magic_t;
 
typedef struct {
    uint32_t filesz;
    uint16_t creator1;
    uint16_t creator2;
    uint32_t bmp_offset;
} bmpfile_header_t;
 
typedef struct {
    uint32_t header_sz;
    int32_t  width;
    int32_t  height;
    uint16_t nplanes;
    uint16_t bitspp;
    uint32_t compress_type;
    uint32_t bmp_bytesz;
    int32_t  hres;
    int32_t  vres;
    uint32_t ncolors;
    uint32_t nimpcolors;
} bitmap_info_header_t;
 
typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
} rgb_t;
 
// Use short int instead `unsigned char' so that we can
// store negative values.
typedef unsigned char pixel_t;


float min3(float a, float b, float c) {
	
	float temp = (a < b) ? a : b;
	float min  = (c < temp) ? c : temp;

	return min;
}
 
rgb_t *load_bmp(const char *filename,
                  bitmap_info_header_t *bitmapInfoHeader)
{
    FILE *filePtr = fopen(filename, "rb");
    if (filePtr == NULL) {
        perror("fopen()");
        return NULL;
    }
 
    bmpfile_magic_t mag;
    if (fread(&mag, sizeof(bmpfile_magic_t), 1, filePtr) != 1) {
        fclose(filePtr);
        return NULL;
    }
 
    // verify that this is a bmp file by check bitmap id
    // warning: dereferencing type-punned pointer will break
    // strict-aliasing rules [-Wstrict-aliasing]
    if (*((uint16_t*)mag.magic) != 0x4D42) {
        fprintf(stderr, "Not a BMP file: magic=%c%c\n",
                mag.magic[0], mag.magic[1]);
        fclose(filePtr);
        return NULL;
    }
 
    bmpfile_header_t bitmapFileHeader; // our bitmap file header
    // read the bitmap file header
    if (fread(&bitmapFileHeader, sizeof(bmpfile_header_t),
              1, filePtr) != 1) {
        fclose(filePtr);
        return NULL;
    }
 
    // read the bitmap info header
    if (fread(bitmapInfoHeader, sizeof(bitmap_info_header_t),
              1, filePtr) != 1) {
        fclose(filePtr);
        return NULL;
    }
 
    if (bitmapInfoHeader->compress_type != 0)
        fprintf(stderr, "Warning, compression is not supported.\n");
 
    printf("[Read] Offset is = %d\n",  bitmapFileHeader.bmp_offset);
    // move file point to the beginning of bitmap data
    if (fseek(filePtr, bitmapFileHeader.bmp_offset, SEEK_SET)) {
        fclose(filePtr);
        return NULL;
    }
 
    // allocate enough memory for the bitmap image data
    rgb_t *bitmapImage = malloc(bitmapInfoHeader->bmp_bytesz *
                                  sizeof(rgb_t));

    if(DEBUG) { 
    	//printf("Total size of bitmapImage in byes is = %d\n",bitmapInfoHeader->bmp_bytesz);
    	//printf("Total size of bitmapImage is=%d\n", bitmapInfoHeader->bmp_bytesz *sizeof(pixel_t));
    }

    // verify memory allocation
    if (bitmapImage == NULL) {
        fclose(filePtr);
        return NULL;
    }
 
    // read in the bitmap image data
    size_t count=0, pad1;
    rgb_t c;
    pad1 = 4*floor((bitmapInfoHeader->bitspp*bitmapInfoHeader->width+31)/32) - (bitmapInfoHeader->width * 3);

    if(DEBUG) {
    	//printf("pad1 = %d\n", pad1);
    	//printf("Bits per pixel=%d\n",bitmapInfoHeader->bitspp);
    	//printf("Pixel per row=%d\n",bitmapInfoHeader->width);
    }

    for(size_t i=0; i<bitmapInfoHeader->height; i++){
	    for(size_t j=0; j<bitmapInfoHeader->width; j++){
		    	if (fread(&c, sizeof(rgb_t), 1, filePtr) != 1) {
				    fclose(filePtr);
				    return NULL;
		    	}
		    	bitmapImage[count++] = c;
	    }
	    fseek(filePtr, pad1, SEEK_CUR);
    }

    if(DEBUG) {
    	//printf("Count is = %d\n", count); 
    }
 
    // close file and return bitmap image data
    fclose(filePtr);
    return bitmapImage;
}
 
// Return: true on error.
bool save_bmp(const char *filename, const bitmap_info_header_t *bmp_ih,
              const rgb_t *data)
{
    FILE* filePtr = fopen(filename, "wb");
    if (filePtr == NULL)
        return true;
 
    bmpfile_magic_t mag = {{0x42, 0x4d}};
    if (fwrite(&mag, sizeof(bmpfile_magic_t), 1, filePtr) != 1) {
        fclose(filePtr);
        return true;
    }

    //Am sters shftarea lui 1 pentru ca nu folosim o tabele de culori intre ultimul header si aray-ul de pixeli 
    const uint32_t offset = sizeof(bmpfile_magic_t) +
                            sizeof(bmpfile_header_t) +
                            sizeof(bitmap_info_header_t);

    if(DEBUG) {
    	printf("Offset is = %d\n", offset); 
    	printf("Size   is = %d\n", bmp_ih->bmp_bytesz);
    }

    const bmpfile_header_t bmp_fh = {
        .filesz = offset + bmp_ih->bmp_bytesz,
        .creator1 = 0,
        .creator2 = 0,
        .bmp_offset = offset
    };
 
    if (fwrite(&bmp_fh, sizeof(bmpfile_header_t), 1, filePtr) != 1) {
        fclose(filePtr);
        return true;
    }
    if (fwrite(bmp_ih, sizeof(bitmap_info_header_t), 1, filePtr) != 1) {
        fclose(filePtr);
        return true;
    }
 
    size_t pad1 = 4*floor((bmp_ih->bitspp*bmp_ih->width+31)/32) - (bmp_ih->width * 3); 
   
    if(DEBUG) { 
	    //printf("Pad1 is = %d\n", pad1);
    }

    unsigned char c;
    int k = 0;
    for(size_t i=0; i < bmp_ih->height; i++) {
	    for(size_t j=0; j < bmp_ih->width; j++) {//Aceeasi idee ca la citire; extragem cate 3 elemente la fiecare iteratie din width
		rgb_t pixel = data[k++];
		if (fwrite(&pixel, sizeof(rgb_t), 1, filePtr) != 1) {
			fclose(filePtr);
			return true;
		}
	    }

	    c = 0;
	    for(size_t j=0; j < pad1; j++)
		    if (fwrite(&c, sizeof(unsigned char), 1, filePtr) != 1) {
			    fclose(filePtr);
			    return true;
		    }
    }
 
    fclose(filePtr);
    return false;
}
 
// if normalize is true, map pixels to range 0..MAX_BRIGHTNESS
void convolution(const rgb_t *in, rgb_t *out, const float *kernel,
                 const int nx, const int ny, const int kn,
                 const bool normalize)
{
    assert(kn % 2 == 1);
    assert(nx > kn && ny > kn);
    const int khalf = kn / 2;
    float min = FLT_MAX, max = -FLT_MAX;
 
    if (normalize) {
        for (int n = 0; n < ny; n++) {
            for (int m = 0; m < nx; m++) {
                	float pixelR = 0.0;
			float pixelG = 0.0;
			float pixelB = 0.0;
                	int c = 0;
			int row = n * nx;
			int col = m;
                	for (int i = -khalf; i <= khalf; i++) {
                    		for (int j = -khalf; j <= khalf; j++) {
					int row_sub = row + i * nx;
					int col_sub = col + j;
					int elem_sub = row_sub + col_sub;
					if (row_sub < 0 || row_sub >= nx*ny || (elem_sub < row_sub || elem_sub > (row_sub +nx))) {
						c++;
						continue;
					}
					pixelR += in[elem_sub].r * kernel[c];
					pixelG += in[elem_sub].g * kernel[c];
					pixelB += in[elem_sub].b * kernel[c];
					c++;
					
                    		}
			}
			float pixel = min3(pixelR, pixelG, pixelB);
                	if (pixel < min) {
                    		min = pixel;
			}
			if (pixel > max) {
                    		max = pixel;
			}
      	     }
	}
     }

	int counter = 0;
	for (int n = 0; n < ny; n++) {
            for (int m = 0; m < nx; m++) {
                	float pixelR = 0.0;
			float pixelG = 0.0;
			float pixelB = 0.0;
                	int c = 0;
			int row = n * nx;
			int col = m;
                	for (int i = -khalf; i <= khalf; i++) {
                    		for (int j = -khalf; j <= khalf; j++) {
					int row_sub = row + i * nx;
					int col_sub = col + j;
					int elem_sub = row_sub + col_sub;
					if (row_sub < 0 || row_sub >= nx*ny || (elem_sub < row_sub || elem_sub > (row_sub +nx))) {
						c++;
						continue;
					}
					pixelR += 1.0*in[elem_sub].r * kernel[c];
					pixelG += 1.0*in[elem_sub].g * kernel[c];
					pixelB += 1.0*in[elem_sub].b * kernel[c];
					c++;
                    		}
			}
			if (normalize) {
				pixelR = MAX_BRIGHTNESS * (pixelR - min) / (max - min);
				pixelG = MAX_BRIGHTNESS * (pixelG - min) / (max - min);
				pixelB = MAX_BRIGHTNESS * (pixelB - min) / (max - min);
			}
			out[counter].r = pixelR;
			out[counter].g = pixelG;
			out[counter].b = pixelB;
			counter++;
	    }
	}
 
    /*for (int n = khalf; n < nx - khalf; n++)
        for (int m = khalf; m < ny - khalf; m++) {
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
        }*/
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
void gaussian_filter(const rgb_t *in, rgb_t *out,
                     const int nx, const int ny, const float sigma)
{
    const int n = 2 * (int)(2 * sigma) + 3;
    const float mean = (float)floor(n / 2.0);
    float kernel[n * n]; // variable length array
 
    fprintf(stderr, "gaussian_filter: kernel size %d, sigma=%g\n",
            n, sigma);
    size_t c = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            kernel[c] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
                                    pow((j - mean) / sigma, 2.0)))
                        / (2 * M_PI * sigma * sigma);
            c++;
        }
    }
 
    convolution(in, out, kernel, nx, ny, n, false);
}

int condition(int index, int row, int width) {

	if ( row * width <= index && (row * width + (width-1)) <= index ) {
		return 1;
	} else {
		return 0;
	}

}


void my_switch(int k, rgb_t* G, rgb_t* nms, int c, int ok) {

	/*
	*	[john] observatii:
	*      -  daca pui 0 in loc de 255 vei observa ca toata poza se face alba mai putin marginile	
	*      -  problema ar putea fi, dar e vaga solutia, ca matrice nms sa nu fie direct matricea out 
			de la blurul gaussian ci o matrice copie ( nu pot justifica de ce - feeling )
	*
	*/
	switch(k) {
		case 0:
			if (ok) {
				
				nms[c].r = G[c].r;
				//nms[c].g = G[c].g;
				//nms[c].b = G[c].b;
				//nms[c].r = 255;
				//nms[c].r = 0;
			} else {
				//nms[c].r = G[c].r;
				nms[c].r = 0;
				//nms[c].g = 0;
				//nms[c].b = 0;
			}
			break;
		case 1:
			if (ok) {
				//nms[c].r = G[c].r;
				nms[c].g = G[c].g;
				//nms[c].b = G[c].b;
				//nms[c].g = 255;
				//nms[c].g = 0;
			} else {
				//nms[c].g = G[c].g;
				//nms[c].r = 0;
				nms[c].g = 0;
				//nms[c].b = 0;
			}
			break;
		case 2:
			if (ok) {
				//nms[c].r = G[c].r;
				//nms[c].g = G[c].g;
				nms[c].b = G[c].b;
				//nms[c].b = 255;
				//nms[c].b = 0;
			} else {
				//nms[c].b = G[c].b;
				//nms[c].r = 0;
				//nms[c].g = 0;
				nms[c].b = 0;
			}
			break;
		default:
			break;
	}
}

void suppression(const float dirR, int k , int i, int nx, rgb_t* G, rgb_t* nms, const int c, const int nn, const int ss,
		 const int ee, const int ww, const int nw, const int se, const int ne, const int sw) {

	/* [john] varianta initiala facea aceaiasi operatie pe toate componentele de culoare comparand cu G[c].r
	*
	*/


	//k == 0 means color is red
	if(k == 0){
	   if (dirR <= 1 || dirR >7) {
		if (condition(ee, i * nx, nx) && condition(ww, i * nx, nx)) {
			if (G[c].r > G[ee].r && G[c].r > G[ww].r) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(ee, i * nx, nx) && (condition(ww, i * nx, nx) && G[c].r > G[ww].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(ww, i * nx, nx) && (condition(ee, i * nx, nx) && G[c].r > G[ee].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 1 && dirR <= 3) {
		if (condition(nw, (i-1) * nx, nx) && condition(se, (i+1) * nx, nx)) {
			if (G[c].r > G[nw].r && G[c].r > G[se].r) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(nw, (i-1) * nx, nx) && (condition(se, (i+1) * nx, nx) && G[c].r > G[se].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(se, (i+1) * nx, nx) && (condition(nw, (i-1) * nx, nx) && G[c].r > G[nw].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 3 && dirR <= 5) {
		if (condition(nn, (i-1) * nx, nx) && condition(ss, (i+1) * nx, nx)) {
			if (G[c].r > G[nn].r && G[c].r > G[ss].r) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(nn, (i-1) * nx, nx) && (condition(ss, (i+1) * nx, nx) && G[c].r > G[ss].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(ss, (i+1) * nx, nx) && (condition(nn, (i-1) * nx, nx) && G[c].r > G[nn].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 5 && dirR <=7) {
		if (condition(ne, (i-1) * nx, nx) && condition(sw, (i+1) * nx, nx)) {
			if (G[c].r > G[ne].r && G[c].r > G[sw].r) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(ne, (i-1) * nx, nx) && (condition(sw, (i+1) * nx, nx) && G[c].r > G[sw].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(sw, (i+1) * nx, nx) && (condition(ne, (i-1) * nx, nx) && G[c].r > G[ne].r)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    }
	}

	//k == 1 means color is green
	if(k == 1){
	    
	   if (dirR <= 1 || dirR >7) {
		if (condition(ee, i * nx, nx) && condition(ww, i * nx, nx)) {
			if (G[c].g > G[ee].g && G[c].g > G[ww].g) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(ee, i * nx, nx) && (condition(ww, i * nx, nx) && G[c].g > G[ww].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(ww, i * nx, nx) && (condition(ee, i * nx, nx) && G[c].g > G[ee].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 1 && dirR <= 3) {
		if (condition(nw, (i-1) * nx, nx) && condition(se, (i+1) * nx, nx)) {
			if (G[c].g > G[nw].g && G[c].g > G[se].g) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(nw, (i-1) * nx, nx) && (condition(se, (i+1) * nx, nx) && G[c].g > G[se].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(se, (i+1) * nx, nx) && (condition(nw, (i-1) * nx, nx) && G[c].g > G[nw].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 3 && dirR <= 5) {
		if (condition(nn, (i-1) * nx, nx) && condition(ss, (i+1) * nx, nx)) {
			if (G[c].g > G[nn].g && G[c].g > G[ss].g) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(nn, (i-1) * nx, nx) && (condition(ss, (i+1) * nx, nx) && G[c].g > G[ss].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(ss, (i+1) * nx, nx) && (condition(nn, (i-1) * nx, nx) && G[c].g > G[nn].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 5 && dirR <=7) {
		if (condition(ne, (i-1) * nx, nx) && condition(sw, (i+1) * nx, nx)) {
			if (G[c].g > G[ne].g && G[c].g > G[sw].g) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(ne, (i-1) * nx, nx) && (condition(sw, (i+1) * nx, nx) && G[c].g > G[sw].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(sw, (i+1) * nx, nx) && (condition(ne, (i-1) * nx, nx) && G[c].g > G[ne].g)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    }
	}

	//k == 2 means color is blue
	if(k == 2){
	    
	   if (dirR <= 1 || dirR >7) {
		if (condition(ee, i * nx, nx) && condition(ww, i * nx, nx)) {
			if (G[c].b > G[ee].b && G[c].b > G[ww].b) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(ee, i * nx, nx) && (condition(ww, i * nx, nx) && G[c].b > G[ww].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(ww, i * nx, nx) && (condition(ee, i * nx, nx) && G[c].b > G[ee].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 1 && dirR <= 3) {
		if (condition(nw, (i-1) * nx, nx) && condition(se, (i+1) * nx, nx)) {
			if (G[c].b > G[nw].b && G[c].b > G[se].b) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(nw, (i-1) * nx, nx) && (condition(se, (i+1) * nx, nx) && G[c].b > G[se].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(se, (i+1) * nx, nx) && (condition(nw, (i-1) * nx, nx) && G[c].b > G[nw].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 3 && dirR <= 5) {
		if (condition(nn, (i-1) * nx, nx) && condition(ss, (i+1) * nx, nx)) {
			if (G[c].b > G[nn].b && G[c].b > G[ss].b) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(nn, (i-1) * nx, nx) && (condition(ss, (i+1) * nx, nx) && G[c].b > G[ss].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(ss, (i+1) * nx, nx) && (condition(nn, (i-1) * nx, nx) && G[c].b > G[nn].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    } else if (dirR > 5 && dirR <=7) {
		if (condition(ne, (i-1) * nx, nx) && condition(sw, (i+1) * nx, nx)) {
			if (G[c].b > G[ne].b && G[c].b > G[sw].b) {
				//nms[c].r = G[c].r;
				my_switch(k, G, nms, c, 1);
			}
		} else if (!condition(ne, (i-1) * nx, nx) && (condition(sw, (i+1) * nx, nx) && G[c].b > G[sw].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else if (!condition(sw, (i+1) * nx, nx) && (condition(ne, (i-1) * nx, nx) && G[c].b > G[ne].b)) {
			//nms[c].r = G[c].r;
			my_switch(k, G, nms, c, 1);
		} else {
			//nms[c].r = 0;
			my_switch(k, G, nms, c, 0);
		}
	    }
	}

}


/*
void hysteresis(rgb_t* out, rgb_t* nms, int row, int width, int c, int* edges, const int tmax, const int tmin) {
 
	int current_row[10000];
	int process = 1;
	switch(type) {
		case 0://R
			if (nms[c].r >= tmax && out[c].r == 0) {
				process = 1;
			} else {
				process = 0;
			}
			break;
		case 1://G
			if (nms[c].g >= tmax && out[c].g == 0) {
				process = 1;
			} else {
				process = 0;
			}
			break;
		case 2://B
			if (nms[c].b >= tmax && out[c].b == 0) {
				process = 1;
			} else {
				process = 0;
			}
			break;
		default:
			break;
	}

	if (process) { // trace edges
                
		if (type == 0) { //R
			out[c].r = MAX_BRIGHTNESS;
		} else if (type == 1) { //G
			out[c].g = MAX_BRIGHTNESS;
		} else {//B
			out[c].b = MAX_BRIGHTNESS;
		}

                int nedges = 1;
                edges[0] = c;
		current_row[0] = row;
 
                do {
                    nedges--;
                    const int t = edges[nedges];
		    int row = current_row[nedges];
 
                    int nbs[8]; // neighbours
                    nbs[0] = t - width;     // nn
                    nbs[1] = t + width;     // ss
                    nbs[2] = t + 1;      // ww
                    nbs[3] = t - 1;      // ee
                    nbs[4] = nbs[0] + 1; // nw
                    nbs[5] = nbs[0] - 1; // ne
                    nbs[6] = nbs[1] + 1; // sw
                    nbs[7] = nbs[1] - 1; // se
 
                    for (int k = 0; k < 8; k++) {
			int ok = 1;
			if (k == 0 || k == 4 || k == 5) {
				ok = condition(nbs[k], row-width, width);//Randul de sus
			} else if (k == 1 || k == 6 || k == 7) {
				ok = condition(nbs[k], row+width, width);//Randul de jos
			} else {
				ok = condition(nbs[k], row, width); //Rand curent
			}
			int process = 1;
			if (ok) {
				switch(type) {
					case 0://R
						if (nms[nbs[k]].r >= tmin && out[nbs[k]].r == 0) {
							process = 1;
						} else {
							process = 0;
						}
						break;
					case 1://G
						if (nms[nbs[k]].g >= tmin && out[nbs[k]].g == 0) {
							process = 1;
						} else {
							process = 0;
						}
						break;
					case 2://B
						if (nms[nbs[k]].b >= tmin && out[nbs[k]].b == 0) {
							process = 1;
						} else {
							process = 0;
						}
						break;
					default:
						break;
				}
			}

                        if (process && ok) {
                           
			   if (type == 0) { //R 
				out[nbs[k]].r = MAX_BRIGHTNESS;
			   } else if (type == 1) {//G
				out[nbs[k]].g = MAX_BRIGHTNESS;
			   } else {//B
				out[nbs[k]].b = MAX_BRIGHTNESS;
			   }

                            edges[nedges] = nbs[k];
			    if (k == 0 || k == 4 || k == 5) {
				    current_row[nedges] = row - width;
			    } else if(k == 1 || k == 6 || k == 7) {
				    current_row[nedges] = row + width;
			    } else {
				    current_row[nedges] = row;
			    }
                            nedges++;
                        }
		    }
                } while (nedges > 0);
            }
}
*/

/*
 * Links:
 * http://en.wikipedia.org/wiki/Canny_edge_detector
 * http://www.tomgibara.com/computer-vision/CannyEdgeDetector.java
 * http://fourier.eng.hmc.edu/e161/lectures/canny/node1.html
 * http://www.songho.ca/dsp/cannyedge/cannyedge.html
 *
 * Note: T1 and T2 are lower and upper thresholds.
 */
rgb_t *canny_edge_detection(const rgb_t *in,
                              const bitmap_info_header_t *bmp_ih,
                              const int tmin, const int tmax,
                              const float sigma)
{
    const int nx = bmp_ih->width;
    const int ny = bmp_ih->height;
 
    rgb_t *G = calloc(nx * ny * sizeof(rgb_t), 1);
    rgb_t *after_Gx = calloc(nx * ny * sizeof(rgb_t), 1);
    rgb_t *after_Gy = calloc(nx * ny * sizeof(rgb_t), 1);
    rgb_t *nms = calloc(nx * ny * sizeof(rgb_t), 1);
    rgb_t *out = malloc(bmp_ih->bmp_bytesz * sizeof(rgb_t));
 
    if (G == NULL || after_Gx == NULL || after_Gy == NULL ||
        nms == NULL || out == NULL) {
        fprintf(stderr, "canny_edge_detection:"
                " Failed memory allocation(s).\n");
        exit(1);
    }
 
    gaussian_filter(in, out, nx, ny, sigma);

    //Comment this part onward if you want to see if gaussian filter is working properly 
    const float Gx[] = {-1, 0, 1,
                        -2, 0, 2,
                        -1, 0, 1};
 
    convolution(out, after_Gx, Gx, nx, ny, 3, true);
 
    const float Gy[] = { 1, 2, 1,
                         0, 0, 0,
                        -1,-2,-1};
 
    convolution(out, after_Gy, Gy, nx, ny, 3, true);
 
    int c = 0;
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            G[c].r = hypot(after_Gx[c].r, after_Gy[c].r);
	    G[c].g = hypot(after_Gx[c].g, after_Gy[c].g);
	    G[c].b = hypot(after_Gx[c].b, after_Gy[c].b);
	    c++;
        }
    }

    /*
     * [john] intrebare: de ce G si nms sunt de dimensiuni diferite fata de matricea mama out
     * - am incercat sa pun la supression
	    suppression(dirR, 0, i, nx, out, nns, c, nn, ss, ee, ww, nw, se, ne, sw); si mi-am luat segfault
	- nu prea stiu de unde vin si de unde pleaca matricele ..asa ca ramane sa te mai intreb 
    */
 
    // Non-maximum suppression, straightforward implementation.
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            const int c = i * nx + j;
            const int nn = c - nx;
            const int ss = c + nx;
            const int ww = c + 1;
            const int ee = c - 1;
            const int nw = nn + 1;
            const int ne = nn - 1;
            const int sw = ss + 1;
            const int se = ss - 1;
 
            const float dirR = (float)(fmod(atan2(after_Gy[c].r,
                                                 after_Gx[c].r) + M_PI,
                                           M_PI) / M_PI) * 8;
 	    const float dirG = (float)(fmod(atan2(after_Gy[c].g,
                                                 after_Gx[c].g) + M_PI,
                                           M_PI) / M_PI) * 8;
	    const float dirB = (float)(fmod(atan2(after_Gy[c].b,
                                                 after_Gx[c].b) + M_PI,
                                           M_PI) / M_PI) * 8;

	   //const float dirR = (float) ((atan2(after_Gy[c].r, after_Gx[c].r) + M_PI) * 180) / M_PI;
	   //const float dirG = (float) ((atan2(after_Gy[c].g, after_Gx[c].g) + M_PI) * 180) / M_PI;
	   //const float dirB = (float) ((atan2(after_Gy[c].b, after_Gx[c].b) + M_PI) * 180) / M_PI;
	   //Al 2-lea argument trimis imi dicteaza mie ce am: r, g sau b
	    //suppression(dirR, 0, i, nx, G, nms, c, nn, ss, ee, ww, nw, se, ne, sw);
	    //suppression(dirG, 1, i, nx, G, nms, c, nn, ss, ee, ww, nw, se, ne, sw);
	    //suppression(dirB, 2, i, nx, G, nms, c, nn, ss, ee, ww, nw, se, ne, sw);

	    suppression(dirR, 0, i, nx, G, out, c, nn, ss, ee, ww, nw, se, ne, sw);
	    suppression(dirG, 1, i, nx, G, out, c, nn, ss, ee, ww, nw, se, ne, sw);
	    suppression(dirB, 2, i, nx, G, out, c, nn, ss, ee, ww, nw, se, ne, sw);
	    /* ce-am facut aici a fost sa iau fiecare directie in parte si sa le tratez pe bucatele
	     * pentru ca baiatul care a facut codul nu a tratat deloc cazurile de la margine
	     * Dupa o sa vezi verificarile sunt simetrice pentru fiecare directie, adica
	     * mai intai verific ca vecinii sa existe, apoi daca nu exista primul iar la final daca nu exista cel de-al 2-lea 
	    if (dirR <= 1 || dirR >7) {
		if (condition(ee, i * nx, nx) && condition(ww, i * nx, nx)) {
			if (G[c].r > G[ee].r && G[c].r > G[ww].r) {
				nms[c].r = G[c].r;
			}
		} else if (!condition(ee, i * nx, nx) && (condition(ww, i * nx, nx) && G[c].r > G[ww].r)) {
			nms[c].r = G[c].r;
		} else if (!condition(ww, i * nx, nx) && (condition(ee, i * nx, nx) && G[c].r > G[ee].r)) {
			nms[c].r = G[c].r;
		} else {
			nms[c].r = 0;
		}
	    } else if (dirR > 1 && dirR <= 3) {
		if (condition(nw, (i-1) * nx, nx) && condition(se, (i+1) * nx, nx)) {
			if (G[c].r > G[nw].r && G[c].r > G[se].r) {
				nms[c].r = G[c].r;
			}
		} else if (!condition(nw, (i-1) * nx, nx) && (condition(se, (i+1) * nx, nx) && G[c].r > G[se].r)) {
			nms[c].r = G[c].r;
		} else if (!condition(se, (i+1) * nx, nx) && (condition(nw, (i-1) * nx, nx) && G[c].r > G[nw].r)) {
			nms[c].r = G[c].r;
		} else {
			nms[c].r = 0;
		}
	    } else if (dirR > 3 && dirR <= 5) {
		if (condition(nn, (i-1) * nx, nx) && condition(ss, (i+1) * nx, nx)) {
			if (G[c].r > G[nn].r && G[c].r > G[ss].r) {
				nms[c].r = G[c].r;
			}
		} else if (!condition(nn, (i-1) * nx, nx) && (condition(ss, (i+1) * nx, nx) && G[c].r > G[ss].r)) {
			nms[c].r = G[c].r;
		} else if (!condition(ss, (i+1) * nx, nx) && (condition(nn, (i-1) * nx, nx) && G[c].r > G[nn].r)) {
			nms[c].r = G[c].r;
		} else {
			nms[c].r = 0;
		}
	    } else if (dirR > 5 && dirR <=7) {
		if (condition(ne, (i-1) * nx, nx) && condition(sw, (i+1) * nx, nx)) {
			if (G[c].r > G[ne].r && G[c].r > G[sw].r) {
				nms[c].r = G[c].r;
			}
		} else if (!condition(ne, (i-1) * nx, nx) && (condition(sw, (i+1) * nx, nx) && G[c].r > G[sw].r)) {
			nms[c].r = G[c].r;
		} else if (!condition(sw, (i+1) * nx, nx) && (condition(ne, (i-1) * nx, nx) && G[c].r > G[ne].r)) {
			nms[c].r = G[c].r;
		} else {
			nms[c].r = 0;
		}
	    }


            if (((dirR <= 1 || dirR > 7) && G[c].r > G[ee].r &&
                 G[c].r > G[ww].r) || // 0 deg
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
	*/ 
	}printf("ajung aici\n");
    // Reuse array
    // used as a stack. nx*ny/2 elements should be enough.
/*
    int *edges = (int*) after_Gy;
    memset(out, 0, sizeof(rgb_t) * bmp_ih->bmp_bytesz);
    memset(edges, 0, sizeof(rgb_t) * nx * ny);
 
    // Tracing edges with hysteresis . Non-recursive implementation.
    size_t c = 0;
    c = 0;
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            if (nms[c] >= tmax && out[c] == 0) { // trace edges
                out[c] = MAX_BRIGHTNESS;
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
                        if (nms[nbs[k]] >= tmin && out[nbs[k]] == 0) {
                            out[nbs[k]] = MAX_BRIGHTNESS;
                            edges[nedges] = nbs[k];
                            nedges++;
                        }
                } while (nedges > 0);
            }
	    hysteresis(out, nms, 0, i, nx, c, edges, tmax, tmin);
            c++;
        }
 */
    //free(after_Gx);
    //free(after_Gy);
    //free(G);
    free(nms);

    return out;
}
 
int main(const int argc, const char ** const argv)
{
    if (argc < 2) {
        printf("Usage: %s image.bmp\n", argv[0]);
        return 1;
    }
 
    static bitmap_info_header_t ih;
    const rgb_t* in_bitmap_data = load_bmp(argv[1], &ih);
    if (in_bitmap_data == NULL) {
        fprintf(stderr, "main: BMP image not loaded.\n");
        return 1;
    }
 
    printf("Info: %d x %d x %d\n", ih.width, ih.height, ih.bitspp);
 
    const rgb_t *out_bitmap_data =
        canny_edge_detection(in_bitmap_data, &ih, 45, 50, 1.0f);
    if (out_bitmap_data == NULL) {
        fprintf(stderr, "main: failed canny_edge_detection.\n");
        return 1;
    }

    //!!!Ca sa testezi blur-ul trebuie sa pui out_bitmap_data
    if (save_bmp("test_out.bmp", &ih, out_bitmap_data)) {
        fprintf(stderr, "main: BMP image not saved.\n");
        return 1;
    }
 
    free((pixel_t*)in_bitmap_data);
    free((pixel_t*)out_bitmap_data);
    return 0;
}
