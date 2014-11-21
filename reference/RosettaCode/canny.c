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
    	printf("Total size of bitmapImage in byes is = %d\n",bitmapInfoHeader->bmp_bytesz);
    	printf("Total size of bitmapImage is=%d\n", bitmapInfoHeader->bmp_bytesz *sizeof(pixel_t));
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
    	printf("pad1 = %d\n", pad1);
    	printf("Bits per pixel=%d\n",bitmapInfoHeader->bitspp);
    	printf("Pixel per row=%d\n",bitmapInfoHeader->width);
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
    	printf("Count is = %d\n", count); 
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
	    printf("Pad1 is = %d\n", pad1);
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
 
    convolution(out, after_Gx, Gx, nx, ny, 3, false);
 
    const float Gy[] = { 1, 2, 1,
                         0, 0, 0,
                        -1,-2,-1};
 
    convolution(out, after_Gy, Gy, nx, ny, 3, false);
 
    int c = 0;
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            G[c].r = hypot(after_Gx[c].r, after_Gy[c].r);
	    G[c].g = hypot(after_Gx[c].g, after_Gy[c].g);
	    G[c].b = hypot(after_Gx[c].b, after_Gy[c].b);
	    c++;
        }
    }
 
    // Non-maximum suppression, straightforward implementation.
    /*for (int i = 0; i < ny; i++)
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
 
    // Reuse array
    // used as a stack. nx*ny/2 elements should be enough.
    int *edges = (int*) after_Gy;
    memset(out, 0, sizeof(pixel_t) * nx * ny);
    memset(edges, 0, sizeof(pixel_t) * nx * ny);
 
    // Tracing edges with hysteresis . Non-recursive implementation.
    size_t c = 1;
    for (int j = 1; j < ny - 1; j++)
        for (int i = 1; i < nx - 1; i++) {
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
            c++;
        }
 
    free(after_Gx);
    free(after_Gy);
    free(G);
    free(nms);
 */
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
