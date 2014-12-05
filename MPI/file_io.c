#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
 
#include "file_io.h"

#define MAX_BRIGHTNESS 255
 
// C99 doesn't define M_PI (GNU-C99 does)
#define M_PI 3.14159265358979323846264338327

static double sRGB_to_linear(double x) {
    if(x < 0.04045) return x/12.92;
    return pow((x+0.055)/1.055, 2.4);
}

static double GrayScaleValue(double R, double G, double B){
    double gray_linear;
    double R_linear, G_linear, B_linear;

    R_linear = sRGB_to_linear(R/255.0);
    G_linear = sRGB_to_linear(G/255.0);
    B_linear = sRGB_to_linear(B/255.0);
    gray_linear = 0.299 * R_linear + 0.587 * G_linear + 0.114 * B_linear;

    return gray_linear;
    
}

static double linear_to_sRGB(double y){
    if (y <= 0.0031308) return 12.92 * y;
    return 1.055 * pow(y, 1/2.4) - 0.055;
}

static double GrayColor(double gray_linear){
    return round(linear_to_sRGB(gray_linear) * 255);
}

static int convert_pixel_to_grayscale(int R, int G, int B){
    
    double gray_linear;
    double gray_color;
    gray_linear = GrayScaleValue((double) R, (double) G, (double) B);  
    gray_color = GrayColor(gray_linear);

    return (int) gray_color;
}

pixel_t *load_bmp(const char *filename,
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
 
    // move file point to the beginning of bitmap data
    if (fseek(filePtr, bitmapFileHeader.bmp_offset, SEEK_SET)) {
        fclose(filePtr);
        return NULL;
    }
 
    // allocate enough memory for the bitmap image data
    pixel_t *bitmapImage = malloc(bitmapInfoHeader->width * bitmapInfoHeader->height * sizeof(pixel_t));
 
    // verify memory allocation
    if (bitmapImage == NULL) {
        fclose(filePtr);
        return NULL;
    }
 
    // read in the bitmap image data
    //size_t pad;
    size_t pad1, count=0;
    //unsigned char c;
    unsigned char c[3];
    
    //pad = 4*ceil(bitmapInfoHeader->bitspp*bitmapInfoHeader->width/32.) - bitmapInfoHeader->width;
    pad1 = 4*floor((bitmapInfoHeader->bitspp*bitmapInfoHeader->width+31)/32) - (bitmapInfoHeader->width * 3);
    for(size_t i=0; i<bitmapInfoHeader->height; i++){
        for(size_t j=0; j<bitmapInfoHeader->width; j++){
            if (fread(c, sizeof(unsigned char) * 3, 1, filePtr) != 1) {
                fclose(filePtr);
                return NULL;
            }
            bitmapImage[count++] = (short int) (convert_pixel_to_grayscale(c[2], c[1], c[0]));
        }
        fseek(filePtr, pad1, SEEK_CUR);
    }
    printf("Total bmp_bytesz=%d from header\n",bitmapInfoHeader->bmp_bytesz);
    printf("Total allocation=%d from header\n",bitmapInfoHeader->bmp_bytesz * sizeof(pixel_t));
    printf("Total bytesz=%d width*height=%d dcount=%d\n", bitmapInfoHeader->height * bitmapInfoHeader->width * 3, 
        bitmapInfoHeader->height * bitmapInfoHeader->width, count);
 
    // If we were using unsigned char as pixel_t, then:
    // fread(bitmapImage, 1, bitmapInfoHeader->bmp_bytesz, filePtr);
 
    // close file and return bitmap image data
    fclose(filePtr);
    return bitmapImage;
}

bool save_bmp(const char *filename, const bitmap_info_header_t *bmp_ih,
              const pixel_t *data)
{
    FILE* filePtr = fopen(filename, "wb");
    if (filePtr == NULL)
        return true;
 
    bmpfile_magic_t mag = {{0x42, 0x4d}};
    if (fwrite(&mag, sizeof(bmpfile_magic_t), 1, filePtr) != 1) {
        fclose(filePtr);
        return true;
    }
 
    const uint32_t offset = sizeof(bmpfile_magic_t) +
                            sizeof(bmpfile_header_t) +
                            sizeof(bitmap_info_header_t) +
                            ((1U << bmp_ih->bitspp) * 4);
 
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
 
    // Palette
    for (size_t i = 0; i < (1U << bmp_ih->bitspp); i++) {
        const rgb_t color = {(uint8_t)i, (uint8_t)i, (uint8_t)i};
        if (fwrite(&color, sizeof(rgb_t), 1, filePtr) != 1) {
            fclose(filePtr);
            return true;
        }
    }
 
    // We use int instead of uchar, so we can't write img
    // in 1 call any more.
    // fwrite(data, 1, bmp_ih->bmp_bytesz, filePtr);
 
    // Padding: http://en.wikipedia.org/wiki/BMP_file_format#Pixel_storage
    // size_t pad = 4*ceil(bmp_ih->bitspp*bmp_ih->width/32.) - bmp_ih->width;
    size_t pad1 = 4*floor((bmp_ih->bitspp*bmp_ih->width+31)/32) - (bmp_ih->width * 3);
    unsigned char c[3];
    for(size_t i=0; i < bmp_ih->height; i++) {
	    for(size_t j=0; j < bmp_ih->width; j++) {
		    c[0] = (unsigned char) data[j + bmp_ih->width*i];
		    c[1] = (unsigned char) data[j + bmp_ih->width*i];
		    c[2] = (unsigned char) data[j + bmp_ih->width*i];
		    if (fwrite(c, sizeof(char) * 3, 1, filePtr) != 1) {
			    fclose(filePtr);
			    return true;
		    }
	    }
	    c[0] = 0;
	    c[1] = 0;
	    c[2] = 0;
	    for(size_t j=0; j<pad1; j++)
		    if (fwrite(c, sizeof(char) * 3, 1, filePtr) != 1) {
			    fclose(filePtr);
			    return true;
		    }
    }
 
    fclose(filePtr);
    return false;
}
