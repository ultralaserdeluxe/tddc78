#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <pthread.h>

#include "ppm.h"
#include "blurfilter.h"
#include "gaussw.h"

#define MAX_RAD 1000
#define ROOT 0

int get_ystart_radius(int ysize, int rank, int world_size, int radius);
int get_yend_radius(int ysize, int rank, int world_size, int radius);
int get_ystart(int ysize, int rank, int world_size);
int get_yend(int ysize, int rank, int world_size);
void check_args(int argc, char** argv, int* radius);
pixel* read_file(char** argv, int* xsize, int* ysize, int* colmax);
pixel* allocate_image(int size);
void write_result(char **argv, int xsize, int ysize, int colmax, pixel* image);
void calc_stuff(int ysize, int xsize, int world_size, int radius, int* counts,
		int* offsets, int* ystarts_rad, int* yends_rad, int* offsets_rad, int* counts_rad);


int main (int argc, char ** argv) {
  /* Check arguments, read file and allocate memory. */
  int radius;
  check_args(argc, argv, &radius);

  int xsize, ysize, colmax;
  pixel* image = read_file(argv, &xsize, &ysize, &colmax);
  printf("Has read the image.\n");
  
  /* Gauss */
  double w[MAX_RAD];
  get_gauss_weights(radius, w);

  /* Filter */
  pixel* intermediary = (pixel*)malloc(sizeof(pixel) * xsize * ysize);
  blurfilter_x(xsize, ysize, image, intermediary, radius, w);
  blurfilter_y(xsize, ysize, intermediary, image, radius, w);
  free(intermediary);

  /* Write result. */    
  printf("Writing output file\n");    
  write_result(argv, xsize, ysize, colmax, image);

  /* Exit */
  return(0);
}

int get_ystart_radius(int ysize, int rank, int world_size, int radius){
  int ystart = (rank * ceil((double)ysize / world_size)) - radius;
  return rank == ROOT  ? 0 : ystart;
}

int get_yend_radius(int ysize, int rank, int world_size, int radius){
  int yend = ((rank + 1) * ceil((double)ysize / world_size)) + radius;
  return yend > ysize ? ysize : yend;
}

int get_ystart(int ysize, int rank, int world_size){
  return rank * ceil((double)ysize / world_size);
}

int get_yend(int ysize, int rank, int world_size){
  int yend = (rank + 1) * ceil((double)ysize / world_size);
  return ysize < yend ? ysize : yend;
}

void check_args(int argc, char** argv, int* radius){
  if (argc != 4) {
    fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
    exit(1);
  }
  *radius = atoi(argv[1]);
  if((*radius > MAX_RAD) || (*radius < 1)) {
    fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", *radius, MAX_RAD);
    exit(1);
  }
}

pixel* read_file(char** argv, int* xsize, int* ysize, int* colmax){
  FILE* infile;
  if (!(infile = fopen(argv[2], "r"))) {
    fprintf(stderr, "Error when opening %s\n", argv[1]);
    exit(1);
  } 

  int magic = ppm_readmagicnumber(infile);
  if (magic != 'P'*256+'6') {
    fprintf(stderr, "Wrong magic number\n");
    exit(1);
  }
  *xsize = ppm_readint(infile);
  *ysize = ppm_readint(infile);
  *colmax = ppm_readint(infile);
  if (*colmax > 255) {
    fprintf(stderr, "Too large maximum color-component value\n");
    exit(1);
  }
    
  pixel* image = allocate_image((*xsize)*(*ysize));
  if (!fread(image, sizeof(pixel), (*xsize)*(*ysize), infile)) {
    fprintf(stderr, "error in fread\n");
    exit(1);
  }

  return image;
}

void write_result(char **argv, int xsize, int ysize, int colmax, pixel* image){
  FILE* outfile;
  if (!(outfile = fopen(argv[3], "w"))) {
    fprintf(stderr, "Error when opening %s\n", argv[2]);
    free(image);
    exit(1);
  }
  fprintf(outfile, "P6 %d %d %d\n", xsize, ysize, colmax);
  if (!fwrite(image, sizeof(pixel), xsize*ysize, outfile)) {
    fprintf(stderr, "error in fwrite");
    free(image);
    exit(1);
  }
}


void calc_stuff(int ysize, int xsize, int world_size, int radius, int* counts, int* offsets, int* ystarts_rad, int* yends_rad, int* offsets_rad, int* counts_rad) {
  int ystart, yend;
  for(int i = 0; i < world_size; i++) {
    ystart = get_ystart(ysize, i, world_size);
    yend = get_yend(ysize, i, world_size);
    counts[i] = (yend - ystart) * xsize;
    offsets[i] = ystart * xsize;
    ystarts_rad[i] = get_ystart_radius(ysize, i, world_size, radius);
    yends_rad[i] = get_yend_radius(ysize, i, world_size, radius);
    offsets_rad[i] = ystarts_rad[i] * xsize;
    counts_rad[i] = (yends_rad[i] - ystarts_rad[i]) * xsize;
  }    
}

pixel* allocate_image(int size){
  pixel* image;
  image = (pixel *)malloc(sizeof(pixel)*size);
  if (!image) {
    fprintf(stderr, "malloc failed");
    exit(1);
  }
  return image;
}
