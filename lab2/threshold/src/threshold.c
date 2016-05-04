#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

#include "ppm.h"

#define ROOT 0
#define MAX_RAD 1000

#define PIXEL(image,x,y) ((image)+((y)*(xsize)+(x)))
typedef struct _pixel {
  unsigned char r,g,b;
} pixel;

typedef struct thread_data_t{
  int tid;
  int* n_workers;
  int* xsize;
  int* ysize;
  int* colmax;
  int own_avg;
  pixel* image;
} thread_data;

int clock_gettime(int clk_id, struct timespec* t);
pixel* allocate_image(int size);
pixel* read_file(char** argv, int* xsize, int* ysize, int* colmax);
void check_args(int argc, char** argv, int* n_workers);
void* threshold_filter(void* param);
void* threshold_average(void* param);
void write_result(char **argv, int xsize, int ysize, int colmax, pixel* image);
int get_yend(int ysize, int rank, int world_size);
int get_ystart(int ysize, int rank, int world_size);

int main(int argc, char** argv)
{
  int n_workers;
  check_args(argc, argv, &n_workers);
  
  int xsize, ysize, colmax;
  pixel* image = read_file(argv, &xsize, &ysize, &colmax);
  printf("Has read the image.\n");

  pthread_t workers[n_workers];
  thread_data t_data[n_workers];

  int i;
  for(i = 0; i < n_workers; i++){
    t_data[i].tid = i;
    t_data[i].n_workers = &n_workers;
    t_data[i].xsize = &xsize;
    t_data[i].ysize = &ysize;
    t_data[i].colmax = &colmax;
    t_data[i].image = image;
  }
  
  /* printf("rank: %d\t ystart: %d\t yend: %d\t diff: %d\t count %d\t displ %d\n", */
  /* 	 rank, ystarts[rank], yends[rank], yends[rank] - ystarts[rank], sendcounts[rank], displs[rank]); */

  /* start timing */
  struct timespec start_time;
  clock_gettime(0, &start_time);
  
  /* calculate average. */
  for(i = 0; i < n_workers; i++){
    pthread_create(&workers[i], NULL, threshold_average, &t_data[i]);
  }

  for(i = 0; i < n_workers; i++){
    pthread_join(workers[i], NULL);
  }

  /* calculate average. */
  int avg = 0;
  for(i = 0; i < n_workers; i++){
    avg += t_data[i].own_avg;
  }
  
  avg /= n_workers;
  printf("average: %d\n", avg);

  
  /* filter */
  printf("Running filter.\n");
  for(i = 0; i < n_workers; i++){
    t_data[i].own_avg = avg;
    pthread_create(&workers[i], NULL, threshold_filter, &t_data[i]);
  }
  
  for(i = 0; i < n_workers; i++){
    pthread_join(workers[i], NULL);
  }

  /* stop timing */
  struct timespec end_time;
  clock_gettime(0, &end_time);
  /* printf("Filtering took %f seconds.\n", end_time - start_time); */
  
  printf("Writing output.\n");
  write_result(argv, xsize, ysize, colmax, image);
  free(image);
  
  exit(0);
}

int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}

pixel* allocate_image(int size)
{
  pixel* image;
  image = (pixel *)malloc(sizeof(pixel)*size);
  if (!image) {
    fprintf(stderr, "malloc failed");
    exit(1);
  }
  return image;
}

pixel* read_file(char** argv, int* xsize, int* ysize, int* colmax){
  FILE* infile;
  if (!(infile = fopen(argv[2], "r"))) {
    fprintf(stderr, "Error when opening %s\n", argv[2]);
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

void check_args(int argc, char** argv, int* n_workers){
  if (argc != 4) {
    fprintf(stderr, "Usage: %s workers infile outfile\n", argv[0]);
    exit(1);
  }
  *n_workers = atoi(argv[1]);
  if(*n_workers <= 0){
    fprintf(stderr, "Workers (%d) must be greater than one\n", *n_workers);
    exit(1);
  }
}

void* threshold_average(void* param) {
  thread_data* t_data = (thread_data*)param;
  pixel* image = t_data->image;
  int xsize = *(t_data->xsize);
  int ysize = *(t_data->ysize);
  int n_workers = *(t_data->n_workers);  
  int ystart = get_ystart(ysize, t_data->tid, n_workers);
  int yend = get_yend(ysize, t_data->tid, n_workers);

  int x, y;
  int sum = 0;
  for (y=ystart; y<yend; y++) {
    for (x=0; x<xsize; x++) {
      sum += PIXEL(image,x,y)->r;
      sum += PIXEL(image,x,y)->g;
      sum += PIXEL(image,x,y)->b;
    }
  }
  t_data->own_avg = sum/(xsize*(yend-ystart));

  return NULL;
}
  
void* threshold_filter(void* param) {
  thread_data* t_data = (thread_data*)param;
  pixel* image = t_data->image;
  int xsize = *(t_data->xsize);
  int ysize = *(t_data->ysize);
  int ystart = get_ystart(ysize, t_data->tid, *(t_data->n_workers));
  int yend = get_yend(ysize, t_data->tid, *(t_data->n_workers));
  int colmax = *(t_data->colmax);
  int avg = t_data->own_avg;

  int x, y;
  for (y=ystart; y<yend; y++) {
    for (x=0; x<xsize; x++) {
      int psum = PIXEL(image,x,y)->r;
      psum += PIXEL(image,x,y)->g;
      psum += PIXEL(image,x,y)->b;
      int pval;
      if (psum > avg)
	pval = colmax;
      else
	pval = 0;
      PIXEL(image,x,y)->r = pval;
      PIXEL(image,x,y)->g = pval;
      PIXEL(image,x,y)->b = pval;
    }
  }
  
  return NULL;
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

int get_ystart(int ysize, int rank, int world_size){
  return rank * ceil(ysize / world_size);
}

int get_yend(int ysize, int rank, int world_size){
  int yend = (rank + 1) * ceil(ysize / world_size);
  return ysize < yend ? ysize : yend;
}
