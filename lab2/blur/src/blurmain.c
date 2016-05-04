#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <pthread.h>

#include "ppm.h"
#include "blurfilter.h"
#include "gaussw.h"

#define MAX_RAD 1000

typedef struct thread_data_t{
  int tid;
  int* n_workers;
  int* xsize;
  int* ysize;
  int* radius;
  int* colmax;
  pixel* image;
  pixel* intermediary;
  double* w;
} thread_data;

int clock_gettime(int clk_id, struct timespec* t);
int get_ystart(int ysize, int rank, int world_size);
int get_yend(int ysize, int rank, int world_size);
void check_args(int argc, char** argv, int* radius, int* n_workers);
pixel* read_file(char** argv, int* xsize, int* ysize, int* colmax);
pixel* allocate_image(int size);
void write_result(char **argv, int xsize, int ysize, int colmax, pixel* image);

void* thread_blur_x(void* param){
  thread_data* t_data = (thread_data*)param;
  
  int xsize = *(t_data->xsize);
  int ysize = *(t_data->ysize);
  int radius = *(t_data->radius);
  int n_workers = *(t_data->n_workers);

  int ystart = get_ystart(ysize, t_data->tid, n_workers);
  int yend = get_yend(ysize, t_data->tid, n_workers);

  /* printf("Running X filter in thread %d.\n", t_data->tid); */
  blurfilter_x(ystart, yend, xsize, ysize, t_data->image, t_data->intermediary, radius, t_data->w);

  return NULL;
}

void* thread_blur_y(void* param){
  thread_data* t_data = (thread_data*)param;
  
  int xsize = *(t_data->xsize);
  int ysize = *(t_data->ysize);
  int radius = *(t_data->radius);
  int n_workers = *(t_data->n_workers);

  int ystart = get_ystart(ysize, t_data->tid, n_workers);
  int yend = get_yend(ysize, t_data->tid, n_workers);

  /* printf("Running Y filter in thread %d.\n", t_data->tid); */
  blurfilter_y(ystart, yend, xsize, ysize, t_data->intermediary, t_data->image, radius, t_data->w);

  return NULL;
}

int main (int argc, char ** argv) {
  /* Check arguments, read file and allocate memory. */
  int radius;
  int n_workers;
  check_args(argc, argv, &radius, &n_workers);

  int xsize, ysize, colmax;
  pixel* image = read_file(argv, &xsize, &ysize, &colmax);
  printf("Has read the image.\n");
  
  /* Gauss */
  double w[MAX_RAD];
  get_gauss_weights(radius, w);

  /* Create threads */
  pthread_t workers[n_workers];
  thread_data t_data[n_workers];
  
  pixel* intermediary = allocate_image(xsize * ysize);
  
  int i;
  for(i = 0; i < n_workers; i++){
    t_data[i].tid = i;
    t_data[i].n_workers = &n_workers;
    t_data[i].xsize = &xsize;
    t_data[i].ysize = &ysize;
    t_data[i].radius = &radius;
    t_data[i].colmax = &colmax;
    t_data[i].w = w;
    t_data[i].image = image;
    t_data[i].intermediary = intermediary;
  }

  /* Start measure time*/
  struct timespec stime;
  clock_gettime(0, &stime);
  
  for(i = 0; i < n_workers; i++){
    pthread_create(&workers[i], NULL, thread_blur_x, &t_data[i]);
  }

  for(i = 0; i < n_workers; i++){
    pthread_join(workers[i], NULL);
  }

  for(i = 0; i < n_workers; i++){
    pthread_create(&workers[i], NULL, thread_blur_y, &t_data[i]);
  }

  for(i = 0; i < n_workers; i++){
    pthread_join(workers[i], NULL);
  }
  
  /* End measure time */
  struct timespec etime;
  clock_gettime(0, &etime);

  free(intermediary);

  /* Write result. */    
  printf("Writing output file\n");    
  write_result(argv, xsize, ysize, colmax, image);
  free(image);

  printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	 1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

  /* Exit */
  return(0);
}

int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}

int get_ystart(int ysize, int rank, int world_size){
  return rank * ceil((double)ysize / world_size);
}

int get_yend(int ysize, int rank, int world_size){
  int yend = (rank + 1) * ceil((double)ysize / world_size);
  return ysize < yend ? ysize : yend;
}

void check_args(int argc, char** argv, int* radius, int* n_workers){
  if (argc != 5) {
    fprintf(stderr, "Usage: %s workers radius infile outfile\n", argv[0]);
    exit(1);
  }
  *radius = atoi(argv[2]);
  if((*radius > MAX_RAD) || (*radius < 1)) {
    fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", *radius, MAX_RAD);
    exit(1);
  }
  *n_workers = atoi(argv[1]);
  if(*n_workers <= 0){
    fprintf(stderr, "Workers (%d) must be greater than one\n", *n_workers);
    exit(1);
  }
}

pixel* read_file(char** argv, int* xsize, int* ysize, int* colmax){
  FILE* infile;
  if (!(infile = fopen(argv[3], "r"))) {
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
  if (!(outfile = fopen(argv[4], "w"))) {
    fprintf(stderr, "Error when opening %s\n", argv[4]);
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

pixel* allocate_image(int size){
  pixel* image;
  image = (pixel *)malloc(sizeof(pixel)*size);
  if (!image) {
    fprintf(stderr, "malloc failed");
    exit(1);
  }
  return image;
}
