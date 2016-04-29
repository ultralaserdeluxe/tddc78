#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <mpi.h>

#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include "mactime.h"

#define MAX_RAD 1000
#define ROOT 0
#define MAC_MEM 700000

MPI_Datatype mpi_pixel_type;
MPI_Status status;

int get_ystart_radius(int ysize, int rank, int world_size, int radius);
int get_yend_radius(int ysize, int rank, int world_size, int radius);
int get_ystart(int ysize, int rank, int world_size);
int get_yend(int ysize, int rank, int world_size);
void create_mpi_pixel_type();
void check_args(int argc, char** argv, int* radius);
void read_file(char** argv, int* xsize, int* ysize, int* colmax, pixel* buf);
void calc_stuff(int ysize, int xsize, int world_size, int radius, int* counts, int* offsets, int* ystarts_rad, int* yends_rad, int* offsets_rad, int* counts_rad);


int main (int argc, char ** argv) {
  int rank, world_size, xsize, ysize, colmax, radius;;
  pixel src[MAX_PIXELS];
  pixel recvbuff[MAC_MEM];
  /* pixel* src; */
  /* pixel* recvbuff;  */
  struct timespec stime, etime;

  /* Set up MPI */
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  create_mpi_pixel_type();

  int ystarts_rad[world_size];
  int yends_rad[world_size];
  int counts[world_size];
  int counts_rad[world_size];
  int offsets[world_size];
  int offsets_rad[world_size];
  
  /* Check arguments and read file */
  if(rank == ROOT) {

    /* src = malloc(sizeof(double) * MAX_PIXELS); */
    /* /\* printf("sizeof(double) * 3: %d\n", sizeof(double * 3)); *\/ */
    
    /* if(src == NULL) { */
    /*   printf("Process %d could not allocate memory, exiting.\n", rank); */
    /*   exit(1); */
    /* } */
    
    check_args(argc, argv, &radius);
    read_file(argv, &xsize, &ysize, &colmax, src);
    printf("Has read the image, generating coefficients\n");
  } 
  
  /* Send size data to processes */
  MPI_Bcast((void*)&xsize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast((void*)&ysize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast((void*)&radius, 1, MPI_INT, ROOT, MPI_COMM_WORLD);  

  /* calculate some stuff */
  calc_stuff(ysize, xsize, world_size, radius, counts, offsets, ystarts_rad, yends_rad, offsets_rad, counts_rad);

  /* if(rank != ROOT) { */

  /*   recvbuff = malloc(sizeof(pixel) * MAX_PIXELS); //counts_rad[rank]); */
  /*   if(src == NULL) { */
  /*     printf("Process %d could not allocate memory, exiting.\n", rank); */
  /*     exit(1); */
  /*   } */
  /* } */
  
  /* Gauss */
  double w[MAX_RAD];
  get_gauss_weights(radius, w);

  /* wait for every task, before getting the time. */
  MPI_Barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_REALTIME, &stime);

  /* send parts to other tasks, and filter */
  if(rank == ROOT) {
    for(int i = 1; i < world_size; i++) {
      MPI_Send((void*)&src[offsets_rad[i]], counts_rad[i], mpi_pixel_type, i, 0, MPI_COMM_WORLD);
    }
    /* printf("Calling filter\n"); */
    printf("Overriding protocols.\n");
    blurfilter(xsize, yends_rad[rank],  src, radius, w);
  } else {
    MPI_Recv((void*)recvbuff, counts_rad[rank], mpi_pixel_type, ROOT, 0, MPI_COMM_WORLD, &status);
    blurfilter(xsize, yends_rad[rank], recvbuff, radius, w);
  }

  /* wait for every task, before getting the time. */
  MPI_Barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_REALTIME, &etime);
  if(rank == ROOT) printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
  	 1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;
  
  /* send stuff back to main task */
  if(rank == ROOT) {
    for(int i = 1; i < world_size; i++) {
      MPI_Recv((void*)&src[offsets[i]], counts[i], mpi_pixel_type, i, 0, MPI_COMM_WORLD, &status);
    }
  } else {
    MPI_Send((void*)&recvbuff[radius*xsize], counts[rank], mpi_pixel_type, ROOT, 0, MPI_COMM_WORLD);
  }
  
  MPI_Finalize();

  if(rank == ROOT) {
    /* write result */    
    printf("Writing output file\n");    
    if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
      exit(1);
  }

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

void create_mpi_pixel_type(){
  const int nitems = 3;
  int blocklengths[3] = {1,1,1};
  MPI_Datatype types[3] = {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};
  MPI_Aint offsets[3];

  offsets[0] = offsetof(pixel, r);
  offsets[1] = offsetof(pixel, g);
  offsets[2] = offsetof(pixel, b);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_pixel_type);
  MPI_Type_commit(&mpi_pixel_type);
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

void read_file(char** argv, int* xsize, int* ysize, int* colmax, pixel* buf){
  if(read_ppm (argv[2], xsize, ysize, colmax, (char *) buf) != 0)
    exit(1);

  if (*colmax > 255) {
    fprintf(stderr, "Too large maximum color-component value\n");
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


 
