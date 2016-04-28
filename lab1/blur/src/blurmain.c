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

MPI_Datatype mpi_pixel_type;

int get_ystart(int ysize, int rank, int world_size);
int get_yend(int ysize, int rank, int world_size);
void create_mpi_pixel_type();
void check_args(int argc, char** argv, int* radius);
void read_file(char** argv, int* xsize, int* ysize, int* colmax, pixel* buf);

int main (int argc, char ** argv) {
  /* Set up MPI */
  int rank, world_size;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  create_mpi_pixel_type();

  /* Check arguments and read file */
  int radius;
  int xsize, ysize, colmax;
  pixel src[MAX_PIXELS];
  if(rank == ROOT) {
    check_args(argc, argv, &radius);
    read_file(argv, &xsize, &ysize, &colmax, src);
  }

  /* Send size data to processes */
  MPI_Bcast((void*)&xsize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast((void*)&ysize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast((void*)&radius, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  /* Allocate and send/recieve memory for own slice of the image */
  int sendcounts[world_size];
  int displs[world_size];
  int ystarts[world_size];
  int yends[world_size];

  int i;
  for(i = 0; i < world_size; i++){
    ystarts[i] = get_ystart(ysize, i, world_size);
    yends[i] = get_yend(ysize, i, world_size);
    int size = (yends[i] - ystarts[i]) * xsize;
    sendcounts[i] = size;
    displs[i] = ystarts[i]*xsize;
  }

  printf("rank: %d\t ystart: %d\t yend: %d\t diff: %d\t count %d\t displ %d\n", rank, ystarts[rank], yends[rank], yends[rank] - ystarts[rank], sendcounts[rank], displs[rank]);

  pixel* data = (pixel*)malloc(sendcounts[rank] * sizeof(pixel));
  MPI_Scatterv(&src, sendcounts, displs, mpi_pixel_type, data, sendcounts[rank], mpi_pixel_type, ROOT, MPI_COMM_WORLD);

  /* TODO: Send shit to other processes */

  /* Gauss */
  double w[MAX_RAD];
  get_gauss_weights(radius, w);

  struct timespec stime, etime;
  clock_gettime(CLOCK_REALTIME, &stime);

  /* Filter */
  blurfilter(xsize, yends[rank]-ystarts[rank], data, radius, w);

  clock_gettime(CLOCK_REALTIME, &etime);

  printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	 1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;


  MPI_Gatherv(data, sendcounts[rank], mpi_pixel_type, &src, sendcounts, displs, mpi_pixel_type, ROOT, MPI_COMM_WORLD);

  /* write result */
  if(rank == ROOT) printf("Writing output file\n");
  if(rank == ROOT && write_ppm (argv[3], xsize, ysize, (char *)src) != 0) exit(1);
  

  /* Exit */
  MPI_Finalize();
  return(0);
}

int get_ystart(int ysize, int rank, int world_size){
  return rank * ceil(ysize / world_size);
}

int get_yend(int ysize, int rank, int world_size){
  int yend = (rank + 1) * ceil(ysize / world_size);
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

 
