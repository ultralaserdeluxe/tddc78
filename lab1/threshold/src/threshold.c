#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "ppm.h"

#define ROOT 0
#define PIXEL(image,x,y) ((image)+((y)*(xsize)+(x)))
typedef struct _pixel {
  unsigned char r,g,b;
} pixel;

MPI_Datatype mpi_pixel_type;

pixel* allocate_image(int size);
pixel* read_file(char** argv, int* xsize, int* ysize, int* colmax);
void check_args(int argc, char** argv);
void threshold_filter(pixel* image, int xsize, int ysize, int colmax, int world_size);
void write_result(char **argv, int xsize, int ysize, int colmax, pixel* image);
void create_mpi_pixel_type();
int get_yend(int ysize, int rank, int world_size);
int get_ystart(int ysize, int rank, int world_size);

int main(int argc, char** argv)
{
  MPI_Init(NULL, NULL);
  int rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  create_mpi_pixel_type();

  if (rank == ROOT) printf("Checking arguments.\n");
  check_args(argc, argv);

  int xsize;
  int ysize;
  int colmax;
  pixel* image = NULL;
  if (rank == ROOT){
    printf("Reading file.\n");
    image = read_file(argv, &xsize, &ysize, &colmax);
  }

  MPI_Bcast(&xsize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&ysize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&colmax, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  /* Calculate displacements and stuff for scatterv and gatherv. */
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

  printf("rank: %d\t ystart: %d\t yend: %d\t diff: %d\t count %d\t displ %d\n",
	 rank, ystarts[rank], yends[rank], yends[rank] - ystarts[rank], sendcounts[rank], displs[rank]);

  pixel* data = (pixel*)malloc(sendcounts[rank] * sizeof(pixel));
  MPI_Scatterv(image, sendcounts, displs, mpi_pixel_type, data, sendcounts[rank],
	       mpi_pixel_type, ROOT, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == ROOT) printf("Running filter.\n");
  double start_time = MPI_Wtime();

  threshold_filter(data, xsize, yends[rank]-ystarts[rank], colmax, world_size);

  MPI_Barrier(MPI_COMM_WORLD);
  double end_time = MPI_Wtime();
  if (rank == ROOT) printf("Filtering took %f seconds.\n", end_time - start_time);

  MPI_Gatherv(data, sendcounts[rank], mpi_pixel_type, image, sendcounts, displs,
	      mpi_pixel_type, ROOT, MPI_COMM_WORLD);
  free(data);

  if (rank == ROOT){
    printf("Writing output.\n");
    write_result(argv, xsize, ysize, colmax, image);
    free(image);
  }
  
  MPI_Finalize();
  exit(0);
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
  if (!(infile = fopen(argv[1], "r"))) {
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

void check_args(int argc, char** argv){
  if (argc != 3) {
    fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
    exit(2);
  }
}

void threshold_filter(pixel* image, int xsize, int ysize, int colmax, int world_size){
  int x, y;
  int sum = 0;
  for (y=0; y<ysize; y++) {
    for (x=0; x<xsize; x++) {
      sum += PIXEL(image,x,y)->r;
      sum += PIXEL(image,x,y)->g;
      sum += PIXEL(image,x,y)->b;
    }
  }
  int own_avg = sum/(xsize*ysize);
  
  int other_avg[world_size];
  MPI_Allgather(&own_avg, 1, MPI_INT, other_avg, 1, MPI_INT, MPI_COMM_WORLD);
  int avg = 0;
  int i; 
  for(i = 0; i < world_size; i++){
    avg += other_avg[i];
  }
  
  avg /= world_size;
  
  for (y=0; y<ysize; y++) {
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
}

void write_result(char **argv, int xsize, int ysize, int colmax, pixel* image){
  FILE* outfile;
  if (!(outfile = fopen(argv[2], "w"))) {
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

int get_ystart(int ysize, int rank, int world_size){
  return rank * ceil(ysize / world_size);
}

int get_yend(int ysize, int rank, int world_size){
  int yend = (rank + 1) * ceil(ysize / world_size);
  return ysize < yend ? ysize : yend;
}
