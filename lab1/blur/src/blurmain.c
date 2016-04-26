#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"

#include "mactime.h"

#include <mpi.h>

#define MAX_RAD 1000
#define ROOT 0

int main (int argc, char ** argv) {
   int radius;
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    double w[MAX_RAD];
    struct timespec stime, etime;
    

    /* Set up MPI */
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int nitems = 3;
    int blocklengths[3] = {1,1,1};
    MPI_Datatype types[3] = {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};
    MPI_Aint offsets[3];

    offsets[0] = offsetof(pixel, r);
    offsets[1] = offsetof(pixel, g);
    offsets[2] = offsetof(pixel, b);

    MPI_Datatype mpi_pixel_type;
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_pixel_type);
    MPI_Type_commit(&mpi_pixel_type);
  
  

    
    /* Take care of the arguments */

    if(rank == ROOT) {
      if (argc != 4) {
	fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
	exit(1);
      }
      radius = atoi(argv[1]);
      if((radius > MAX_RAD) || (radius < 1)) {
	fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
	exit(1);
      }

      /* read file */
      if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

      if (colmax > 255) {
	fprintf(stderr, "Too large maximum color-component value\n");
	exit(1);
      }
    
    }

    MPI_Bcast((void*)&xsize,
	      1,
	      MPI_INT,
	      ROOT,
	      MPI_COMM_WORLD);

    MPI_Bcast((void*)&ysize,
	      1,
	      MPI_INT,
	      ROOT,
	      MPI_COMM_WORLD);

    MPI_Bcast((void*)&radius,
	      1,
	      MPI_INT,
	      ROOT,
	      MPI_COMM_WORLD);

    MPI_Bcast((void*)src,
	      ysize*xsize,
	      mpi_pixel_type,
	      ROOT,
	      MPI_COMM_WORLD);

    printf("Done waiting for Bcast\n");
    

    printf("Has read the image, generating coefficients\n");

    /* filter */
    get_gauss_weights(radius, w);

    printf("Calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);

    blurfilter(xsize, ysize, src, radius, w);

    clock_gettime(CLOCK_REALTIME, &etime);

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(rank == ROOT && write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
      exit(1);

    MPI_Finalize();

    return(0);
}
