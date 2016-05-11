#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 10
#define MAX_ITER 60

void print_T(double T[][N+1]) {
  for(int y = 0; y < N+1; y++) {
    printf("  ");
    for(int x = 0; x < N+1; x++) {
      printf("%f  ", T[y][x]);
    }
    printf("\n");
  }
}

int main() {
  double tol = pow(10, -3);
  double T[N+1][N+1];
  double tmp1[N],tmp2[N];
  double error;
  struct timespec stime, etime;
  int x,y;

  for(y = 0; y < N+1; y++) {
    for(x = 0; x < N+1; x++) {
      if(y == N) {
	T[y][x] = 2;
      } else if(x == 0 || x == N) {
	T[y][x] = 1;
      } else {
	T[y][x] = 0;
      }
    }
  }

  print_T(T);
  printf("\n");
  
  /* /\* start time. *\/ */
  /* clock_gettime(0, &stime); */

  int k;
  for(k = 0; k < MAX_ITER; k++) {
    memcpy(tmp1, T[0], N*sizeof(double));
    error = 0;
    int i;
    for(i = 1; i < N; i++) {
      memcpy(tmp2, T[i], N*sizeof(double));

      int j;
      for(j = 1; j < N; j++){
	T[i][j] = (tmp1[j] + tmp2[j-1] + tmp2[j+1] + T[i+1][j]) / 4.0;
      }
      
      /* TODO: Fix calc error */
      int l;
      for(l = 0; l < N; l++){
	float new_error = fabs(tmp2[l] - T[i][l]);
	error = new_error > error ? new_error : error;
      }
      if(error < tol) break;


      memcpy(tmp1, tmp2, N*sizeof(double));
    }
  }

  print_T(T);
  printf("Number of Iterations: %d\n", k);
  printf("The temperature of element T(5,5): %Lf\n", T[5][5]);
  
  /* /\* end time.  *\/ */
  /* clock_gettime(0, &etime); */

  /* /\* print time. *\/ */
  /* printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) + */
  /* 	 1e-9*(etime.tv_nsec  - stime.tv_nsec)) ; */

  return 0;
}


