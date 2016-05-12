#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#ifdef __MACH__
#include <sys/time.h>

#define CLOCK_REALTIME 1234

int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}

#endif /* __MACH__ */

#define MAX(A,B) (A > B ? A : B)

#define N 10
#define MAX_ITER 1000

void print_T(double T[][N+2]) {
  for(int y = 0; y < N+2; y++) {
    printf("  ");
    for(int x = 0; x < N+2; x++) {
      printf("%f  ", T[y][x]);
    }
    printf("\n");
  }
}

void init_T(double T[][N+2]){
  for(int y = 0; y < N+2; y++) {
    for(int x = 0; x < N+2; x++) {
      if(y == N+1) {
	T[y][x] = 2;
      } else if(x == 0 || x == N+1) {
	T[y][x] = 1;
      } else {
	T[y][x] = 0;
      }
    }
  }
}

double calc_error(double tmp2[], double T[][N+2], int y, double old_error){
  double max_error = 0;
  for(int x = 1; x <= N; x++){
    double new_error = fabs(tmp2[x] - T[y][x]);
    max_error = MAX(new_error, max_error);
  }
  return MAX(max_error, old_error);
}

int main() {
  double T[N+2][N+2];
  init_T(T);
  double tol = pow(10, -3);
  double global_error = 0.0;

  struct timespec start;
  clock_gettime(CLOCK_REALTIME, &start);


  int k=0;
#pragma omp parallel private(k) num_threads(2)
  {
    double tmp1[N+2], tmp2[N+2];

    for(k = 1; k <= MAX_ITER; k++) {
      memcpy(tmp1, T[0], (N+2)*sizeof(double));
      double local_error = 0.0;

#pragma omp for
      for(int y = 1; y <= N; y++) {
	memcpy(tmp2, T[y], (N+2)*sizeof(double));
      
	for(int x = 1; x <= N; x++){
	  T[y][x] = (tmp1[x] + tmp2[x-1] + tmp2[x+1] + T[y+1][x])/4.0;
	}

	local_error = calc_error(tmp2, T, y, local_error);

#pragma omp critical
	{
	  global_error = MAX(local_error, global_error);
	}

	memcpy(tmp1, tmp2, (N+2)*sizeof(double));
      }

#pragma omp barrier
#pragma omp single
      {
	if(global_error < tol){
	  k = MAX_ITER;
	}
      }
    }
  }

  struct timespec end;
  clock_gettime(CLOCK_REALTIME, &end);

  printf("Time: %f Number of Iterations: %d\n", (float)(end.tv_sec - start.tv_sec), k);
  printf("The temperature of element T(5,5): %f\n", T[5][5]);

  /* print_T(T); */

  return 0;
}


