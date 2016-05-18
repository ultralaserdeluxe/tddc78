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

#define N 1000
#define MAX_ITER 1000

/* void print_T(double T[][N+2]) { */
/*   for(int y = 0; y < N+2; y++) { */
/*     printf("  "); */
/*     for(int x = 0; x < N+2; x++) { */
/*       printf("%f  ", T[y][x]); */
/*     } */
/*     printf("\n"); */
/*   } */
/* } */

void init_T(double* T){
  for(int y = 0; y < N+2; y++) {
    for(int x = 0; x < N+2; x++) {
      if(y == N+1) {
	*(T + (N+2)*y + x) = 2;
	/* T[y][x] = 2; */
      } else if(x == 0 || x == N+1) {
	*(T + (N+2)*y + x) = 1;
	/* T[y][x] = 1; */
      } else {
	*(T + (N+2)*y + x) = 0;
	/* T[y][x] = 0; */
      }
    }
  }
}

int main() {
  /* double T[N+2][N+2]; */
  double* T = (double*)malloc((N+2)*(N+2)*sizeof(double));
  init_T(T);
  double tol = pow(10, -3);
  double error;
  int stop = 0;
  /* double tmp1[N+2], tmp2[N+2]; */
  double* tmp1 = (double*)malloc((N+2)*sizeof(double));
  double* tmp2 = (double*)malloc((N+2)*sizeof(double));
  int iterations = 0;
  double max_error = 0;
  double new_error = 0;

  double starttime = omp_get_wtime();

# pragma omp parallel firstprivate(iterations, max_error, new_error) shared(tmp1, tmp2, T, error, tol, stop) num_threads(4)
  {
    while(!stop){
#     pragma omp single
      {
	memcpy(tmp1, T, (N+2)*sizeof(double));
      }

      for(int y = 1; y <= N; y++) {
#       pragma omp single
	{
	  memcpy(tmp2, T + (N+2)*y, (N+2)*sizeof(double));
	}

#       pragma omp for schedule(static, 1)
	for(int x = 1; x <= N; x++){
	  *(T + (N+2)*y + x) = (tmp1[x] + tmp2[x-1] + tmp2[x+1] + *(T + (N+2)*(y+1) + x))/4.0;
	  new_error = fabs(tmp2[x] - *(T + (N+2)*y + x));
	  max_error = MAX(new_error, max_error);
	}

#       pragma omp critical
	{
	  error = MAX(max_error, error);
	  max_error = 0.0;
	}
#       pragma omp single
	{
	  memcpy(tmp1, tmp2, (N+2)*sizeof(double));
	}
      /* printf("thread: %d iterations: %d error: %f\n", omp_get_thread_num(), iterations, error); */
      }

      iterations++;


#     pragma omp single
      {
	if(error < tol || iterations >= MAX_ITER){
	  stop = 1;
	}else{
	  error = 0.0; 
	}
      }
    }
    printf("EXIT thread: %d iterations: %d\n", omp_get_thread_num(), iterations);
  }

  double endtime = omp_get_wtime();
  printf("Time: %f\n", endtime-starttime);
  printf("The temperature of element T(5,5): %f\n", *(T + (N+2)*5 + 5));

  /* print_T(T); */

  free(tmp1);
  free(tmp2);
  free(T);

  return 0;
}


