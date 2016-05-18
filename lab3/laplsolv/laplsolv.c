#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>


#define MAX(A,B) (A > B ? A : B)
#define N 1000
#define MAX_ITER 1000

void print_T(double* T) {
  for(int y = 0; y < N+2; y++) {
    printf("  ");
    for(int x = 0; x < N+2; x++) {
      printf("%f  ", T[(N+2)*y + x]);
    }
    printf("\n");
  }
}

void init_T(double* T){
  for(int y = 0; y < N+2; y++) {
    for(int x = 0; x < N+2; x++) {
      if(y == N+1) {
	T[(N+2)*y + x] = 2;
      } else if(x == 0 || x == N+1) {
	T[(N+2)*y + x] = 1;
      } else {
	T[(N+2)*y + x] = 0;
      }
    }
  }
}

int get_yend(int ysize, int rank, int world_size){
  int yend = (rank + 1) * ceil((double)ysize / world_size);
  return ysize < yend ? ysize : yend;
}

int main() {
  double* T = (double*)malloc((N+2)*(N+2)*sizeof(double));
  init_T(T);
  double tol = pow(10, -3);
  double error;
  int stop = 0;
  int iterations = 0;
  double max_error = 0;
  double new_error = 0;

  double starttime = omp_get_wtime();

# pragma omp parallel					\
  firstprivate(iterations, max_error, new_error)	\
  shared(T, error, tol, stop)				\
  num_threads(7)
  {
    double* tmp1 = (double*)malloc((N+2)*sizeof(double));
    double* tmp2 = (double*)malloc((N+2)*sizeof(double));
    double* tmp3 = (double*)malloc((N+2)*sizeof(double));

    int yend = get_yend(N, omp_get_thread_num(), omp_get_num_threads());
    int ystart = get_yend(N, omp_get_thread_num()-1, omp_get_num_threads()) + 1;

    while(!stop){
      if(omp_get_thread_num() == 0){
       	memcpy(tmp1, T, (N+2)*sizeof(double));
      }else{	
	memcpy(tmp1, T + (N+2)*yend, (N+2)*sizeof(double));
      }

      int yend = get_yend(N, omp_get_thread_num(), omp_get_num_threads());
      memcpy(tmp3, T + (N+2)*(yend+1), (N+2)*sizeof(double));

#     pragma omp barrier
#     pragma omp for schedule(static, (int)ceil((double)N/omp_get_num_threads()))
      for(int y = 1; y <= N; y++) {
	memcpy(tmp2, T + (N+2)*y, (N+2)*sizeof(double));
	for(int x = 1; x <= N; x++){
	  if(y == yend){
	    T[(N+2)*y + x] = (tmp1[x] + tmp2[x-1] + tmp2[x+1] + tmp3[x])/4.0;
	  }else{
	    T[(N+2)*y + x] = (tmp1[x] + tmp2[x-1] + tmp2[x+1] + T[(N+2)*(y+1) + x])/4.0;
	  }
	  new_error = fabs(tmp2[x] - T[ + (N+2)*y + x]);
	  max_error = MAX(new_error, max_error);
	}

	memcpy(tmp1, tmp2, (N+2)*sizeof(double));
      }

      iterations++;
#     pragma omp critical
      {
	error = MAX(max_error, error);
	max_error = 0.0;
      }
#     pragma omp barrier

#     pragma omp single
      {
	if(error < tol || iterations >= MAX_ITER){
	  stop = 1;
	}else{
	  error = 0.0; 
	}
      }
#     pragma omp barrier
    }
    printf("EXIT thread: %d iterations: %d\n", omp_get_thread_num(), iterations);

    free(tmp1);
    free(tmp2);
    free(tmp3);
  }

  double endtime = omp_get_wtime();
  printf("Time: %f\n", endtime-starttime);
  printf("The temperature of element T(5,5): %f\n", *(T + (N+2)*5 + 5));

  /* print_T(T); */

  free(T);

  return 0;
}


