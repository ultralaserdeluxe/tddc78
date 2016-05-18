#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define MAX(A,B) (A > B ? A : B)

#define N 2000
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

int get_ystart(int ysize, int rank, int world_size){
  return rank * ceil((double)ysize / world_size);
}


int main() {
  double* T = (double*)malloc((N+2)*(N+2)*sizeof(double));
  init_T(T);
  double tol = pow(10, -3);
  double global_error;
  double starttime = omp_get_wtime();

# pragma omp parallel
  {
    double* tmp1 = (double*)malloc((N+2)*sizeof(double));
    double* tmp2 = (double*)malloc((N+2)*sizeof(double));
    double* tmp3 = (double*)malloc((N+2)*sizeof(double));

    int ystart = get_ystart(N, omp_get_thread_num(), omp_get_num_threads());
    int yend = get_yend(N, omp_get_thread_num(), omp_get_num_threads());

    int iterations = 0;
    double local_error = 0;

    do{
#     pragma omp barrier
      global_error = 0.0;
      local_error = 0.0;
      
      memcpy(tmp1, T + (N+2)*ystart, (N+2)*sizeof(double));
      memcpy(tmp3, T + (N+2)*(yend+1), (N+2)*sizeof(double));

#     pragma omp for schedule(static, (int)ceil((double)N/omp_get_num_threads()))
      for(int y = 1; y <= N; y++) {
	memcpy(tmp2, T + (N+2)*y, (N+2)*sizeof(double));

	for(int x = 1; x <= N; x++){
	  double below;
	  if(y == yend) below = tmp3[x];
	  else below = T[(N+2)*(y+1) + x];

	  T[(N+2)*y + x] = (tmp1[x] + tmp2[x-1] + tmp2[x+1] + below)/4.0;
	  local_error = MAX(fabs(tmp2[x] - T[(N+2)*y + x]), local_error);
	}

	double* swap = tmp1;
	tmp1 = tmp2;
	tmp2 = swap;
      }

      iterations++;

#     pragma omp critical
      {
	global_error = MAX(local_error, global_error);
      }
#     pragma omp barrier
    }while(global_error >= tol && iterations < MAX_ITER);
    
    printf("EXIT thread: %d iterations: %d\n", omp_get_thread_num(), iterations);

    free(tmp1);
    free(tmp2);
    free(tmp3);
  }

  double endtime = omp_get_wtime();
  printf("Time: %f\n", endtime-starttime);
  printf("The temperature of element T(5,5): %f\n", T[(N+2)*5 + 5]);

  /* print_T(T); */

  free(T);

  return 0;
}


