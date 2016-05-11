#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

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

int main() {
  double T[N+2][N+2];

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

  double tol = pow(10, -3);
  double tmp1[N+2],tmp2[N+2];

  int k;
  for(k = 1; k <= MAX_ITER; k++) {
    memcpy(tmp1, T[0], (N+2)*sizeof(double));
    double error = 0.0;

    for(int y = 1; y <= N; y++) {
      memcpy(tmp2, T[y], (N+2)*sizeof(double));
      
      for(int x = 1; x <= N; x++){
	T[y][x] = (tmp1[x] + tmp2[x-1] + tmp2[x+1] + T[y+1][x]) / 4.0;
      }
      
      double max_error = 0;
      for(int x = 1; x <= N; x++){
	double new_error = fabs(tmp2[x] - T[y][x]);
	max_error = new_error > max_error ? new_error : max_error;
      }
      error = max_error > error ? max_error : error;

      memcpy(tmp1, tmp2, (N+2)*sizeof(double));
    }

    if(error < tol) {
      break;
    }
  }

  printf("Number of Iterations: %d\n", k);
  printf("The temperature of element T(5,5): %f\n", T[5][5]);

  print_T(T);

  return 0;
}


