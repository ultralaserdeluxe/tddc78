#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "coordinate.h"
#include "physics.h"
#include "definitions.h"

#define N_PARTICLES 500
#define TIME 1000

void init_particles(pcord_t particles[]){
  for(int i = 0; i < N_PARTICLES; i++){
    particles[i].x = BOX_HORIZ_SIZE*((double)rand()/(double)RAND_MAX);
    particles[i].y = BOX_VERT_SIZE*((double)rand()/(double)RAND_MAX);

    float r = MAX_INITIAL_VELOCITY*((double)rand()/(double)RAND_MAX);
    float theta = 2*PI*((double)rand()/(double)RAND_MAX);

    particles[i].vx = r*cos(theta);
    particles[i].vy = r*sin(theta);
  }
}

void print_particles(pcord_t particles[]){
  for(int i = 0; i < N_PARTICLES; i++){
    printf("i=%d x=%f y=%f vx=%f vy=%f\n", i, particles[i].x, particles[i].y, particles[i].vx, particles[i].vy);
  }
}

float calc_circumference(cord_t walls){
  return 2*(walls.x1 - walls.x0 + walls.y1 - walls.y0);
}

float calc_pressure(float momentum, cord_t walls){
  return momentum / (TIME*calc_circumference(walls));
}

int main(){
  srand(time(NULL));

  pcord_t particles[N_PARTICLES];  
  init_particles(particles);
  /* print_particles(particles); */

  cord_t walls = {0, BOX_HORIZ_SIZE, 0, BOX_VERT_SIZE};

  float total_momentum = 0;
  float collision_times[N_PARTICLES];

  for(int t = 0; t < TIME; t++){
    for(int i = 0; i < N_PARTICLES; i++){
      for(int j = 0; j < N_PARTICLES; j++){
	float collision_time = collide(&particles[i], &particles[j]);
	interact(&particles[i], &particles[j], collision_time);
	collision_times[i] = collision_times[j] = collision_time;
      }
    }

    	
    for(int i = 0; i < N_PARTICLES; i++){
      if(collision_times[i] == -1) feuler(&particles[i], t);
    }

    for(int i = 0; i < N_PARTICLES; i++){
      total_momentum += wall_collide(&particles[i], walls);
    }

    printf("t=%d ", t);
    printf("Total momentum=%f ", total_momentum);
    printf("Pressure=%f\n", calc_pressure(total_momentum, walls));
  }

  printf("Circumference=%f \n", calc_circumference(walls));

  return 0;
}
