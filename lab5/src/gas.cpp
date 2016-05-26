#include <vector>
#include <math.h>
#include <time.h>

#include <mpi.h>

#include "coordinate.h"
#include "physics.h"
#include "definitions.h"

#define N_PARTICLES 10
#define TIME 1000

#define rand01 ((double)rand()/(double)RAND_MAX)

using namespace std;

void init_particles(vector<pcord_t>& particles, cord_t local_walls){
  for(int i = 0; i < N_PARTICLES; i++){
    pcord_t particle;
    float r = MAX_INITIAL_VELOCITY*rand01;
    float theta = 2*PI*rand01;
    particle.x = (local_walls.x1 - local_walls.x0)*rand01;
    particle.y = (local_walls.y1 - local_walls.y0)*rand01;
    particle.vx = r*cos(theta);
    particle.vy = r*sin(theta);

    particles.push_back(particle);
  }
}

void print_particle(pcord_t particle){
  cout << "x=" << particle.x << " y=" <<  particle.y << " vx=" << particle.vx << " vy=" << particle.vy << endl;
}

void print_particles(vector<pcord_t> particles){
  for(int i = 0; i < N_PARTICLES; i++){
    cout << "i=" << i << " ";
    print_particle(particles.at(i));
  }
}

float calc_circumference(cord_t walls){
  return 2*(walls.x1 - walls.x0 + walls.y1 - walls.y0);
}

float calc_pressure(float momentum, cord_t walls){
  return momentum / (TIME*calc_circumference(walls));
}

int main(int argc, char** argv){
  srand(time(NULL));

  MPI::Init(argc, argv);
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();

  cord_t global_walls = {0, BOX_HORIZ_SIZE, 0, BOX_VERT_SIZE};

  vector<pcord_t> particles;
  init_particles(particles, global_walls);

  print_particles(particles);


  
  float local_momentum = 0;

  /* for(int t = 0; t < TIME; t += STEP_SIZE){ */
  /*   for(int i = 0; i < N_PARTICLES-1; i++){ */
  /*     float collision_time; */
  /*     for(int j = i+1; j < N_PARTICLES; j++){ */
  /* 	collision_time = collide(&particles[i], &particles[j]); */
  /* 	interact(&particles[i], &particles[j], collision_time); */
  /*     } */
  /*     if(collision_time == -1) feuler(&particles[i], t); */
  /*     total_momentum += wall_collide(&particles[i], walls); */
  /*   } */

  /*   MPI::COMM_WORLD::Barrier(); */
  /* } */
  MPI::Finalize();
  return 0;
}
