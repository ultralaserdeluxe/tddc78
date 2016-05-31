#include <vector>
#include <math.h>
#include <time.h>

#include <mpi.h>

#include "coordinate.h"
#include "physics.h"
#include "definitions.h"

#define N_PARTICLES 1000
#define TIME 1000

#define rand01 ((double)rand()/(double)RAND_MAX)

using namespace std;

MPI::Datatype mpi_pcord_t;

void init_particles(vector<pcord_t>& particles, int n_particles, cord_t local_walls){
  srand(time(NULL));

  for(int i = 0; i < n_particles; i++){
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

int get_ystart(int ysize, int rank, int world_size){
  return rank * ceil((double)ysize / world_size);
}

int get_yend(int ysize, int rank, int world_size){
  int yend = (rank + 1) * ceil((double)ysize / world_size);
  return ysize < yend ? ysize : yend;
}

cord_t get_local_walls(cord_t global_walls, int rank, int world_size){
  int ysize = global_walls.y1 - global_walls.y0;
  int ystart = get_ystart(ysize, rank, world_size);
  int yend = get_yend(ysize, rank, world_size);
  cord_t walls = {global_walls.x0, global_walls.x1, ystart, yend};
  return walls;
}

int main(int argc, char** argv){
  MPI::Init(argc, argv);
  int rank = MPI::COMM_WORLD.Get_rank();
  int world_size = MPI::COMM_WORLD.Get_size();

  cord_t global_walls = {0, 10, 0, 1000000};
  cord_t local_walls = get_local_walls(global_walls, rank, world_size);

  int local_particles = N_PARTICLES/world_size;
  vector<pcord_t> particles;
  init_particles(particles, local_particles, local_walls);

#ifdef DEBUG
  if(rank == 0) cout << "global_walls.y0=" << global_walls.y0 << " global_walls.y1=" << global_walls.y1 << endl;
  cout << "rank=" << rank << " local_walls.y0=" << local_walls.y0 << " local_walls.y1=" << local_walls.y1 << " local_particles=" << local_particles << endl;
  MPI::COMM_WORLD.Barrier();
#endif

  vector<pcord_t> up;
  vector<pcord_t> down;
  
  float local_momentum = 0;
  for(unsigned t = 0; t < TIME; t += STEP_SIZE){
    for(unsigned i = 0; i < particles.size()-1; i++){
      float collision_time;
      for(unsigned j = i+1; j < particles.size(); j++){
  	collision_time = collide(&particles.at(i), &particles.at(j));
  	interact(&particles.at(i), &particles.at(j), collision_time);
      }
      if(collision_time == -1) feuler(&particles.at(i), t);

      if(below_wall(&particles.at(i), local_walls) && rank != (world_size - 1)){
	down.push_back(particles.at(i));
	particles.erase(particles.begin() + i);
      }else if(above_wall(&particles.at(i), local_walls) && rank != 0){
	up.push_back(particles.at(i));
	particles.erase(particles.begin() + i);
      }

      local_momentum += wall_collide(&particles.at(i), global_walls);
    }

    // Send stuff upwards
    int up_sends = (int)up.size();
    if(rank != 0){
      MPI::COMM_WORLD.Send(&up_sends, 1, MPI::INTEGER, rank - 1, 0);
      for(int i = 0; i < up_sends; i++) MPI::COMM_WORLD.Send(&up.at(i), 4, MPI::FLOAT, rank - 1, 0);
    }
    
    // Recv stuff from below
  int up_recvs = 0;
    if(rank != world_size - 1){
      MPI::COMM_WORLD.Recv(&up_recvs, 1, MPI::INTEGER, rank + 1, 0);
      pcord_t up_particle;
      for(int i = 0; i < up_recvs; i++){
	MPI::COMM_WORLD.Recv(&up_particle, 4, MPI::FLOAT, rank + 1, 0);
	particles.push_back(up_particle);
      }
    }

    // Send stuff downwards
    int down_sends = (int)down.size();
    if(rank != world_size - 1){
      MPI::COMM_WORLD.Send(&down_sends, 1, MPI::INTEGER, rank + 1, 0);
      for(int i = 0; i < down_sends; i++) MPI::COMM_WORLD.Send(&down.at(i), 4, MPI::FLOAT, rank + 1, 0);
    }

    // Recv stuff from above
    int down_recvs = 0;
    if(rank != 0){
      MPI::COMM_WORLD.Recv(&down_recvs, 1, MPI::INTEGER, rank - 1, 0);
      pcord_t down_particle;
      for(int i = 0; i < down_recvs; i++){
	MPI::COMM_WORLD.Recv(&down_particle, 4, MPI::FLOAT, rank - 1, 0);
	particles.push_back(down_particle);
      }
    }

    // Reset transfer vectors
    up.clear();
    down.clear();      

#ifdef DEBUG
    cout << "t=" << t << "\trank=" << rank << "\tup_sends=" << up_sends << "\tup_recvs=" << up_recvs << "\tdown_sends=" << down_sends <<  "\tdown_recvs=" << down_recvs << "\tpart_count=" << particles.size() << endl;
#endif

    MPI::COMM_WORLD.Barrier();
  }

  float global_momentum;
  MPI::COMM_WORLD.Reduce(&local_momentum, &global_momentum, 1, MPI::FLOAT, MPI::SUM, 0);

  if(rank == 0){
    printf("global_momentum=%f", global_momentum);
    printf(" pressure=%f\n", calc_pressure(global_momentum, global_walls));
  }


  MPI::Finalize();
  return 0;
}
