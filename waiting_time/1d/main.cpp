#include "parameter.hpp"
#include "population.hpp"
#include <fstream>
#include <chrono>

int main(int argc, char* argv[]){
  int pop_density, n_div, para_set1, para_set2, run_index;
  double vm, s1, s2, r, rs, dist;

  if(argc == 12){
    sscanf(argv[1], "%d", &pop_density);
    sscanf(argv[2], "%d", &n_div);
    sscanf(argv[3], "%lf", &vm);
    sscanf(argv[4], "%lf", &s1);
    sscanf(argv[5], "%lf", &s2);
    sscanf(argv[6], "%lf", &r);
    sscanf(argv[7], "%lf", &rs);
    sscanf(argv[8], "%d", &para_set1);
    sscanf(argv[9], "%d", &para_set2);
    sscanf(argv[10], "%lf", &dist);
    sscanf(argv[11], "%d", &run_index);
  }else{
    std::cerr << "Error!" << std::endl;
    std::exit(1);
  }

  Parameter para;

  para.pop_density = pop_density;
  para.n_div = n_div;
  para.vm = vm;
  para.s1 = s1;
  para.s2 = s2;
  para.r = r;

  // population size for each grid
  para.pop_size = para.pop_density / para.n_div;

  // twice of migration rate between adjacent grids 
  para.m = para.vm * para.n_div * para.n_div;

  para.dist = dist;
  para.rs = rs;

  // maximum length of each replication
  int max_gen = 2e+7;
  int generation = 0;

  Population pop(para, generation);
  int state = 0;

  std::chrono::system_clock::time_point last_output;
  last_output = std::chrono::system_clock::now();

  // while alleles A and a coexist and generation is smaller than max_gen
  while(state == 0 && generation < max_gen){
    pop.migration();
    state = pop.reproduction();
    generation++;

    // make an intermediate file every five minutes
    if(generation % 1000 == 0){
      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      auto time = now - last_output;
      double minute = std::chrono::duration_cast<std::chrono::minutes>(time).count();

      if(minute > 5){
        pop.regi_state();
        last_output = now;
      }
    }
  }

  std::ofstream ofs("wait_time.txt", std::ios::app);
  ofs << dist << "\t" << para.vm << "\t" << para.s1 << "\t" << para.s2
    << "\t" << rs << "\t" << state << "\t" << generation << "\t" 
    << para_set1 << "\t" << para_set2 << "\t" << run_index << std::endl;
  ofs.close();

  return(0);
}
