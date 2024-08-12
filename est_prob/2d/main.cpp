#include "parameter.hpp"
#include "population.hpp"
#include <fstream>
#include <chrono>

int main(int argc, char* argv[]){
  int pop_density, n_div, rep, para_set1, para_set2, ini, done, count;
  double vm, s1, s2, r;

  if(argc == 13){
    sscanf(argv[1], "%d", &pop_density);
    sscanf(argv[2], "%d", &n_div);
    sscanf(argv[3], "%lf", &vm);
    sscanf(argv[4], "%lf", &s1);
    sscanf(argv[5], "%lf", &s2);
    sscanf(argv[6], "%lf", &r);
    sscanf(argv[7], "%d", &rep);
    sscanf(argv[8], "%d", &para_set1);
    sscanf(argv[9], "%d", &para_set2);
    sscanf(argv[10], "%d", &ini);
    sscanf(argv[11], "%d", &done);
    sscanf(argv[12], "%d", &count);
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
  para.pop_size = para.pop_density / para.n_div / para.n_div;

  // twice of migration rate between adjacent grids 
  para.m = para.vm * para.n_div * para.n_div;

  // distance in the absolute scale
  double pos = 1.0 * ini / para.n_div;

  std::chrono::system_clock::time_point last_output;
  last_output = std::chrono::system_clock::now();

  for(int i = done + 1; i <= rep; i++){
    Population pop(para, ini, 0);
    int state = 0;

    // while alleles A and a coexist
    while(state == 0){
      pop.migration();
      state = pop.reproduction();
    }

    // established
    if(state == 2){
      count++;
    }

    // make an intermediate file every five minutes
    if(i % 1000 == 0){
      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      auto time = now - last_output;
      double minute = std::chrono::duration_cast<std::chrono::minutes>(time).count();

      if(minute > 5){
        std::ofstream regi("regi_state.txt");
        regi << i << "\t" << count << std::endl;
        last_output = now;
      }
    }
  }

  std::ofstream ofs("est_prob.txt", std::ios::app);
  ofs << pos << "\t" << para.vm << "\t" << para.s1 << "\t" << para.s2
    << "\t" << count << "\t" << rep << "\t" << para_set1 <<
    "\t" << para_set2 << std::endl;
  ofs.close();

  return(0);
}