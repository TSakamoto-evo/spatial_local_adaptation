#ifndef POPULATION
#define POPULATION

#include <unordered_map>
#include <random>
#include <boost/random.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "parameter.hpp"

class Population{
private:
  Parameter para;
  std::unordered_map<int, std::unordered_map<int, Grid>> metapop;
  std::vector<std::vector<int>> active_list;
  std::vector<std::vector<int>> erase_list;
  std::vector<double> center1;
  boost::random::mt19937 mt;
  
  double threshold;
  int gen;

  double calculate_area(const std::vector<int>& pos, 
    const std::vector<double>& center, const double rs) const;

public:
  Population(const Parameter input_para, const int ini_pos1, const int ini_pos2);
  void migration();
  int reproduction();
};

#endif