#include "population.hpp"

Population::Population(const Parameter input_para, const int ini_pos){
  para = input_para;
  metapop.clear();

  std::random_device seed;
  mt.seed(seed());

  // calculate the area in the focal grid that belong to the target patch
  double area1 = calculate_area(ini_pos, 0.0, para.r);
  Grid ini(para.pop_size, 1, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
  metapop.emplace(ini_pos, ini);

  // make a list of grids where allele A exists
  active_list.clear();
  active_list.push_back(ini_pos);

  threshold = 2 * para.r * para.pop_density * 0.5;
  gen = 0;
}

double Population::calculate_area(const int pos, const double center, 
  const double rs) const{

  double x_min = (pos - 0.5) / para.n_div;
  double x_max = (pos + 0.5) / para.n_div;

  double overlap = (std::min(center + rs, x_max) - 
    std::max(center - rs, x_min)) / (x_max - x_min);
  if(overlap < 0.0){
  	return(0.0);
  }else{
  	return(overlap);
  }
}

void Population::migration(){
  for(const auto& i: active_list){
  	metapop.at(i).emigration(para.m);
  }

  for(const auto& i: active_list){
  	double mig_add = metapop.at(i).return_num_allele() * 1.0 / para.pop_size * para.m / 2.0;

    // if a destination grid does not exist in the map, prepare that grid
  	if(metapop.count(i - 1) == 0){
  	  double area1 = calculate_area(i - 1, 0, para.r);
  	  Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
  	  metapop.emplace(i - 1, ini);
  	}
  	if(metapop.count(i + 1) == 0){
  	  double area1 = calculate_area(i + 1, 0, para.r);
  	  Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
  	  metapop.emplace(i + 1, ini);
  	}

  	metapop.at(i - 1).immigration(mig_add);
  	metapop.at(i + 1).immigration(mig_add);
  }
}

int Population::reproduction(){
  active_list.clear();
  erase_list.clear();

  for(auto& i: metapop){
  	bool active = i.second.reproduction(mt);
  	if(active){
  	  active_list.push_back(i.first);
  	}else if(gen % 10000 == 0){
      erase_list.push_back(i.first);
    }
  }

  if(erase_list.size() > 0){
    for(const auto& i: erase_list){
      metapop.erase(i);
    }
  }

  double sum = 0.0;
  for(const auto& i: active_list){
    sum += metapop.at(i).return_area_num();
  }

  gen++;

  if(active_list.size() == 0){
  	return(1);
  }else if(sum > threshold){
    return(2);
  }else{
  	return(0);
  }
}
