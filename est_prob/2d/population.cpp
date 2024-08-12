#include "population.hpp"

Population::Population(const Parameter input_para, const int ini_pos1, 
  const int ini_pos2){

  para = input_para;
  metapop.clear();

  std::random_device seed;
  mt.seed(seed());

  std::vector<int> tmp_vec = {ini_pos1, ini_pos2};
  
  // position of the center of the target patch
  center1.clear();
  center1.push_back(0.0);
  center1.push_back(0.0);

  // calculate the area in the focal grid that belong to the target patch
  double area1 = calculate_area(tmp_vec, center1, para.r);
  Grid ini(para.pop_size, 1, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
  std::unordered_map<int, Grid> tmp_map;
  tmp_map.emplace(ini_pos2, ini);
  metapop.emplace(ini_pos1, tmp_map);

  // make a list of grids where allele A exists
  active_list.clear();
  active_list.push_back(tmp_vec);

  threshold = std::acos(-1.0) * para.r * para.r * para.pop_density * 0.5;
  gen = 0;
}

double Population::calculate_area(const std::vector<int>& pos, 
  const std::vector<double>& center, const double rs) const{
  
  int sep = 1000;

  double x_min = (pos[0] - 0.5) / para.n_div;
  double x_max = (pos[0] + 0.5) / para.n_div;

  double y_min = (pos[1] - 0.5) / para.n_div;
  double y_max = (pos[1] + 0.5) / para.n_div;

  double ori_x = center[0];
  double ori_y = center[1];

  double x0, x1, y0, y1;
  if((x_min - ori_x) * (x_max - ori_x) < 0.0){
    x0 = 0.0;
  }else{
    x0 = std::min(std::abs(x_min - ori_x), std::abs(x_max - ori_x));
  }
  x1 = std::max(std::abs(x_min - ori_x), std::abs(x_max - ori_x));

  if((y_min - ori_y) * (y_max - ori_y) < 0.0){
    y0 = 0.0;
  }else{
    y0 = std::min(std::abs((y_min - ori_y)), std::abs((y_max - ori_y)));
  }
  y1 = std::max(std::abs(y_min - ori_y), std::abs(y_max - ori_y));

  double r2_min = x0 * x0 + y0 * y0;
  double r2_max = x1 * x1 + y1 * y1;

  if(r2_min > rs * rs){
    return(0.0);
  }else if(r2_max < rs * rs){
    return(1.0);
  }else{
    double sum = 0.0;
    double h = (x_max - x_min) / sep;

    for(int i = 0; i <= sep; i++){
      double xx = x_min + (x_max - x_min) / sep * i;
      double yy;

      if(rs * rs > (xx - ori_x) * (xx - ori_x)){
        double r_min = -std::sqrt(rs * rs - (xx - ori_x) * (xx - ori_x));
        double r_max = std::sqrt(rs * rs - (xx - ori_x) * (xx - ori_x));

        yy = std::min(ori_y + r_max, y_max) - std::max(ori_y + r_min, y_min);

        if(yy < 0.0){
          yy = 0.0;
        }
      }else{
        yy = 0.0;
      }

      if(i == 0 || i == sep){
        sum += yy * h / 2.0;
      }else{
        sum += yy * h;
      }
    }

    double ret = sum / (x_max - x_min) / (y_max - y_min);
    return(ret);
  }
}

void Population::migration(){
  for(const auto& i: active_list){
  	metapop.at(i[0]).at(i[1]).emigration(2.0 * para.m);
  }

  for(const auto& i: active_list){
  	double mig_add = metapop.at(i[0]).at(i[1]).return_num_allele() * 1.0 / para.pop_size * para.m / 2.0;

    // if a destination grid does not exist in the map, prepare that grid
  	if(metapop.count(i[0] - 1) == 0){
      std::vector<int> tmp_vec = {i[0] - 1, i[1]};
      double area1 = calculate_area(tmp_vec, center1, para.r);
      Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
      std::unordered_map<int, Grid> tmp_map;
      tmp_map.emplace(i[1], ini);
      metapop.emplace(i[0] - 1, tmp_map);
    }else if(metapop.at(i[0] - 1).count(i[1]) == 0){
      std::vector<int> tmp_vec = {i[0] - 1, i[1]};
      double area1 = calculate_area(tmp_vec, center1, para.r);
      Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
      metapop.at(i[0] - 1).emplace(i[1], ini);
    }

    if(metapop.count(i[0] + 1) == 0){
      std::vector<int> tmp_vec = {i[0] + 1, i[1]};
      double area1 = calculate_area(tmp_vec, center1, para.r);
      Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
      std::unordered_map<int, Grid> tmp_map;
      tmp_map.emplace(i[1], ini);
      metapop.emplace(i[0] + 1, tmp_map);
    }else if(metapop.at(i[0] + 1).count(i[1]) == 0){
      std::vector<int> tmp_vec = {i[0] + 1, i[1]};
      double area1 = calculate_area(tmp_vec, center1, para.r);
      Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
      metapop.at(i[0] + 1).emplace(i[1], ini);
    }

    if(metapop.at(i[0]).count(i[1] - 1) == 0){
      std::vector<int> tmp_vec = {i[0], i[1] - 1};
      double area1 = calculate_area(tmp_vec, center1, para.r);
      Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
      metapop.at(i[0]).emplace(i[1] - 1, ini);
    }
    if(metapop.at(i[0]).count(i[1] + 1) == 0){
      std::vector<int> tmp_vec = {i[0], i[1] + 1};
      double area1 = calculate_area(tmp_vec, center1, para.r);
      Grid ini(para.pop_size, 0, area1, para.s1 * area1 - para.s2 * (1.0 - area1));
      metapop.at(i[0]).emplace(i[1] + 1, ini);
    }
  	metapop.at(i[0] - 1).at(i[1]).immigration(mig_add);
    metapop.at(i[0] + 1).at(i[1]).immigration(mig_add);
    metapop.at(i[0]).at(i[1] - 1).immigration(mig_add);
    metapop.at(i[0]).at(i[1] + 1).immigration(mig_add);
  }
}

int Population::reproduction(){
  active_list.clear();
  erase_list.clear();

  for(auto& i: metapop){
    for(auto& j: i.second){
      bool active = j.second.reproduction(mt);
      if(active){
        std::vector<int> tmp = {i.first, j.first};
        active_list.push_back(tmp);
      }else if(gen % 10000 == 0){
        std::vector<int> tmp = {i.first, j.first};
        erase_list.push_back(tmp);
      }
    }
  }

  if(erase_list.size() > 0){
    for(const auto& i: erase_list){
      metapop.at(i[0]).erase(i[1]);
    }
  }

  double sum = 0.0;
  for(const auto& i: active_list){
    sum += metapop.at(i[0]).at(i[1]).return_area_num();
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
