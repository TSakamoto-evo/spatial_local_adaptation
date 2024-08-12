#include "population.hpp"

Population::Population(const Parameter input_para, int& generation){
  para = input_para;
  metapop.clear();

  std::random_device seed;
  mt.seed(seed());

  // if there is a suspended file, restart with that state
  if(std::filesystem::is_regular_file("regi_state.txt")){
    initialize_from_regi();
  }else{
    // add grids belonging to the source patch
    int x_min = (para.dist - para.rs) * para.n_div + 0.5;
    int x_max = (para.dist + para.rs) * para.n_div + 0.5;

    active_list.clear();

    for(int i = x_min; i <= x_max; i++){
      double area1 = calculate_area(i, 0.0, para.r);
      double area2 = calculate_area(i, para.dist, para.rs);
      double area = area1 + area2;
      boost::random::binomial_distribution<> det_ini(para.pop_size, area2);
      int ini_num = det_ini(mt);
      Grid ini(para.pop_size, ini_num, area1, area2,
        para.s1 * area - para.s2 * (1.0 - area));
      metapop.emplace(i, ini);

      active_list.push_back(i);
    }
    gen = 0;
  }

  threshold = 2 * para.r * para.pop_density * 0.5;
  generation = gen;
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
  	double mig_add = metapop.at(i).return_num_allele() *
      1.0 / para.pop_size * para.m / 2.0;

    // if a destination grid does not exist in the map, prepare that grid
  	if(metapop.count(i - 1) == 0){
  	  double area1 = calculate_area(i - 1, 0, para.r);
      double area2 = calculate_area(i - 1, para.dist, para.rs);
      double area = area1 + area2;
      Grid ini(para.pop_size, 0, area1, area2,
        para.s1 * area - para.s2 * (1.0 - area));
  	  metapop.emplace(i - 1, ini);
  	}
  	if(metapop.count(i + 1) == 0){
      double area1 = calculate_area(i + 1, 0, para.r);
      double area2 = calculate_area(i + 1, para.dist, para.rs);
      double area = area1 + area2;
      Grid ini(para.pop_size, 0, area1, area2,
        para.s1 * area - para.s2 * (1.0 - area));
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
    sum += metapop.at(i).return_area1_num();
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

void Population::regi_state(){
  std::ofstream ofs("regi_state.txt");
  ofs << gen << std::endl;

  std::vector<int> ret1;
  std::vector<double> ret2;

  for(const auto& i: active_list){
    metapop.at(i).return_all_state(ret1, ret2);
    ofs << i << "\t" << ret1.at(0) << "\t" << ret1.at(1) << "\t"
      << ret2.at(0) << "\t" << ret2.at(1) << "\t" << ret2.at(2) << std::endl;
  }
}

void Population::initialize_from_regi(){
  std::ifstream ifs("regi_state.txt");
  if(!ifs){
    std::cerr << "Fail to open the file!" << std::endl;
    std::exit(1);
  }

  std::vector<std::string> lines;
  std::string line;
  int num_line = 0;

  active_list.clear();

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<std::string> list;

    while(getline(iss, tmp_list, '\t')){
      list.push_back(tmp_list);
    }

    if(num_line == 0){
      gen = std::stoi(list[0]);
    }else{
      int i0 = std::stoi(list[0]);

      int pop_size = std::stoi(list[1]);
      int num_allele = std::stoi(list[2]);

      double area1 = std::stod(list[3]);
      double area2 = std::stod(list[4]);
      double s = std::stod(list[5]);

      Grid ini(pop_size, num_allele, area1, area2, s);
      metapop.emplace(i0, ini);
      
      active_list.push_back(i0);
    }

    num_line++;
  }
}
