#ifndef GRID
#define GRID

#include <boost/random.hpp>
#include <cmath>
#include <vector>

class Grid{
private:
  int pop_size;
  int num_allele;
  double expect;
  double area1;
  double area2;
  double s;
  bool active;

public:
  Grid(const int input_pop_size, const int input_num_allele,
  		const double input_area1, const double input_area2,
      const double input_s){

    pop_size = input_pop_size;
    num_allele = input_num_allele;
    expect = 0.0;
    area1 = input_area1;
    area2 = input_area2;
    s = input_s;
    active = 1;
  }
  void emigration(const double m){
  	double p = 1.0 * num_allele / pop_size;
  	expect = p / (p + std::exp(-s) * (1.0 - p)) - m * p;
  }
  void immigration(const double add){
  	if(active){
  	  expect += add;
  	}else{
  	  expect = add;
  	  active = 1;
  	}
  }
  bool reproduction(boost::random::mt19937& mt){
  	if(active){
	    boost::random::binomial_distribution<> det_num(pop_size, expect);
  	  num_allele = det_num(mt);

  	  if(num_allele == 0){
  	    active = 0;
  	  }
  	}
  	return(active);
  }
  bool return_active() const{ return(active); };
  int return_num_allele() const{ return(num_allele); };
  double return_area1_num() const{ return(num_allele * area1); };
  double return_area2_num() const{ return(num_allele * area2); };

  void return_all_state(std::vector<int>& ret1, std::vector<double>& ret2) const{
    ret1.clear();
    ret1.push_back(pop_size);
    ret1.push_back(num_allele);

    ret2.clear();
    ret2.push_back(area1);
    ret2.push_back(area2);
    ret2.push_back(s);
  }
};

#endif
