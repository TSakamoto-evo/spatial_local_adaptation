#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <boost/math/special_functions/bessel.hpp>

void initialize_1d(std::vector<double>& r1, std::vector<double>& u1){
  std::ifstream ifs("output_int_1d.txt");
  if(!ifs){
    std::cerr << "Fail to open the file!" << std::endl;
    std::exit(1);
  }

  r1.clear();
  u1.clear();

  std::vector<std::string> lines;
  std::string line;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<std::string> list;

    while(getline(iss, tmp_list, '\t')){
      list.push_back(tmp_list);
    }

    r1.push_back(std::stod(list.at(0)));
    u1.push_back(std::stod(list.at(1)));
  }
}

void initialize_2d(std::vector<double>& r2, std::vector<double>& u2){
  std::ifstream ifs("output_int_2d.txt");
  if(!ifs){
    std::cerr << "Fail to open the file!" << std::endl;
    std::exit(1);
  }

  r2.clear();
  u2.clear();

  std::vector<std::string> lines;
  std::string line;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<std::string> list;

    while(getline(iss, tmp_list, '\t')){
      list.push_back(tmp_list);
    }

    r2.push_back(std::stod(list.at(0)));
    u2.push_back(std::stod(list.at(1)));
  }
}

void initialize_3d(std::vector<double>& r3, std::vector<double>& u3){
  std::ifstream ifs("output_int_3d.txt");
  if(!ifs){
    std::cerr << "Fail to open the file!" << std::endl;
    std::exit(1);
  }

  r3.clear();
  u3.clear();

  std::vector<std::string> lines;
  std::string line;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<std::string> list;

    while(getline(iss, tmp_list, '\t')){
      list.push_back(tmp_list);
    }

    r3.push_back(std::stod(list.at(0)));
    u3.push_back(std::stod(list.at(1)));
  }
}

double calculate_1d(const std::vector<double>& r1, 
  const std::vector<double>& u1, const double d, const double rs){
  
  int index = 0;
  double sum = 0.0;
  double r_low, r_high, u_low, u_high, u;

  double dist = d - rs;
  while(r1.at(index) < dist){
    index++;
  }
  while(r1.at(index) > dist){
    index--;
  }

  r_low = r1.at(index);
  r_high = r1.at(index + 1);
  u_low = u1.at(index);
  u_high = u1.at(index + 1);

  u = u_low + (u_high - u_low) * (dist - r_low) / (r_high - r_low);
  sum += u;


  dist = d + rs;
  while(r1.at(index) < dist){
    index++;
  }
  while(r1.at(index) > dist){
    index--;
  }

  r_low = r1.at(index);
  r_high = r1.at(index + 1);
  u_low = u1.at(index);
  u_high = u1.at(index + 1);

  u = u_low + (u_high - u_low) * (dist - r_low) / (r_high - r_low);
  sum += u;

  return(sum);
}

double calculate_2d(const std::vector<double>& r2, 
  const std::vector<double>& u2, const double d, const double rs){
  
  int index = 0;
  double sum = 0.0;
  double r_low, r_high, u_low, u_high;
  
  int sep = 1000;
  double pi = std::acos(-1.0);
  double h = pi / sep;

  for(int i = 0; i < sep; i++){
    double theta0 = 1.0 * i * pi / sep;
    double theta1 = 1.0 * (i + 0.5) * pi / sep;
    double theta2 = 1.0 * (i + 1.0) * pi / sep;

    double dist0 = std::sqrt((d + rs * std::cos(theta0)) * (d + rs * std::cos(theta0)) +
      rs * std::sin(theta0) * rs * std::sin(theta0));
    while(r2.at(index) < dist0){
      index++;
    }
    while(r2.at(index) > dist0){
      index--;
    }

    r_low = r2.at(index);
    r_high = r2.at(index + 1);
    u_low = u2.at(index);
    u_high = u2.at(index + 1);

    double uval0 = u_low + (u_high - u_low) * (dist0 - r_low) / (r_high - r_low);

    double dist1 = std::sqrt((d + rs * std::cos(theta1)) * (d + rs * std::cos(theta1)) +
      rs * std::sin(theta1) * rs * std::sin(theta1));
    while(r2.at(index) < dist1){
      index++;
    }
    while(r2.at(index) > dist1){
      index--;
    }

    r_low = r2.at(index);
    r_high = r2.at(index + 1);
    u_low = u2.at(index);
    u_high = u2.at(index + 1);

    double uval1 = u_low + (u_high - u_low) * (dist1 - r_low) / (r_high - r_low);

    double dist2 = std::sqrt((d + rs * std::cos(theta2)) * (d + rs * std::cos(theta2)) +
      rs * std::sin(theta2) * rs * std::sin(theta2));
    while(r2.at(index) < dist2){
      index++;
    }
    while(r2.at(index) > dist2){
      index--;
    }

    r_low = r2.at(index);
    r_high = r2.at(index + 1);
    u_low = u2.at(index);
    u_high = u2.at(index + 1);

    double uval2 = u_low + (u_high - u_low) * (dist2 - r_low) / (r_high - r_low);

    sum += h / 6.0 * (uval0 + 4.0 * uval1 + uval2);
  }

  return(sum);
}

double calculate_3d(const std::vector<double>& r3, 
  const std::vector<double>& u3, const double d, const double rs){
  
  int index = 0;
  double sum = 0.0;
  double r_low, r_high, u_low, u_high;
  
  int sep = 1000;
  double pi = std::acos(-1.0);
  double h = pi / sep;

  for(int i = 0; i < sep; i++){
    double theta0 = 1.0 * i * pi / sep;
    double theta1 = 1.0 * (i + 0.5) * pi / sep;
    double theta2 = 1.0 * (i + 1.0) * pi / sep;

    double dist0 = std::sqrt((d + rs * std::cos(theta0)) * (d + rs * std::cos(theta0)) +
      rs * std::sin(theta0) * rs * std::sin(theta0));
    while(r3.at(index) < dist0){
      index++;
    }
    while(r3.at(index) > dist0){
      index--;
    }

    r_low = r3.at(index);
    r_high = r3.at(index + 1);
    u_low = u3.at(index);
    u_high = u3.at(index + 1);

    double uval0 = u_low + (u_high - u_low) * (dist0 - r_low) / (r_high - r_low);

    double dist1 = std::sqrt((d + rs * std::cos(theta1)) * (d + rs * std::cos(theta1)) +
      rs * std::sin(theta1) * rs * std::sin(theta1));
    while(r3.at(index) < dist1){
      index++;
    }
    while(r3.at(index) > dist1){
      index--;
    }

    r_low = r3.at(index);
    r_high = r3.at(index + 1);
    u_low = u3.at(index);
    u_high = u3.at(index + 1);

    double uval1 = u_low + (u_high - u_low) * (dist1 - r_low) / (r_high - r_low);

    double dist2 = std::sqrt((d + rs * std::cos(theta2)) * (d + rs * std::cos(theta2)) +
      rs * std::sin(theta2) * rs * std::sin(theta2));
    while(r3.at(index) < dist2){
      index++;
    }
    while(r3.at(index) > dist2){
      index--;
    }

    r_low = r3.at(index);
    r_high = r3.at(index + 1);
    u_low = u3.at(index);
    u_high = u3.at(index + 1);

    double uval2 = u_low + (u_high - u_low) * (dist2 - r_low) / (r_high - r_low);

    sum += h / 6.0 * (std::sin(theta0) * uval0 + 
      4.0 * std::sin(theta1) * uval1 + std::sin(theta2) * uval2);
  }

  return(sum);
}

int main(){
  double sigma_sq = 0.005;
  double s2 = 0.005;

  double r = 1.0;
  double rs = 1.0;
  double beta2 = std::sqrt(2.0 * s2 / sigma_sq);
  double pi = std::acos(-1.0);
  double rho = 65536.0;

  double d_max = 16.0;
  int sep = 100;

  std::vector<double> r1, u1, r2, u2, r3, u3;

  initialize_1d(r1, u1);
  initialize_2d(r2, u2);
  initialize_3d(r3, u3);

  std::ofstream ofs("lambda.txt");

  for(int i = 0; i <= sep; i++){
    double d = (r + rs) + (d_max - r - rs) * i / sep;

    // 1d
    double l1 = calculate_1d(r1, u1, d, rs);
    l1 *= beta2 * sigma_sq * rho / (1.0 + std::exp(-2.0 * beta2 * rs));
    ofs << d << "\t" << 1.0 / l1 << "\t" << 1.0 / l1 + (d - r - rs) / beta2 / sigma_sq;

    // 2d
    double l2 = calculate_2d(r2, u2, d, rs);
    l2 *= sigma_sq * rho / boost::math::cyl_bessel_k(0.0, beta2 * rs) / 
      boost::math::cyl_bessel_i(0.0, beta2 * rs);
    ofs << "\t" << 1.0 / l2 << "\t" << 1.0 / l2 + (d - r - rs) / beta2 / sigma_sq;

    // 3d
    double l3 = calculate_3d(r3, u3, d, rs);
    l3 *= 2.0 * beta2 * sigma_sq * pi * rs * rs * rho / (1.0 - std::exp(-2.0 * beta2 * rs));
    ofs << "\t" << 1.0 / l3 << "\t" << 1.0 / l3 + (d - r - rs) / beta2 / sigma_sq << std::endl;;
  }
}