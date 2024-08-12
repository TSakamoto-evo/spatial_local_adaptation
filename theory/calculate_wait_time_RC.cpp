#include <iostream>
#include <fstream>
#include <cmath>

double calculate_a_one(const double r, const double d, const double a){
  if(a > d + r){
    return(2.0 * r);
  }else{
    return(a - (d - r));
  }
}

double calculate_a_two(const double r, const double d, const double a){
  double sum = 0.0;

  int sep = 1000;
  double h = 2.0 * r / sep;

  for(int i = 0; i < sep; i++){
    double x0 = -r + 2.0 * r * i / sep;
    double x1 = -r + 2.0 * r * (i + 0.5) / sep;
    double x2 = -r + 2.0 * r * (i + 1) / sep;

    double y0, y1, y2;
    if(a > d - x0){
      double tmp = std::min(std::sqrt(r * r - x0 * x0), std::sqrt(a * a - (d - x0) * (d - x0)));
      y0 = 2.0 * tmp;
    }else{
      y0 = 0.0;
    }
    if(a > d - x1){
      double tmp = std::min(std::sqrt(r * r - x1 * x1), std::sqrt(a * a - (d - x1) * (d - x1)));
      y1 = 2.0 * tmp;
    }else{
      y1 = 0.0;
    }
    if(a > d - x2){
      double tmp = std::min(std::sqrt(r * r - x2 * x2), std::sqrt(a * a - (d - x2) * (d - x2)));
      y2 = 2.0 * tmp;
    }else{
      y2 = 0.0;
    }

    sum += h / 6.0 * (y0 + 4.0 * y1 + y2);
  }
  return(sum);
}

double calculate_a_three(const double r, const double d, const double a){
  double sum = 0.0;

  int sep = 1000;
  double h = 2.0 * r / sep;

  for(int i = 0; i < sep; i++){
    double x0 = -r + 2.0 * r * i / sep;
    double x1 = -r + 2.0 * r * (i + 0.5) / sep;
    double x2 = -r + 2.0 * r * (i + 1) / sep;

    double y0, y1, y2;
    if(a > d - x0){
      double tmp = std::min(std::sqrt(r * r - x0 * x0), std::sqrt(a * a - (d - x0) * (d - x0)));
      y0 = std::acos(-1.0) * tmp * tmp;
    }else{
      y0 = 0.0;
    }
    if(a > d - x1){
      double tmp = std::min(std::sqrt(r * r - x1 * x1), std::sqrt(a * a - (d - x1) * (d - x1)));
      y1 = std::acos(-1.0) * tmp * tmp;
    }else{
      y1 = 0.0;
    }
    if(a > d - x2){
      double tmp = std::min(std::sqrt(r * r - x2 * x2), std::sqrt(a * a - (d - x2) * (d - x2)));
      y2 = std::acos(-1.0) * tmp * tmp;
    }else{
      y2 = 0.0;
    }

    sum += h / 6.0 * (y0 + 4.0 * y1 + y2);
  }
  return(sum);
}

int main(){
  double s1 = 0.05;
  double s2 = 0.005;
  double sigma_sq = 0.001;
  double rho = 65536.0;

  double beta2 = std::sqrt(2.0 * s2 / sigma_sq);

  double r = 1.0;
  double rs = 1.0;

  double min_d = r + rs;
  double max_d = 16.0;

  int sep = 100;

  std::ofstream ofs("predict_RC.txt");

  for(int i = 0; i <= sep; i++){
    double d = min_d + (max_d - min_d) * i / sep;
    double rr = d - r - rs;

    if(rr * beta2 > 1.0){
      double a_one = calculate_a_one(r, d, rr + 1.0 / beta2 + rs);
      double a_two = calculate_a_two(r, d, rr + 1.0 / beta2 + rs);
      double a_three = calculate_a_three(r, d, rr + 1.0 / beta2 + rs);

      double l1 = a_one * rho * s2 * std::min(s2, 2.0 * s1) * 
        std::pow(beta2 * rr, -(1.0 - 1.0) / 2.0) * std::exp(-rr * beta2);
      double l2 = a_two * rho * s2 * std::min(s2, 2.0 * s1) * 
        std::pow(beta2 * rr, -(2.0 - 1.0) / 2.0) * std::exp(-rr * beta2);
      double l3 = a_three * rho * s2 * std::min(s2, 2.0 * s1) * 
        std::pow(beta2 * rr, -(3.0 - 1.0) / 2.0) * std::exp(-rr * beta2);

      ofs << d << "\t" << 1.0 / l1 << "\t" << 1.0 / l2 << "\t" << 1.0 / l3 << std::endl;
    }else{
      ofs << d << "\tNA\tNA\tNA" << std::endl;
    }
  }
}