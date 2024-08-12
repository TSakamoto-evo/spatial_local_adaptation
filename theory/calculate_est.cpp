// use with options g++ calculate_est.cpp -Wall -Wextra -O3 -std=c++17 -lgmp -lgmpxx -o test.out

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <gmp.h>
#include <gmpxx.h>

namespace mp = boost::multiprecision;
typedef mp::number<mp::gmp_float<50> > mpf_dec_float_prec;

std::vector<mpf_dec_float_prec> pertubation(
  const mpf_dec_float_prec u0, const mpf_dec_float_prec s1,
  const mpf_dec_float_prec sigma_sq, const mpf_dec_float_prec d,
  const mpf_dec_float_prec r0, const int sep, const int term,
  std::vector<mpf_dec_float_prec>& rs, std::vector<mpf_dec_float_prec>& us){

  // calculation of u(r) for small r by perturbation in r
  rs.clear();
  us.clear();

  // u0 should be specified in advance
  std::vector<mpf_dec_float_prec> cs = {u0};

  // determine polynomial coefficients
  for(int i = 0; i <= term; i++){
    mpf_dec_float_prec sum(0);
    mpf_dec_float_prec coef(0);

    for(int j = 0; j <= i; j++){
      sum += cs.at(j) * cs.at(i - j) / mpf_dec_float_prec(2);
    }
    sum -= s1 * cs.at(i);

    coef += (mpf_dec_float_prec(i) + mpf_dec_float_prec(1)) *
      (d - mpf_dec_float_prec(1)) * sigma_sq +
      (mpf_dec_float_prec(i) + mpf_dec_float_prec(1)) *
      (mpf_dec_float_prec(2) * mpf_dec_float_prec(i) + mpf_dec_float_prec(1)) *
      sigma_sq;

    cs.push_back(sum / coef);
  }

  mpf_dec_float_prec r_max;

  if(r0 > mpf_dec_float_prec("0.01")){
    r_max = mpf_dec_float_prec("0.01");
  }else{
    r_max = r0;
  }

  // calculate u(d) as a function of r
  for(int i = 0; i <= sep; i++){
    mpf_dec_float_prec r = mpf_dec_float_prec(i) * r_max / mpf_dec_float_prec(sep);

    mpf_dec_float_prec val(0);
    mpf_dec_float_prec pow(1);
    for(int j = 0; j <= term; j++){
      val += cs.at(j) * pow;
      pow *= r * r;
    }

    rs.push_back(r);
    us.push_back(val);
  }

  mpf_dec_float_prec u0_end(0);
  mpf_dec_float_prec u1_end(0);
  mpf_dec_float_prec pow(1);

  for(int i = 0; i <= term; i++){
    u0_end += cs.at(i) * pow;

    if(i > 0){
      u1_end += (mpf_dec_float_prec(2) * mpf_dec_float_prec(i)) * cs.at(i) *
        pow / r_max;
    }

    pow *= r_max * r_max;
  }

  // return a set of r_max, u(r_max), and du(r_max)/dr
  std::vector<mpf_dec_float_prec> ret = {r_max, u0_end, u1_end};

  return(ret);
}

int runge_kutta(const mpf_dec_float_prec initial_u0,
  const mpf_dec_float_prec initial_u1, const mpf_dec_float_prec initial_r,
  const mpf_dec_float_prec s1, const mpf_dec_float_prec s2,
  const mpf_dec_float_prec sigma_sq, const mpf_dec_float_prec d,
  const mpf_dec_float_prec r0, const int sep,
  const mpf_dec_float_prec r_max,
  std::vector<mpf_dec_float_prec>& rs, std::vector<mpf_dec_float_prec>& us,
  double& ur){

  // use a result of perturbation calculation as an input
  mpf_dec_float_prec r = initial_r;
  mpf_dec_float_prec u0 = initial_u0;
  mpf_dec_float_prec u1 = initial_u1;
  mpf_dec_float_prec h = (r0 - initial_r) / mpf_dec_float_prec(sep);

  int result = 0;
  ur = 0.0;

  // runge-kutta integration (for r < R)
  for(int i = 0; i < sep; i++){
    mpf_dec_float_prec tmp_u0 = u0;
    mpf_dec_float_prec tmp_u1 = u1;
    mpf_dec_float_prec tmp_r = r;

    mpf_dec_float_prec k1_u0 = tmp_u1;
    mpf_dec_float_prec k1_u1 = tmp_u0 * tmp_u0 / sigma_sq -
      s1 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    tmp_u0 = u0 + h / mpf_dec_float_prec(2) * k1_u0;
    tmp_u1 = u1 + h / mpf_dec_float_prec(2) * k1_u1;
    tmp_r = r + h / mpf_dec_float_prec(2);

    mpf_dec_float_prec k2_u0 = tmp_u1;
    mpf_dec_float_prec k2_u1 = tmp_u0 * tmp_u0 / sigma_sq -
      s1 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    tmp_u0 = u0 + h / mpf_dec_float_prec(2) * k2_u0;
    tmp_u1 = u1 + h / mpf_dec_float_prec(2) * k2_u1;
    tmp_r = r + h / mpf_dec_float_prec(2);

    mpf_dec_float_prec k3_u0 = tmp_u1;
    mpf_dec_float_prec k3_u1 = tmp_u0 * tmp_u0 / sigma_sq -
      s1 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    tmp_u0 = u0 + h * k3_u0;
    tmp_u1 = u1 + h * k3_u1;
    tmp_r = r + h;

    mpf_dec_float_prec k4_u0 = tmp_u1;
    mpf_dec_float_prec k4_u1 = tmp_u0 * tmp_u0 / sigma_sq -
      s1 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    u0 += h / mpf_dec_float_prec(6) * (k1_u0 + mpf_dec_float_prec(2) * k2_u0 +
      mpf_dec_float_prec(2) * k3_u0 + k4_u0);
    u1 += h / mpf_dec_float_prec(6) * (k1_u1 + mpf_dec_float_prec(2) * k2_u1 +
      mpf_dec_float_prec(2) * k3_u1 + k4_u1);
    r = initial_r + mpf_dec_float_prec(i + 1) *
      (r0 - initial_r) / mpf_dec_float_prec(sep);

    rs.push_back(r);
    us.push_back(u0);

    if(u0 < 0.0){
      result = -1;
      return(result);
    }else if(u0 > 1.0){
      result = 1;
      return(result);
    }
  }

  ur = static_cast<double>(us.back());

  // runge-kutta integration (for r > R)
  while(r < r_max){
    mpf_dec_float_prec tmp_u0 = u0;
    mpf_dec_float_prec tmp_u1 = u1;
    mpf_dec_float_prec tmp_r = r;

    mpf_dec_float_prec k1_u0 = tmp_u1;
    mpf_dec_float_prec k1_u1 = tmp_u0 * tmp_u0 / sigma_sq +
      s2 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    tmp_u0 = u0 + h / mpf_dec_float_prec(2) * k1_u0;
    tmp_u1 = u1 + h / mpf_dec_float_prec(2) * k1_u1;
    tmp_r = r + h / mpf_dec_float_prec(2);

    mpf_dec_float_prec k2_u0 = tmp_u1;
    mpf_dec_float_prec k2_u1 = tmp_u0 * tmp_u0 / sigma_sq +
      s2 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    tmp_u0 = u0 + h / mpf_dec_float_prec(2) * k2_u0;
    tmp_u1 = u1 + h / mpf_dec_float_prec(2) * k2_u1;
    tmp_r = r + h / mpf_dec_float_prec(2);

    mpf_dec_float_prec k3_u0 = tmp_u1;
    mpf_dec_float_prec k3_u1 = tmp_u0 * tmp_u0 / sigma_sq +
      s2 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    tmp_u0 = u0 + h * k3_u0;
    tmp_u1 = u1 + h * k3_u1;
    tmp_r = r + h;

    mpf_dec_float_prec k4_u0 = tmp_u1;
    mpf_dec_float_prec k4_u1 = tmp_u0 * tmp_u0 / sigma_sq +
      s2 * mpf_dec_float_prec(2) / sigma_sq * tmp_u0 -
      (d - mpf_dec_float_prec(1)) / tmp_r * tmp_u1;

    u0 += h / mpf_dec_float_prec(6) * (k1_u0 + mpf_dec_float_prec(2) * k2_u0 +
      mpf_dec_float_prec(2) * k3_u0 + k4_u0);
    u1 += h / mpf_dec_float_prec(6) * (k1_u1 + mpf_dec_float_prec(2) * k2_u1 +
      mpf_dec_float_prec(2) * k3_u1 + k4_u1);
    r += h;

    rs.push_back(r);
    us.push_back(u0);

    if(u0 < 0.0){
      result = -1;
      return(result);
    }else if(u0 > 1.0){
      result = 1;
      return(result);
    }
  }

  return(result);
}

mpf_dec_float_prec search(const mpf_dec_float_prec s1,
  const mpf_dec_float_prec s2, const mpf_dec_float_prec sigma_sq,
  const mpf_dec_float_prec d, const mpf_dec_float_prec r0, const int sep,
  const int term, const mpf_dec_float_prec r_max,
  std::vector<mpf_dec_float_prec>& r, std::vector<mpf_dec_float_prec>& u,
  double& ur){

  // determine u0 that gives u(r) -> 0 as r -> infty

  // maximum possible value (that is u0 under global positive selection)
  mpf_dec_float_prec ini_u0 = mpf_dec_float_prec(2) * s1;

  mpf_dec_float_prec h = mpf_dec_float_prec("0.1") * s1;

  std::cout << "start search: " << std::flush;

  for(int i = 0; i < 50; i++){
    while(1){
      std::vector<mpf_dec_float_prec> ret = pertubation(ini_u0, s1, sigma_sq, d, r0, sep, term, r, u);
      int result = runge_kutta(ret.at(1), ret.at(2), ret.at(0), s1, s2, sigma_sq, d, r0, sep, r_max, r, u, ur);

      if(result >= 0){
        ini_u0 -= h;
      }else{
        ini_u0 += h;
        break;
      }
    }
    h /= mpf_dec_float_prec(10);
    std::cout << "*" << std::flush;
  }
  std::cout << "|" << std::endl;

  std::vector<mpf_dec_float_prec> ret = pertubation(ini_u0, s1, sigma_sq, d, r0, sep, term, r, u);
  runge_kutta(ret.at(1), ret.at(2), ret.at(0), s1, s2, sigma_sq, d, r0, sep, r_max, r, u, ur);

  return(ini_u0);
}

double calculate_sum(const std::vector<mpf_dec_float_prec>& r,
  const std::vector<mpf_dec_float_prec>& u, const int d, const double r_max){

  double sum = 0.0;

  int index = 0;
  while(r.at(index + 1) < r_max && u.at(index + 1) <= u.at(index)){
    double r0 = static_cast<double>(r.at(index));
    double r1 = static_cast<double>(r.at(index + 1));

    double val0 = static_cast<double>(u.at(index));
    if(d == 1){
      val0 *= 2;
    }else if(d == 2){
      val0 *= 2 * std::acos(-1.0) * r0;
    }else if(d == 3){
      val0 *= 4 * std::acos(-1.0) * r0 * r0;
    }

    double val1 = static_cast<double>(u.at(index + 1));
    if(d == 1){
      val1 *= 2;
    }else if(d == 2){
      val1 *= 2 * std::acos(-1.0) * r1;
    }else if(d == 3){
      val1 *= 4 * std::acos(-1.0) * r1 * r1;
    }

    sum += (val0 + val1) / 2.0 * (r1 - r0);

    index++;
  }

  std::cout << "OK" << std::endl;
  return(sum);
}

int main(){
  mpf_dec_float_prec sigma_sq("0.005");
  mpf_dec_float_prec s1("0.05");
  mpf_dec_float_prec s2("0.005");

  mpf_dec_float_prec r1("1.0");
  mpf_dec_float_prec r2("1.0");
  mpf_dec_float_prec r3("1.0");

  int sep = 1000;

  {
    mpf_dec_float_prec r0 = r1;
    mpf_dec_float_prec d(1);

    std::vector<mpf_dec_float_prec> r;
    std::vector<mpf_dec_float_prec> u;

    double ur;

    // determine u0 such that u(r) -> 0 as r -> infty
    mpf_dec_float_prec u0 = search(s1, s2, sigma_sq, d, r0, sep, 100, 100.0, r, u, ur);
    std::cout << u0 << std::endl;

    double sum = calculate_sum(r, u, 1, 50.0);

    std::ofstream ofs("output_int_1d.txt");
    for(int i = 0; i < static_cast<int>(r.size()); i++){
      ofs << r.at(i) << "\t" << u.at(i) << "\t"
        << 2 * u.at(i) / sum << "\t" << sum << std::endl;
    }
    ofs.close();
  }

  {
    mpf_dec_float_prec r0 = r2;
    mpf_dec_float_prec d(2);

    std::vector<mpf_dec_float_prec> r;
    std::vector<mpf_dec_float_prec> u;

    double ur;

    // determine u0 such that u(r) -> 0 as r -> infty
    mpf_dec_float_prec u0 = search(s1, s2, sigma_sq, d, r0, sep, 100, 100.0, r, u, ur);
    std::cout << u0 << std::endl;

    double sum = calculate_sum(r, u, 2, 50.0);

    std::ofstream ofs("output_int_2d.txt");
    for(int i = 0; i < static_cast<int>(r.size()); i++){
      ofs << r.at(i) << "\t" << u.at(i) << "\t" <<
        2 * std::acos(-1.0) * r.at(i) * u.at(i) / sum << "\t" << sum << std::endl;
    }
    ofs.close();
  }

  {
    mpf_dec_float_prec r0 = r3;
    mpf_dec_float_prec d(3);

    std::vector<mpf_dec_float_prec> r;
    std::vector<mpf_dec_float_prec> u;

    double ur;

    // determine u0 such that u(r) -> 0 as r -> infty
    mpf_dec_float_prec u0 = search(s1, s2, sigma_sq, d, r0, sep, 100, 100.0, r, u, ur);
    std::cout << u0 << std::endl;

    double sum = calculate_sum(r, u, 3, 50.0);

    std::ofstream ofs("output_int_3d.txt");
    for(int i = 0; i < static_cast<int>(r.size()); i++){
      ofs << r.at(i) << "\t" << u.at(i) << "\t" <<
        4 * std::acos(-1.0) * r.at(i) * r.at(i) * u.at(i) / sum << "\t" << sum << std::endl;
    }
    ofs.close();
  }
}
