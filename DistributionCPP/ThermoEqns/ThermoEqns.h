#ifndef __THERMOEQNS_H__
#define __THERMOEQNS_H__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <Airfoil.h>

class ThermoEqns {
 public:
  ThermoEqns(const char* filenameCF,Airfoil& airfoil);
  ~ThermoEqns();
  void NewtonKrylovIteration(void (*f)(std::vector<double>& X),std::vector<double>& X,std::vector<double>& u0);

 private:
  // Functions to read in data files
  Eigen::MatrixXd readCHCF(const char* filenameCHCF);
  Eigen::MatrixXd readBetaXY(const char* filenameBeta);
  // Size of grid
  int NPts_;
  // Grid (s) and unknowns (film height, temperature, ice rate)
  std::vector<double> s_upper_;
  std::vector<double> hf_upper_;
  std::vector<double> ts_upper_;
  std::vector<double> mice_upper_;
  // Auxiliary parameters
  std::vector<double> cF_upper_;
  std::vector<double> cH_;

};

#endif
