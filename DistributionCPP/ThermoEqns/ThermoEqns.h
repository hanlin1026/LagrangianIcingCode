#ifndef __THERMOEQNS_H__
#define __THERMOEQNS_H__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <Airfoil/Airfoil.h>
#include <Grid/FluidScalars.h>

class ThermoEqns {
 public:
  ThermoEqns(const char* filenameCF,const char* filenameBETA,Airfoil& airfoil,FluidScalars& fluid);
  ~ThermoEqns();
  std::vector<double> NewtonKrylovIteration(const char* balance,std::vector<double>& u0,double globaltol);
  std::vector<double> trapz(std::vector<double>& X, std::vector<double>& Y);
  // Mass/Energy balance equations
  std::vector<double> JX(int func, std::vector<double>& X, std::vector<double>& u0);
  std::vector<double> massBalanceUpper(std::vector<double>& X);
  std::vector<double> energyBalanceUpper(std::vector<double>& Y);
  std::vector<double> testBalance(std::vector<double>& X);
  std::vector<double> SolveThermoForIceRate(std::vector<double>& X, std::vector<double>& Y, const char* surface);
  std::vector<double> integrateMassEqnUpper();
  void SolveIcingEqns();
  // Set/get routines
  void setHF_upper(std::vector<double>& hf);
  void setTS_upper(std::vector<double>& ts);
  void setMICE_upper(std::vector<double>& mice);
  std::vector<double> getS_upper();

 private:
  // Functions to read in data files
  Eigen::MatrixXd readCHCF(const char* filenameCHCF);
  Eigen::MatrixXd readBetaXY(const char* filenameBeta);
  void interpUpperSurface(const char* filenameBeta, Airfoil& airfoil, const char* parameter);
  // Size of grid
  int NPts_;
  // Grid (s) and unknowns (film height, temperature, ice rate, impinging water mass)
  std::vector<double> s_upper_;
  std::vector<double> hf_upper_;
  std::vector<double> ts_upper_;
  std::vector<double> mice_upper_;
  // Auxiliary parameters
  std::vector<double> cF_upper_;
  std::vector<double> cH_upper_;
  std::vector<double> beta_upper_;
  double rhoL_, muL_, LWC_, Uinf_;
  double Td_, cW_, ud_, cICE_, Lfus_;
  double rhoINF_, pINF_, TINF_;

};

#endif
