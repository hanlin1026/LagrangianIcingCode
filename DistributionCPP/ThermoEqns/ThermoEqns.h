#ifndef __THERMOEQNS_H__
#define __THERMOEQNS_H__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <Airfoil/Airfoil.h>
#include <Grid/FluidScalars.h>
#include <Cloud/Cloud.h>
#include <Grid/PLOT3D.h>

class ThermoEqns {
 public:
  ThermoEqns(const char* filenameCF,const char* filenameBETA,Airfoil& airfoil,FluidScalars& fluid,Cloud& cloud,PLOT3D& p3d,const char* strSurf);
  ~ThermoEqns();
  std::vector<double> NewtonKrylovIteration(const char* balance,std::vector<double>& u0,double globaltol);
  std::vector<double> trapz(std::vector<double>& X, std::vector<double>& Y);
  // Mass/Energy balance equations
  std::vector<double> JX(int func, std::vector<double>& X, std::vector<double>& u0);
  std::vector<double> massBalance(std::vector<double>& X);
  std::vector<double> energyBalance(std::vector<double>& Y);
  std::vector<double> testBalance(std::vector<double>& X);
  std::vector<double> SolveThermoForIceRate(std::vector<double>& X, std::vector<double>& Y);
  std::vector<double> integrateMassEqn(bool& C_filmHeight);
  std::vector<double> explicitSolver(const char* balance, std::vector<double>& y0, double eps, double tol);
  void SolveIcingEqns();
  void computeMevap();
  void computePstat(PLOT3D& p3d);
  // Set/get routines
  void setHF(std::vector<double>& hf);
  void setTS(std::vector<double>& ts);
  void setMICE(std::vector<double>& mice);
  std::vector<double> getS();
  std::vector<double> getMICE();

 private:
  // Functions to read in data files
  Eigen::MatrixXd readCHCF(const char* filenameCHCF);
  Eigen::MatrixXd readBetaXY(const char* filenameBeta);
  void interpUpperSurface(const char* filenameBeta, Airfoil& airfoil, const char* parameter);
  // Size of grid
  int NPts_;
  // Upper or lower surface string
  const char* strSurf_;
  // Grid (s), unknowns (film height, temperature, ice rate, impinging water mass), and other parameter functions of s
  std::vector<double> s_;
  std::vector<double> hf_;
  std::vector<double> ts_;
  std::vector<double> mice_;
  std::vector<double> mevap_;
  std::vector<double> pstat_;
  // Auxiliary parameters
  std::vector<double> cF_;
  std::vector<double> cH_;
  std::vector<double> beta_;
  double rhoL_, muL_, LWC_, Uinf_;
  double Td_, cW_, ud_, cICE_, Lfus_;
  double rhoINF_, pINF_, TINF_;
  double chord_;
  double cpAir_;
  // PLOT3D parameters (to assist in interpolation calculations)
  int indFirst_,indLast_; // First/last indices for where interpolation starts/ends for parameters taken from PLOT3D grid
  std::vector<double> sP3D_; // S-coordinates of airfoil wrap in P3D variables 

};

#endif
