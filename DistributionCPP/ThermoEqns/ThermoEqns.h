#ifndef __THERMOEQNS_H__
#define __THERMOEQNS_H__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <eigen3/Eigen/Dense>

class ThermoEqns {
 public:
  ThermoEqns(const char* filenameCF);
  ~ThermoEqns();

 private:
  // Functions to read in data files
  Eigen::MatrixXd readSkinFrictionCoeff(const char* fileSkinFriction);
  Eigen::MatrixXd readHeatFlux(const char* fileHeatFlux);
  // Size of grid
  int NPTS_;
  // Grid (s) and unknowns (film height, temperature, ice rate)
  std::vector<double> s_;
  std::vector<double> filmHeight_;
  std::vector<double> temperature_;
  std::vector<double> iceRate_;
  // Auxiliary parameters
  std::vector<double> cF_;
  std::vector<double> cH_;

};

#endif
