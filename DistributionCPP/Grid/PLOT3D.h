#ifndef PLOT3D_H_
#define PLOT3D_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <eigen3/Eigen/Dense>
#include "FluidScalars.h"

using namespace std;

class PLOT3D {
 public:
  // Constructor: read in mesh/soln
  PLOT3D(const char *meshfname, const char *solnfname, FluidScalars* scalars);
  ~PLOT3D();
  // Get methods
  Eigen::MatrixXd getX();
  Eigen::MatrixXd getY();
  Eigen::MatrixXf getRHO();
  Eigen::MatrixXf getRHOU();
  Eigen::MatrixXf getRHOV();
  Eigen::MatrixXf getE();
  void getPROPS(float* PROPS);
  // Cell metrics methods
  void computeCellAreas();
  void computeCellCenters();
  void computeGridMetrics();
  void transformXYtoIJ(int ind, Eigen::VectorXd* xq, Eigen::VectorXd* yq, Eigen::VectorXd* Iq, Eigen::VectorXd* Jq);

 private:
  // Grid coordinates/solution
  Eigen::MatrixXd x_;
  Eigen::MatrixXd y_;
  Eigen::MatrixXf rho_; 
  Eigen::MatrixXf rhou_; 
  Eigen::MatrixXf rhov_; 
  Eigen::MatrixXf E_;
  // Properties of the grid/soln
  int nx_, ny_;
  float mach_, alpha_, reynolds_, time_;
  double pinf_, R_, Tinf_, rhoinf_, Ubar_, rhol_;
  double Uinf_;
  // Properties for the dual grid (grid of centroid locations)
  Eigen::MatrixXd xCENT_;
  Eigen::MatrixXd yCENT_;
  Eigen::MatrixXf rhoCENT_; 
  Eigen::MatrixXf rhouCENT_; 
  Eigen::MatrixXf rhovCENT_; 
  Eigen::MatrixXf ECENT_;
  Eigen::MatrixXd cellArea_;
  Eigen::MatrixXd Jxx_;
  Eigen::MatrixXd Jxy_;
  Eigen::MatrixXd Jyx_;
  Eigen::MatrixXd Jyy_;
  Eigen::MatrixXd Lmin_;
};
#endif // PLOT3D_H
