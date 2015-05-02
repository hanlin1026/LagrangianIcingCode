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
  Eigen::MatrixXd getX();      double getX(int ind);
  Eigen::MatrixXd getY();      double getY(int ind);
  Eigen::MatrixXf getRHO();    float  getRHO(int ind);
  Eigen::MatrixXf getU();      float  getU(int ind);
  Eigen::MatrixXf getV();      float  getV(int ind);
  Eigen::MatrixXf getE();      float  getE(int ind);
  Eigen::MatrixXd getXCENT();  double getXCENT(int ind);
  Eigen::MatrixXd getYCENT();  double getYCENT(int ind);
  Eigen::MatrixXd getLMIN();   double getLMIN(int ind);
  void getPROPS(FluidScalars& PROPS);
  int getNX();
  int getNY();
  // Cell metrics methods
  void computeCellAreas();
  void computeCellCenters();
  void computeGridMetrics();
  void transformXYtoIJ(int ind, Eigen::MatrixXd& xq, Eigen::MatrixXd& yq, Eigen::MatrixXd& Iq, Eigen::MatrixXd& Jq);
  void transformXYtoIJ(int ind, double xq, double yq, double Iq, double Jq);
  
 private:
  // Grid coordinates/solution
  Eigen::MatrixXd x_;
  Eigen::MatrixXd y_;
  Eigen::MatrixXf rho_; 
  Eigen::MatrixXf u_; 
  Eigen::MatrixXf v_; 
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
  Eigen::MatrixXf uCENT_; 
  Eigen::MatrixXf vCENT_; 
  Eigen::MatrixXf ECENT_;
  Eigen::MatrixXd cellArea_;
  Eigen::MatrixXd Jxx_;
  Eigen::MatrixXd Jxy_;
  Eigen::MatrixXd Jyx_;
  Eigen::MatrixXd Jyy_;
  Eigen::MatrixXd Lmin_;
  // Method to return 
};
#endif // PLOT3D_H
