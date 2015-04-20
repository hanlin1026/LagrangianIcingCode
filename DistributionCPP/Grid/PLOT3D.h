#ifndef PLOT3D_H_
#define PLOT3D_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "FluidScalars.h"

using namespace std;

class PLOT3D {
 public:
  // Constructor: read in mesh/soln
  PLOT3D(const char *meshfname, const char *solnfname, FluidScalars* scalars);
  ~PLOT3D();
  // Get methods
  void getXY(double* X, double* Y);
  void getRHO(float* RHO);
  void getRHOU(float* RHOU);
  void getRHOV(float* RHOV);
  void getE(float* E);
  void getPROPS(float* PROPS);

 private:
  // Grid coordinates/solution
  double* x_;
  double* y_;
  float* rho_;
  float* rhou_;
  float* rhov_;
  float* E_;
  // Size of the grid
  int nx_, ny_;
  float mach_, alpha_, reynolds_, time_;
  double pinf_, R_, Tinf_, rhoinf_, Ubar_, rhol_;
  double Uinf_;
  double* cellArea_;
  // Functions for grid metrics
  void computeCellAreas();
};
#endif // PLOT3D_H
