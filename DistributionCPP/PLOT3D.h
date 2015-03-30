#ifndef PLOT3D_H_
#define PLOT3D_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using namespace std;

class PLOT3D {
 public:
  PLOT3D(const char *meshfname, const char *solnfname);
  ~PLOT3D();
  double* xy_;
  int nx_, ny_;
  float mach_, alpha_, reynolds_, time_;
  float* rho_;
  float* rhou_;
  float* rhov_;
  float* E_;
  //void getXY(int** XY);

 private:
  //int[] xy_;
  //int** X_;
  //int** Y_;
  //double rho_;
  //double rhoU_;
  //double rhoV_;
  //double E_;
  //double mach_;
  //double alpha_;
  //double reynolds_;
  //double time_;

};
#endif // PLOT3D_H
