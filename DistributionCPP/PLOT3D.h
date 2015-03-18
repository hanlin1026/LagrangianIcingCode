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
  PLOT3D(const std::string& meshfile, const std::string& solnfile);
  ~PLOT3D();
  int xy_[];
  //void getXY(int** XY);

 private:
  //int[] xy_;
  //int** X_;
  //int** Y_;
  double rho_;
  double rhoU_;
  double rhoV_;
  double E_;
  double mach_;
  double alpha_;
  double reynolds_;
  double time_;

};
#endif // PLOT3D_H