#ifndef PLOT3D_H_
#define PLOT3D_H_

#include <stdlib.h>

using namespace std;

class PLOT3D {
 public:
  PLOT3D(const std::string& meshfile, const std::string& solnfile);
  ~PLOT3D();
  void getXY(int** XY);

 private:
  int** X;
  int** Y;
  double rho;
  double rhoU;
  double rhoV;
  double E;
  double mach;
  double alpha;
  double reynolds;
  double time;

};
#endif // PLOT3D_H
