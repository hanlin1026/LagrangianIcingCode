#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iterator>
#include <fstream>
#include "RBFInterpolant.h"

// Test driver program for reading input 

using namespace std;

int main(int argc, const char *argv[]) {
  int N = 5;
  int d = 2;
  Eigen::MatrixXd X(d,N);
  Eigen::VectorXd Y(N);
  X << 0,0,1,1,0.5,
    0,1,0,1,0.5;
  Y << 0,1,1,2,0;
  RBFInterpolant spline = RBFInterpolant(X,Y);
  spline.calcInterpolationWeights();
  Eigen::VectorXd Xq(2);
  // Interpolate on a grid
  int Xres = 10; double dx = 1.0/(Xres-1);
  int Yres = 10; double dy = 1.0/(Yres-1);
  std::vector<double> Xgrid(Xres*Yres);
  std::vector<double> Ygrid(Xres*Yres);
  std::vector<double> Zgrid(Xres*Yres);
  double Xtmp = 0;
  double Ytmp = 0;
  double zq;
  for (int i=0; i<Xres; i++) {
    Ytmp = 0;
    for (int j=0; j<Yres; j++) {
      Xq(0) = Xtmp; Xq(1) = Ytmp;
      Xgrid.push_back(Xtmp);
      Ygrid.push_back(Ytmp);
      zq = spline.evaluateInterpolant(Xq);
      Zgrid.push_back(zq);
      Ytmp += dy;
    }
    Xtmp += dx;
  }
  // Output to file
  ofstream foutX("Xgrid.out");
  ofstream foutY("Ygrid.out");
  ofstream foutZ("Zgrid.out");
  ostream_iterator<double> out_itX (foutX,"\n");
  copy ( Xgrid.begin(), Xgrid.end(), out_itX );
  ostream_iterator<double> out_itY (foutY,"\n");
  copy ( Ygrid.begin(), Ygrid.end(), out_itY );
  ostream_iterator<double> out_itZ (foutZ,"\n");
  copy ( Zgrid.begin(), Zgrid.end(), out_itZ );
}
