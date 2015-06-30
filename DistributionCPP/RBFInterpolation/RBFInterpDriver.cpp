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
  int N = 20;
  int d = 2;
  Eigen::MatrixXd X(d,N);
  Eigen::VectorXd Y(N);
  X << 0,0,0,0,.25,.25,.25,.25,.5,.5,.5,.5,.75,.75,.75,.75,1,1,1,1,
    0,.25,.5,.75,1,0,.25,.5,.75,1,0,.25,.5,.75,1,0,.25,.5,.75,1;
  Eigen::Vector2d Mean;
  Mean << 0.7,0.6;
  for (int i=0; i<N; i++) {
    Y(i) = exp(-((X.col(i)-Mean).squaredNorm()));
  }
  RBFInterpolant spline = RBFInterpolant(X,Y);
  spline.calcInterpolationWeights();
  Eigen::VectorXd Xq(2);
  // Interpolate on a grid
  int Xres = 10; double dx = 1.0/(Xres-1);
  int Yres = 10; double dy = 1.0/(Yres-1);
  std::vector<double> Xgrid;
  std::vector<double> Ygrid;
  std::vector<double> Zgrid;
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
