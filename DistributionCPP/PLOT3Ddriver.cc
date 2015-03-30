#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "PLOT3D.h"

// Driver program to test PLOT3D class

int main(int argc, const char *argv[]) {
  // Initialize plot3D object
  PLOT3D* p3d = new PLOT3D("MESH.P3D", "q103.0.50E+01.bin");
  // Output mach,alpha,reynolds,time as a test
  printf("mach = %f\nalpha = %f\nreynolds = %f\ntime = %f\n", p3d->mach_, p3d->alpha_, p3d->reynolds_, p3d->time_);
  // Output solution file flow variable data
  ofstream foutSoln;
  foutSoln.open("outputSOLN.dat");
  float* tmp = p3d->rho_;
  int nx = p3d->nx_;
  int ny = p3d->ny_;
  int n = nx*ny;
  for (int i=0; i<n; i++) {
    foutSoln << tmp[i] << "\n";
  }
  /*
  // Output xy coordinates as a test
  ofstream fout;
  fout.open("outputXY.dat");
  double tmp;
  int nx = p3d->nx_;
  int ny = p3d->ny_;
  int sizeXY = 2*nx*ny;
  for (int i=0; i<sizeXY; i++) {
    tmp = p3d->xy_[i];
    fout << tmp << "\n";
  }
  // Output mach,alpha,reynolds,time from solnfile as a test
  ofstream foutSoln;
  foutSoln.open("outputSOLN.dat");
  foutSoln << p3d->mach_ << p3d->alpha_ << p3d->reynolds_ << p3d->time_;
  */
  /**
  for (int i=0; i<10; i++) {
    foutSoln << p3d->soln_[i];
  }
  fout.close(); 
  **/
  delete p3d;
  foutSoln.close();
}
