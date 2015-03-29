#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "PLOT3D.h"

// Driver program to test PLOT3D class

int main(int argc, const char *argv[]) {
  ifstream meshfile, solnfile;
  meshfile.open("MESH.P3D");
  solnfile.open("q103.0.50E+01.bin");
  // Initialize plot3D object
  PLOT3D* p3d;
  p3d = new PLOT3D(meshfile,solnfile);
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
  delete p3d;
  fout.close();
}
