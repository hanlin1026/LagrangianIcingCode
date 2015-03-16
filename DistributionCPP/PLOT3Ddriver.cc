#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "PLOT3D.h"

// Driver program to test PLOT3D class

int main(int argc, const char *argv[]) {
  string meshfile, solnfile;
  meshfile = "MESH.P3D";
  solnfile = "q103.0.50E+01.bin";
  // Initialize plot3D object
  PLOT3D* p3d;
  p3d = new PLOT3D(meshfile,solnfile);
  // Output xy coordinates as a test
  FILE* fout;
  string outfile = "outputXY.dat";
  fout = fopen(outfile.c_str(),"w");
  double tmp;
  int i=0;
  while (!fout.eof()) {
    tmp = p3d.xy_[i];
    fprintf(fout,"%d \n",tmp);
    i++;
  }
  fclose(fout);
}
