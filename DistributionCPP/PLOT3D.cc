#include "PLOT3D.h"
#include <stdio.h>
#include <assert.h>

using namespace std;

PLOT3D::PLOT3D(const char *meshfname, const char *solnfname) {
  // Read in size of mesh
  int nx, ny;
  ifstream meshfile;
  meshfile.open(meshfname);
  meshfile >> nx; meshfile >> ny;
  nx_ = nx; ny_ = ny;
  // Read in mesh coordinates
  xy_ = new double[2*nx*ny];
  double tmp;
  for (int i=0; i<2*nx*ny; i++) {
    meshfile >> tmp;
    xy_[i] = tmp;
  }
  meshfile.close();
  // Read in mach,alpha,reynolds,time
  FILE *solnfile = fopen(solnfname, "r");
  assert(solnfile != NULL);
  char buff[BUFSIZ];
  fread(buff, sizeof(int), 2, solnfile); // First 2 ints are just nx,ny
  fread(&mach_, sizeof(float), 1, solnfile);
  fread(&alpha_, sizeof(float), 1, solnfile);
  fread(&reynolds_, sizeof(float), 1, solnfile);
  fread(&time_, sizeof(float), 1, solnfile);
  // Read in solution data
  int n = nx*ny;
  rho_ = new float[n];
  rhou_ = new float[n];
  rhov_ = new float[n];
  E_ = new float[n];
  fread(rho_, sizeof(float), n, solnfile);
  fread(rhou_, sizeof(float), n, solnfile);
  fread(rhov_, sizeof(float), n, solnfile);
  fread(E_, sizeof(float), n, solnfile);
  /**
  for (int i=0; i<10; i++) {
    solnfile.read(buff,4);
    soln_[i] = atof(buff);
  }
  **/
  // Close input file
  fclose(solnfile);
}

PLOT3D::~PLOT3D() {
  delete[] xy_;
  delete[] rho_, rhou_, rhov_, E_;
}
