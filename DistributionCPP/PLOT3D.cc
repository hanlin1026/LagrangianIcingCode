#include "PLOT3D.h"
#include <assert.h>

using namespace std;

PLOT3D::PLOT3D(ifstream& meshfile, ifstream& solnfile) {
  // Read in size of mesh
  int nx, ny; 
  meshfile >> nx; meshfile >> ny;
  nx_ = nx; ny_ = ny;
  // Read in mesh coordinates
  xy_ = new double[2*nx*ny];
  double tmp;
  for (int i=0; i<2*nx*ny; i++) {
    meshfile >> tmp;
    xy_[i] = tmp;
  }
  // Read in soln data
  soln_ = new double[10];
  char buff[4];
  solnfile.read(buff,2);
  solnfile.read(buff,4); mach_ = atof(buff);
  solnfile.read(buff,4); alpha_ = atof(buff);
  solnfile.read(buff,4); reynolds_ = atof(buff);
  solnfile.read(buff,4); time_ = atof(buff);
  /**
  for (int i=0; i<10; i++) {
    solnfile.read(buff,4);
    soln_[i] = atof(buff);
  }
  **/
  // Close input file
  meshfile.close(); solnfile.close();
}

PLOT3D::~PLOT3D() {
  delete[] xy_;
  delete[] soln_;
}
