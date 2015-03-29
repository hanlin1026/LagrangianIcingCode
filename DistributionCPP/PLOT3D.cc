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

  // Close input file
  meshfile.close(); solnfile.close();
}

PLOT3D::~PLOT3D() {
  delete[] xy_;

}
