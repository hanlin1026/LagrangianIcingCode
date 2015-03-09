#include "PLOT3D.h"
#include <assert.h>

using namespace std;

PLOT3D::PLOT3D(const std::string& meshfile, const std::string& solnfile) {
  // Initialize input file stream to read in meshfile
  FILE* finput;
  finput = fopen(meshfile.c_str(),"r");
  if (finput == NULL) {
    cout << "Error, didn't load meshfile" << endl;
    exit(1);
  }
  // Read in size of mesh
  int nx, ny; 
  fscan(finput,"%i",&nx);
  fscan(finput,"%i",&ny);
  // Read in mesh coordinates
  xy = int[nx*ny];
  int i = 0;
  while (!finput.eof()) {
    fscan(finput,"%i",&tmp);
    xy[i] = tmp;
    i++;
  }
  xy_ = xy;

  // Close input file
  fclose(finput);
}

PLOT3D::~PLOT3D() {
  delete[] xy_;

}
