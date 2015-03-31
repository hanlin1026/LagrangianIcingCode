#include "PLOT3D.h"
#include <stdio.h>
#include <assert.h>

using namespace std;

PLOT3D::PLOT3D(const char *meshfname, const char *solnfname, double* scalars) {
  // Constructor: read in mesh/soln file data

  // Read in size of mesh
  int nx, ny, n;
  ifstream meshfile;
  meshfile.open(meshfname);
  meshfile >> nx; meshfile >> ny;
  nx_ = nx; ny_ = ny; n = nx*ny;
  // Read in mesh coordinates
  double* xy = new double[2*n];
  x_ = new double[n];
  y_ = new double[n];
  for (int i=0; i<2*n; i++) {
    meshfile >> xy[i];
  }
  for (int i=0; i<n; i++) {
      x_[i] = xy[i];
      y_[i] = xy[i+n];
  }
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
  rho_ = new float[n];
  rhou_ = new float[n];
  rhov_ = new float[n];
  E_ = new float[n];
  fread(rho_, sizeof(float), n, solnfile);
  fread(rhou_, sizeof(float), n, solnfile);
  fread(rhov_, sizeof(float), n, solnfile);
  fread(E_, sizeof(float), n, solnfile);
  // Read in scalars
  pinf_ =   scalars[0];
  R_ =      scalars[1];
  Tinf_ =   scalars[2];
  rhoinf_ = scalars[3];
  Ubar_ =   scalars[4];
  rhol_ =   scalars[5];
  Uinf_ = mach_*340;
  // Close input streams/files, free any allocated memory
  delete[] xy;
  meshfile.close();
  fclose(solnfile);
}

PLOT3D::~PLOT3D() {
  // Destructor: free any allocated memory
  
  delete[] x_, y_, rho_, rhou_, rhov_, E_;
}

void PLOT3D::getXY(double* X, double* Y) {
  // Function to return grid coordinates x,y
  
  int nx = this->nx_;
  int ny = this->ny_;
  double* x = this->x_;
  double* y = this->y_;
  int n = nx*ny;
  for (int i=0; i<n; i++) {
    X[i] = x[i];
    Y[i] = y[i];
  }

}

void PLOT3D::getRHO(float* RHO) {
  // Function to return rho
  
  int nx = this->nx_;
  int ny = this->ny_;
  float* rho = this->rho_;
  int n = nx*ny;
  for (int i=0; i<n; i++) {
    RHO[i] = rho[i];
  }

}

void PLOT3D::getRHOU(float* RHOU) {
  // Function to return rhou
  
  int nx = this->nx_;
  int ny = this->ny_;
  float* rhou = this->rhou_;
  int n = nx*ny;
  for (int i=0; i<n; i++) {
    RHOU[i] = rhou[i];
  }

}

void PLOT3D::getRHOV(float* RHOV) {
  // Function to return rhov
  
  int nx = this->nx_;
  int ny = this->ny_;
  float* rhov = this->rhov_;
  int n = nx*ny;
  for (int i=0; i<n; i++) {
    RHOV[i] = rhov[i];
  }

}

void PLOT3D::getE(float* E) {
  // Function to return E
  
  int nx = this->nx_;
  int ny = this->ny_;
  float* e = this->E_;
  int n = nx*ny;
  for (int i=0; i<n; i++) {
    E[i] = e[i];
  }

}

void PLOT3D::getPROPS(float* PROPS) {
  // Function to return [nx,ny,mach,alpha,reynolds,time]
  
  int nx = this->nx_;
  int ny = this->ny_;
  float mach = this->mach_;
  float alpha = this->alpha_;
  float reynolds = this->reynolds_;
  float time = this->time_;

  PROPS[0] = (float)nx;
  PROPS[1] = (float)ny;
  PROPS[2] = mach;
  PROPS[3] = alpha;
  PROPS[4] = reynolds;
  PROPS[5] = time;
  
}
