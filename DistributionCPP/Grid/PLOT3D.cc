#include "PLOT3D.h"
#include <stdio.h>
#include <assert.h>
#include <cmath>

using namespace std;

PLOT3D::PLOT3D(const char *meshfname, const char *solnfname, FluidScalars* scalars) {
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
  pinf_ =   scalars->pinf_;
  R_ =      scalars->R_;
  Tinf_ =   scalars->Tinf_;
  rhoinf_ = scalars->rhoinf_;
  Ubar_ =   scalars->Ubar_;
  rhol_ =   scalars->rhol_;
  Uinf_ = mach_*340;
  // Close input streams/files, free any allocated memory
  delete[] xy;
  meshfile.close();
  fclose(solnfile);
}

PLOT3D::~PLOT3D() {
  // Destructor: free any allocated memory
  
  delete[] x_, y_, rho_, rhou_, rhov_, E_;
  delete[] cellArea_;
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

void PLOT3D::computeCellAreas() {
  // Function to compute cell areas

  double* x = this->x_;
  double* y = this->y_;
  int nx = this->nx_; 
  int ny = this->ny_;
  double XX[nx][ny];
  double YY[nx][ny];
  double DI_x[nx-1][ny-1];
  double DI_y[nx-1][ny-1];
  double DJ_x[nx-1][ny-1];
  double DJ_y[nx-1][ny-1];
  double DI, DJ;
  int iter = 0;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      XX[i][j] = x[iter];
      YY[i][j] = y[iter];
      iter++;
    }
  }
  for (int i=0; i<nx-1; i++) {
    for (int j=0; j<ny-1; j++) {
      DI_x[i][j] = XX[i+1][j] - XX[i][j];
      DI_y[i][j] = YY[i+1][j] - YY[i][j];
      DJ_x[i][j] = XX[i][j+1] - XX[i][j];
      DJ_y[i][j] = YY[i][j+1] - YY[i][j];
    }
  }
  this->cellArea_ = new double[(nx-1)*(ny-1)];
  iter = 0;
  for (int i=0; i<nx-1; i++) {
    for (int j=0; j<ny-1; j++) {
      DI = sqrt(pow(DI_x[i][j],2) + pow(DI_y[i][j],2));
      DJ = sqrt(pow(DJ_x[i][j],2) + pow(DJ_y[i][j],2));
      cellArea_[iter] = DI*DJ;
      iter++;
    }
  }
  

}
