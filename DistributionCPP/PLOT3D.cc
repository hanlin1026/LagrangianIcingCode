#include "PLOT3D.h"
#include <stdio.h>
#include <assert.h>

using namespace std;

PLOT3D::PLOT3D(const char *meshfname, const char *solnfname) {
  // Constructor: read in mesh/soln file data

  // Read in size of mesh
  int nx, ny, n;
  ifstream meshfile;
  meshfile.open(meshfname);
  meshfile >> nx; meshfile >> ny;
  nx_ = nx; ny_ = ny; n = nx*ny;
  // Read in mesh coordinates
  double* xy = new double[2*nx*ny];
  x_ = new double*[nx];
  y_ = new double*[nx];
  for (int i=0; i<nx; i++) {
    x_[i] = new double[ny];
    y_[i] = new double[ny];
  }
  for (int i=0; i<2*nx*ny; i++) {
    meshfile >> xy[i];
  }
  int iter = 0;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      x_[i][j] = xy[iter];
      y_[i][j] = xy[iter+n];
      iter++; 
    }
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
  float* rho = new float[n];
  float* rhou = new float[n];
  float* rhov = new float[n];
  float* E = new float[n];
  fread(rho, sizeof(float), n, solnfile);
  fread(rhou, sizeof(float), n, solnfile);
  fread(rhov, sizeof(float), n, solnfile);
  fread(E, sizeof(float), n, solnfile);
  // Reshape solution data vectors
  rho_ = new float*[nx];
  rhou_ = new float*[nx];
  rhov_ = new float*[nx];
  E_ = new float*[nx];
  for (int i=0; i<nx; i++) {
    rho_[i] = new float[ny];
    rhou_[i] = new float[ny];
    rhov_[i] = new float[ny];
    E_[i] = new float[ny];
  }
  iter = 0;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      rho_[i][j] = rho[iter];
      rhou_[i][j] = rhou[iter];
      rhov_[i][j] = rhov[iter];
      E_[i][j] = E[iter];
      iter++;
    }
  }
  // Close input streams/files, free any allocated memory
  delete[] xy, rho, rhou, rhov, E;
  meshfile.close();
  fclose(solnfile);
}

PLOT3D::~PLOT3D() {
  // Destructor: free any allocated memory

  for (int i=0; i<this->nx_; i++) {
    delete[] x_[i];
    delete[] y_[i];
    delete[] rho_[i];
    delete[] rhou_[i];
    delete[] rhov_[i];
    delete[] E_[i];
  }
  delete[] x_, y_, rho_, rhou_, rhov_, E_;
}

void PLOT3D::getXY(double**X, double** Y) {
  // Function to return grid coordinates x,y
  
  int nx = this->nx_;
  int ny = this->ny_;
  double** x = this->x_;
  double** y = this->y_;
  int n = nx*ny;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      X[i][j] = x[i][j];
      Y[i][j] = y[i][j];
    }
  }

}

void PLOT3D::getRHO(float** RHO) {
  // Function to return rho
  
  int nx = this->nx_;
  int ny = this->ny_;
  float** rho = this->rho_;
  int n = nx*ny;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      RHO[i][j] = rho[i][j];
    }
  }

}

void PLOT3D::getRHOU(float** RHOU) {
  // Function to return rhou
  
  int nx = this->nx_;
  int ny = this->ny_;
  float** rhou = this->rhou_;
  int n = nx*ny;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      RHOU[i][j] = rhou[i][j];
    }
  }

}

void PLOT3D::getRHOV(float** RHOV) {
  // Function to return rhov
  
  int nx = this->nx_;
  int ny = this->ny_;
  float** rhov = this->rhov_;
  int n = nx*ny;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      RHOV[i][j] = rhov[i][j];
    }
  }

}

void PLOT3D::getE(float** E) {
  // Function to return E
  
  int nx = this->nx_;
  int ny = this->ny_;
  float** e = this->E_;
  int n = nx*ny;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      E[i][j] = e[i][j];
    }
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
