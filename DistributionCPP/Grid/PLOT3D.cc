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
  delete[] xCENT_, yCENT_;
  delete[] rhoCENT_, rhouCENT_, rhovCENT_, ECENT_;
  delete[] cellArea_;
}

void PLOT3D::getXY(double* X, double* Y) {
  // Function to return grid coordinates x,y
  
  int n = nx_*ny_;
  for (int i=0; i<n; i++) {
    X[i] = x_[i];
    Y[i] = y_[i];
  }

}

void PLOT3D::getRHO(float* RHO) {
  // Function to return rho
  
  int n = nx_*ny_;
  for (int i=0; i<n; i++) {
    RHO[i] = rho_[i];
  }

}

void PLOT3D::getRHOU(float* RHOU) {
  // Function to return rhou
  
  int n = nx_*ny_;
  for (int i=0; i<n; i++) {
    RHOU[i] = rhou_[i];
  }

}

void PLOT3D::getRHOV(float* RHOV) {
  // Function to return rhov
  
  int n = nx_*ny_;
  for (int i=0; i<n; i++) {
    RHOV[i] = rhov_[i];
  }

}

void PLOT3D::getE(float* E) {
  // Function to return E
  
  int n = nx_*ny_;
  for (int i=0; i<n; i++) {
    E[i] = E_[i];
  }

}

void PLOT3D::getPROPS(float* PROPS) {
  // Function to return [nx,ny,mach,alpha,reynolds,time]

  PROPS[0] = (float)nx_;
  PROPS[1] = (float)ny_;
  PROPS[2] = mach_;
  PROPS[3] = alpha_;
  PROPS[4] = reynolds_;
  PROPS[5] = time_;
  
}

void PLOT3D::computeCellAreas() {
  // Function to compute cell areas

  double XX[nx_][ny_];
  double YY[nx_][ny_];
  double DI_x[nx_-1][ny_-1];
  double DI_y[nx_-1][ny_-1];
  double DJ_x[nx_-1][ny_-1];
  double DJ_y[nx_-1][ny_-1];
  double DI, DJ;
  // Convert vectorized grid into 2D array
  int iter = 0;
  for (int i=0; i<nx_; i++) {
    for (int j=0; j<ny_; j++) {
      XX[i][j] = x_[iter];
      YY[i][j] = y_[iter];
      iter++;
    }
  }
  for (int i=0; i<nx_-1; i++) {
    for (int j=0; j<ny_-1; j++) {
      DI_x[i][j] = XX[i+1][j] - XX[i][j];
      DI_y[i][j] = YY[i+1][j] - YY[i][j];
      DJ_x[i][j] = XX[i][j+1] - XX[i][j];
      DJ_y[i][j] = YY[i][j+1] - YY[i][j];
    }
  }
  cellArea_ = new double[(nx_-1)*(ny_-1)];
  iter = 0;
  for (int i=0; i<nx_-1; i++) {
    for (int j=0; j<ny_-1; j++) {
      DI = sqrt(pow(DI_x[i][j],2) + pow(DI_y[i][j],2));
      DJ = sqrt(pow(DJ_x[i][j],2) + pow(DJ_y[i][j],2));
      cellArea_[iter] = DI*DJ;
      iter++;
    }
  }
}

void PLOT3D::computeCellCenters() {
  // Function to compute cell centroids of the grid
  
  double XX[nx_][ny_];
  double YY[nx_][ny_];
  float RHO[nx_][ny_];
  float RHOU[nx_][ny_];
  float RHOV[nx_][ny_];
  float E[nx_][ny_];
  xCENT_ = new double[(nx_-1)*(ny_-1)];
  yCENT_ = new double[(nx_-1)*(ny_-1)];
  rhoCENT_ = new float[(nx_-1)*(ny_-1)];
  rhouCENT_ = new float[(nx_-1)*(ny_-1)];
  rhovCENT_ = new float[(nx_-1)*(ny_-1)];
  ECENT_ = new float[(nx_-1)*(ny_-1)];
  
  // Convert vectorized grid/soln into 2D array
  int iter = 0;
  for (int i=0; i<nx_; i++) {
    for (int j=0; j<ny_; j++) {
      XX[i][j] = x_[iter];
      YY[i][j] = y_[iter];
      RHO[i][j] = rho_[iter];
      RHOU[i][j] = rhou_[iter];
      RHOV[i][j] = rhov_[iter];
      E[i][j] = E_[iter];
      iter++;
    }
  }
  // Compute centroid locations and flow variables
  iter = 0;
  for (int i=0; i<nx_-1; i++) {
    for (int j=0; j<ny_-1; j++) {
      xCENT_[iter] =    0.25*(XX[i][j] + XX[i+1][j] + XX[i][j+1] + XX[i+1][j+1]);
      yCENT_[iter] =    0.25*(YY[i][j] + YY[i+1][j] + YY[i][j+1] + YY[i+1][j+1]);
      rhoCENT_[iter] =  0.25*(RHO[i][j] + RHO[i+1][j] + RHO[i][j+1] + RHO[i+1][j+1]);
      rhouCENT_[iter] = 0.25*(RHOU[i][j] + RHOU[i+1][j] + RHOU[i][j+1] + RHOU[i+1][j+1]);
      rhovCENT_[iter] = 0.25*(RHOV[i][j] + RHOV[i+1][j] + RHOV[i][j+1] + RHOV[i+1][j+1]);
      ECENT_[iter] =    0.25*(E[i][j] + E[i+1][j] + E[i][j+1] + E[i+1][j+1]);
      iter++;
    }
  }
}
