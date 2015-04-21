#include "PLOT3D.h"
#include <stdio.h>
#include <assert.h>
#include <cmath>

using namespace std;
using namespace Eigen;

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
  for (int i=0; i<2*n; i++) {
    meshfile >> xy[i];
  }
  x_.resize(nx_,ny_);
  y_.resize(nx_,ny_);
  int iter = 0;
  for (int i=0; i<nx_; i++) {
    for (int j=0; j<ny_; j++) {
      x_(i,j) = xy[iter];
      y_(i,j) = xy[iter+n];
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
  rho_.resize(nx_,ny_);
  rhou_.resize(nx_,ny_);
  rhov_.resize(nx_,ny_);
  E_.resize(nx_,ny_);
  double rhoTMP[nx_*ny_];
  double rhouTMP[nx_*ny_];
  double rhovTMP[nx_*ny_];
  double ETMP[nx_*ny_];
  fread(&rhoTMP, sizeof(float), nx_*ny_, solnfile);
  fread(&rhouTMP, sizeof(float), nx_*ny_, solnfile);
  fread(&rhovTMP, sizeof(float), nx_*ny_, solnfile);
  fread(&ETMP, sizeof(float), nx_*ny_, solnfile);
  iter = 0;
  for (int i=0; i<nx_; i++) {
    for (int j=0; j<ny_; j++) {
      rho_(i,j) = rhoTMP[iter];
      rhou_(i,j) = rhouTMP[iter];
      rhov_(i,j) = rhovTMP[iter];
      E_(i,j) = ETMP[iter];
      iter++;
    }
  }
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
  
}

MatrixXd PLOT3D::getX() {
  return x_;
}
MatrixXd PLOT3D::getY() {
  return y_;
}
MatrixXf PLOT3D::getRHO() {
  return rho_;
}
MatrixXf PLOT3D::getRHOU() {
  return rhou_;
}
MatrixXf PLOT3D::getRHOV() {
  return rhov_;
}
MatrixXf PLOT3D::getE() {
  return E_;
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

  double DI_x[nx_-1][ny_-1];
  double DI_y[nx_-1][ny_-1];
  double DJ_x[nx_-1][ny_-1];
  double DJ_y[nx_-1][ny_-1];
  double DI, DJ;
  for (int i=0; i<nx_-1; i++) {
    for (int j=0; j<ny_-1; j++) {
      DI_x[i][j] = x_(i+1,j) - x_(i,j);
      DI_y[i][j] = y_(i+1,j) - y_(i,j);
      DJ_x[i][j] = x_(i,j+1) - x_(i,j);
      DJ_y[i][j] = y_(i,j+1) - y_(i,j);
    }
  }
  for (int i=0; i<nx_-1; i++) {
    for (int j=0; j<ny_-1; j++) {
      DI = sqrt(pow(DI_x[i][j],2) + pow(DI_y[i][j],2));
      DJ = sqrt(pow(DJ_x[i][j],2) + pow(DJ_y[i][j],2));
      cellArea_(i,j) = DI*DJ;
    }
  }
}

void PLOT3D::computeCellCenters() {
  // Function to compute cell centroids of the grid
  
  xCENT_.resize(nx_-1,ny_-1);
  yCENT_.resize(nx_-1,ny_-1);
  rhoCENT_.resize(nx_-1,ny_-1);
  rhouCENT_.resize(nx_-1,ny_-1);
  rhovCENT_.resize(nx_-1,ny_-1);
  ECENT_.resize(nx_-1,ny_-1);
  
  // Compute centroid locations and flow variables
  for (int i=0; i<nx_-1; i++) {
    for (int j=0; j<ny_-1; j++) {
      xCENT_(i,j) =    0.25*(x_(i,j) + x_(i+1,j) + x_(i,j+1) + x_(i+1,j+1));
      yCENT_(i,j) =    0.25*(y_(i,j) + y_(i+1,j) + y_(i,j+1) + y_(i+1,j+1));
      rhoCENT_(i,j) =  0.25*(rho_(i,j) + rho_(i+1,j) + rho_(i,j+1) + rho_(i+1,j+1));
      rhouCENT_(i,j) = 0.25*(rhou_(i,j) + rhou_(i+1,j) + rhou_(i,j+1) + rhou_(i+1,j+1));
      rhovCENT_(i,j) = 0.25*(rhov_(i,j) + rhov_(i+1,j) + rhov_(i,j+1) + rhov_(i+1,j+1));
      ECENT_(i,j) =    0.25*(E_(i,j) + E_(i+1,j) + E_(i,j+1) + E_(i+1,j+1));
    }
  }
}

void PLOT3D::computeGridMetrics() {
  // Function to compute grid metrics (Jacobians)

  // Compute Jacobians and minimum cell lengths (for use in CFL condition)
  double LI, LJ;
  for (int i=0; i<nx_-1; i++) {
    for (int j=0; j<ny_-1; j++) {
      Jxx_(i,j) = 0.5*( (x_(i+1,j+1)-x_(i,j+1)) + (x_(i+1,j)-x_(i,j)) );
      Jxy_(i,j) = 0.5*( (x_(i,j+1)-x_(i,j)) + (x_(i+1,j+1)-x_(i+1,j)) );
      Jyx_(i,j) = 0.5*( (y_(i+1,j+1)-y_(i,j+1)) + (y_(i+1,j)-y_(i,j)) );
      Jyy_(i,j) = 0.5*( (y_(i,j+1)-y_(i,j)) + (y_(i+1,j+1)-y_(i+1,j)) );
      LI = sqrt( pow(Jxx_(i,j),2) + pow(Jyx_(i,j),2) );
      LJ = sqrt( pow(Jxy_(i,j),2) + pow(Jyy_(i,j),2) );
      Lmin_(i,j) = std::min(LI,LJ);
    }
  }

}

void PLOT3D::transformXYtoIJ(int ind, Eigen::VectorXd* xq, Eigen::VectorXd* yq, Eigen::VectorXd* Iq, Eigen::VectorXd* Jq) {
  // Function to transform physical domain coordinates (x,y) to computational 
  // domain coordinates (I,J), centered at cell center 'ind'

  double xC = xCENT_(ind);
  double yC = yCENT_(ind);
  double area = cellArea_(ind);
  double Jxx = Jxx_(ind);
  double Jxy = Jxy_(ind);
  double Jyx = Jyx_(ind);
  double Jyy = Jyy_(ind);

  // Inverse Jacobian transformation
  double X, Y;
  for (int i=0; i<xq->size(); i++) {
    X = (*xq)(i) - xC;
    Y = (*yq)(i) - yC;
    (*Iq)(i) = (X*Jyy - Y*Jxy)/area;
    (*Jq)(i) = (-X*Jyx + Y*Jxx)/area;
  }

}
