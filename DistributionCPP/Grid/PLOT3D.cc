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
  for (int j=0; j<ny_; j++) {
    for (int i=0; i<nx_; i++) {
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
  // Read in scalars
  pinf_ =   scalars->pinf_;
  R_ =      scalars->R_;
  Tinf_ =   scalars->Tinf_;
  rhoinf_ = scalars->rhoinf_;
  Ubar_ =   scalars->Ubar_;
  rhol_ =   scalars->rhol_;
  Uinf_ = mach_*340;
  // Read in solution data
  rho_.resize(nx_,ny_);
  u_.resize(nx_,ny_);
  v_.resize(nx_,ny_);
  E_.resize(nx_,ny_);
  float rhoTMP[nx_*ny_];
  float uTMP[nx_*ny_];
  float vTMP[nx_*ny_];
  float ETMP[nx_*ny_];
  fread(&rhoTMP, sizeof(float), nx_*ny_, solnfile);
  fread(&uTMP, sizeof(float), nx_*ny_, solnfile);
  fread(&vTMP, sizeof(float), nx_*ny_, solnfile);
  fread(&ETMP, sizeof(float), nx_*ny_, solnfile);
  iter = 0;
  // Normalize solution data
  for (int j=0; j<ny_; j++) {
    for (int i=0; i<nx_; i++) {
      rho_(i,j) = rhoinf_*rhoTMP[iter];
      u_(i,j) = Ubar_*rhoinf_*uTMP[iter]/rho_(i,j);
      v_(i,j) = Ubar_*rhoinf_*vTMP[iter]/rho_(i,j);
      E_(i,j) = ETMP[iter];
      iter++;
    }
  }
  // Initialize searcher object
  this->computeCellCenters();
  this->computeCellAreas();
  this->computeGridMetrics();
  this->createQuadTree();
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
MatrixXf PLOT3D::getU() {
  return u_;
}
MatrixXf PLOT3D::getV() {
  return v_;
}
MatrixXf PLOT3D::getE() {
  return E_;
}
MatrixXd PLOT3D::getXCENT() {
  return xCENT_;
}
MatrixXd PLOT3D::getYCENT() {
  return yCENT_;
}
MatrixXd PLOT3D::getLMIN() {
  return Lmin_;
}
double PLOT3D::getX(int ind) {
  return x_(ind);
}
double PLOT3D::getY(int ind) {
  return y_(ind);
}
float PLOT3D::getRHO(int ind) {
  return rho_(ind);
}
float PLOT3D::getU(int ind) {
  return u_(ind);
}
float PLOT3D::getV(int ind) {
  return v_(ind);
}
float PLOT3D::getE(int ind) {
  return E_(ind);
}
double PLOT3D::getXCENT(int ind) {
  return xCENT_(ind);
}
double PLOT3D::getYCENT(int ind) {
  return yCENT_(ind);
}
double PLOT3D::getRHOCENT(int ind) {
  return rhoCENT_(ind);
}
double PLOT3D::getUCENT(int ind) {
  return uCENT_(ind);
}
double PLOT3D::getVCENT(int ind) {
  return vCENT_(ind);
}
double PLOT3D::getLMIN(int ind) {
  return Lmin_(ind);
}

void PLOT3D::getPROPS(FluidScalars& PROPS) {
  // Function to return [nx,ny,mach,alpha,reynolds,time]

  PROPS.nx_ = (float)nx_;
  PROPS.ny_ = (float)ny_;
  PROPS.mach_ = mach_;
  PROPS.alpha_ = alpha_;
  PROPS.reynolds_ = reynolds_;
  PROPS.time_ = time_;
  
}

int PLOT3D::getNX() {
  return nx_;
}
int PLOT3D::getNY() {
  return ny_;
}
double PLOT3D::getTINF() {
  return Tinf_;
}

void PLOT3D::computeCellAreas() {
  // Function to compute cell areas

  MatrixXd DI_x(nx_-1,ny_-1);
  MatrixXd DI_y(nx_-1,ny_-1);
  MatrixXd DJ_x(nx_-1,ny_-1);
  MatrixXd DJ_y(nx_-1,ny_-1);
  double DI, DJ;
  for (int j=0; j<ny_-1; j++) {
    for (int i=0; i<nx_-1; i++) {
      DI_x(i,j) = x_(i+1,j) - x_(i,j);
      DI_y(i,j) = y_(i+1,j) - y_(i,j);
      DJ_x(i,j) = x_(i,j+1) - x_(i,j);
      DJ_y(i,j) = y_(i,j+1) - y_(i,j);
    }
  }
  cellArea_.resize(nx_-1,ny_-1);
  for (int j=0; j<ny_-1; j++) {
    for (int i=0; i<nx_-1; i++) {
      DI = sqrt(pow(DI_x(i,j),2) + pow(DI_y(i,j),2));
      DJ = sqrt(pow(DJ_x(i,j),2) + pow(DJ_y(i,j),2));
      cellArea_(i,j) = DI*DJ;
    }
  }
}

void PLOT3D::computeCellCenters() {
  // Function to compute cell centroids of the grid
  
  xCENT_.resize(nx_-1,ny_-1);
  yCENT_.resize(nx_-1,ny_-1);
  rhoCENT_.resize(nx_-1,ny_-1);
  uCENT_.resize(nx_-1,ny_-1);
  vCENT_.resize(nx_-1,ny_-1);
  ECENT_.resize(nx_-1,ny_-1);
  
  // Compute centroid locations and flow variables
  for (int j=0; j<ny_-1; j++) {
    for (int i=0; i<nx_-1; i++) {
      xCENT_(i,j) =    0.25*(x_(i,j) + x_(i+1,j) + x_(i,j+1) + x_(i+1,j+1));
      yCENT_(i,j) =    0.25*(y_(i,j) + y_(i+1,j) + y_(i,j+1) + y_(i+1,j+1));
      rhoCENT_(i,j) =  0.25*(rho_(i,j) + rho_(i+1,j) + rho_(i,j+1) + rho_(i+1,j+1));
      uCENT_(i,j) = 0.25*(u_(i,j) + u_(i+1,j) + u_(i,j+1) + u_(i+1,j+1));
      vCENT_(i,j) = 0.25*(v_(i,j) + v_(i+1,j) + v_(i,j+1) + v_(i+1,j+1));
      ECENT_(i,j) =    0.25*(E_(i,j) + E_(i+1,j) + E_(i,j+1) + E_(i+1,j+1));
    }
  }
}

void PLOT3D::computeGridMetrics() {
  // Function to compute grid metrics (Jacobians)

  // Compute Jacobians and minimum cell lengths (for use in CFL condition)
  double LI, LJ;
  Jxx_.resize(nx_-1,ny_-1);
  Jxy_.resize(nx_-1,ny_-1);
  Jyx_.resize(nx_-1,ny_-1);
  Jyy_.resize(nx_-1,ny_-1);
  Lmin_.resize(nx_-1,ny_-1);
  for (int j=0; j<ny_-1; j++) {
    for (int i=0; i<nx_-1; i++) {
      Jxx_(i,j) = 0.5*( (x_(i+1,j+1) - x_(i,j+1)) + (x_(i+1,j)   - x_(i,j))   );
      Jxy_(i,j) = 0.5*( (x_(i,j+1)   - x_(i,j))   + (x_(i+1,j+1) - x_(i+1,j)) );
      Jyx_(i,j) = 0.5*( (y_(i+1,j+1) - y_(i,j+1)) + (y_(i+1,j)   - y_(i,j))   );
      Jyy_(i,j) = 0.5*( (y_(i,j+1)   - y_(i,j))   + (y_(i+1,j+1) - y_(i+1,j)) );
      LI = sqrt( pow(Jxx_(i,j),2) + pow(Jyx_(i,j),2) );
      LJ = sqrt( pow(Jxy_(i,j),2) + pow(Jyy_(i,j),2) );
      Lmin_(i,j) = std::min(LI,LJ);
    }
  }

}

void PLOT3D::transformXYtoIJ(int ind, Eigen::MatrixXd& xq, Eigen::MatrixXd& yq, Eigen::MatrixXd& Iq, Eigen::MatrixXd& Jq) {
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
  int n = xq.rows()*xq.cols();
  for (int i=0; i<n; i++) {
    X = xq(i) - xC;
    Y = yq(i) - yC;
    Iq(i) = (X*Jyy - Y*Jxy)/area;
    Jq(i) = (-X*Jyx + Y*Jxx)/area;
  }

}

void PLOT3D::transformXYtoIJ(int ind, double xq, double yq, double& Iq, double& Jq) {
  // Function to transform a single query point from physical domain coordinates (x,y)
  // to computational domain coordinates (I,J), centered at cell center 'ind'
  
  double xC = xCENT_(ind);
  double yC = yCENT_(ind);
  double area = cellArea_(ind);
  double Jxx = Jxx_(ind);
  double Jxy = Jxy_(ind);
  double Jyx = Jyx_(ind);
  double Jyy = Jyy_(ind);

  // Inverse Jacobian transformation
  double X = xq - xC;
  double Y = yq - yC;
  Iq = (X*Jyy - Y*Jxy)/area;
  Jq = (-X*Jyx + Y*Jxx)/area;

}

void PLOT3D::createQuadTree() {
  // Function to create a quadtree searcher of the grid cell centers
  
  // Set bounds of quadtree object
  double minX,minY,maxX,maxY;
  minX = xCENT_.minCoeff(); maxX = xCENT_.maxCoeff();
  minY = yCENT_.minCoeff(); maxY = yCENT_.maxCoeff();
  double SW[2] = {minX, minY};
  double SE[2] = {maxX, minY};
  double NW[2] = {minX, maxY};
  double NE[2] = {maxX, maxY};
  QT_.setBounds(&SW[0],&SE[0],&NW[0],&NE[0]);
  // Create quadtree search object
  QT_.calcQuadTree(xCENT_.data(),yCENT_.data(),(nx_-1)*(ny_-1));

}

void PLOT3D::pointSearch(double xq, double yq, double& xnn, double& ynn, int& indnn) {
  // Function to search the quadtree for a query point
  
  QT_.knnSearch(&xq,&yq,&xnn,&ynn,&indnn);

}
