#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <random>
#include "Grid/PLOT3D.h"
#include "QuadTree/Bucket.h"
#include "Cloud/Cloud.h"
#include "Airfoil/Airfoil.h"

// Airfoil icing code driver program

int main(int argc, const char *argv[]) {
  // Initialize scalars
  FluidScalars scalars;
  scalars.pinf_ = 1.01325e5;
  scalars.R_ = 287.058;
  scalars.Tinf_ = 300;
  scalars.rhoinf_ = scalars.pinf_/scalars.R_/scalars.Tinf_;
  scalars.Ubar_ = sqrt(1.4*scalars.pinf_/scalars.rhoinf_);
  scalars.rhol_ = 1000;
  // Initialize plot3D object, read in basic problem data
  PLOT3D p3d = PLOT3D("Grid/MESH.P3D", "Grid/q103.0.50E+01.bin", &scalars);
  FluidScalars PROPS;
  p3d.getPROPS(PROPS);
  int nx = PROPS.nx_; int ny = PROPS.ny_;
  float mach = PROPS.mach_; float alpha = PROPS.alpha_;
  float reynolds = PROPS.reynolds_; float time = PROPS.time_;
  printf("mach = %f\nalpha = %f\nreynolds = %f\ntime = %f\n", mach, alpha, reynolds, time);
  // Initialize quadtree search object
  double SW[2] = {-22.5,-32.5};
  double SE[2] = {18.5,-32.5};
  double NW[2] = {-22.5,33.0};
  double NE[2] = {18.5,33.0};
  int n = nx*ny;
  p3d.computeCellCenters();
  p3d.computeCellAreas();
  p3d.computeGridMetrics();
  Eigen::MatrixXd xCENT = p3d.getXCENT();
  Eigen::MatrixXd yCENT = p3d.getYCENT();
  Bucket* QT = new Bucket(&SW[0],&SE[0],&NW[0],&NE[0]);
  QT->calcQuadTree(xCENT.data(),yCENT.data(),n);
  // Initialize cloud of particles
  int particles = 10;
  State state = State(particles);
  double Rmean = 100e-6;
  double Tmean = 273.15;
  double pg, ug, vg, Xnn, Ynn, indnn;
  // Search for a query point
  default_random_engine generator;
  uniform_real_distribution<double> distX(-5.1,-5);
  uniform_real_distribution<double> distY(-0.1,0.1);
  for (int i=0; i<particles; i++) {
    state.x_(i) = distX(generator);
    state.y_(i) = distY(generator);
    QT->knnSearch(&state.x_(i),&state.y_(i),&Xnn,&Ynn,&indnn);
    state.u_(i) = p3d.getU(indnn);
    state.v_(i) = p3d.getV(indnn);
    state.r_(i) = Rmean;
    state.temp_(i) = Tmean;
    state.time_(i) = 0;
    state.numDrop_(i) = 1;
  }
  Cloud cloud(state,*QT,scalars.rhol_);
  printf("Cloud initialization successful.\n");
  // Intialize airfoil object
  Eigen::MatrixXd Xgrid = p3d.getX();
  Eigen::MatrixXd Ygrid = p3d.getY();
  Eigen::VectorXd X(Xgrid.rows());
  Eigen::VectorXd Y(Ygrid.rows());
  int iter = 0;
  for (int i=0; i<Xgrid.rows(); i++) {
    if (Xgrid(i,0) <= 1) {
      X(iter) = Xgrid(i,0);
      Y(iter) = Ygrid(i,0);
      iter++;
    }
  }
  X = X.block(0,0,iter,1);
  Y = Y.block(0,0,iter,1);
  Airfoil airfoil = Airfoil(X,Y);
  
  // Clear any allocated memory, close files/streams
  delete QT;
}
