#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <random>
#include "Grid/PLOT3D.h"
#include "QuadTree/Bucket.h"
#include "Cloud/Cloud.h"

// Driver program to test PLOT3D class

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
  PLOT3D* p3d = new PLOT3D("Grid/MESH.P3D", "Grid/q103.0.50E+01.bin", &scalars);
  float* PROPS = new float[6];
  p3d->getPROPS(PROPS);
  int nx = (int)PROPS[0]; int ny = (int)PROPS[1];
  float mach = PROPS[2]; float alpha = PROPS[3];
  float reynolds = PROPS[4]; float time = PROPS[5];
  printf("mach = %f\nalpha = %f\nreynolds = %f\ntime = %f\n", mach, alpha, reynolds, time);
  // Output solution file flow variable data
  int n = nx*ny;
  FILE* foutSoln = fopen("outputSOLN.dat","w");
  Eigen::MatrixXf RHO  = p3d->getRHO();
  Eigen::MatrixXf RHOU = p3d->getRHOU();
  Eigen::MatrixXf RHOV = p3d->getRHOV();
  Eigen::MatrixXf E    = p3d->getE();
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      fprintf(foutSoln,"%f \t %f \t %f \t %f \t \n", RHO(i,j), RHOU(i,j), RHOV(i,j), E(i,j));
    }
  }
  // Get grid x,y coordinates, output to file
  Eigen::MatrixXd X = p3d->getX();
  Eigen::MatrixXd Y = p3d->getY();
  ofstream foutXY; foutXY.open("outputGRID.dat");
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      foutXY << X(i,j) << "\t" << Y(i,j) << "\n";
    }
  }
  // Quadtree test
  double SW[2] = {-22.5,-32.5};
  double SE[2] = {18.5,-32.5};
  double NW[2] = {-22.5,33.0};
  double NE[2] = {18.5,33.0};
  Bucket* QT = new Bucket(&SW[0],&SE[0],&NW[0],&NE[0]);
  QT->calcQuadTree(X.data(),Y.data(),n);
  // Search for a query point
  FILE* foutKNN = fopen("outputKNNSEARCH.dat","w");
  default_random_engine generator;
  uniform_real_distribution<double> distX(-0.2,-0.1);
  uniform_real_distribution<double> distY(-0.1,0.1);
  double Xq, Yq, Xnn, Ynn, indnn;
  for (int i=0; i<1000; i++) {
    Xq = distX(generator);
    Yq = distY(generator);
    QT->knnSearch(&Xq,&Yq,&Xnn,&Ynn,&indnn);
    fprintf(foutKNN,"%f\t%f\t%f\t%f\t%f\n",Xq,Yq,Xnn,Ynn,indnn);
  }
  printf("KNN search successful. Data written to outputKNNSEARCH.dat.\n");

  // Test grid metric calculations
  FILE* foutTRANS = fopen("outputIJCOORDS.dat","w");
  double ind = (nx-1)*10 + 256;
  Eigen::MatrixXd I(nx,ny);
  Eigen::MatrixXd J(nx,ny);
  p3d->computeCellCenters();
  p3d->computeCellAreas();
  p3d->computeGridMetrics();
  p3d->transformXYtoIJ(ind,X,Y,I,J);
  for (int i=0; i<nx*ny; i++) {
    fprintf(foutTRANS,"%f\t%f\n",I(i),J(i));
  }
  printf("Jacobian transformation successful. Data written to outputIJCOORDS.dat.\n");

  // Initialize cloud of particles
  int particles = 1000;
  State state;
  state.x.resize(particles);
  state.y.resize(particles);
  state.u.resize(particles);
  state.v.resize(particles);
  state.r.resize(particles);
  state.temp.resize(particles);
  state.time.resize(particles);
  state.numDrop.resize(particles);
  double Rmean = 100e-6;
  double Tmean = 0;
  double pg, ug, vg;
  for (int i=0; i<particles; i++) {
    state.x(i) = distX(generator);
    state.y(i) = distY(generator);
    QT->knnSearch(&state.x(i),&state.y(i),&Xnn,&Ynn,&indnn);
    state.u(i) = p3d->getRHOU(indnn);
    state.v(i) = p3d->getRHOV(indnn);
    state.r(i) = Rmean;
    state.temp(i) = Tmean;
    state.time(i) = 0;
    state.numDrop(i) = 1;
  }
  Cloud cloud(state,*QT,scalars.rhol_);
  printf("Cloud initialization successful.\n");
  
  // Clear any allocated memory, close files/streams
  delete p3d;
  foutXY.close();
  fclose(foutSoln); fclose(foutKNN);
  delete[] PROPS, scalars;
  delete QT;
}
