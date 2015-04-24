#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <random>
#include "Grid/PLOT3D.h"
#include "QuadTree/Bucket.h"

// Driver program to test PLOT3D class

int main(int argc, const char *argv[]) {
  // Initialize scalars
  FluidScalars scalars;
  scalars.pinf_ = 0;
  scalars.R_ = 1;
  scalars.Tinf_ = 2;
  scalars.rhoinf_ = 3;
  scalars.Ubar_ = 4;
  scalars.rhol_ = 5;
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
  
  // Clear any allocated memory, close files/streams
  delete p3d;
  foutXY.close();
  fclose(foutSoln); fclose(foutKNN);
  delete[] PROPS, scalars;
  delete QT;
}
