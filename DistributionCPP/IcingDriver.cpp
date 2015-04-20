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
  // Initialize plot3D object
  PLOT3D* p3d = new PLOT3D("Grid/MESH.P3D", "Grid/q103.0.50E+01.bin", &scalars);
  // Output mach,alpha,reynolds,time as a test
  float* PROPS = new float[6];
  p3d->getPROPS(PROPS);
  int nx = (int)PROPS[0]; int ny = (int)PROPS[1];
  float mach = PROPS[2]; float alpha = PROPS[3];
  float reynolds = PROPS[4]; float time = PROPS[5];
  printf("mach = %f\nalpha = %f\nreynolds = %f\ntime = %f\n", mach, alpha, reynolds, time);
  // Output solution file flow variable data
  int n = nx*ny;
  FILE* foutSoln = fopen("outputSOLN.dat","w");
  float* RHO = new float[n];
  float* RHOU = new float[n];
  float* RHOV = new float[n];
  float* E = new float[n];
  p3d->getRHO(RHO);
  p3d->getRHOU(RHOU);
  p3d->getRHOV(RHOV);
  p3d->getE(E);
  for (int i=0; i<n; i++) {
    fprintf(foutSoln,"%f \t %f \t %f \t %f \t \n", RHO[i], RHOU[i], RHOV[i], E[i]);
  }
  // Get grid x,y coordinates
  double* X = new double[n];
  double* Y = new double[n];
  p3d->getXY(X,Y);
  // Output grid x,y coordinates to file
  ofstream foutXY;
  foutXY.open("outputXY.dat");
  for (int i=0; i<n; i++) {
    foutXY << X[i] << "\t" << Y[i] << "\n";
  }

  // Quadtree test
  double SW[2] = {-22.5,-32.5};
  double SE[2] = {18.5,-32.5};
  double NW[2] = {-22.5,33.0};
  double NE[2] = {18.5,33.0};
  Bucket* QT = new Bucket(&SW[0],&SE[0],&NW[0],&NE[0]);
  QT->calcQuadTree(&X[0],&Y[0],n);

  // Search for a query point
  default_random_engine generator;
  uniform_real_distribution<double> distX(-0.2,-0.1);
  uniform_real_distribution<double> distY(-0.1,0.1);
  double Xq, Yq, Xnn, Ynn;
  for (int i=0; i<1000; i++) {
    Xq = distX(generator);
    Yq = distY(generator);
    QT->knnSearch(&Xq,&Yq,&Xnn,&Ynn);
    //printf("Xq = %f, Yq = %f\nXnn = %f, Ynn = %f\n",Xq,Yq,Xnn,Ynn);
  }
  
  // Clear any allocated memory, close files/streams
  delete p3d;
  foutXY.close();
  fclose(foutSoln);
  delete[] X,Y,RHO,RHOU,RHOV,E;
  delete[] PROPS, scalars;
  delete QT;
}
