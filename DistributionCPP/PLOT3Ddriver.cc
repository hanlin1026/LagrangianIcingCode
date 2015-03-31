#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "PLOT3D.h"

// Driver program to test PLOT3D class

int main(int argc, const char *argv[]) {
  // Initialize scalars
  double* scalars = new double[6];
  scalars[0] = 0; scalars[1] = 1; scalars[2] = 2;
  scalars[3] = 3; scalars[4] = 4; scalars[5] = 5;
  // Initialize plot3D object
  PLOT3D* p3d = new PLOT3D("MESH.P3D", "q103.0.50E+01.bin", scalars);
  // Output mach,alpha,reynolds,time as a test
  float* PROPS = new float[6];
  p3d->getPROPS(PROPS);
  int nx = (int)PROPS[0]; int ny = (int)PROPS[1];
  float mach = PROPS[2]; float alpha = PROPS[3]; float reynolds = PROPS[4]; float time = PROPS[5];
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
  
  // Clear any allocated memory, close files/streams
  delete p3d;
  foutXY.close();
  fclose(foutSoln);
  delete[] X,Y,RHO,RHOU,RHOV,E;
  delete[] PROPS, scalars;
}
