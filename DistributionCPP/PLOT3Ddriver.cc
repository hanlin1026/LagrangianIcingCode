#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "PLOT3D.h"

// Driver program to test PLOT3D class

int main(int argc, const char *argv[]) {
  // Initialize plot3D object
  PLOT3D* p3d = new PLOT3D("MESH.P3D", "q103.0.50E+01.bin");
  // Output mach,alpha,reynolds,time as a test
  float* PROPS = new float[6];
  p3d->getPROPS(PROPS);
  int nx = (int)PROPS[0]; int ny = (int)PROPS[1];
  float mach = PROPS[2]; float alpha = PROPS[3]; float reynolds = PROPS[4]; float time = PROPS[5];
  printf("mach = %f\nalpha = %f\nreynolds = %f\ntime = %f\n", mach, alpha, reynolds, time);
  // Output solution file flow variable data
  int n = nx*ny;
  FILE* foutSoln = fopen("outputSOLN.dat","w");
  float** RHO = new float*[nx];
  float** RHOU = new float*[nx];
  float** RHOV = new float*[nx];
  float** E = new float*[nx];
  for (int i=0; i<nx; i++) {
    RHO[i] = new float[ny];
    RHOU[i] = new float[ny];
    RHOV[i] = new float[ny];
    E[i] = new float[ny];
  }
  p3d->getRHO(RHO);
  p3d->getRHOU(RHOU);
  p3d->getRHOV(RHOV);
  p3d->getE(E);
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      fprintf(foutSoln,"%f \t %f \t %f \t %f \t \n", RHO[i][j], RHOU[i][j], RHOV[i][j], E[i][j]);
    }
  }
  // Get grid x,y coordinates
  double** X = new double*[nx];
  double** Y = new double*[nx];
  for (int i=0; i<nx; i++) {
    X[i] = new double[ny];
    Y[i] = new double[ny];
  }
  p3d->getXY(X,Y);
  // Output grid x,y coordinates to file
  ofstream foutXY;
  foutXY.open("outputXY.dat");
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      foutXY << X[i][j] << "\t" << Y[i][j] << "\n";
    }
  }
  
  // Clear any allocated memory, close files/streams
  delete p3d;
  foutXY.close();
  fclose(foutSoln);
  for (int i=0; i<nx; i++) {
    delete[] X[i];
    delete[] Y[i];
    delete[] RHO[i]; delete[] RHOU[i]; delete[] RHOV[i]; delete[] E[i];
  }
  delete[] X,Y,RHO,RHOU,RHOV,E;
  delete[] PROPS;
}
