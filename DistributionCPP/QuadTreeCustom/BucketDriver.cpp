#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <random>
#include "Bucket.h"

using namespace std;

// Driver program to test Bucket class

int main(int argc, const char *argv[]) {
  // Initialize first bucket
  double SW[2] = {0.,0.};
  double SE[2] = {1.,0.};
  double NW[2] = {0.,1.};
  double NE[2] = {1.,1.};
  Bucket* QT = new Bucket(&SW[0],&SE[0],&NW[0],&NE[0]);

  // Test set/get points
  const int nrolls=1000;  // number of experiments

  default_random_engine generator;
  normal_distribution<double> distX(0.3,.1);
  normal_distribution<double> distY(0.4,0.2);
  FILE* fout = fopen("DataSet.dat","w");
  
  double sampsX[nrolls];
  double sampsY[nrolls];
  for (int i=0; i<nrolls; i++) {
    double nx = distX(generator);
    double ny = distY(generator);
    sampsX[i] = nx;
    sampsY[i] = ny;   
    fprintf(fout,"%f\t%f\n",nx,ny);
  }

  // Divide buckets
  QT->calcQuadTree(&sampsX[0],&sampsY[0],nrolls);
  
  // Delete allocated memory
  delete QT;
  fclose(fout);

}
