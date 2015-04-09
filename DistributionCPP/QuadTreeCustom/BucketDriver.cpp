#include <stdio.h>
#include <stdlib.h>
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

  // Print four corners of first bucket
  printf("xSW = %f, ySW = %f \n",QT->SW_[0],QT->SW_[1]);
  printf("xSE = %f, ySE = %f \n",QT->SE_[0],QT->SE_[1]);
  printf("xNW = %f, yNW = %f \n",QT->NW_[0],QT->NW_[1]);
  printf("xNE = %f, yNE = %f \n",QT->NE_[0],QT->NE_[1]);

  // Test set/get points
  const int nrolls=10000;  // number of experiments

  default_random_engine generator;
  uniform_real_distribution<double> distX(0.,1);
  uniform_real_distribution<double> distY(0.,1.);
  
  double sampsX[nrolls];
  double sampsY[nrolls];
  for (int i=0; i<nrolls; i++) {
    double nx = distX(generator);
    double ny = distY(generator);
    sampsX[i] = nx;
    sampsY[i] = ny;    
  }
  
  /**
  vector<double> PX; vector<double> PY;
  QT->getPoints(&PX,&PY);
  int N = QT->getNPts();
  printf("N = %d \n", N);
  **/

  // Divide buckets
  QT->calcQuadTree(&sampsX[0],&sampsY[0],nrolls);
  
  // Delete allocated memory
  delete QT;

}
