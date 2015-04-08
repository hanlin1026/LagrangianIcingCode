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

  // Add 4 subbuckets
  double SW2[2] = {0,0};
  double SE2[2] = {0.5,0};
  double NW2[2] = {0,0.5};
  double NE2[2] = {0.5,0.5};
  QT->addNode(0,&SW2[0],&SE2[0],&NW2[0],&NE2[0]);
  
  // Print four corners of child node
  printf("xSW = %f, ySW = %f \n",QT->buckets_[0]->SW_[0],QT->buckets_[0]->SW_[1]);
  printf("xSE = %f, ySE = %f \n",QT->buckets_[0]->SE_[0],QT->buckets_[0]->SE_[1]);
  printf("xNW = %f, yNW = %f \n",QT->buckets_[0]->NW_[0],QT->buckets_[0]->NW_[1]);
  printf("xNE = %f, yNE = %f \n",QT->buckets_[0]->NE_[0],QT->buckets_[0]->NE_[1]);

  // Test set/get points
  const int nrolls=100;  // number of experiments

  default_random_engine generator;
  uniform_real_distribution<double> distX(-0.5,1.5);
  uniform_real_distribution<double> distY(-0.5,1.5);
  
  double sampsX[nrolls];
  double sampsY[nrolls];
  for (int i=0; i<nrolls; i++) {
    double nx = distX(generator);
    double ny = distY(generator);
    sampsX[i] = nx;
    sampsY[i] = ny;    
  }
  
  QT->setPoints(&sampsX[0],&sampsY[0],nrolls);
  vector<double> PX; vector<double> PY;
  QT->getPoints(&PX,&PY);
  int N = QT->getNPts();
  printf("N = %d \n", N);
  for (int i=0; i<N; i++) {
    printf("%f\t%f\n",PX[i],PY[i]);
  }

  // Divide buckets
  Bucket* bPoint = QT->buckets_[1];
  if (bPoint) {
    printf("Bucket has children %f \n",bPoint->SE_[1]);
  }
  else {
    printf("Empty Pointer \n");
  }

  QT->divideBucket();
  for (int i=0; i<4; i++) {
    printf("CHILD NODE %d:\n", i);
    printf("xSW = %f, ySW = %f \n",QT->buckets_[i]->SW_[0],QT->buckets_[i]->SW_[1]);
    printf("xSE = %f, ySE = %f \n",QT->buckets_[i]->SE_[0],QT->buckets_[i]->SE_[1]);
    printf("xNW = %f, yNW = %f \n",QT->buckets_[i]->NW_[0],QT->buckets_[i]->NW_[1]);
    printf("xNE = %f, yNE = %f \n",QT->buckets_[i]->NE_[0],QT->buckets_[i]->NE_[1]);
  }
  
  // Delete allocated memory
  delete QT;

}
