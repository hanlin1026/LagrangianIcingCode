#include <stdio.h>
#include <stdlib.h>
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
  SW[0] = 0.; SW[1] = 0.; SE[0] = 0.5; SE[1] = 0.; 
  NW[0] = 0.; NW[1] = 0.5; NE[0] = 0.5; NE[1] = 0.5;
  QT->addNode(0,&SW[0],&SE[0],&NW[0],&NE[0]);
  
  // Print four corners of child node
  printf("xSW = %f, ySW = %f \n",QT->buckets_[0]->SW_[0],QT->buckets_[0]->SW_[1]);
  printf("xSE = %f, ySE = %f \n",QT->buckets_[0]->SE_[0],QT->buckets_[0]->SE_[1]);
  printf("xNW = %f, yNW = %f \n",QT->buckets_[0]->NW_[0],QT->buckets_[0]->NW_[1]);
  printf("xNE = %f, yNE = %f \n",QT->buckets_[0]->NE_[0],QT->buckets_[0]->NE_[1]);
  
  delete[] QT;

}
