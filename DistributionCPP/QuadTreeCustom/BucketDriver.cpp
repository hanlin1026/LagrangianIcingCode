#include <stdio.h>
#include <stdlib.h>
#include "Bucket.h"

using namespace std;

// Driver program to test Bucket class

int main(int argc, const char *argv[]) {
  // Initialize first bucket
  double SW[2] = {0,0};
  double SE[2] = {1,0};
  double NW[2] = {0,1};
  double NE[2] = {1,1};
  Bucket* QT = new Bucket(&SW[2],&SE[2],&NW[2],&NE[2]);

  printf("xSW = %d, ySW = %d \n",QT->SW_[1],QT->SW_[2]);

}
