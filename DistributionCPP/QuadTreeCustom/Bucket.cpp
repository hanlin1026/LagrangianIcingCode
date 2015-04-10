#include "Bucket.h"
#include <stdlib.h>
#include <stdio.h>

using namespace std;

Bucket::Bucket(double* SW, double* SE, double* NW, double* NE) {
  SW_[0] = SW[0]; SW_[1] = SW[1];
  SE_[0] = SE[0]; SE_[1] = SE[1];
  NW_[0] = NW[0]; NW_[1] = NW[1];
  NE_[0] = NE[0]; NE_[1] = NE[1];
  buckets_ = new Bucket*[4];
  for (int i=0; i<4; i++) {
    buckets_[i] = NULL;
  }
  // Set default bucket size
  BucketSize_ = 20;
  // Set level to be zero (default)
  level_ = 0;
}

Bucket::Bucket(double* SW, double* SE, double* NW, double* NE, int level) {
  SW_[0] = SW[0]; SW_[1] = SW[1];
  SE_[0] = SE[0]; SE_[1] = SE[1];
  NW_[0] = NW[0]; NW_[1] = NW[1];
  NE_[0] = NE[0]; NE_[1] = NE[1];
  buckets_ = new Bucket*[4];
  for (int i=0; i<4; i++) {
    buckets_[i] = NULL;
  }
  // Set default bucket size
  BucketSize_ = 20;
  // Set level
  level_ = level;
}

Bucket::~Bucket() {
  for (int i=0; i<4; i++) {
    delete buckets_[i];
  }
  delete buckets_;
}

void Bucket::addNode(int ind, double* sw, double* se, double* nw, double* ne) {
  // Add a node to the bucket
  // 'ind' specifies which bucket is being added
  // [0,1,2,3] = [NE,NW,SW,SE]
  
  buckets_[ind] = new Bucket(sw,se,nw,ne,level_+1);
}

void Bucket::setPoints(double* dataX, double* dataY, int NumPts) {
  // Function to pass in a data set and determine the subset contained in bucket
  
  int count = 0;
  PX_.reserve(NumPts);
  PY_.reserve(NumPts);
  bool flag1, flag2, flag3, flag4, flag;
  for (int i=0; i<NumPts; i++) {
    flag1 = (dataX[i] > SW_[0]);
    flag2 = (dataX[i] < SE_[0]);
    flag3 = (dataY[i] > SW_[1]);
    flag4 = (dataY[i] < NW_[1]);

    flag = flag1 && flag2 && flag3 && flag4;
    if (flag==true) {
	PX_[count] = dataX[i];
	PY_[count] = dataY[i];
	count++;
    }
  }
  NumPts_ = count;

}

void Bucket::getPoints(std::vector<double>* PX, std::vector<double>* PY) {
  // Function to return points in the data set
  
  PX->reserve(NumPts_); PY->reserve(NumPts_);
  for (int i=0; i<NumPts_; i++) {
    PX->push_back(PX_[i]);
    PY->push_back(PY_[i]);
  }
}

int Bucket::getNPts() {
  // Return NumPts_

  return NumPts_;
}

int Bucket::getLevel() {
  // Return level

  return level_;
}

void Bucket::setBucketSize(int BS) {
  // Set bucket size

  BucketSize_ = BS;
}

void Bucket::divideBucket() {
  // Function to divide a bucket if the number of points inside of it exceeds threshold
  
  if (NumPts_ > BucketSize_) {
    // Calculate centroids
    double S[2]; S[0] = 0.5*(SW_[0]+SE_[0]); S[1] = 0.5*(SW_[1]+SE_[1]);
    double N[2]; N[0] = 0.5*(NW_[0]+NE_[0]); N[1] = 0.5*(NW_[1]+NE_[1]);
    double E[2]; E[0] = 0.5*(SE_[0]+NE_[0]); E[1] = 0.5*(SE_[1]+NE_[1]);
    double W[2]; W[0] = 0.5*(SW_[0]+NW_[0]); W[1] = 0.5*(SW_[1]+NW_[1]);
    double C[2]; C[0] = 0.5*(SW_[0]+NE_[0]); C[1] = 0.5*(SW_[1]+NE_[1]);

    // Add nodes as [0,1,2,3] = [NE,NW,SW,SE] buckets
    this->addNode(0,&C[0],&E[0],&N[0],&NE_[0]);
    this->addNode(1,&W[0],&C[0],&NW_[0],&N[0]);
    this->addNode(2,&SW_[0],&S[0],&W[0],&C[0]);
    this->addNode(3,&S[0],&SE_[0],&C[0],&E[0]);

    // Calculate data inside each of the new child buckets
    this->buckets_[0]->setPoints(&PX_[0],&PY_[0],NumPts_);
    this->buckets_[1]->setPoints(&PX_[0],&PY_[0],NumPts_);
    this->buckets_[2]->setPoints(&PX_[0],&PY_[0],NumPts_);
    this->buckets_[3]->setPoints(&PX_[0],&PY_[0],NumPts_);
  }
}

void Bucket::calcQuadTree(double* dataX, double* dataY, int NumPts) {
  // Function to handle the entire construction of the quadtree

  // Set points from dataset given
  this->setPoints(dataX,dataY,NumPts);
  // Recursive division of the quadtree until finest levels 
  // of resolution agree with bucket size
  vector<Bucket*> current;
  current.push_back(this);
  vector<Bucket*> next;
  int sizeCurrent = 1;
  int numNext = 0;
  bool flag = false;
  FILE* fout = fopen("QuadTreeXY.dat","w");
  Bucket* child;
  // Print initial bucket to file
  fprintf(fout,"%f\t%f\n", this->SW_[0], this->SW_[1]);
  fprintf(fout,"%f\t%f\n", this->SE_[0], this->SE_[1]);
  fprintf(fout,"%f\t%f\n", this->NE_[0], this->NE_[1]);
  fprintf(fout,"%f\t%f\n", this->NW_[0], this->NW_[1]);

  while(flag==false) {
    // Divide current
    for (int i=0; i<sizeCurrent; i++) {
      current[i]->divideBucket();
      // Check 4 children
      for (int j=0; j<4; j++) {
	child = current[i]->buckets_[j];
	if (child->getNPts() > BucketSize_) {
	  next.push_back(child);
	  numNext++;
	}
	// Print to file
	fprintf(fout,"%f\t%f\n", child->SW_[0], child->SW_[1]);
	fprintf(fout,"%f\t%f\n", child->SE_[0], child->SE_[1]);
	fprintf(fout,"%f\t%f\n", child->NE_[0], child->NE_[1]);
	fprintf(fout,"%f\t%f\n", child->NW_[0], child->NW_[1]);
      }
    }
    if (numNext>0) {
      // Reset current nodes
      current.swap(next);
      next.clear();
      sizeCurrent = numNext;
      numNext = 0;
    }
    else {
      // Exit; no more divisions needed
      printf("Quadtree successfully created. Bucket coordinates written to QuadTreeXY.dat. \n");
      flag = true;
      break;
    }
  
  }
  fclose(fout);
  
}
