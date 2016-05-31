#include "Bucket.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>

using namespace std;

Bucket::Bucket() {
  // Default constructor
  
  buckets_ = new Bucket*[4];
  for (int i=0; i<4; i++) {
    buckets_[i] = NULL;
  }
  // Set default bucket size
  BucketSize_ = 20;
  // Set level to be zero (default)
  level_ = 0;
}

Bucket::Bucket(const std::string workDir) {
  // Constructor + set output working directory
  
  workDir_ = workDir;
  buckets_ = new Bucket*[4];
  for (int i=0; i<4; i++) {
    buckets_[i] = NULL;
  }
  // Set default bucket size
  BucketSize_ = 20;
  // Set level to be zero (default)
  level_ = 0;
}

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

void Bucket::setBounds(double* SW, double* SE, double* NW, double* NE) {
  // Function to set bounds of bucket (assuming it has not been done in constructor)

  SW_[0] = SW[0]; SW_[1] = SW[1];
  SE_[0] = SE[0]; SE_[1] = SE[1];
  NW_[0] = NW[0]; NW_[1] = NW[1];
  NE_[0] = NE[0]; NE_[1] = NE[1];

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

void Bucket::setPoints(double* dataX, double* dataY, vector<int>& indData, int NumPts) {
  // Function to pass in a data set and determine the subset contained in bucket
  
  int count = 0;
  PX_.reserve(NumPts);
  PY_.reserve(NumPts);
  indData_.reserve(NumPts);
  bool flag1, flag2, flag3, flag4, flag;
  for (int i=0; i<NumPts; i++) {
    flag1 = (dataX[i] >= SW_[0]);
    flag2 = (dataX[i] <= SE_[0]);
    flag3 = (dataY[i] >= SW_[1]);
    flag4 = (dataY[i] <= NW_[1]);

    flag = flag1 && flag2 && flag3 && flag4;
    if (flag==true) {
	PX_[count] = dataX[i];
	PY_[count] = dataY[i];
        indData_.push_back(indData[i]);
	count++;
    }
  }
  NumPts_ = count;

}

void Bucket::getPoints(vector<double>* PX, vector<double>* PY, vector<int>* indXY) {
  // Function to return points in the data set and their indices
  
  PX->reserve(NumPts_); PY->reserve(NumPts_); indXY->reserve(NumPts_);
  for (int i=0; i<NumPts_; i++) {
    PX->push_back(PX_[i]);
    PY->push_back(PY_[i]);
    indXY->push_back(indData_[i]);
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
    this->buckets_[0]->setPoints(&PX_[0],&PY_[0],indData_,NumPts_);
    this->buckets_[1]->setPoints(&PX_[0],&PY_[0],indData_,NumPts_);
    this->buckets_[2]->setPoints(&PX_[0],&PY_[0],indData_,NumPts_);
    this->buckets_[3]->setPoints(&PX_[0],&PY_[0],indData_,NumPts_);
  }
}

void Bucket::calcQuadTree(double* dataX, double* dataY, int NumPts) {
  // Function to handle the entire construction of the quadtree
  
  // Initialize indices of mother bucket
  vector<int> indData;
  indData.reserve(NumPts);
  for (int i=0; i<NumPts; i++) {
    indData.push_back(i);
  }
  // Set points from dataset given
  this->setPoints(dataX,dataY,indData,NumPts);
  // Recursive division of the quadtree until finest levels 
  // of resolution agree with bucket size
  vector<Bucket*> current;
  current.push_back(this);
  vector<Bucket*> next;
  int sizeCurrent = 1;
  int numNext = 0;
  bool flag = false;
  FILE* fout;
  if (!workDir_.empty()) {
    const std::string outFile = workDir_ + "/QuadTreeXY.dat";
    fout = fopen(outFile.c_str(),"w");
  }
  Bucket* child;
  // Print initial bucket to file
  if (!workDir_.empty()) {
    fprintf(fout,"%f\t%f\n", this->SW_[0], this->SW_[1]);
    fprintf(fout,"%f\t%f\n", this->SE_[0], this->SE_[1]);
    fprintf(fout,"%f\t%f\n", this->NE_[0], this->NE_[1]);
    fprintf(fout,"%f\t%f\n", this->NW_[0], this->NW_[1]);
  }
  
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
        if (!workDir_.empty()) {
	  fprintf(fout,"%f\t%f\n", child->SW_[0], child->SW_[1]);
	  fprintf(fout,"%f\t%f\n", child->SE_[0], child->SE_[1]);
	  fprintf(fout,"%f\t%f\n", child->NE_[0], child->NE_[1]);
	  fprintf(fout,"%f\t%f\n", child->NW_[0], child->NW_[1]);
        }
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
      //printf("Quadtree successfully created. Bucket coordinates written to QuadTreeXY.dat. \n");
      flag = true;
      break;
    }
  
  }
  if (!workDir_.empty()) {
    fclose(fout);
  }
  
}

bool Bucket::calcInBucket(double* Xq, double* Yq) {
  // Function to determine whether a query point is inside a bucket

  bool flag1, flag2, flag3, flag4, flag;
  flag1 = (*Xq >= this->SW_[0]);
  flag2 = (*Xq <= this->SE_[0]);
  flag3 = (*Yq >= this->SW_[1]);
  flag4 = (*Yq <= this->NW_[1]);
  flag = flag1 && flag2 && flag3 && flag4;
  
  return flag;
}

void Bucket::knnSearch(double* Xq, double* Yq, double* Xnn, double* Ynn, int* indnn) {
  // Function that takes a query point and finds the nearest 
  // neighbor in the data set of the quadtree

  bool flagFinal = false;
  bool flagChild = false;
  bool flagHasPts;
  Bucket* current = this;
  int iter;
  // Iteratively search down the tree
  while(flagFinal==false) {
    // Determine if current bucket has child iter
    iter = -1;
    while((flagChild==false) && (iter<3)) {
      iter++;
      if (current->buckets_[iter] != NULL) {
        // If child iter exists, search it and make sure there are pts inside it
        flagChild = current->buckets_[iter]->calcInBucket(Xq,Yq);
        flagHasPts = current->buckets_[iter]->NumPts_ > 0;
      }
    }
    // If we have found a valid child containing (Xq,Yq), reset current
    if ((flagChild==true) && (flagHasPts==true)) {
      current = current->buckets_[iter];
      flagChild = false;
    }
    // Otherwise, we are done
    else {
      flagFinal = true;
    }

  }

  // Calculate nearest neighbor inside bucket
  int BS = BucketSize_;
  double* dist = new double[BS];
  vector<double> x;
  vector<double> y;
  vector<int> ind;
  current->getPoints(&x,&y,&ind);
  for (int i=0; i<BS; i++) {
    dist[i] = pow(x[i]-*Xq,2) + pow(y[i]-*Yq,2);
  }
  int indMin = 0;
  double distMin = dist[0];
  for (int i=1; i<current->NumPts_; i++) {
    if (dist[i] < distMin) {
      distMin = dist[i];
      indMin = i;
    }
  }
  *Xnn = x[indMin];
  *Ynn = y[indMin];
  *indnn = ind[indMin];

  delete[] dist;

}

void Bucket::setOutDir(const std::string workDir) {
  workDir_ = workDir;
}
