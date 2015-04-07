#include "Bucket.h"

using namespace std;

Bucket::Bucket(double* SW, double* SE, double* NW, double* NE)
{
  SW_[0] = SW[0]; SW_[1] = SW[1];
  SE_[0] = SE[0]; SE_[1] = SE[1];
  NW_[0] = NW[0]; NW_[1] = NW[1];
  NE_[0] = NE[0]; NE_[1] = NE[1];
  buckets_ = new Bucket*[4];
}

Bucket::~Bucket() {
  for (int i=0; i<4; i++) {
    delete[] buckets_[i];
  }
  delete[] buckets_;
}

void Bucket::addNode(int ind, double* sw, double* se, double* nw, double* ne) {
  // Add a node to the bucket
  // 'ind' specifies which bucket is being added
  // [1,2,3,4] = [SW,SE,NW,NE]

  buckets_[ind] = new Bucket(sw,se,nw,ne);
}
