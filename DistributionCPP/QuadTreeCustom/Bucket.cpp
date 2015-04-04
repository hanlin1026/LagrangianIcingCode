#include "Bucket.h"

using namespace std;

Bucket::Bucket(double* SW, double* SE, double* NW, double* NE):
  SW_[1](SW[1]), SW_[2](SW[2]),
  SE_[1](SE[1]), SE_[2](SE[2]),
  NW_[1](NW[1]), NW_[2](NW[2]),
  NE_[1](NE[1]), NE_[2](NE[2])

{
  buckets_** = new Bucket*[4];
}

Bucket::~Bucket() {
  for (int i=0; i<4; i++) {
    delete[] buckets_[i];
  }
  delete[] buckets_;
}

void Bucket::addNode(int ind, double sw, double se, double nw, double ne) {
  // Add a node to the bucket
  // 'ind' specifies which bucket is being added
  // [1,2,3,4] = [SW,SE,NW,NE]

  buckets_[ind] = Bucket(sw,se,nw,ne);
}
