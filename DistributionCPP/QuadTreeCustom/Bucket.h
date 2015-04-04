#ifndef __BUCKET_H__
#define __BUCKET_H__

#include <vector>

using namespace std;

class Bucket {
 public:
  Bucket(double* SW, double* SE, double* NW, double* NE);
  ~Bucket();
  void addNode(int ind, double sw, double se, double nw, double ne);
  double SW_[2], SE_[2], NW_[2], NE_[2];
  Bucket** buckets_;

};


#endif
