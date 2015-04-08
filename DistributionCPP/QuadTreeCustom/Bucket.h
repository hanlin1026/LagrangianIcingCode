#ifndef __BUCKET_H__
#define __BUCKET_H__

#include <vector>

class Bucket {
 public:
  Bucket(double* SW, double* SE, double* NW, double* NE);
  ~Bucket();
  double SW_[2], SE_[2], NW_[2], NE_[2];
  Bucket** buckets_;
  void addNode(int ind, double* sw, double* se, double* nw, double* ne);
  void setPoints(double* dataX, double* dataY, int NumPts);
  void getPoints(std::vector<double>* PX, std::vector<double>* PY);
  int getNPts();
  void setBucketSize(int BS);
  void divideBucket();

 private:
  std::vector<double> PX_;
  std::vector<double> PY_;
  int NumPts_;
  int BucketSize_;
};


#endif
