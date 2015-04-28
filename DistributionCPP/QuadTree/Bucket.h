#ifndef __BUCKET_H__
#define __BUCKET_H__

#include <vector>

class Bucket {
 public:
  Bucket();
  Bucket(double* SW, double* SE, double* NW, double* NE);
  Bucket(double* SW, double* SE, double* NW, double* NE, int level);
  ~Bucket();
  double SW_[2], SE_[2], NW_[2], NE_[2];
  Bucket** buckets_;
  void setBounds(double* SW, double* SE, double* NW, double* NE);
  void addNode(int ind, double* sw, double* se, double* nw, double* ne);
  void setPoints(double* dataX, double* dataY, std::vector<double>& indData, int NumPts);
  void getPoints(std::vector<double>* PX, std::vector<double>* PY, std::vector<double>* indXY);
  int getNPts();
  int getLevel();
  void setBucketSize(int BS);
  void divideBucket();
  void calcQuadTree(double* dataX, double* dataY, int NumPts);
  void knnSearch(double* Xq, double* Yq, double* Xnn, double* Ynn, double* indnn);

 private:
  std::vector<double> PX_;
  std::vector<double> PY_;
  std::vector<double> indData_;
  int NumPts_;
  int BucketSize_;
  int level_;
  bool calcInBucket(double* Xq, double* Yq);
};


#endif
