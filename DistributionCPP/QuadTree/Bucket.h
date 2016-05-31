#ifndef __BUCKET_H__
#define __BUCKET_H__

#include <vector>
#include <string>

class Bucket {
 public:
  Bucket();
  Bucket(const std::string workDir);
  Bucket(double* SW, double* SE, double* NW, double* NE);
  Bucket(double* SW, double* SE, double* NW, double* NE, int level);
  ~Bucket();
  double SW_[2], SE_[2], NW_[2], NE_[2];
  Bucket** buckets_;
  void setBounds(double* SW, double* SE, double* NW, double* NE);
  void addNode(int ind, double* sw, double* se, double* nw, double* ne);
  void setPoints(double* dataX, double* dataY, std::vector<int>& indData, int NumPts);
  void getPoints(std::vector<double>* PX, std::vector<double>* PY, std::vector<int>* indXY);
  int getNPts();
  int getLevel();
  void setBucketSize(int BS);
  void divideBucket();
  void calcQuadTree(double* dataX, double* dataY, int NumPts);
  void knnSearch(double* Xq, double* Yq, double* Xnn, double* Ynn, int* indnn);
  void setOutDir(const std::string workDir);

 private:
  std::string workDir_;
  std::vector<double> PX_;
  std::vector<double> PY_;
  std::vector<int> indData_;
  int NumPts_;
  int BucketSize_;
  int level_;
  bool calcInBucket(double* Xq, double* Yq);
};


#endif
