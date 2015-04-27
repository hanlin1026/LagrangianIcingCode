#ifndef CLOUD_H_
#define CLOUD_H_

#include <stdio.h>
#include "State.h"
#include "../QuadTree/Bucket.h"
#include <eigen3/Eigen/Dense>

class Cloud {
 public:
  Cloud(State& state, Bucket& gridQT, double rhol);
  ~Cloud();

 private:
  State state_;
  double rhol_;
  double particles_;
  double sigma_;
  Eigen::VectorXd indCell_;

};

#endif
