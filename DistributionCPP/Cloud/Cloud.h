#ifndef CLOUD_H_
#define CLOUD_H_

#include <stdio.h>
#include "State.h"
#include <QuadTree/Bucket.h>
#include <Airfoil/Airfoil.h>
#include <Grid/PLOT3D.h>
#include <eigen3/Eigen/Dense>

class Cloud {
 public:
  Cloud(State& state, Bucket& gridQT, double rhol);
  ~Cloud();
  void addParticle(State& state, Bucket& gridQT);
  State getState();
  // Methods for SLD dynamics
  void calcDtandImpinge(Airfoil& airfoil, PLOT3D& grid);

 private:
  State state_;
  double rhol_;
  int particles_;
  double sigma_;
  Eigen::VectorXd indCell_;

};

#endif
