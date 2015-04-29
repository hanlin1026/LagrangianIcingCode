#ifndef CLOUD_H_
#define CLOUD_H_

#include <stdio.h>
#include "State.h"
#include <QuadTree/Bucket.h>
#include <Airfoil/Airfoil.h>
#include <Grid/PLOT3D.h>
#include <findAll.h>
#include <eigen3/Eigen/Dense>

class Cloud {
 public:
  Cloud(State& state, Bucket& gridQT, double rhol);
  ~Cloud();
  void addParticle(State& state, Bucket& gridQT);
  // Methods for SLD dynamics
  void findInSimulation();
  void computeNewCellLocations(PLOT3D& grid);
  void calcDtandImpinge(Airfoil& airfoil, PLOT3D& grid);
  // Set/get methods
  State getState();
  void setIndAdv(std::vector<int>& indAdv);
  std::vector<int> getIndAdv();

 private:
  State state_;
  double rhol_;
  int particles_;
  double sigma_;
  std::vector<int> impingeTotal_;
  std::vector<int> indCell_;
  std::vector<int> indAdv_;

};

#endif
