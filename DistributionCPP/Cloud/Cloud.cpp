#include "Cloud.h"

using namespace std;
using namespace Eigen;

Cloud::Cloud(State& state, Bucket& gridQT, double rhol) {
  // Set initial state of particles
  state_ = state;
  rhol_ = rhol;
  particles_ = state.size_;
  sigma_ = 75.64e-3;
  // Search grid QT for initial cell indices
  indCell_.resize(particles_);
  double xq, yq, Xnn, Ynn, indCell;
  for (int i=0; i<particles_; i++) {
    xq = state.x_(i);
    yq = state.y_(i);
    gridQT.knnSearch(&xq,&yq,&Xnn,&Ynn,&indCell);
    indCell_(i) = indCell;
  }

}

Cloud::~Cloud() {

}

void Cloud::addParticle(State& state, Bucket& gridQT) {
  // Function to add new particles to the cloud

  // Append new state elements 
  state_.appendState(state);
  // Search grid for new state cell indices
  indCell_.conservativeResize(particles_ + state.size_);
  double xq, yq, Xnn, Ynn, indCell;
  for (int i=0; i<state.size_; i++) {
    xq = state.x_(i);
    yq = state.y_(i);
    gridQT.knnSearch(&xq,&yq,&Xnn,&Ynn,&indCell);
    indCell_(particles_+i) = indCell;
  }
}

State Cloud::getState() {

  return state_;
}

void Cloud::calcDtandImpinge(Airfoil& airfoil, PLOT3D& grid) {
  // Function which does the following:
  // (1) calculate which particles are currently being advected
  // (2) set local time steps based on CFL condition
  // (3) determine which droplets are currently impinging

  double indT[particles_];
  for (int i=0; i<particles_; i++) {
    indT[i] = i;
  }
  

}
