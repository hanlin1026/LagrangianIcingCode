#include "Cloud.h"

using namespace std;
using namespace Eigen;

Cloud::Cloud(State& state, Bucket& gridQT, double rhol) {
  // Set initial state of particles
  state_ = state;
  rhol_ = rhol;
  particles_ = state.x.rows();
  sigma_ = 75.64e-3;
  // Search grid QT for initial cell indices
  indCell_.resize(particles_);
  double xq, yq, Xnn, Ynn, indCell;
  for (int i=0; i<particles_; i++) {
    xq = state.x(i);
    yq = state.y(i);
    gridQT.knnSearch(&xq,&yq,&Xnn,&Ynn,&indCell);
    indCell_(i) = indCell;
  }

}

Cloud::~Cloud() {

}
