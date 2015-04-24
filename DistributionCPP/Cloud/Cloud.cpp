#include "Cloud.h"

using namespace std;
using namespace Eigen;

Cloud::Cloud(State& state, Bucket& gridQT, double rhol) {
  // Set initial state of particles
  state_ = state;
  rhol_ = rhol;
  particles_ = state.x.rows();
  // Search grid QT for initial cell indices
  VectorXd Xnn(particles_);
  VectorXd Ynn(particles_);
  VectorXd indCell(particles_);
  for (int i=0; i<particles_; i++) {
    double xq = state.x;
    double yq = state.y;
    gridQT.knnSearch(&xq,&yq,);
  }

}
