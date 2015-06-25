#include <math.h>
#include "RBFInterpolant.h"

using namespace Eigen;

RBFInterpolant::RBFInterpolant(Eigen::MatrixXd& X, Eigen::VectorXd& Y) {
  // Constructor for RBF interpolant

  N_ = X.cols();
  d_ = X.rows();
  X_ = X;
  Y_ = Y;

}

RBFInterpolant::~RBFInterpolant() {

}

double RBFInterpolant::evaluateRBF(double r) {
  // Evaluate RBF interpolant at specified radius

  double rsq,phi;
  if (r < 1) {
    phi = r*log(pow(r,r));
  }
  else {
    rsq = pow(r,2);
    phi = rsq*log(r);
  }
  return phi;

}

void RBFInterpolant::calcInterpolationWeights() {
  // Calculate interpolation weights

  // Form interpolation matrix
  MatrixXd A(N_,N_);
  MatrixXd V(d_+1,N_);
  MatrixXd Ainterp(N_+d_+1,N_+d_+1);
  MatrixXd Z = MatrixXd::Zero(d_+1,d_+1);
  double r;
  for (int j=0; j<N_; j++) {
    for (int i=0; i<N_; i++) {
      r = (X_.col(j)-X_.col(i)).norm();
      A(i,j) = this->evaluateRBF(r);
    }
  }
  
  for (int j=0; j<N_; j++) {
    V(0,j) = 1;
  }
  for (int j=0; j<N_; j++) {
    for (int i=1; i<d_+1; i++) {
      V(i,j) = X_(i-1,j);
    }
  }
  Ainterp << A,V.transpose(),
             V,Z;
  // Assemble RHS vector
  VectorXd b(N_+d_+1);
  VectorXd Z2 = VectorXd::Zero(d_+1);
  b << Y_,
       Z2;
  // Invert interpolation matrix to get weights
  VectorXd x(N_+d_+1);
  w_.resize(N_);
  v_.resize(d_+1);
  x = Ainterp.colPivHouseholderQr().solve(b);
  for (int i=0; i<N_; i++) {
    w_(i) = x(i);
  }
  for (int i=0; i<d_+1; i++) {
    v_(i) = x(i+N_);
  }

}

double RBFInterpolant::evaluateInterpolant(VectorXd& Xq) {
  // Evaluate RBF interpolant at query point  

  // Sum RBFs with weights given by w_
  double interp = 0;
  double r, phi;
  for (int i=0; i<N_; i++) {
    r = (Xq - X_.col(i)).norm();
    phi = this->evaluateRBF(r);
    interp += w_(i)*phi;
  }
  // Sum in contribution of constant and d monomials
  interp += v_(0);
  for (int i=1; i<d_+1; i++) {
    interp += v_(i)*Xq(i-1);
  }

  return interp;
}
