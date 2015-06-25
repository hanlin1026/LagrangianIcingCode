#ifndef __RBFINTERPOLANT_H__
#define __RBFINTERPOLANT_H__

#include <eigen3/Eigen/Dense>

class RBFInterpolant {
  public:
    RBFInterpolant(Eigen::MatrixXd& X,Eigen::VectorXd& Y);
    ~RBFInterpolant();
    double evaluateRBF(double r);
    void calcInterpolationWeights();
    double evaluateInterpolant(Eigen::VectorXd& Xq);

  private:
    int d_,N_;
    Eigen::MatrixXd X_;
    Eigen::VectorXd Y_;
    Eigen::VectorXd w_;
    Eigen::VectorXd v_;
    
};



#endif
