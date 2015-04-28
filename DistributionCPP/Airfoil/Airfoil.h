#ifndef __AIRFOIL_H__
#define __AIRFOIL_H__

#include <eigen3/Eigen/Dense>
#include <QuadTree/Bucket.h>

class Airfoil {
  public:
    Airfoil(Eigen::VectorXd& X, Eigen::VectorXd& Y);
    ~Airfoil();

  private:
    Eigen::VectorXd panelX_;
    Eigen::VectorXd panelY_;
    Eigen::VectorXd panelS_;
    Eigen::MatrixXd tangent_;
    Eigen::MatrixXd normal_;
    Eigen::MatrixXd FilmSplash_;
    Eigen::MatrixXd FilmSpread_;
    Eigen::MatrixXd Film_;
    double impingeLimitUP_;
    double impingeLimitDOWN_;
    double stagPt_;
    Bucket panelSearcher_;
    
};



#endif
