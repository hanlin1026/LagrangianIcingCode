#ifndef __STATE_H__
#define __STATE_H__

#include <eigen3/Eigen/Dense>

struct State {
  // Scalars defining particle state
  Eigen::VectorXd x; 
  Eigen::VectorXd y;
  Eigen::VectorXd u; 
  Eigen::VectorXd v;
  Eigen::VectorXd r;
  Eigen::VectorXd temp;
  Eigen::VectorXd time;
  Eigen::VectorXd numDrop;

};

#endif
