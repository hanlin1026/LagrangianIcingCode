#ifndef __STATE_H__
#define __STATE_H__

#include <eigen3/Eigen/Dense>

class State {
  public:
    State();
    State(int size);
    ~State();
    void appendState(State& addition);
    int size_;
    // Scalars defining particle state
    Eigen::VectorXd x_; 
    Eigen::VectorXd y_;
    Eigen::VectorXd u_; 
    Eigen::VectorXd v_;
    Eigen::VectorXd r_;
    Eigen::VectorXd temp_;
    Eigen::VectorXd time_;
    Eigen::VectorXd numDrop_;
 private:

};

#endif
