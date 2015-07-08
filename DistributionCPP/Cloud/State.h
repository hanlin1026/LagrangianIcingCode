#ifndef __STATE_H__
#define __STATE_H__

#include <eigen3/Eigen/Dense>
#include <Grid/PLOT3D.h>
#include <Cloud/ParcelScalars.h>

class State {
  public:
    State(const char* distributionType, ParcelScalars& scalarsParcel, PLOT3D& p3d);
    State(int size);
    State();
    ~State();
    // Append method
    void appendState(State& addition);
    // Scalars defining particle state
    int size_;
    Eigen::VectorXd x_; 
    Eigen::VectorXd y_;
    Eigen::VectorXd u_; 
    Eigen::VectorXd v_;
    Eigen::VectorXd r_;
    Eigen::VectorXd temp_;
    Eigen::VectorXd time_;
    Eigen::VectorXd numDrop_;
 private:
    double Xmin_,Xmax_,Ymin_,Ymax_;

};

#endif
