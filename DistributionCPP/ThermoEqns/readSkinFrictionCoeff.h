#ifndef __READSKINFRICTIONCOEFF_H__
#define __READSKINFRICTIONCOEFF_H__

#include <stdio.h>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>

void readSkinFrictionCoeff(const char* fileSkinFriction,Eigen::MatrixXd& cF);

#endif
