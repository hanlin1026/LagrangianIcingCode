#ifndef __READTHERMOFILES_H__
#define __READTHERMOFILES_H__

#include <stdio.h>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>

void readSkinFrictionCoeff(const char* fileSkinFriction,Eigen::MatrixXd& cF);
void readHeatFlux(const char* fileHeatFlux,Eigen::MatrixXd& cH);

#endif
