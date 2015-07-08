#ifndef __CALCIMPINGEMENTLIMITS_H__
#define __CALCIMPINGEMENTLIMITS_H__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Cloud/Cloud.h>
#include <Grid/PLOT3D.h>
#include <Airfoil/Airfoil.h>

std::vector<double> calcImpingementLimits(double Xloc,double R,double T,double rhoL,PLOT3D& p3d);
void resetCloud(Cloud& cloud,PLOT3D& p3d,double X,double Ylower,double Yupper);
void calcHitMissLower(double& Yhit,double& Ymiss,Cloud& cloud,PLOT3D& p3d,Airfoil& airfoil);
void calcHitMissUpper(double& Yhit,double& Ymiss,Cloud& cloud,PLOT3D& p3d,Airfoil& airfoil);
void findInitialHit(Cloud& cloud, PLOT3D& p3d, Airfoil& airfoil, double& Yhit);

#endif
