#include "calcImpingementLimits.h"
#include "ParcelScalars.h"
#include <eigen3/Eigen/Dense>
#include <Cloud/Cloud.h>
#include <Grid/PLOT3D.h>

using namespace std;
int iter = 0;

std::vector<double> calcImpingementLimits(double Xloc,double R,double T,double rhoL,PLOT3D& p3d) {
  // Function to calculate impingement limits for a particular droplet location/size/temperature/rhoL

  // Initial guesses at impingement limits
  double Ylower = -1.0;
  double Yupper = 0.0;
  // Initialize test cloud of particles at X1
  int numParticles = 1000;
  int indnn;
  double Xnn,Ynn;
  State state(numParticles);
  double dY = (Yupper-Ylower)/(numParticles-1);
  for (int i=0; i<numParticles; i++) {
    state.x_(i) = Xloc;
    state.y_(i) = Ylower + i*dY;
    p3d.pointSearch(state.x_(i),state.y_(i),Xnn,Ynn,indnn);
    state.u_(i) = p3d.getUCENT(indnn);
    state.v_(i) = p3d.getVCENT(indnn);
    state.r_(i) = R;
    state.temp_(i) = T;
    state.time_(i) = 0;
    state.numDrop_(i) = 1; 
  }
  Cloud cloud(state,p3d,rhoL);
  // Intialize airfoil object
  Eigen::MatrixXd Xgrid = p3d.getX();
  Eigen::MatrixXd Ygrid = p3d.getY();
  Eigen::VectorXd X(Xgrid.rows());
  Eigen::VectorXd Y(Ygrid.rows());
  iter = 0;
  for (int i=0; i<Xgrid.rows(); i++) {
    if (Xgrid(i,0) <= 1) {
      X(iter) = Xgrid(i,0);
      Y(iter) = Ygrid(i,0);
      iter++;
    }
  }
  X = X.block(0,0,iter,1);
  Y = Y.block(0,0,iter,1);
  Airfoil airfoil = Airfoil(X,Y);
  // First advection to find lowest trajectory that hits
  double Yhit,Ymiss;
  findInitialHit(cloud,p3d,airfoil,Yhit);
  Ymiss = Ylower;
  // Reset cloud state so particles are between missing and hitting lower trajectories
  resetCloud(cloud,p3d,Xloc,Ylower,Yhit);
  double TOL = 1e-5;
  iter = 0;
  int indHit;
  // Iterative procedure to determine lower limit
  while (abs(Ymiss-Yhit)>TOL) {
    calcHitMissLower(Yhit,Ymiss,cloud,p3d,airfoil);
  }
  vector<double> limits(2);
  limits[0] = Yhit;
  // Reset cloud state so particles are between missing and hitting upper trajectories
  resetCloud(cloud,p3d,Xloc,Yhit,Yupper);
  // Iterative procedure to determine upper limit
  Ymiss = Yupper;
  while (abs(Ymiss-Yhit)>TOL) {
    calcHitMissUpper(Yhit,Ymiss,cloud,p3d,airfoil);
  }
  limits[1] = Yhit;
  
  return limits;

}

void resetCloud(Cloud& cloud,PLOT3D& p3d,double X,double Ylower,double Yupper) {
  // Function to re-initialize cloud to screen of particles between Ylower/Yupper

  // Get old properties that will not change
  State oldState = cloud.getState();
  double R = oldState.r_(0);
  double T = oldState.temp_(0);
  // Create new screen of particles
  int numParticles = 1000;
  int indnn;
  double Xnn,Ynn;
  State state(numParticles);
  double dY = (Yupper-Ylower)/(numParticles-1);
  for (int i=0; i<numParticles; i++) {
    state.x_(i) = X;
    state.y_(i) = Ylower + i*dY;
    p3d.pointSearch(state.x_(i),state.y_(i),Xnn,Ynn,indnn);
    state.u_(i) = p3d.getUCENT(indnn);
    state.v_(i) = p3d.getVCENT(indnn);
    state.r_(i) = R;
    state.temp_(i) = T;
    state.time_(i) = 0;
    state.numDrop_(i) = 1; 
  }
  // Update cloud
  cloud.setState(state,p3d);

}

void calcHitMissLower(double& Yhit,double& Ymiss,Cloud& cloud,PLOT3D& p3d,Airfoil& airfoil) {
  // Function to calculate hit and miss y-locations for a cloud

  // Advect screen of particles
  int maxiter = 2000;
  vector<int> impinge;
  for (int i=0; i<maxiter; i++) {
    cloud.calcDtandImpinge(airfoil,p3d);
    cloud.transportSLD(p3d);
  }
  // Find hit and miss
  impinge = cloud.getIMPINGETOTAL();
  State state = cloud.getState();
  int indHit;
  if (!impinge.empty()) {
    indHit = *min_element(impinge.begin(),impinge.end());
  }
  else {
    indHit = state.size_;
  }
  Yhit = state.y_(indHit);
  Ymiss = state.y_(indHit-1);

}

void calcHitMissUpper(double& Yhit,double& Ymiss,Cloud& cloud,PLOT3D& p3d,Airfoil& airfoil) {
  // Function to calculate hit and miss y-locations for a cloud

  // Advect screen of particles
  int maxiter = 2000;
  vector<int> impinge;
  for (int i=0; i<maxiter; i++) {
    cloud.calcDtandImpinge(airfoil,p3d);
    cloud.transportSLD(p3d);
  }
  // Find hit and miss
  impinge = cloud.getIMPINGETOTAL();
  State state = cloud.getState();
  int indHit;
  if (!impinge.empty()) {
    indHit = *min_element(impinge.begin(),impinge.end());
  }
  else {
    indHit = 0;
  }
  Yhit = state.y_(indHit);
  Ymiss = state.y_(indHit+1);

}


void findInitialHit(Cloud& cloud, PLOT3D& p3d, Airfoil& airfoil, double& Yhit) {
  // Function to calculate a location in y that does hit the airfoil

  State initialState = cloud.getState();
  double X = initialState.x_(0);
  // Advect screen of particles
  int maxiter = 2000;
  vector<int> impinge;
  for (int i=0; i<maxiter; i++) {
    cloud.calcDtandImpinge(airfoil,p3d);
    cloud.transportSLD(p3d);
  }
  impinge = cloud.getIMPINGETOTAL();
  State state = cloud.getState();
  double Ylower = state.y_(0);
  double Yupper = state.y_(state.size_);
  while (impinge.empty()) {
    // Reset the screen limits
    Yupper = Yupper - (Yupper-Ylower)/4;
    Ylower = Ylower + (Yupper-Ylower)/4;
    resetCloud(cloud,p3d,X,Ylower,Yupper);
    // Re-advect particles
    for (int i=0; i<maxiter; i++) {
      cloud.calcDtandImpinge(airfoil,p3d);
      cloud.transportSLD(p3d);
    }
    impinge = cloud.getIMPINGETOTAL();
  }
  // Record hit
  int indHit = *min_element(impinge.begin(),impinge.end());
  Yhit = state.y_(indHit);

}
