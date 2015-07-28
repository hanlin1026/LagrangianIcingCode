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
  double Ylower = -0.7;
  double Yupper = -0.1;
  // Initialize test cloud of particles at X1
  int numParticles = 100;
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
  std::vector<double> X;
  std::vector<double> Y;
  int iter = 0;
  for (int i=0; i<Xgrid.rows(); i++) {
    if (Xgrid(i,0) <= 1) {
      X.push_back(Xgrid(i,0));
      Y.push_back(Ygrid(i,0));
      iter++;
    }
  }
  Airfoil airfoil = Airfoil(X,Y);
  // First advection to find lowest trajectory that hits
  double Yhit,Ymiss;
  printf("Determining initial impinging trajectory...\n");
  findInitialHit(cloud,p3d,airfoil,Yhit);
  Ymiss = Ylower;
  // Reset cloud state so particles are between missing and hitting lower trajectories
  resetCloud(cloud,p3d,Xloc,Ylower,Yhit,numParticles);
  double TOL = 1e-5;
  iter = 0;
  int indHit;
  // Iterative procedure to determine lower limit
  printf("Determining lower impingement limit...\n");
  while (abs(Ymiss-Yhit)>TOL) {
    resetCloud(cloud,p3d,Xloc,Ymiss,Yhit,numParticles);
    calcHitMissLower(Yhit,Ymiss,cloud,p3d,airfoil);
  }
  printf("Ylower = %f\n",Yhit);
  vector<double> limits(2);
  limits[0] = Yhit;
  // Iterative procedure to determine upper limit
  Ymiss = Yupper;
  printf("Determining upper impingement limit...\n");
  while (abs(Ymiss-Yhit)>TOL) {
    resetCloud(cloud,p3d,Xloc,Yhit,Ymiss,numParticles);
    calcHitMissUpper(Yhit,Ymiss,cloud,p3d,airfoil);
  }
  printf("Yupper = %f\n",Yhit);
  limits[1] = Yhit;
  // Padding on limits
  dY = limits[1]-limits[0];
  limits[0] = limits[0] - 0.05*dY;
  limits[1] = limits[1] + 0.05*dY;
  
  return limits;

}

void resetCloud(Cloud& cloud,PLOT3D& p3d,double X,double Ylower,double Yupper,int numParticles) {
  // Function to re-initialize cloud to screen of particles between Ylower/Yupper

  // Get old properties that will not change
  State oldState = cloud.getState();
  double R = oldState.r_(0);
  double T = oldState.temp_(0);
  // Create new screen of particles
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
  // Delete old cloud properties
  cloud.clearData();
  // Update cloud
  cloud.setState(state,p3d);

}

void calcHitMissLower(double& Yhit,double& Ymiss,Cloud& cloud,PLOT3D& p3d,Airfoil& airfoil) {
  // Function to calculate hit and miss y-locations for a cloud

  // Save initial state
  State state = cloud.getState();
  // Advect screen of particles
  int maxiter = 1500;
  vector<int> impinge;
  vector<int> impingeTotal;
  for (int i=0; i<maxiter; i++) {
    cloud.calcDtandImpinge(airfoil,p3d);
    cloud.transportSLD(p3d);
    impinge = cloud.getIMPINGE();
    if (!impinge.empty()) {
      cloud.computeImpingementRegimes(airfoil);
    }
  }
  // Find hit and miss
  impingeTotal = cloud.getIMPINGETOTAL();
  int indHit;
  if (!impingeTotal.empty()) {
    indHit = *min_element(impingeTotal.begin(),impingeTotal.end());
  }
  else {
    indHit = state.size_-1;
  }
  Yhit = state.y_(indHit);
  if (indHit != 0) {
    Ymiss = state.y_(indHit-1);
  }
  else {
    Ymiss = Yhit;
  }

  printf("IndHit = %d, Yhit = %f, Ymiss = %f\n",indHit,Yhit,Ymiss);

}

void calcHitMissUpper(double& Yhit,double& Ymiss,Cloud& cloud,PLOT3D& p3d,Airfoil& airfoil) {
  // Function to calculate hit and miss y-locations for a cloud

  // Save initial state
  State state = cloud.getState();
  // Advect screen of particles
  int maxiter = 1500;
  vector<int> impinge;
  vector<int> impingeTotal;
  for (int i=0; i<maxiter; i++) {
    cloud.calcDtandImpinge(airfoil,p3d);
    cloud.transportSLD(p3d);
    impinge = cloud.getIMPINGE();
    if (!impinge.empty()) {
      cloud.computeImpingementRegimes(airfoil);
    }
  }
  // Find hit and miss
  impingeTotal = cloud.getIMPINGETOTAL();
  int indHit;
  if (!impingeTotal.empty()) {
    indHit = *max_element(impingeTotal.begin(),impingeTotal.end());
  }
  else {
    indHit = 0;
  }
  Yhit = state.y_(indHit);
  if (indHit != state.size_-1) {
    Ymiss = state.y_(indHit+1);
  }
  else {
    Ymiss = Yhit;
  }

  printf("IndHit = %d, Yhit = %f, Ymiss = %f\n",indHit,Yhit,Ymiss);

}


void findInitialHit(Cloud& cloud, PLOT3D& p3d, Airfoil& airfoil, double& Yhit) {
  // Function to calculate a location in y that does hit the airfoil

  State initialState = cloud.getState();
  double X = initialState.x_(0);
  double Ylower,Yupper,indHit,dY;
  // Advect screen of particles
  int maxiter = 1500;
  vector<int> impinge;
  vector<int> impingeTotal;
  for (int i=0; i<maxiter; i++) {
    cloud.calcDtandImpinge(airfoil,p3d);
    cloud.transportSLD(p3d);
    impinge = cloud.getIMPINGE();
    if (!impinge.empty()) {
      cloud.computeImpingementRegimes(airfoil);
    }
  }
  impingeTotal = cloud.getIMPINGETOTAL();
  State state = cloud.getState();
  int iterations = 0;
  Ylower = initialState.y_(0);
  Yupper = initialState.y_(initialState.size_-1);
  int numParticles = initialState.size_;
  while (impingeTotal.empty()) {
    // Increase resolution, keep screen limits the same
    numParticles = 5*numParticles;
    resetCloud(cloud,p3d,X,Ylower,Yupper,numParticles);
    printf("No hit; increasing screen resolution to %d\n",numParticles);
    // Save initial state
    initialState = cloud.getState();
    // Re-advect particles
    for (int i=0; i<maxiter; i++) {
      cloud.calcDtandImpinge(airfoil,p3d);
      cloud.transportSLD(p3d);
      impinge = cloud.getIMPINGE();
      if (!impinge.empty()) {
	cloud.computeImpingementRegimes(airfoil);
      }
    }
    // Check for impingements
    impingeTotal = cloud.getIMPINGETOTAL();
    iterations++;
  }
  // Record hit
  indHit = *min_element(impingeTotal.begin(),impingeTotal.end());
  int indHitMax = *max_element(impingeTotal.begin(),impingeTotal.end());
  Yhit = initialState.y_(indHit);
  double YhitMax = initialState.y_(indHitMax);
  printf("Yhit = %f\n",Yhit);

}
