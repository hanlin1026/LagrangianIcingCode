#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <random>
#include "Grid/PLOT3D.h"
#include "QuadTree/Bucket.h"
#include "Cloud/Cloud.h"
#include "Cloud/ParcelScalars.h"
#include "Airfoil/Airfoil.h"
#include "InputData/readInputParams.h"
#include <iterator>
#include <findAll.h>

// Airfoil icing code driver program

int main(int argc, const char *argv[]) {
  // Specify initialization files
  const char *inFileName = "InputData/Input.dat";
  const char *meshFileName = "Grid/MESH.P3D";
  const char *solnFileName = "Grid/q103.0.50E+01.bin";
  // Read in initialization scalars from input file
  FluidScalars scalars;
  ParcelScalars scalarsParcel;
  readInputParams(scalars,scalarsParcel,inFileName);
  // Initialize plot3D object, read in basic problem data
  PLOT3D p3d = PLOT3D(meshFileName, solnFileName, &scalars);
  // Initialize cloud of particles
  int particles = scalarsParcel.particles_;
  double Rmean  = scalarsParcel.Rmean_;
  double Tmean  = scalarsParcel.Tmean_;
  double Xmin   = scalarsParcel.Xmin_;
  double Xmax   = scalarsParcel.Xmax_;
  double Ymin   = scalarsParcel.Ymin_;
  double Ymax   = scalarsParcel.Ymax_;
  int maxiter   = scalarsParcel.maxiter_;
  State state   = State(particles);
  double pg, ug, vg, Xnn, Ynn;
  int indnn;
  // Select particle locations randomly and set state
  default_random_engine generator;
  uniform_real_distribution<double> distX(Xmin,Xmax);
  uniform_real_distribution<double> distY(Ymin,Ymax);
  for (int i=0; i<particles; i++) {
    state.x_(i) = distX(generator);
    state.y_(i) = distY(generator);
    p3d.pointSearch(state.x_(i),state.y_(i),Xnn,Ynn,indnn);
    state.u_(i) = p3d.getUCENT(indnn);
    state.v_(i) = p3d.getVCENT(indnn);
    state.r_(i) = Rmean;
    state.temp_(i) = Tmean;
    state.time_(i) = 0;
    state.numDrop_(i) = 1;
  }
  Cloud cloud(state,p3d,scalars.rhol_);
  // Intialize airfoil object
  Eigen::MatrixXd Xgrid = p3d.getX();
  Eigen::MatrixXd Ygrid = p3d.getY();
  Eigen::VectorXd X(Xgrid.rows());
  Eigen::VectorXd Y(Ygrid.rows());
  int iter = 0;
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
  // Advect (no splashing/fracture)
  ofstream foutX("CloudX.out");
  ofstream foutY("CloudY.out");
  ofstream foutCELLX("CloudCELLX.out");
  ofstream foutCELLY("CloudCELLY.out");
  ofstream foutBINS("BetaBins.out");
  ofstream foutMASS("Beta.out");
  State stateCloud;
  iter = 0;
  int totalImpinge = 0;
  vector<double> x;
  vector<double> y;
  vector<int> impinge;
  vector<int> totalImpingeInd;
  vector<int> indAdv;
  vector<int> indCell;
  vector<int> indSplash;
  double xCENT,yCENT;
  vector<double> XCENT;
  vector<double> YCENT;
  int indtmp = 0;
  int numSplash = 0;
  printf("maxiter = %d\n",maxiter);
  while ((totalImpinge < particles) && (iter < maxiter)) {
    cloud.calcDtandImpinge(airfoil,p3d);
    cloud.transportSLD(p3d);
    impinge = cloud.getIMPINGE();
    if (!impinge.empty()) {
      cloud.computeImpingementRegimes(airfoil);
      cloud.bounceDynamics(airfoil);
      cloud.spreadDynamics(airfoil);
      cloud.splashDynamics(airfoil);
    }
    totalImpingeInd = cloud.getIMPINGETOTAL();
    totalImpinge = totalImpingeInd.size();
    stateCloud = cloud.getState();
    particles = stateCloud.size_;
    // Save states
    if (iter % 25==0) {
      indCell = cloud.getINDCELL();
      for (int i=0; i<particles; i++) {
        x.push_back(stateCloud.x_(i));
        y.push_back(stateCloud.y_(i));
      }
      for (int i=0; i<indCell.size(); i++) {
        xCENT = p3d.getXCENT(indCell[i]);
        yCENT = p3d.getYCENT(indCell[i]);
        XCENT.push_back(xCENT);
        YCENT.push_back(yCENT);
      }
    }
    indAdv = cloud.getIndAdv();
    printf("ITER = %d\t%d\t%d\n",iter,particles,indAdv.size());
    iter++;

  }
  // Get collection efficiency and output to file
  int numBins = 35;
  gsl_histogram *h = airfoil.calcCollectionEfficiency(numBins);
  std::vector<double> bins(numBins+1);
  std::vector<double> mass(numBins);
  double upper,lower;
  for (int i=0; i<numBins; i++) {
    gsl_histogram_get_range(h,i,&lower,&upper);
    mass[i] = gsl_histogram_get(h,i);
    bins[i] = lower;
    bins[i+1] = upper;
  }
  ostream_iterator<double> out_itBINS(foutBINS,"\n");
  copy ( bins.begin(), bins.end(), out_itBINS );
  ostream_iterator<double> out_itMASS(foutMASS,"\n");
  copy ( mass.begin(), mass.end(), out_itMASS );
  // Output particle state history to file
  ostream_iterator<double> out_itX (foutX,"\n");
  copy ( x.begin(), x.end(), out_itX );
  ostream_iterator<double> out_itY (foutY,"\n");
  copy ( y.begin(), y.end(), out_itY );
  ostream_iterator<double> out_itCELLX (foutCELLX,"\n");
  copy ( XCENT.begin(), XCENT.end(), out_itCELLX );
  ostream_iterator<double> out_itCELLY (foutCELLY,"\n");
  copy ( YCENT.begin(), YCENT.end(), out_itCELLY );
  // Clear any allocated memory, close files/streams
  
}
