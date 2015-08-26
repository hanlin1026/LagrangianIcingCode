#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <random>
#include "Grid/PLOT3D.h"
#include "QuadTree/Bucket.h"
#include "Cloud/Cloud.h"
#include "Cloud/ParcelScalars.h"
#include "Airfoil/Airfoil.h"
#include "InputData/readInputParams.h"
#include "Cloud/calcImpingementLimits.h"
#include <iterator>
#include <findAll.h>

// Airfoil icing code driver program

int main(int argc, const char *argv[]) {
  // Check that user has specified an input filepath
  if (argc < 2) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " <InputFilePath>" << std::endl;
    return 1;
  }
  // Specify initialization files
  const char *inFileName = argv[1];
  const char *meshFileName = "/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/MESH.P3D";
  const char *solnFileName = "/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/q103.0.25E+01.bin";
  // Read in initialization scalars from input file
  FluidScalars scalarsFluid;
  ParcelScalars scalarsParcel;
  readInputParams(scalarsFluid,scalarsParcel,inFileName);
  // Initialize plot3D object, read in basic problem data
  PLOT3D p3d = PLOT3D(meshFileName, solnFileName, &scalarsFluid);
  double dY;
  if (scalarsFluid.calcImpingementLimits_ == 1) { 
    // Over-ride input screen and determine impingement limits
    std::vector<double> Ylimits(2);
    Ylimits = calcImpingementLimits(scalarsParcel.Xmax_,scalarsParcel.Rmean_,scalarsParcel.Tmean_,scalarsFluid.rhol_,p3d);
    scalarsParcel.Ymin_ = Ylimits[0];
    scalarsParcel.Ymax_ = Ylimits[1];
    dY = Ylimits[1]-Ylimits[0];
  }
  else {
    // Use input screen provided in input file
    dY = scalarsParcel.Ymax_ - scalarsParcel.Ymin_;
  }
  // Initialize cloud of particles
  State state = State("MonoDispersed",scalarsParcel,p3d);
  Cloud cloud(state,p3d,scalarsFluid.rhol_,scalarsParcel);
  // Calculate initial total droplet mass in cloud
  double massTotal = cloud.calcTotalMass();
  double fluxFreeStream = massTotal/dY;
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
  airfoil.calcStagnationPt(p3d);
  //airfoil.setStagPt(1.0238);
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
  int maxiter = scalarsParcel.maxiter_;
  int refreshRate = scalarsParcel.refreshRate_;
  int particles = scalarsParcel.particles_;
  printf("maxiter = %d\n",maxiter);

  // TEMPORARY CODE TO CONVERT (I,0) TO S-COORDS *****
  const char *File_CF = "/home/adegenna/LagrangianIcingCode/DistributionCPP/ThermoEqns/SkinFrictionCoeff.dat";
  std::ifstream inFile_CF;
  inFile_CF.open(File_CF);
  int size_CF = 384;
  vector<int> I_cf(size_CF);
  std::string line; std::istringstream lin;
  vector<double> x_cf(size_CF);
  vector<double> y_cf(size_CF);
  vector<double> xy_cf(2);
  vector<double> s_cf(size_CF);
  vector<double> cf(size_CF);
  for (int i=0; i<size_CF; i++) {
    std::getline(inFile_CF, line);
    lin.clear();
    lin.str(line);
    lin >> I_cf[i]; printf("%d\n",I_cf[i]);
    lin >> cf[i];
    x_cf[i] = Xgrid(I_cf[i],0);
    y_cf[i] = Ygrid(I_cf[i],0);
    xy_cf[0] = x_cf[i]; xy_cf[1] = y_cf[i];
    s_cf[i] = airfoil.interpXYtoS(xy_cf);
  }
  ofstream foutS_CF("CF_Scoords.out");
  ostream_iterator<double> out_itS_CF (foutS_CF,"\n");
  copy ( s_cf.begin(), s_cf.end(), out_itS_CF );

  inFile_CF.close();
  
  
  // *************************************************



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
    if (iter % refreshRate==0) {
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
  // Output particle state history to file
  ostream_iterator<double> out_itX (foutX,"\n");
  copy ( x.begin(), x.end(), out_itX );
  ostream_iterator<double> out_itY (foutY,"\n");
  copy ( y.begin(), y.end(), out_itY );
  ostream_iterator<double> out_itCELLX (foutCELLX,"\n");
  copy ( XCENT.begin(), XCENT.end(), out_itCELLX );
  ostream_iterator<double> out_itCELLY (foutCELLY,"\n");
  copy ( YCENT.begin(), YCENT.end(), out_itCELLY );
  // Get collection efficiency and output to file
  double dS = 0.001;
  airfoil.calcCollectionEfficiency(fluxFreeStream,dS);
  std::vector<double> BetaBins = airfoil.getBetaBins();
  std::vector<double> Beta = airfoil.getBeta();
  ostream_iterator<double> out_itBINS(foutBINS,"\n");
  copy ( BetaBins.begin(), BetaBins.end(), out_itBINS );
  ostream_iterator<double> out_itMASS(foutMASS,"\n");
  copy ( Beta.begin(), Beta.end(), out_itMASS );
  // Clear any allocated memory, close files/streams
  
}
