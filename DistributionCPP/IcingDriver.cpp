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
#include "ThermoEqns/ThermoEqns.h"
#include <iterator>
#include <findAll.h>

// *******************************************************
// AIRFOIL ICING CODE DRIVER PROGRAM
// *******************************************************

int main(int argc, const char *argv[]) {
  // Check that user has specified an input filepath
  if (argc < 2) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " <InputFilePath>" << std::endl;
    return 1;
  }
  // Specify initialization files
  const char *inFileName = argv[1];
  //const char *meshFileName = "/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/NACA0012/MESH.P3D";
  //const char *solnFileName = "/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/NACA0012/q103.0.40E+01.bin";
  // Read in initialization scalars from input file
  FluidScalars scalarsFluid;
  ParcelScalars scalarsParcel;
  readInputParams(scalarsFluid,scalarsParcel,inFileName);
  // Read in grid/flow solution files
  const char *meshFileName = scalarsFluid.gridfile_.c_str();
  const char *solnFileName = scalarsFluid.solnfile_.c_str();
  const char *filenameCHCF = scalarsFluid.heatfile_.c_str();
  const char *filenameBETA = scalarsFluid.betafile_.c_str();
  // Initialize plot3D object, read in basic problem data
  double chord = scalarsFluid.chord_;
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
    if (Xgrid(i,0) <= chord) {
      X.push_back(Xgrid(i,0));
      Y.push_back(Ygrid(i,0));
      iter++;
    }
  }
  Airfoil airfoil = Airfoil(X,Y);
  airfoil.calcStagnationPt(p3d);
  //airfoil.setStagPt(1.0238);
  // Advect (no splashing/fracture)
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

  // *******************************************************
  // DROPLET ADVECTION MODULE
  // *******************************************************
  /**
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
  // Get collection efficiency and output to file
  double dS = 0.0025;
  airfoil.calcCollectionEfficiency(fluxFreeStream,dS);
  std::vector<double> BetaBins = airfoil.getBetaBins();
  std::vector<double> Beta = airfoil.getBeta();
  // Output particle state history to file
  FILE* outfileDROP;
  FILE* outfileBETA;
  outfileDROP = fopen("DropletXY.out","w");
  outfileBETA = fopen("BETA.out","w");
  for (int i=0; i<x.size(); i++)
    fprintf(outfileDROP,"%lf\t%lf\n",x[i],y[i]);
  for (int i=0; i<Beta.size(); i++) 
    fprintf(outfileBETA,"%lf\t%lf\n",BetaBins[i],Beta[i]*.74/.83);
  fclose(outfileDROP);
  fclose(outfileBETA);
  **/
  // *******************************************************
  // THERMO EQUATIONS
  // *******************************************************
  
  // Initialize thermo eqns solver
  //const char *filenameCHCF = "/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/NACA0012/heatflux";
  //const char *filenameBETA = "/home/adegenna/LagrangianIcingCode/DistributionCPP/BETATMP.out";
  // Solve upper surface
  printf("SOLVING UPPER SURFACE...\n\n");
  ThermoEqns thermoUPPER = ThermoEqns(filenameCHCF,filenameBETA,airfoil,scalarsFluid,cloud,p3d,"UPPER");
  thermoUPPER.SolveLEWICEformulation();
  printf("...DONE\n\n");
  // Solve lower surface
  printf("SOLVING LOWER SURFACE...\n\n");
  ThermoEqns thermoLOWER = ThermoEqns(filenameCHCF,filenameBETA,airfoil,scalarsFluid,cloud,p3d,"LOWER");
  thermoLOWER.SolveLEWICEformulation();
  printf("...DONE\n\n");
  // Get old grid XY coordinates
  vector<double> XOLD = airfoil.getX();
  vector<double> YOLD = airfoil.getY();
  // Concatenate upper/lower surface ice growth rates
  vector<double> sUP        = thermoUPPER.getS(); sUP[0] = 0.0;
  vector<double> miceUP     = thermoUPPER.getMICE();
  vector<double> sLOW       = thermoLOWER.getS(); sLOW[sLOW.size()-1] = 0.0;
  vector<double> miceLOW    = thermoLOWER.getMICE();
  miceLOW.insert( miceLOW.end(), miceUP.begin(), miceUP.end() );
  sLOW.insert( sLOW.end(), sUP.begin(), sUP.end() );
  vector<double> mice = miceLOW;
  vector<double> s    = sLOW; 
  // Update grid (grow ice)
  double DT = 60.0*1.0;
  printf("GROWING ICE FOR DT = %lf SECONDS...\n\n",DT);
  airfoil.growIce(s,mice,DT,chord,"ENTIRE");
  printf("...DONE\n\n");
  // Output new grid coordinates to file
  vector<double> XNEW = airfoil.getX();
  vector<double> YNEW = airfoil.getY();
  FILE* outfileXYOLDNEW; FILE* outfileXYNEW;
  outfileXYOLDNEW = fopen("XY_OLD_NEW.out","w");
  outfileXYNEW = fopen("XY_NEW.out","w");
  for (int i=0; i<XNEW.size(); i++) {
    fprintf(outfileXYOLDNEW,"%lf\t%lf\t%lf\t%lf\n",XOLD[i]/chord,YOLD[i]/chord,XNEW[i]/chord,YNEW[i]/chord);
    fprintf(outfileXYNEW,"%lf\t%lf\n",XNEW[i]/chord,YNEW[i]/chord);
  }
  fclose(outfileXYOLDNEW);
  fclose(outfileXYNEW);
  
}
