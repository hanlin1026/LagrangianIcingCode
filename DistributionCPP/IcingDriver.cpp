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
#include "AutoGridGen/autoGridGen.h"
#include <iterator>
#include <findAll.h>

// ***********************************************************
// COUPLED AERO-THERMODYNAMICS FOR ICING SIMULATIONS (CATFISh)
// ***********************************************************

int main(int argc, const char *argv[]) {

  const std::string banner =
R"(#                  -+/:-.                     `                  `-://:`                      `:o:   #)" "\n"
R"(#                /o/--:/oyyy+-          `+ymNMMMNmhssso/:.   .+ymNhhhmMMd-                      `dm- #)" "\n"
R"(#              +s-   `.-.` `:ohy:      -NMMs//+ymd-   .-/shhydd+.     .yMN.                      `mM:#)" "\n"
R"(#            :y-        `:/+.-/smh+++++omM:  +hdo-oyo-      .-   .sdmd+ yMo                       /Mm#)" "\n"
R"(#          `y+    `--:/++soo/-`         `ss-dMMMMh  .+s+`      `oNMMMMM/:Mo                       :MM#)" "\n"
R"(#         -h.  `..`-oys/.                 .omMMMMM::/+oyo`     `/ymMMMM-yM:                       yMy#)" "\n"
R"(#        /h`  ``:hd+.                                              .odMyMm    `-+syys+/-`       .hMm`#)" "\n"
R"(#       /y `.-sds.                         -ohNMMMMMmhs/`             .+hNo:sdMMMMmdmNMMMNdhyyhmMNs` #)" "\n"
R"(#      -y `.od+                          -dMMMMmdhhdmMMMMmhhhh-          `+NMMmo-     .:+shhdhhs/`   #)" "\n"
R"(#     `d``om+                           oMMNo-        .:+oyyhh-           `/NMmmmds`                 #)" "\n"
R"(#     s:.dy`                           sMMy`  oddmNMMMMNNmddhyso+/::/osydNMMMMMMMMMy                 #)" "\n"
R"(#    `h/m:                            /MMs      +NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMs                 #)" "\n"
R"(#    /dm.                             mMh        `:---:/+shdNMMMMMMMMMMMMMMmddNMMMd`                 #)" "\n"
R"(#    hd`                             /MM.                    `.-:///////-.    hMNo                   #)" "\n"
R"(#   om`                              dMs                                      Ny`                    #)" "\n"
R"(#  -M-                              /MN`                     `       -`   /d+sN`                     #)" "\n"
R"(#  ds                              .NM:   `                 oMN.    NMd   .mMm-                      #)" "\n"
R"(# /N`                             /NM+    `+yyo/-.`        `MMy     sMM-`/ydhM+                      #)" "\n"
R"(# ds           -+:.            -+mMN/        -+ydNMMMNNNmdhhMMh++++osMMMy/.  /Nh.                    #)" "\n"
R"(#-M.            `/hNdhsoooosydMMMd+`              .mMMMh+-``mMN:---.` :mm.     +my-                  #)" "\n"
R"(#od +y+/.           -+shddmdhyo:`                 :NMMd     .NMo        :y+      -oo:                #)" "\n"
R"(#ds oy `:++++/:-....`                         `/hmMs`No      .hM/         `-`        .               #)" "\n"
R"(#N/ .N.   `` `....` `.                      -yhN+.M` /d        /d+                                   #)" "\n"
R"(#M:  -d.   `--:////:.-.                   -hM+ d. d-  :y.        -:                                  #)" "\n"
R"(#M:   .h+    `.-:////:- `-               sN/m: y: .d:  `+o+.                                         #)" "\n"
R"(#m+     :yo.   `-:/++--o/               hd` /o .h` `+` `+s/                                          #)" "\n"
R"(#yy       .oss/:--:+yh+`               sN`   y` .. `:+o+.                                            #)" "\n"
R"(#/N           .:/+/:`                  Mo    `ooooo/-                                                #)" "\n"
R"(# m+                                  .M+                                                            #)" "\n"
R"(# :N`                                  Nd                                         `:oyhhhhhhhdhs.    #)" "\n"
R"(#  sd                                  +My                                     .ohh+:`     :yh/      #)" "\n"
R"(#   hh                                  yMh-                                `/hy:``:+ssyhNNs.        #)" "\n"
R"(#    yd`                                 yyyh:                            -sh+`.+oo:`  +Ny.          #)" "\n"
R"(#     +m/                                 +o.sds:`                     -ohs-`/+:`    `hm-            #)" "\n"
R"(#      .hh.                                .+` -+yhyo/:.`       `.-/oyhs:`.--      `/md`             #)" "\n"
R"(#        :dy-                                .`    `-:+osyyyyyysso/:.           .++/sN.              #)" "\n"
R"(#          -yd+`                                                   .`      `.://:`  dy               #)" "\n"
R"(#            `/hho-                                                :d      `        hy               #)" "\n"
R"(#               `:shho:`                                          :Nm       `.-:/+/-/N.              #)" "\n"
R"(#                   `-+yhhs+/-`                               .:smNy.             `:/dm`             #)" "\n"
R"(#                         .-/oshhhyyssoo+++//////+++oosssyyddmmhs/.         ``.--`    oN:            #)" "\n"
R"(#                                   `.--:://////////::--------::/+ossssss+:.    `:/+/` :No           #)" "\n"
R"(#                                                                        .:+syyo:`  .+o:.my          #)" "\n"
R"(#                                                                              -+yhs:  -o/Ns         #)" "\n"
R"(#                                                                                  -smy- -yM-        #)" "\n"
R"(#                                                                                     :yd/`my        #)" "\n"
R"(#                                                                                       .sdhh        #)" "\n"
R"(#                                                                                         .h+        #)" "\n"
R"(#                           ______ ___   ______ ______ ____ _____  __                                #)" "\n"
R"(#                          / ____//   | /_  __// ____//  _// ___/ / /_                               #)" "\n"
R"(#                         / /    / /| |  / /  / /_    / /  \__ \ / __ \                              #)" "\n"
R"(#                        / /___ / ___ | / /  / __/  _/ /  ___/ // / / /                              #)" "\n"
R"(#                        \____//_/  |_|/_/  /_/    /___/ /____//_/ /_/                               #)" "\n";

  printf("\n\n%s\n\n\n",banner.c_str());

  // Check that user has specified an input filepath
  if (argc < 4) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << "<IcingInputFile> " << "<InputDirectory> " << "<OutputDirectory>" <<  std::endl;
    return 1;
  }

  // Specify initialization files
  const std::string s_inFileName(argv[1]);
  const std::string s_inDir(argv[2]);
  const std::string s_outDir(argv[3]);

  // Read in initialization scalars from input file
  FluidScalars scalarsFluid;
  ParcelScalars scalarsParcel;
  readInputParams(scalarsFluid,scalarsParcel,s_inFileName.c_str());

  // Read in grid/flow solution files
  const std::string s_meshFileName = s_inDir + "/MESH.P3D";
  const std::string s_solnFileName = s_inDir + "/q103.0.40E+01.bin";
  const std::string s_filenameCHCF = s_inDir + "/heatflux";
  const std::string s_filenameBETA = s_inDir + "/BETA.out";

  // Initialize plot3D object, read in basic problem data
  double chord = scalarsFluid.chord_;
  PLOT3D p3d = PLOT3D(s_meshFileName.c_str(), s_solnFileName.c_str(), &scalarsFluid);  
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
  Airfoil airfoil = Airfoil(s_inDir,X,Y);
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
  int numIndAdv = 0;
  int maxiter = scalarsParcel.maxiter_;
  int refreshRate = scalarsParcel.refreshRate_;
  int particles = scalarsParcel.particles_;
  printf("maxiter = %d\n",maxiter);
  
  // *******************************************************
  // DROPLET ADVECTION MODULE
  // *******************************************************
  
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
    numIndAdv = indAdv.size();
    printf("ITER = %d\t%d\t%d\n",iter,particles,numIndAdv);
    iter++;

  }
  // Get collection efficiency and output to file
  const std::string s_dropName = s_inDir + "/DropletXY.out";
  double dS = 0.0025;
  airfoil.calcCollectionEfficiency(fluxFreeStream,dS);
  std::vector<double> BetaBins = airfoil.getBetaBins();
  std::vector<double> Beta = airfoil.getBeta();
  // Output particle state history to file
  FILE* outfileDROP;
  FILE* outfileBETA;
  outfileDROP = fopen(s_dropName.c_str(),"w");
  outfileBETA = fopen(s_filenameBETA.c_str(),"w");
  for (int i=0; i<x.size(); i++)
    fprintf(outfileDROP,"%lf\t%lf\n",x[i],y[i]);
  for (int i=0; i<Beta.size(); i++) 
    fprintf(outfileBETA,"%lf\t%lf\n",BetaBins[i],Beta[i]*.74/.83);
  fclose(outfileDROP);
  fclose(outfileBETA);
  
  // *******************************************************
  // THERMO EQUATIONS
  // *******************************************************
  
  // Solve upper surface
  printf("SOLVING UPPER SURFACE...\n\n");
  ThermoEqns thermoUPPER = ThermoEqns(s_inDir,s_filenameCHCF.c_str(),s_filenameBETA.c_str(),airfoil,scalarsFluid,cloud,p3d,"UPPER","MULTISHOT");
  thermoUPPER.SolveLEWICEformulation();
  //thermoUPPER.SolveIcingEqns();
  //thermoUPPER.explicitSolverSimultaneous(5.0e-1,1.0e-4);
  printf("...DONE\n\n");

  // Solve lower surface
  printf("SOLVING LOWER SURFACE...\n\n");
  ThermoEqns thermoLOWER = ThermoEqns(s_inDir,s_filenameCHCF.c_str(),s_filenameBETA.c_str(),airfoil,scalarsFluid,cloud,p3d,"LOWER","MULTISHOT");
  thermoLOWER.SolveLEWICEformulation();
  //thermoLOWER.SolveIcingEqns();
  //thermoLOWER.explicitSolverSimultaneous(5.0e-1,1.0e-4);
  printf("...DONE\n\n");

  // Get old grid XY coordinates
  vector<double> XOLD = airfoil.getX();
  vector<double> YOLD = airfoil.getY();
  // Concatenate upper/lower surface ice growth rates
  vector<double> sUP        = thermoUPPER.getS(); //sUP[0] = 0.0;
  vector<double> miceUP     = thermoUPPER.getMICE();
  vector<double> sLOW       = thermoLOWER.getS(); //sLOW[sLOW.size()-1] = 0.0;
  vector<double> miceLOW    = thermoLOWER.getMICE();
  miceLOW.insert( miceLOW.end(), miceUP.begin(), miceUP.end() );
  sLOW.insert( sLOW.end(), sUP.begin(), sUP.end() );
  vector<double> mice = miceLOW;
  vector<double> s    = sLOW; 

  // Update grid (grow ice)
  double DT = scalarsFluid.DT_;
  printf("GROWING ICE FOR DT = %lf SECONDS...\n\n",DT);
  airfoil.growIce(s,mice,DT,chord,"ENTIRE");
  printf("...DONE\n\n");

  // Output new grid coordinates to file
  vector<double> XNEW = airfoil.getX();
  vector<double> YNEW = airfoil.getY();
  FILE* outfileXYOLDNEW; FILE* outfileXYNEW;
  const std::string s_xyOldNew = s_outDir + "/XY_OLD_NEW.out";
  const std::string s_xyNew    = s_outDir + "/XY_NEW.out";
  outfileXYOLDNEW = fopen(s_xyOldNew.c_str(),"w");
  outfileXYNEW = fopen(s_xyNew.c_str(),"w");
  for (int i=0; i<XNEW.size(); i++) {
    fprintf(outfileXYOLDNEW,"%lf\t%lf\t%lf\t%lf\n",XOLD[i]/chord,YOLD[i]/chord,XNEW[i]/chord,YNEW[i]/chord);
    fprintf(outfileXYNEW,"%lf\t%lf\n",XNEW[i]/chord,YNEW[i]/chord);
  }
  fclose(outfileXYOLDNEW);
  fclose(outfileXYNEW);
  
  // *******************************************************
  // INPUT FILE CREATION FOR NEW GRID GENERATION
  // *******************************************************

  // Generate input files for GAIR/HYPERG mesh generation
  autoGridGen(s_xyNew.c_str(),s_outDir.c_str());
  
}
