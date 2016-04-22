#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <istream>
#include "readInputParams.h"

void readInputParams(FluidScalars& PROPS, ParcelScalars& PARCEL, const char *inFileName) {
  // Function to read in simulation parameters from specified input
  // file and return them in a property struct
  
  // **********************************
  // INITIALIZE INPUT FILE STREAM
  // **********************************

  std::string line = "";
  std::ifstream inFile;
  inFile.open(inFileName);

  // **********************************
  // GRID/FLOW SOLUTION FILENAMES
  // **********************************

  // Scan through headers
  for (int i=0; i<3; i++) {
    std::getline(inFile,line);
  }
  // Get names of files
  std::getline(inFile,line,'/'); inFile.putback('/');
  std::getline(inFile,line,'\n');
  PROPS.gridfile_.assign(line);
  printf("GRID: %s\n",PROPS.gridfile_.c_str());
  std::getline(inFile,line,'/'); inFile.putback('/');
  std::getline(inFile,line,'\n');
  PROPS.solnfile_.assign(line);
  printf("SOLN: %s\n",PROPS.solnfile_.c_str());
  std::getline(inFile,line,'/'); inFile.putback('/');
  std::getline(inFile,line,'\n');
  PROPS.heatfile_.assign(line);
  printf("HEAT: %s\n",PROPS.heatfile_.c_str());
  std::getline(inFile,line,'/'); inFile.putback('/');
  std::getline(inFile,line,'\n');
  PROPS.betafile_.assign(line);
  printf("BETA: %s\n",PROPS.betafile_.c_str());
  std::getline(inFile,line,'/'); inFile.putback('/');
  std::getline(inFile,line,'\n');
  PROPS.outfile_.assign(line);
  printf("OUT: %s\n",PROPS.outfile_.c_str());

  // **********************************
  // DROPLET ADVECTION PARAMETERS
  // **********************************
  
  // Scan through headers
  for (int i=0; i<4; i++)
    std::getline(inFile,line);
  // Physical parameters (pinf,R,Tinf,rhol)
  inFile >> PROPS.pinf_;
  inFile >> PROPS.R_;
  inFile >> PROPS.Tinf_;
  inFile >> PROPS.rhol_;
  // Parcel cloud properties (particles,Rmean,Tmean)
  std::getline(inFile,line);
  std::getline(inFile,line);
  inFile >> PARCEL.particles_;
  inFile >> PARCEL.Rmean_;
  inFile >> PARCEL.Tmean_;
  // Domain box properties (Xmin,Xmax,Ymin,Ymax)
  std::getline(inFile,line);
  std::getline(inFile,line);
  inFile >> PARCEL.Xmin_;
  inFile >> PARCEL.Xmax_;
  inFile >> PARCEL.Ymin_;
  inFile >> PARCEL.Ymax_;
  // Simulation properties (maxiter)
  std::getline(inFile,line);
  std::getline(inFile,line);
  inFile >> PARCEL.maxiter_;
  // Driver refresh rate (refreshRate)
  inFile >> PARCEL.refreshRate_;
  // Calculate impingement limits (calcImpingementLimits)
  inFile >> PROPS.calcImpingementLimits_;
  // Splashing flags
  std::getline(inFile,line);
  std::getline(inFile,line);
  inFile >> PARCEL.SplashFlag_;
  inFile >> PARCEL.TrackSplashFlag_;
  // Compute derived parameters
  PROPS.rhoinf_ = PROPS.pinf_/PROPS.R_/PROPS.Tinf_;
  PROPS.Ubar_ = sqrt(1.4*PROPS.pinf_/PROPS.rhoinf_);

  // **********************************
  // THERMO PARAMETERS
  // **********************************
  
  // Read through headers
  for (int i=0; i<5; i++)
    std::getline(inFile,line);
  // First line (NPts,Uinf,LWC,Td)
  double tmp;
  inFile >> tmp;
  inFile >> tmp;
  inFile >> tmp;
  inFile >> PROPS.Td_;
  inFile >> PROPS.chord_;
  // Close file stream
  inFile.close();

}
