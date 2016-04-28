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
  inFile >> PROPS.NPts_;
  inFile >> PROPS.Uinf_;
  inFile >> PROPS.LWC_;
  inFile >> PROPS.Td_;
  inFile >> PROPS.chord_;
  inFile >> PROPS.mach_;
  inFile >> PROPS.DT_;
  // Close file stream
  inFile.close();

}
