#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "readInputParams.h"

void readInputParams(FluidScalars& PROPS, ParcelScalars& PARCEL, const char *inFileName) {
  // Function to read in simulation parameters from specified input
  // file and return them in a property struct

  // Initialize input file stream
  std:string line = "";
  std::ifstream inFile;
  inFile.open(inFileName);
  std::stringstream theLine(line);
  for (int i=0; i<4; i++) {
    getline(inFile,line);
  }
  // Physical parameters (pinf,R,Tinf,rhol)
  std::getline(inFile,line,'\t');
  inFile >> PROPS.pinf_;
  std::getline(inFile,line,'\t');
  inFile >> PROPS.R_;
  std::getline(inFile,line,'\t');
  inFile >> PROPS.Tinf_;
  std::getline(inFile,line,'\t');
  inFile >> PROPS.rhol_;
  // Parcel cloud properties (particles,Rmean,Tmean,Xmin,Xmax,Ymin,Ymax,maxiter)
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.particles_;
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.Rmean_;
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.Tmean_;
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.Xmin_;
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.Xmax_;
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.Ymin_;
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.Ymax_;
  std::getline(inFile,line,'\t');
  inFile >> PARCEL.maxiter_;
  // Compute derived parameters
  PROPS.rhoinf_ = PROPS.pinf_/PROPS.R_/PROPS.Tinf_;
  PROPS.Ubar_ = sqrt(1.4*PROPS.pinf_/PROPS.rhoinf_);
  // Close file stream
  inFile.close();

}
