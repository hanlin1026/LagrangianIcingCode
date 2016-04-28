#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <istream>
#include "autoGridGen.h"

void autoGridGen(const char* XYfile, const char *outDir) {
  // Function to read in airfoil XY coordinates and collect/write GAIR/HYPERG/FLO103 input files
  
  // **********************************
  // INITIALIZE INPUT/OUTPUT FILE STREAMS
  // **********************************

  // Declarations
  std::string line = "";
  std::ifstream inFile;
  std::ifstream header1; std::ifstream header2;
  std::ofstream outFile;
  char buf[256]; 
  // Open filestreams used for GAIR/HYPERG (fort.30, XYfile)
  inFile.open(XYfile);
  strcpy(buf,outDir); strcat(buf,"/fort.30");
  outFile.open(buf);

  // **********************************
  // INPUT NEW GRID COORDINATES
  // **********************************
  
  // Read in XY data
  std::vector<double> X;
  std::vector<double> Y;
  double x,y;
  while (inFile >> x) {
    inFile >> y;
    X.push_back(x); Y.push_back(y);
    std::getline(inFile,line,'\n');
  }
  int N = X.size();
  X[0] = 1.0;   Y[0] = 0.0;
  X[N-1] = 1.0; Y[N-1] = 0.0;

  // **********************************
  // GENERATE INPUT FILES FOR GAIR
  // **********************************
  
  // Read/write header1
  header1.open("header1");
  while (std::getline(header1,line,'\n')) {
    outFile << line; outFile << '\n';
  }
  // Write new grid coordinates
  for (int i=0; i<N; i++) {
    outFile << X[i]; outFile << '\t'; outFile << Y[i];
    outFile << '\n';
  }
  // Read/write header2
  header2.open("header2");
  while (std::getline(header2,line,'\n')) {
    outFile << line; outFile << '\n';
  }
  // Close file streams
  inFile.close(); outFile.close();
  header1.close(); header2.close();
  
}
