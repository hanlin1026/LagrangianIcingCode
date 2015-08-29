#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <iterator>
#include <algorithm>
#include "readSkinFrictionCoeff.h"

void readSkinFrictionCoeff(const char* fileSkinFriction,Eigen::MatrixXd& cF) {
  // Function to read in skin friction coefficient from file

  // Initialize file stream
  std::ifstream inFile_CF;
  inFile_CF.open(fileSkinFriction);
  std::string line; std::istringstream lin;
  // Determine size of input file
  int sizeCF = std::count(std::istreambuf_iterator<char>(inFile_CF), 
			  std::istreambuf_iterator<char>(), '\n') + 1;
  // Resize cF matrix
  cF.resize(sizeCF,2);
  double a,b;
  // Read from file (two column format = (s,cf))
  for (int i=0; i<sizeCF; i++) {
    std::getline(inFile_CF, line);
    lin.clear();
    lin.str(line);
    lin >> a; printf("%f\n",a); // First number is s-coordinate
    lin >> b; // Second number is cF
    cF(i,0) = a; cF(i,1) = b;
  }
  // Close input filestream
  inFile_CF.close();

}
