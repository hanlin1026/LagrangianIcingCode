#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <iterator>
#include <algorithm>
#include "readThermoFiles.h"

void readSkinFrictionCoeff(const char* fileSkinFriction,Eigen::MatrixXd& cF) {
  // Function to read in skin friction coefficient from file

  // Initialize file stream
  FILE* cfFile = fopen(fileSkinFriction,"r");
  assert(cfFile != NULL);
  // Determine size of input file
  int c;
  int sizeCF = 1;
  while ( (c=fgetc(cfFile)) != EOF ) {
    if ( c == '\n' )
      sizeCF++;
  }
  rewind(cfFile);
  // Resize cF matrix
  cF.resize(sizeCF,2);
  double a,b;
  for (int i=0; i<sizeCF; i++) {
    fscanf(cfFile,"%lf %lf",&a,&b);
    cF(i,0) = a; cF(i,1) = b;
    printf("%lf %lf\n",a,b);
  }
  // Output to file
  std::ofstream foutS_CF("CF_Scoords.out");
  foutS_CF << cF;
  // Close file streams
  fclose(cfFile);
  foutS_CF.close();  

}

void readHeatFlux(const char* fileHeatFlux,Eigen::MatrixXd& cH) {
  // Function to read in convective heat transfer coefficient from file

}

