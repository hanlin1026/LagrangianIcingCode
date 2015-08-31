#include "ThermoEqns.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <iterator>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;

ThermoEqns::ThermoEqns(const char* filenameCF) {
  // Constructor to read in input files and initialize thermo eqns

  // Import data from file
  MatrixXd CF = this->readSkinFrictionCoeff(filenameCF);
  VectorXd s = CF.col(0);
  VectorXd cf = CF.col(1);
  // Establish grid
  
  

}











ThermoEqns::~ThermoEqns() {

}

MatrixXd ThermoEqns::readSkinFrictionCoeff(const char* fileSkinFriction) {
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
  MatrixXd cF(sizeCF,2);
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

  return cF;

}

MatrixXd ThermoEqns::readHeatFlux(const char* fileHeatFlux) {


}
