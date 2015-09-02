#include "ThermoEqns.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <iterator>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
using namespace Eigen;

ThermoEqns::ThermoEqns(const char* filenameCF,Airfoil& airfoil) {
  // Constructor to read in input files and initialize thermo eqns

  // Import data from file
  MatrixXd data = this->readCHCF(filenameCF);
  VectorXd s = data.col(0);
  VectorXd ch = data.col(1);
  VectorXd cf = data.col(2);
  double stagPt = airfoil.getStagPt();
  for (int i=0; i<s.size(); i++) {
    s(i) -= stagPt;
  }
  // Find stagnation point
  double zero = 1.0; double zero_new = 1.0;
  int zero_ind = 0;
  for (int i=0; i<s.size(); i++) {
    zero_new = std::abs(s(i));
    if (zero_new < zero) {
      zero = zero_new;
      zero_ind = i;
    }
  }
  // Split into upper/lower surface grids
  int NPts_orig = s.size() - zero_ind;
  vector<double> s_upper_orig(NPts_orig);
  vector<double> cf_upper_orig(NPts_orig);
  for (int i=0; i<NPts_orig; i++) {
    s_upper_orig[i] = s(zero_ind+i);
    cf_upper_orig[i] = cf(zero_ind+i);
    printf("%lf %lf\n",s_upper_orig[i],cf_upper_orig[i]);
  }
  // Interpolate upper/lower surface grids
  NPts_ = 1000;
  s_upper_.resize(NPts_);
  cF_upper_.resize(NPts_);
  double s_max = s.maxCoeff();
  double ds = (s_max-zero)/NPts_;
  for (int i=0; i<NPts_; i++) {
    s_upper_[i] = zero + i*ds;
  }
  // Interpolate values on grids
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *splineCF = gsl_spline_alloc(gsl_interp_linear, NPts_orig);
  gsl_spline_init(splineCF, &s_upper_orig[0], &cf_upper_orig[0], NPts_orig);
  for (int i=0; i<NPts_; i++) {
    cF_upper_[i] = gsl_spline_eval(splineCF, s_upper_[i], acc);
  }

}











ThermoEqns::~ThermoEqns() {

}

MatrixXd ThermoEqns::readCHCF(const char* filenameCHCF) {
  // Function to read in skin friction coefficient from file

  // Initialize file stream
  FILE* filept = fopen(filenameCHCF,"r");
  assert(filept != NULL);
  // Determine size of input file
  int c;
  int sizeCHCF = 0;
  while ( (c=fgetc(filept)) != EOF ) {
    if ( c == '\n' )
      sizeCHCF++;
  }
  rewind(filept);
  // Resize cF matrix
  MatrixXd s_cH_cF(sizeCHCF,3);
  double a,b,d;
  for (int i=0; i<sizeCHCF; i++) {
    fscanf(filept,"%lf %lf %lf",&a,&b,&d);
    s_cH_cF(i,0) = a; s_cH_cF(i,1) = b; s_cH_cF(i,2) = d;
  }
  // Close file streams
  fclose(filept);

  return s_cH_cF;

}

MatrixXd ThermoEqns::readBetaXY(const char* filenameBeta) {
  // Function to read in skin friction coefficient from file

  // Initialize file stream
  FILE* filept = fopen(filenameBeta,"r");
  assert(filept != NULL);
  // Determine size of input file
  int c;
  int sizeBeta = 0;
  while ( (c=fgetc(filept)) != EOF ) {
    if ( c == '\n' )
      sizeBeta++;
  }
  rewind(filept);
  // Resize cF matrix
  MatrixXd s_Beta(sizeBeta,3);
  double a,b,d;
  for (int i=0; i<sizeBeta; i++) {
    fscanf(filept,"%lf %lf %lf",&a,&b,&d);
    s_Beta(i,0) = a; s_Beta(i,1) = b;
  }
  // Close file streams
  fclose(filept);

  return s_Beta;

}
