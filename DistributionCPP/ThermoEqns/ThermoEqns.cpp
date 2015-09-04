#include "ThermoEqns.h"
#include <Eigen/Dense>
#include <gsl_errno.h>
#include <gsl_spline.h>
#include "compcol_double.h"                 // Compressed column matrix header
#include "iohb_double.h"                    // Harwell-Boeing matrix I/O header
#include "mvblasd.h"                        // MV_Vector level 1 BLAS
#include "diagpre_double.h"                 // Diagonal preconditioner
#include MATRIX_H                           // dense matrix header
#include <GMRES/include/gmres.h>            // IML++ GMRES template

using namespace std;
using namespace Eigen;

ThermoEqns::ThermoEqns(const char* filenameCF,Airfoil& airfoil) {
  // Constructor to read in input files and initialize thermo eqns

  // Import (s,beta)
  MatrixXd datatmp = this->readBetaXY(filenameBETA);
  VectorXd sBeta = datatmp.col(0);
  VectorXd Beta = datatmp.col(1);
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
  gsl_spline *splineBETA = gsl_spline_alloc(gsl_interp_linear, NPts_orig);
  gsl_spline_init(splineCF, &s_upper_orig[0], &cf_upper_orig[0], NPts_orig);
  gsl_spline_init(splineBETA, &sBeta_upper[0], &Beta_upper[0], NPts_orig_beta);
  for (int i=0; i<NPts_; i++) {
    cF_upper_[i] = gsl_spline_eval(splineCF, s_upper_[i], acc);
  }
  // Set rhoL_,muL_;
  rhoL_ = 1000;
  muL_ = 1.787e-3;

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

// Define action of Jacobian on vector
inline vector<double> JX(vector<double> (*f)(vector<double>& Xq), vector<double>& X, vector<double>& u0) {
  vector<double> jx(u0.size());
  double eps = 1.e-6;
  vector<double> X2(u0.size());
  for (int i=0; i<u0.size(); i++) {
    X2[i] = u0[i] + eps*X[i];
  }
  vector<double> f2 = f(X2);
  vector<double> f1 = f(u0);
  for (int i=0; i<u0.size(); i++) {
    jx[i] = (1./eps)*(f2[i]-f1[i]);
  }
  return jx;
}

vector<double> ThermoEqns::massBalanceUpper(vector<double>& x) {
  // Function to compute mass balance

  vector<double> err(x.size());
  vector<double> F(x.size());
  vector<double> f(x.size());
  vector<double> xFACE(x.size()-1);
  vector<double> cfFACE(x.size()-1);
  vector<double> DF(x.size()-1);
  vector<double> D_flux(x.size()-2);
  vector<double> I_sources(x.size()-2);
  // Calculate body centered fluxes
  for (int i=0; i<x.size(); i++) {
    F[i] = (0.5/muL_)*pow(x[i],2)*cF_upper_[i];
  }
  // Calculate fluxes at cell faces (Roe scheme upwinding)
  for (int i=0; i<xFACE.size(); i++) {
    xFACE[i] = 0.5*(x[i]+x[i+1]);
    cfFACE[i] = 0.5*(cF_upper_[i]+cF_upper_[i+1]);
    DF[i] = (1/muL_)*(xFACE[i]*cfFACE[i]);
    f[i] = 0.5*(F[i]+F[i+1]) - 0.5*std::abs(DF[i])*(x[i+1]-x[i]);
  }
  // Calculate error for internal cells
  double ds;
  for (int i=1; i<x.size()-1; i++) {
    ds = s_upper_[i+1]-s_upper_[i];
    D_flux[i-1] = f[i]-f[i-1];
    I_sources[i-1] = (1./rhoL_)*ds*(mimp_[i]-mice_upper_[i]);
    err[i] = D_flux[i] - I_sources[i];
  }
  // Boundary conditions
  err[0] = x[0]-0;
  err[x.size()-1] = err[x.size()-2];
}




void ThermoEqns::NewtonKrylovIteration(vector<double> (*f)(vector<double>& X), vector<double>& u0) {
  // Function to take a balance of form f(X) = 0 and do Newton-Krylov iteration

  double tol = 1.e-3;                       // Convergence tolerance
  int result, maxit = 100, restart = 10;    // GMRES Maximum, restart iterations

  // Initialize Jacobian and RHS, solution vectors
  int stateSize = u0.size();
  VECTOR_double b(stateSize, 0.0); vector<double> bvec(stateSize);
  vector<double> x = u0;
  vector<double> x0;
  vector<double> dx0;
  vector<double> jx;
  // Storage for upper Hessenberg H
  MATRIX_double H(restart+1, restart, 0.0);

  // Begin iteration
  int nitermax = 30; double eps = 1.e-6;
  vector<double> globalerr; double globaltol = 1.0e-5;
  vector<double> r(stateSize);
  double normR,normGlob;
  for (int i=0; i<nitermax; i++) {
    x0 = x;
    // Compute RHS
    bvec = f(x0);
    for (int j=0; j<stateSize; j++) {
      b(i) = -bvec[i];
    }
    // Compute approximate Jacobian
    result = GMRES(f, x, x0, b, H, restart, maxit, tol);  // Solve system
    // Compute global error
    globalerr = f(x);
    for (int j=0; j<x.size(); j++) {
      dx0[j] = x[j]-x0[j];
    }
    jx = JX(f,dx0,x0);
    for (int j=0; j<dx0.size(); j++) {
      r[j] = -globalerr[j] - jx[j];
    }
    normR=0.0; normGlob=0.0;
    for (int j=0; j<stateSize; j++) {
      normR = normR + pow(r[j],2);
      normGlob = normGlob + pow(globalerr[j],2);
    }
    normR = pow(normR,0.5);
    normGlob = pow(normGlob,0.5);
    // Test to see if converged
    if (normGlob < globaltol) {
      return;
    }
  }
}
