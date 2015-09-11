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

ThermoEqns::ThermoEqns(const char* filenameCHCF, const char* filenameBETA, Airfoil& airfoil) {
  // Constructor to read in input files and initialize thermo eqns

  this->interpUpperSurface(filenameCHCF,airfoil,"CHCF");
  this->interpUpperSurface(filenameBETA,airfoil,"BETA");
  // Set rhoL_,muL_
  rhoL_ = 1000.0;
  muL_ = 1.787e-3;
  // Set LWC_,Uinf_
  LWC_ = 1.0;
  Uinf_ = 100;
  // Initial guess for ice accretion
  mice_upper_.resize(NPts_);
  for (int i=0; i<NPts_; i++)
    mice_upper_[i] = 0.0;

}

void ThermoEqns::interpUpperSurface(const char* filename, Airfoil& airfoil, const char* parameter) {
  // Function to interpolate upper surface

  MatrixXd data; VectorXd s; 
  VectorXd beta;
  VectorXd ch; VectorXd cf;
  if (strcmp(parameter,"BETA") == 0) {
    // Import (s,beta)
    data = this->readBetaXY(filename);
    s = data.col(0);
    beta = data.col(1);
  }
  else if (strcmp(parameter,"CHCF") == 0) {
    // Import (s,ch,cf)
    data = this->readCHCF(filename);
    s = data.col(0);
    ch = data.col(1);
    cf = data.col(2);
    // Set stagPt at s=0
    double stagPt = airfoil.getStagPt();
    for (int i=0; i<s.size(); i++) {
      s(i) -= stagPt;
    }
  }
  // Find stagnation point
  double zero = 1.0; double zero_new = 1.0;
  int zero_ind = 0;
  for (int i=0; i<s.size(); i++) {
    zero_new = abs(s(i));
    if (zero_new < zero) {
      zero = zero_new;
      zero_ind = i;
    }
  }
  // Split into upper/lower surface grids
  int NPts_orig = s.size() - zero_ind;
  vector<double> s_upper_orig(NPts_orig);
  for (int i=0; i<NPts_orig; i++) {
    s_upper_orig[i] = s(zero_ind+i);
  }
  // Refine upper surface grid
  NPts_ = 1000;
  s_upper_.resize(NPts_);
  //double s_max = s.maxCoeff();
  double s_max = 0.4;
  double ds = (s_max-zero)/NPts_;
  for (int i=0; i<NPts_; i++) {
    s_upper_[i] = zero + i*ds;
  }
  // Interpolate parameter values on grids
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  if (strcmp(parameter,"CHCF") == 0) {    
    cH_upper_.resize(NPts_);
    cF_upper_.resize(NPts_);
    gsl_spline *splineCH = gsl_spline_alloc(gsl_interp_linear, NPts_orig);
    gsl_spline *splineCF = gsl_spline_alloc(gsl_interp_linear, NPts_orig);
    gsl_spline_init(splineCH, &s_upper_orig[0], &ch[zero_ind], NPts_orig);
    gsl_spline_init(splineCF, &s_upper_orig[0], &cf[zero_ind], NPts_orig);
    for (int i=0; i<NPts_; i++) {
      if (s_upper_[i] <= s.maxCoeff()) { 
	cH_upper_[i] = gsl_spline_eval(splineCH, s_upper_[i], acc);
        cF_upper_[i] = gsl_spline_eval(splineCF, s_upper_[i], acc);
      }
      else {
	cH_upper_[i] = 0.0;
        cF_upper_[i] = 0.0;
      }
    }
  }
  else if (strcmp(parameter,"BETA") == 0) {
    beta_upper_.resize(NPts_);
    gsl_spline *splineBETA = gsl_spline_alloc(gsl_interp_linear, NPts_orig);
    gsl_spline_init(splineBETA, &s_upper_orig[0], &beta[zero_ind], NPts_orig);
    for (int i=0; i<NPts_; i++) {
      if (s_upper_[i] <= s.maxCoeff())
	beta_upper_[i] = gsl_spline_eval(splineBETA, s_upper_[i], acc);
      else
	beta_upper_[i] = 0.0;
    }
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
  MatrixXd s_Beta(sizeBeta,2);
  double a,b;
  for (int i=0; i<sizeBeta; i++) {
    fscanf(filept,"%lf,%lf",&a,&b);
    s_Beta(i,0) = a; s_Beta(i,1) = b;
  }
  // Close file streams
  fclose(filept);

  return s_Beta;

}

// Define action of Jacobian on vector
vector<double> ThermoEqns::JX(int func, vector<double>& X, vector<double>& u0) {
  vector<double> jx(u0.size());
  double eps = 1.e-4;
  vector<double> X2(u0.size());
  for (int i=0; i<u0.size(); i++) {
    X2[i] = u0[i] + eps*X[i];
  }
  // Determine which mass/energy balance to use
  vector<double> f1;
  vector<double> f2;
  if (func==0) {
    f2 = massBalanceUpper(X2);
    f1 = massBalanceUpper(u0);
  }
  else if (func==2) {
    f2 = testBalance(X2);
    f1 = testBalance(u0);
  }
  for (int i=0; i<u0.size(); i++) {
    jx[i] = (1./eps)*(f2[i]-f1[i]);
  }
  return jx;
}

vector<double> ThermoEqns::massBalanceUpper(vector<double>& x) {
  // Function to compute mass balance

  vector<double> err(x.size());
  vector<double> F(x.size());
  vector<double> f(x.size()-1);
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
  for (int i=0; i<x.size()-1; i++) {
    xFACE[i] = 0.5*(x[i]+x[i+1]);
    cfFACE[i] = 0.5*(cF_upper_[i]+cF_upper_[i+1]);
    DF[i] = (1/muL_)*(xFACE[i]*cfFACE[i]);
    f[i] = 0.5*(F[i]+F[i+1]) - 0.5*std::abs(DF[i])*(x[i+1]-x[i]);
  }
  // Calculate error for internal cells
  double ds,mimp;
  for (int i=1; i<x.size()-1; i++) {
    ds = s_upper_[i+1]-s_upper_[i];
    mimp = beta_upper_[i]*LWC_*Uinf_;
    D_flux[i-1] = f[i]-f[i-1];
    I_sources[i-1] = (1./rhoL_)*ds*(mimp-mice_upper_[i]);
    err[i] = D_flux[i-1] - I_sources[i-1];
  }
  // Boundary conditions
  err[0] = x[0]-0;
  err[x.size()-1] = err[x.size()-2];

  return err;
}

vector<double> ThermoEqns::testBalance(vector<double>& x) {
  // Test function RHS

  vector<double> RHS(3);
  RHS[0]=1.0; RHS[1]=-1.0; RHS[2]=2.0;
  vector<double> err(3);
  Eigen::MatrixXd A(3,3);
  A << 8,1,6,3,5,7,4,9,2;
  err[0] = A(0,0)*x[0] + A(0,1)*x[1] + A(0,2)*x[2] - RHS[0];
  err[1] = A(1,0)*x[0] + A(1,1)*x[1] + A(1,2)*x[2] - RHS[1];
  err[2] = A(2,0)*x[0] + A(2,1)*x[1] + A(2,2)*x[2] - RHS[2];

  return err;
}




void ThermoEqns::NewtonKrylovIteration(const char* balance, vector<double>& u0) {
  // Function to take a balance of form f(X) = 0 and do Newton-Krylov iteration

  // Set balance flag
  int balFlag;
  if (strcmp(balance,"MASS")==0)
    balFlag = 0;
  else if (strcmp(balance,"TEST")==0)
    balFlag = 2;

  double tol = 1.e-3;                       // Convergence tolerance
  int result, maxit = 1000, restart = 5;    // GMRES Maximum, restart iterations

  // Initialize Jacobian and RHS, solution vectors
  int stateSize = u0.size();
  vector<double> b(stateSize);
  vector<double> x(stateSize);
  vector<double> x0(stateSize);
  vector<double> dx0(stateSize);
  vector<double> jx;
  // Storage for upper Hessenberg H
  MatrixXd H(restart+1, restart);
  // Begin iteration
  int nitermax = 20;
  vector<double> globalerr; double globaltol = 2.0e-5;
  vector<double> r(stateSize);
  double normR,normGlob;
  // Initialize linearization point
  vector<double> un = u0;
  
  for (int i=0; i<nitermax; i++) {
    // Reset linearization point
    u0 = un;
    x = u0;
    // Compute RHS
    if (balFlag==0)
      b = massBalanceUpper(u0);
    else if (balFlag==2)
      b = testBalance(u0);
    b = b*-1.0;
    // Compute approximate Jacobian
    result = GMRES(this, balFlag, x, u0, b, H, restart, maxit, tol);  // Solve system
    un = u0 + x;
    // Compute global error
    if (balFlag==0)
      globalerr = massBalanceUpper(un);
    else if (balFlag==2)
      globalerr = testBalance(un);
    jx = JX(balFlag,x,u0);
    r = globalerr*-1.0 + jx*-1.0;
    normR = NORM(r);
    normGlob = NORM(globalerr);
    printf("JFNK ERROR = %lf\tGLOBAL ERROR = %lf\n",normR,normGlob);
    // Test to see if converged
    if (normGlob < globaltol) {
      break;
    }
  }
  for (int ii=0; ii<stateSize; ii++)
    printf("%e\t%e\n",s_upper_[ii],un[ii]);
}
