//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <eigen3/Eigen/Dense>

double ABS(double x)
{
  return (x > 0 ? x : -x);
}

double NORM(std::vector<double>& x) {
  double L2norm = 0.0;
  for (int i=0; i<x.size(); i++) {
    L2norm += pow(x[i],2);
  }
  return sqrt(L2norm);
}

double DOT(std::vector<double>& a, std::vector<double>& b) {
  double dotProd = 0.0;
  for (int i=0; i<a.size(); i++) {
    dotProd += a[i]*b[i];
  }
  return dotProd;

}

std::vector<double> operator*(std::vector<double> vec, double scal) {
  for (int i=0; i<vec.size(); i++) {
    vec[i] *= scal;
  }
  return vec;
}

std::vector<double> operator/(std::vector<double> vec, double scal) {
  for (int i=0; i<vec.size(); i++) {
    vec[i] *= (1.0/scal);
  }
  return vec;
}

std::vector<double> operator+(std::vector<double> a, std::vector<double> b) {
  std::vector<double> c(a.size());
  for (int i=0; i<a.size(); i++)
    c[i] = a[i] + b[i];
  return c;
}

std::vector<double> operator-(std::vector<double> a, std::vector<double> b) {
  std::vector<double> c(a.size());
  for (int i=0; i<a.size(); i++)
    c[i] = a[i] - b[i];
  return c;
}

void SETEQ(std::vector<double> vec, double num) {
  for (int i; i<vec.size(); i++)
    vec[i] = num;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int scalarsign(d) {
  int sign = sgn(d);
  if (sign == 0)
    sign = 1;
}


int GMRES(ThermoEqns* thermo, int balFlag,
	  std::vector<double>& x, std::vector<double>& u0, std::vector<double>& b,
	  Eigen::MatrixXd& H, int& restart, int& max_iter, double& tol) {
  // Main routine: GMRES implemented using Householder transformations

  int stateSize = x.size();
  //Set up for the method
  int flag = 1;
  Eigen::VectorXd xmin = x;            // Iterate which has minimal residual so far
  int imin = 0;                        // "Outer" iteration at which xmin was computed
  int jmin = 0;                        // "Inner" iteration at which xmin was computed
  int evalxm = 0;
  int stag = 0;
  int moresteps = 0;
  int maxmsteps = std::min([std::floor(stateSize/50),5,stateSize-maxit]);
  int maxstagsteps = 3;
  int minupdated = 0;
  int inner = restart;
  int outer = max_iter;
  double beta;
  std::vector<double> jx;
  // Compute initial residual
  jx = thermo->JX(balFlag,x,u0);
  r = b - jx;
  double normR = NORM(r);
  double normB = NORM(b);
  if (normR == 0.0)
    normR = 1;
  // If initial guess is already good enough, we are done
  if ((resid = normR / normB) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  // Preallocations for the method
  std::vector<double> resvec(inner*outer+1,1);          // Norm of residuals
  resvec[0] = normR;                                    // resvec[0] = norm(b-A*x0)
  double normrmin = normR;                              // Norm of residual from xmin
  // Preallocate J to hold the Given's rotation constants.
  Eigen::MatrixXd J(2,inner);
  Eigen::MatrixXd U(stateSize,inner);
  Eigen::MatrixXd R(inner,inner);
  std::vector<double> w(inner+1);
  std::vector<double> u(stateSize);
  std::vector<double> v(stateSize);
  double tmpv;
  double Uv = 0.0;
  double alpha;
  double rho;
  double normr_act;
  // Begin main GMRES routine
  for (int outiter=0; outiter<outer; outiter++) {
    // Construct u for Householder reflector
    // u = r + sign(r[0])*||r||*e1
    u = r;
    normR = NORM(r);
    beta = scalarsign(r[0])*normR;
    u[0] += beta;
    u = u / norm(u);
    for (int i=0; i<stateSize; i++)
      U(i,0) = u[i];
    // Apply Householder projection to r.
    // w = r - 2*u*u'*r;
    w[0] = -beta;
    for (int initer=0; initer<inner; initer++) {
      // Form P1*P2*...*Pj*ej
      // v = Pj*ej = ej - 2*u*u'*ej
      v = u*(-2*u[initer]);
      v[initer] += 1;
      // v = P1*P2*...Pjm1*(Pj*ej)
      for (int k=initer-1; k>0; k--) {
	Uv = 0.0;
	for (int kk = 0; kk<stateSize; kk++)
	  Uv += U(kk,k)*v[kk];
	for (int kk=0; kk<stateSize; kk++)
	  v[kk] -= U(kk,k)*2*Uv;
      }
      // Explicitly normalize v to reduce the effects of round-off
      v = v/NORM(v);
      // Apply Jacobian to v
      jx.clear();
      jx = thermo->JX(balFlag,v,u0);
      // Form Pj*Pj-1*...P1*Av.
      for (int k=0; k<initer; k++) {
	Uv = 0.0;
	for (int kk=0; kk<stateSize; kk++)
	  Uv += U(kk,k)*v[kk];
	for (int kk=0; kk<stateSize; kk++)
	  v[kk] -= U(kk,k)*2*Uv;
      }
      // Determine Pj+1
      if (initer != v.size()-1) {
	// Construct u for Householder reflector Pj+1
	for (int i=0; i<initer; i++)
	  u[i] = 0;
	for (int i=initer; i<stateSize; i++)
	  u[i] = v[i];
	alpha = NORM(u);
	if (alpha != 0) {
	  alpha = scalarsign(v[initer])*alpha;
	  // u = v(initer+1:end) +
          //     sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
	  u[initer] += alpha;
	  u = u/NORM(u);
	  for (int i=0; i<stateSize; i++)
	    U(i,initer) = u[i];
	  // Apply Pj+1 to v
	  for (int i=initer+1; i<stateSize; i++)
	    v[i] = 0.0;
	  v[initer] = -alpha;
	}
      }
      // Apply Given's rotations to the newly formed v
      for (int colJ=0; colJ<initer-1; colJ++) {
	tmpv = v[colJ];
        v[colJ]   = J(1,colJ)*v[colJ] + J(2,colJ)*v[colJ+1];
        v[colJ+1] = -J(2,colJ)*tmpv + J(1,colJ)*v[colJ+1];
      }
      // Compute Given's rotation Jm.
      if (initer != v.size()-1) {
	rho = sqrt(pow(v[initer-1],2) + pow(v[initer],2));
	J(0,initer-1) = v[initer-1]/rho;
	J(1,initer-1) = v[initer]/rho;
	w[initer] = -2*J(1,initer-1)*w[initer-1];
	w[initer-1] = J(0,initer-1)*w[initer-1];
	v[initer-1] = rho;
	v[initer] = 0;
      }
      for (int i=0; i<inner; i++)
	R(i,initer-1) = v[i];
      normR = ABS(w[initer]);
      resvec[(outiter-1)*inner+initer] = normR;
      normr_act = normR;




    }

  }



}
