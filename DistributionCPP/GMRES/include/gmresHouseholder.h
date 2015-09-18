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

std::vector<double> operator+(std::vector<double> a, Eigen::VectorXd b) {
  std::vector<double> c(a.size());
  for (int i=0; i<a.size(); i++)
    c[i] = a[i] + b(i);
  return c;

}

void SETEQ(std::vector<double> vec, double num) {
  for (int i; i<vec.size(); i++)
    vec[i] = num;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int scalarsign(double d) {
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
  std::vector<double> xmin = x;        // Iterate which has minimal residual so far
  int imin = 0;                        // "Outer" iteration at which xmin was computed
  int jmin = 0;                        // "Inner" iteration at which xmin was computed
  int evalxm = 0;
  int stag = 0;
  int moresteps = 0;
  int stateSizeFloor = std::floor(stateSize/50);
  int min1 = std::min(stateSizeFloor,5);
  int maxmsteps = std::min(min1,stateSize-max_iter);
  int maxstagsteps = 3;
  int minupdated = 0;
  int inner;
  if (restart < stateSize)
    inner = restart;
  else
    inner = stateSize;
  int outer = max_iter;
  double beta;
  int initerFINAL = 0;
  int outiterFINAL = 0;
  std::vector<double> jx;
  double eps = 1.0e-6;
  double n2b = NORM(b);
  double tolb = tol*n2b;
  if (tol < eps)
    tol = eps;
  // Compute initial residual
  jx = thermo->JX(balFlag,x,u0);
  std::vector<double> r = b - jx;
  double normR = NORM(r);
  double normB = NORM(b);
  if (normR == 0.0)
    normR = 1;
  // If initial guess is already good enough, we are done
  double resid = normR/normB;
  if (resid <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  // Preallocations for the method
  std::vector<double> resvec(inner*outer+1);            // Norm of residuals
  resvec[0] = normR;                                    // resvec[0] = norm(b-A*x0)
  double normrmin = normR;                              // Norm of residual from xmin
  // Preallocate J to hold the Given's rotation constants.
  Eigen::MatrixXd J(2,inner);
  for (int i=0; i<2; i++) {
    for (int j=0; j<inner; j++)
      J(i,j) = 0.0;
  }
  Eigen::MatrixXd U(stateSize,inner);
  for (int i=0; i<stateSize; i++) {
    for (int j=0; j<inner; j++)
      U(i,j) = 0.0;
  }
  Eigen::MatrixXd R(inner,inner);
  for (int i=0; i<inner; i++) {
    for (int j=0; j<inner; j++)
      R(j,i) = 0.0;
  }
  std::vector<double> w(inner+1);
  std::vector<double> u(stateSize);
  std::vector<double> v(stateSize);
  std::vector<double> xm(stateSize);
  std::vector<double> additive(stateSize);
  std::vector<double> addvc;
  Eigen::VectorXd ytmp;
  Eigen::VectorXd wtmp;
  Eigen::VectorXd addvc1;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> linearSolver;
  Eigen::MatrixXd RTMP;
  double tmpv;
  double Uv = 0.0;
  double alpha;
  double rho;
  double normr_act;
  int idx;
  // Begin main GMRES routine
  for (int outiter=0; outiter<outer; outiter++) {
    // Construct u for Householder reflector
    // u = r + sign(r[0])*||r||*e1
    u = r;
    normR = NORM(r);
    beta = scalarsign(r[0])*normR;
    u[0] += beta;
    u = u / NORM(u);
    for (int i=0; i<stateSize; i++) {
      U(i,0) = u[i];
    }
    // Apply Householder projection to r.
    // w = r - 2*u*u'*r;
    w[0] = -beta;
    for (int initer=0; initer<inner; initer++) {
      // Form P1*P2*...*Pj*ej
      // v = Pj*ej = ej - 2*u*u'*ej
      v = u*(-2*u[initer]);
      v[initer] += 1;
      // v = P1*P2*...Pjm1*(Pj*ej)
      for (int k=initer-1; k>-1; k--) {
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
      v = jx;
      // Form Pj*Pj-1*...P1*Av.
      for (int k=0; k<initer+1; k++) {
	Uv = 0.0;
	for (int kk=0; kk<stateSize; kk++)
	  Uv += U(kk,k)*v[kk];
	for (int kk=0; kk<stateSize; kk++)
	  v[kk] -= U(kk,k)*2*Uv;
      }
      // Determine Pj+1
      if (initer != v.size()-1) {
	// Construct u for Householder reflector Pj+1
	for (int i=0; i<initer+1; i++)
	  u[i] = 0;
	for (int i=initer+1; i<stateSize; i++)
	  u[i] = v[i];
	alpha = NORM(u);
	if (alpha != 0) {
	  alpha = scalarsign(v[initer+1])*alpha;
	  // u = v(initer+1:end) +
          //     sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
	  u[initer+1] += alpha;
	  u = u/NORM(u);
	  for (int i=0; i<stateSize; i++)
	    U(i,initer+1) = u[i];
	  // Apply Pj+1 to v
	  for (int i=initer+2; i<stateSize; i++)
	    v[i] = 0.0;
	  v[initer+1] = -alpha;
	}
      }
      // Apply Given's rotations to the newly formed v
      for (int colJ=0; colJ<initer-1; colJ++) {
	tmpv = v[colJ];
        v[colJ]   = J(0,colJ)*v[colJ] + J(1,colJ)*v[colJ+1];
        v[colJ+1] = -J(1,colJ)*tmpv + J(0,colJ)*v[colJ+1];
      }
      // Compute Given's rotation Jm.
      if (initer != v.size()-1) {
	rho = sqrt(pow(v[initer],2) + pow(v[initer+1],2));
	J(0,initer) = v[initer]/rho;
	J(1,initer) = v[initer+1]/rho;
	w[initer+1] = -J(1,initer)*w[initer];
	w[initer] = J(0,initer)*w[initer];
	v[initer] = rho;
	v[initer+1] = 0;
      }
      for (int i=0; i<inner; i++)
	R(i,initer) = v[i];
      normR = ABS(w[initer+1]);
      resvec[(outiter)*inner+initer+1] = normR;
      normr_act = normR;
      // MORE DEBUGGING STARTING HERE
      if ((normR <= tolb) || (stag>= maxstagsteps) || moresteps) {
	if (evalxm == 0) {
	  wtmp.resize(initer+1);
	  for (int i=0; i<initer+1; i++)
	    wtmp(i) = w[i];
	  linearSolver.compute(R.block(0,0,initer+1,initer+1));
	  ytmp = linearSolver.solve(wtmp);
	  for (int i=0; i<stateSize; i++) {
	    additive[i] = U(i,initer)*(-2*ytmp(initer))*U(initer,initer);
	  }
	  additive[initer] += ytmp(initer);
	  for (int k=initer-1; k>-1; k--) {
	    additive[k] += ytmp(k);
	    Uv = 0.0;
	    for (int kk=0; kk<stateSize; kk++)
	      Uv += U(kk,k)*additive[kk];
	    for (int kk=0; kk<stateSize; kk++)
	      additive[kk] += -U(kk,k)*(2*Uv);
	  }
	  if (NORM(additive) < eps*NORM(x))
	    stag += 1;
	  else
	    stag = 0;
	  xm = x + additive;
	  evalxm = 1;
	}
	else if (evalxm == 1) {
	  linearSolver.compute(R.block(0,0,initer,initer));
	  addvc1.resize(initer-1);
	  addvc.resize(initer);
	  addvc1 = -1*linearSolver.solve(R.block(0,initer,initer,initer+1));
	  addvc1 *= w[initer]/R(initer,initer);
	  for (int i=0; i<initer-1; i++)
	    addvc[i] = addvc1(i);
	  addvc[initer] = w[initer]/R(initer,initer);
	  if (NORM(addvc) < eps*NORM(xm))
	    stag += 1;
	  else
	    stag = 0;
	  for (int i=0; i<stateSize; i++)
	    additive[i] = U(i,initer)*(-2*addvc[initer]*U(initer,initer));
	  additive[initer] += addvc[initer];
	  for (int k=initer-1; k>-1; k--) {
	    additive[k] += addvc[k];
	    Uv = 0.0;
	    for (int kk=0; kk<stateSize; kk++)
	      Uv += U(kk,k)*additive[kk];
	    for (int kk=0; kk<stateSize; kk++)
	      additive[kk] += -U(kk,k)*(2*Uv);
	  }
	  xm = xm + additive;
	}
	jx.clear();
	jx = thermo->JX(balFlag,xm,u0);
	r = b - jx;
	if (NORM(r) <= tol*n2b) {
	  x = xm;
	  flag = 0;
	  break;
	}
	normr_act = NORM(r);
	resvec[(outiter)*inner+initer+1] = normr_act;

	if (normr_act <= normrmin) {
	  normrmin = normr_act;
	  imin = outiter;
	  jmin = initer;
	  xmin = xm;
	  minupdated = 1;
	}
	if (normr_act <= tolb) {
	  x = xm;
	  flag = 0;
	  break;
	}
	else {
	  if ((stag >= maxstagsteps) && (moresteps == 0))
	    stag = 0;
	  moresteps += 1;
	  if (moresteps >= maxmsteps) {
	    flag = 3;
	    break;
	  }
	}
      }
      if (normr_act <= normrmin) {
	normrmin = normr_act;
	imin = outiter;
	jmin = initer;
	minupdated = 1;
      }
      if (stag >= maxstagsteps) {
	flag = 3;
	break;
      }
      initerFINAL = initer;
    } // End inner loop
    evalxm = 0;
    if (flag != 0) {
      if (minupdated == 1)
	idx = jmin;
      else
	idx = initerFINAL;
      wtmp.resize(idx+1);
      for (int i=0; i<idx+1; i++)
	wtmp(i) = w[i];
      linearSolver.compute(R.block(0,0,idx+1,idx+1));
      ytmp = linearSolver.solve(wtmp);
      for (int i=0; i<stateSize; i++)
	additive[i] = U(i,idx)*(-2*ytmp(idx)*U(idx,idx));
      additive[idx] += ytmp(idx);
      for (int k=idx-1; k>-1; k--) {
	additive[k] += ytmp(k);
	Uv = 0.0;
	for (int kk=0; kk<stateSize; kk++)
	  Uv += U(kk,k)*additive[kk];
	for (int kk=0; kk<stateSize; kk++)
	  additive[kk] += -U(kk,k)*(-2*Uv);
      }
      x = x + additive;
      xmin = x;
      jx.clear();
      jx = thermo->JX(balFlag,x,u0);
      r = b - jx;
      normr_act = NORM(r);      
    }
    if (normr_act <= normrmin) {
      xmin = x;
      normrmin = normr_act;
      imin = outiter;
      jmin = initerFINAL;
    }
    if (flag == 3)
      break;
    if (normr_act <= tolb) {
      flag = 0;
      break;
    }
    minupdated = 0;
    outiterFINAL = outiter;
  } // Ends outer loop

  // Returned solution is that with minimum residual
  if (flag != 0) {
    x = xmin;
  }
  // Truncate the zeros from resvec
  std::vector<double> resvecTMP;
  if ((flag <= 1) || (flag == 3)) {
    resvecTMP.resize((outiterFINAL)*inner+initerFINAL+1);
    for (int i=0; i<resvecTMP.size(); i++) {
      if (resvec[i] != 0)
	resvecTMP.push_back(resvec[i]);
    }
  }
  else {
    if (initerFINAL == 0) {
      resvecTMP.resize((outiterFINAL)*inner+initerFINAL+1);
      for (int i=0; i<resvecTMP.size(); i++)
	resvecTMP[i] = resvec[i];
    }
    else {
      resvecTMP.resize((outiterFINAL)*inner+initerFINAL+initerFINAL);
      for (int i=0; i<resvecTMP.size(); i++)
	resvecTMP[i] = resvec[i];
    }
  }
  resvec.clear();
  resvec = resvecTMP;

}
