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

void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (ABS(dy) > ABS(dx)) {
    double temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    double temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}
 
void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
  double temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

void Update(std::vector<double>& x, int k, Eigen::MatrixXd& h, std::vector<double> &s, std::vector<std::vector<double>>& v)
{
  std::vector<double> y(s.size());
  for (int i=0; i<s.size(); i++)
    y[i] = s[i];
  std::vector<double> dx(x.size());
  for (int i=0; i<x.size(); i++) {
    dx[i] = 0;
  }

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[i] /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y[j] -= h(j,i) * y[i];
  }

  std::vector<double> vtmp;
  for (int j = 0; j <= k; j++) {
    vtmp = v[j];
    for (int jj = 0; jj<vtmp.size(); jj++) {
      dx[jj] += vtmp[jj] * y[j];
    }
  }
  for (int i=0; i<x.size(); i++) {
    x[i] += dx[i];
  }
}



std::vector<double> operator*(std::vector<double> vec, double scal) {
  for (int i=0; i<vec.size(); i++) {
    vec[i] *= scal;
  }
  return vec;
}



int GMRES(ThermoEqns* thermo, int balFlag,
	  std::vector<double>& x, std::vector<double>& u0, std::vector<double>& b,
	  Eigen::MatrixXd& H, int& m, int& max_iter, double& tol)
{
  double resid;
  int i, j = 1, k;
  std::vector<double> s(m+1), cs(m+1), sn(m+1);
  
  double normb = NORM(b);
  std::vector<double> jx;
  jx = thermo->JX(balFlag,x,u0);
  std::vector<double> r(b.size());
  for (int ii=0; ii<jx.size(); ii++) {
    r[ii] = b[ii] - jx[ii];
  }
  double beta = NORM(r);
  
  if (normb == 0.0)
    normb = 1;
  
  if ((resid = NORM(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  std::vector<std::vector<double>> v(m+1,vector<double>(r.size()));
  std::vector<double> vtmp;
  std::vector<double> vv(m+1);
  std::vector<double> w(b.size());
  jx.clear();
  while (j <= max_iter) {
    printf("GMRESiter = %d\n",j);
    v[0] = r * (1.0 / beta);    // ??? r / beta
    for (int jj=0; jj<s.size(); jj++)
      s[jj] = 0.0;
    s[0] = beta;
    
    for (i = 0; i < m && j <= max_iter; i++, j++) {
      printf("I = %d\n",i);
      vtmp = v[i];
      for (int ii=0; ii<vtmp.size(); ii++) {
	vv[ii] = vtmp[ii];
      }
      jx = thermo->JX(balFlag,vv,u0);
      for (int ii=0; ii<jx.size(); ii++) {
	w[ii] = jx[ii];
      }
      for (k = 0; k <= i; k++) {
	vtmp = v[k];
        H(k, i) = DOT(w, vtmp);
	for (int kk=0; kk<vtmp.size(); kk++) {
	  w[kk] -= H(k,i)*vtmp[kk];
	}
      }
      H(i+1, i) = NORM(w);
      v[i+1] = w * (1.0 / H(i+1, i)); // ??? w / H(i+1, i)

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      
      if ((resid = ABS(s[i+1]) / normb) < tol) {
        Update(x, i, H, s, v);
        tol = resid;
        max_iter = j;
        return 0;
      }
    }
    Update(x, m - 1, H, s, v);
    jx.clear();
    jx = thermo->JX(balFlag,x,u0);
    for (int ii=0; ii<jx.size(); ii++) {
      r[ii] = b[ii] - jx[ii];
    }
    beta = NORM(r);
    if ((resid = beta / normb) < tol) {
      tol = resid;
      max_iter = j;
      return 0;
    }
  }
  
  tol = resid;
  return 1;
}



