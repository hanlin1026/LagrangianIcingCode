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

template<class Real> 
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (abs(dy) > abs(dx)) {
    Real temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    Real temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}


template<class Real> 
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  Real temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}



template < class Matrix, class Vector >
void 
Update(std::vector<double>& x, int k, Matrix &h, Vector &s, Vector v[])
{
  Vector y(s);
  Vector dx(x.size());
  for (int i=0; i<x.size(); i++) {
    dx(i) = 0;
  }

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y(i) /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y(j) -= h(j,i) * y(i);
  }

  for (int j = 0; j <= k; j++)
    dx += v[j] * y(j);
  for (int i=0; i<x.size(); i++) {
    x[i] += dx(i);
  }
}


template < class Real >
Real 
abs(Real x)
{
  return (x > 0 ? x : -x);
}

// Define action of Jacobian on vector
inline vector<double> Jx(vector<double> (*f)(vector<double>& Xq), vector<double>& X, vector<double>& u0) {
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

template < class Vector, class Matrix, class Real >
int 
  GMRES(std::vector<double> (*f)(std::vector<double>& xq),
	std::vector<double>& x, std::vector<double>& u0, const Vector &b,
	Matrix &H, int &m, int &max_iter, Real &tol)
{
  Real resid;
  int i, j = 1, k;
  Vector s(m+1), cs(m+1), sn(m+1), w;
  
  Real normb = norm(b);
  std::vector<double> jx;
  jx = Jx(f,x,u0);
  Vector r;
  for (int ii=0; ii<jx.size(); ii++) {
    r(ii) = b(ii) - jx[ii];
  }
  Real beta = norm(r);
  
  if (normb == 0.0)
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  Vector *v = new Vector[m+1];
  Vector vtmp;
  std::vector<double> vv(m+1);
  while (j <= max_iter) {
    v[0] = r * (1.0 / beta);    // ??? r / beta
    s = 0.0;
    s(0) = beta;
    
    for (i = 0; i < m && j <= max_iter; i++, j++) {
      vtmp = v[i];
      for (int ii=0; ii<vtmp.size(); ii++) {
	vv[ii] = vtmp(ii);
      }
      jx = Jx(f,vv,u0);
      for (int ii=0; ii<jx.size(); ii++) {
	w(ii) = jx[ii];
      }
      for (k = 0; k <= i; k++) {
        H(k, i) = dot(w, v[k]);
        w -= H(k, i) * v[k];
      }
      H(i+1, i) = norm(w);
      v[i+1] = w * (1.0 / H(i+1, i)); // ??? w / H(i+1, i)

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
      
      if ((resid = abs(s(i+1)) / normb) < tol) {
        Update(x, i, H, s, v);
        tol = resid;
        max_iter = j;
        delete [] v;
        return 0;
      }
    }
    Update(x, m - 1, H, s, v);
    jx = Jx(f,x,u0);
    for (int ii=0; ii<jx.size(); ii++) {
      r(ii) = b(ii) - jx[ii];
    }
    beta = norm(r);
    if ((resid = beta / normb) < tol) {
      tol = resid;
      max_iter = j;
      delete [] v;
      return 0;
    }
  }
  
  tol = resid;
  delete [] v;
  return 1;
}



