#include "Airfoil.h"
#include <cmath>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>
#include <VectorOperations/VectorOperations.h>

using namespace std;

Airfoil::Airfoil(std::vector<double>& X, std::vector<double>& Y) {
  // Constructor: takes grid points X,Y
  
  // Set panel center points, tangent/normal vectors
  int gridPts = X.size();
  panelX_.resize(gridPts-1);
  panelY_.resize(gridPts-1);
  tangent_.resize(gridPts-1,2);
  normal_.resize(gridPts-1,2);
  double ds_x, ds_y;
  for (int i=0; i<gridPts-1; i++) {
    ds_x = X[i+1]-X[i];
    ds_y = Y[i+1]-Y[i];
    panelX_(i) = X[i] + 0.5*ds_x;
    panelY_(i) = Y[i] + 0.5*ds_y;
    tangent_(i,0) = ds_x/sqrt( pow(ds_x,2) + pow(ds_y,2) );
    tangent_(i,1) = ds_y/sqrt( pow(ds_x,2) + pow(ds_y,2) );
    //normal_(i,0) = -tangent_(i,1);
    //normal_(i,1) = tangent_(i,0);
  }
  // 3-point moving average normal vectors
  normal_(0,0) = -tangent_(0,1); normal_(0,1) = tangent_(0,0);
  normal_(gridPts-2,0) = -tangent_(gridPts-2,1); normal_(gridPts-2,1) = tangent_(gridPts-2,0);
  for (int i=1; i<gridPts-3; i++) {
    normal_(i,0) = -(tangent_(i-1,1)+tangent_(i,1)+tangent_(i+1,1))/3.0;
    normal_(i,1) = (tangent_(i-1,0)+tangent_(i,0)+tangent_(i+1,0))/3.0;
  }
  // Set quadtree search object bounds
  double minX, minY, maxX, maxY;
  minX = panelX_.minCoeff() - 0.1; maxX = panelX_.maxCoeff() + 0.1;
  minY = panelY_.minCoeff() - 0.1; maxY = panelY_.maxCoeff() + 0.1;
  double SW[2] = {minX, minY};
  double SE[2] = {maxX, minY};
  double NW[2] = {minX, maxY};
  double NE[2] = {maxX, maxY};
  panelSearcher_.setBounds(&SW[0],&SE[0],&NW[0],&NE[0]);
  // Create quadtree search object
  panelSearcher_.calcQuadTree(panelX_.data(),panelY_.data(),panelX_.rows());
  // Calculate s-coordinates of panel points
  this->calcSCoords();
  // Output to file
  FILE* fout = fopen("AirfoilXY.out","w");
  FILE* foutT = fopen("AirfoilTxTy.out","w");
  FILE* foutN = fopen("AirfoilNxNy.out","w");
  for (int i=0; i<gridPts-1; i++) {
    fprintf(fout,"%f\t%f\n",panelX_[i],panelY_[i]);
    fprintf(foutT,"%f\t%f\n",tangent_(i,0),tangent_(i,1));
    fprintf(foutN,"%f\t%f\n",normal_(i,0),normal_(i,1));
  }
  fclose(fout); fclose(foutT); fclose(foutN);
}

Airfoil::~Airfoil() {
  
}

void Airfoil::findPanel(std::vector<double>& XYq, std::vector<double>& XYnn, std::vector<double>& NxNy, std::vector<double>& TxTy) {
  // Function to return (x,y) coordinates and normal/tangential 
  // vectors of the point on the airfoil closest to the query

  double xq = XYq[0];
  double yq = XYq[1];
  double xnn,ynn;
  int indnn;
  panelSearcher_.knnSearch(&xq,&yq,&xnn,&ynn,&indnn);
  XYnn[0] = xnn; XYnn[1] = ynn;
  NxNy[0] = normal_(indnn,0); 
  NxNy[1] = normal_(indnn,1);
  TxTy[0] = tangent_(indnn,0);
  TxTy[1] = tangent_(indnn,1);

}

void Airfoil::findPanel(std::vector<double>& XYq, std::vector<double>& XYnn, std::vector<double>& NxNy, std::vector<double>& TxTy, int& indexNN) {
  // Function to return (x,y) coordinates and normal/tangential 
  // vectors of the point on the airfoil closest to the query
  // AS WELL AS the index of the closest panel point

  double xq = XYq[0];
  double yq = XYq[1];
  double xnn,ynn;
  int indnn;
  panelSearcher_.knnSearch(&xq,&yq,&xnn,&ynn,&indnn);
  XYnn[0] = xnn; XYnn[1] = ynn;
  NxNy[0] = normal_(indnn,0); 
  NxNy[1] = normal_(indnn,1);
  TxTy[0] = tangent_(indnn,0);
  TxTy[1] = tangent_(indnn,1);
  indexNN = indnn;

}

double Airfoil::calcIncidenceAngle(std::vector<double>& XYq, std::vector<double>& UVq, int indNN) {
  // Function to calculate the incidence angle of a droplet impinging
  // on the airfoil surface

  // Find closest panel
  std::vector<double> XYa(2);
  std::vector<double> NxNy(2);
  std::vector<double> TxTy(2);
  std::vector<double> velUnitNorm(2);
  NxNy[0] = normal_(indNN,0);
  NxNy[1] = normal_(indNN,1);
  // Find angle between airfoil surface normal vector and query velocity
  double velNorm = sqrt(pow(UVq[0],2) + pow(UVq[1],2));
  velUnitNorm[0] = UVq[0]/velNorm;
  velUnitNorm[1] = UVq[1]/velNorm;
  double projection = velUnitNorm[0]*(-NxNy[0]) + velUnitNorm[1]*(-NxNy[1]);
  double theta = acos(projection);
  // Incidence angle is defined as the angle between the surface and
  // the velocity
  theta = abs(M_PI/2.0) - abs(theta);

  return theta;

}

void Airfoil::calcSCoords() {
  // Function to calculate s-coordinates of panel center points

  double dx,dy,ds;
  int numPts = panelX_.size();
  panelS_.resize(numPts);
  panelS_(0) = 0;
  for (int i=0; i<numPts-1; i++) {
    dx = panelX_(i+1) - panelX_(i);
    dy = panelY_(i+1) - panelY_(i);
    ds = sqrt(pow(dx,2) + pow(dy,2));
    panelS_(i+1) = panelS_(i) + ds;
  }

}

double Airfoil::interpXYtoS(std::vector<double>& XYq) {
  // Function to approximate s-coords of query pt (x,y) on airfoil surface
  // NOTE: assumes directionality of tangent vectors!

  double xq = XYq[0]; double yq = XYq[1];
  // Find closest panel point
  std::vector<double> XYa(2);
  std::vector<double> NxNy(2);
  std::vector<double> TxTy(2);
  std::vector<double> velUnitNorm(2);
  int indNN;
  this->findPanel(XYq,XYa,NxNy,TxTy,indNN);
  // Find coordinates of query point in panel frame
  double dx = xq - XYa[0];
  double dy = yq - XYa[1];
  // Project displacement vector onto tangent vector
  double tangDisplacement = dx*TxTy[0] + dy*TxTy[1];
  double sCoord;
  if (tangDisplacement > 0) {
    sCoord = panelS_(indNN) + sqrt(pow(dx,2) + pow(dy,2));
  }
  else {
    sCoord = panelS_(indNN) - sqrt(pow(dx,2) + pow(dy,2));
  }
  // Estimate s-coordinate of query point by adding tangent displacement
  // sCoord = panelS_(indNN) + tangDisplacement;

  return sCoord;

}

void Airfoil::appendFilm(double sCoord, double mass) {
  FilmScoords_.push_back(sCoord);
  FilmMass_.push_back(mass);

}

void Airfoil::calcCollectionEfficiency(double fluxFreeStream,double dS) {
  // Function to calculate collection efficiency of airfoil

  if (!FilmScoords_.empty()) {
    // Create bins
    double minS = *min_element(FilmScoords_.begin(),FilmScoords_.end());
    double maxS = *max_element(FilmScoords_.begin(),FilmScoords_.end());
    // double dS = (maxS-minS)/(numBins-1);
    int numBins = (maxS-minS)/dS + 1;
    gsl_histogram *h = gsl_histogram_alloc(numBins);
    // Histogram count of mass in each bin
    gsl_histogram_set_ranges_uniform(h,minS,maxS);
    for (int i=0; i<FilmMass_.size(); i++) {
      gsl_histogram_accumulate(h,FilmScoords_[i],FilmMass_[i]);
    }
    BetaBins_.resize(numBins);
    Beta_.resize(numBins);
    double upper,lower,cent,mass,fluxLocal;
    for (int i=0; i<numBins; i++) {
      gsl_histogram_get_range(h,i,&lower,&upper);
      cent = 0.5*(upper+lower);
      BetaBins_[i] = cent-stagPt_;
      mass = gsl_histogram_get(h,i);
      fluxLocal = mass/dS;
      Beta_[i] = fluxLocal/fluxFreeStream;
    }

  }
  
}

void Airfoil::calcStagnationPt(PLOT3D& grid) {
  // Function to calculate the stagnation point

  int I = grid.getNX()-1;
  Eigen::MatrixXf U = grid.getUCENT();
  Eigen::MatrixXf V = grid.getVCENT();
  Eigen::MatrixXd X = grid.getXCENT();
  Eigen::MatrixXd Y = grid.getYCENT();
  double u; double v;
  int i1 = floor(0.45*I);
  int i2 = floor(0.55*I);
  vector<double> VelMagSq(i2-i1);
  // Get first wrap of (u,v)
  for (int i=i1; i<i2; i++) {
    u = U(i,0); v = V(i,0);
    VelMagSq[i-i1] = pow(u,2) + pow(v,2);
  }
  // Minimize velMag to find stagnation point
  int indMin = min_element(VelMagSq.begin(),VelMagSq.end()) - VelMagSq.begin();
  stagPtX_ = X(indMin+i1,0);
  stagPtY_ = Y(indMin+i1,0);
  std::vector<double> XYstag(2);
  XYstag[0] = stagPtX_; XYstag[1] = stagPtY_;
  stagPt_ = this->interpXYtoS(XYstag);
  printf("stagPtS = %f, stagPtX = %f, stagPtY = %f\n",stagPt_,stagPtX_,stagPtY_);

}

void Airfoil::correctJagged(int id1, int id2, int id3, int id4) {
  // Function to correct jaggedness of 4 points by using an
  // area-preserving trapezoid

  // Calculate rotation of the p1-p4 line segment
  double segX = panelX_(id4)-panelX_(id1);
  double segY = panelY_(id4)-panelY_(id1);
  double d3 = sqrt(pow(segX,2) + pow(segY,2));
  segX = segX/d3;
  segY = segY/d3;
  double theta = atan2(segY,segX);
  // Calculations for area-preserving trapezoid
  double d1 = sqrt(pow(panelX_(id2)-panelX_(id1),2) + pow(panelY_(id2)-panelY_(id1),2));
  double d2 = sqrt(pow(panelX_(id3)-panelX_(id2),2) + pow(panelY_(id3)-panelY_(id2),2));
  double A = 0.5*((panelX_(id3)-panelX_(id1))*(panelY_(id2)-panelY_(id4)) + (panelX_(id4)-panelX_(id2))*(panelY_(id3)-panelY_(id1)));
  double B = (panelX_(id3)-panelX_(id1))*(panelY_(id2)-panelY_(id4)) + (panelX_(id4)-panelX_(id2))*(panelY_(id3)-panelY_(id1));
  double G = pow(sqrt(4*pow(B,6)+pow(B,4))-pow(B,2),0.3333);
  double r;
  if (A != 0) {
    r = (1./3./sqrt(2.))*sqrt(-6.*pow(2.,0.3333)*pow(B,2)/G + (3.*pow(2.,0.6667))*G + 2.) + 
      0.5*sqrt(4.*pow(2.,0.3333)*pow(B,2.)/(3*G) - 2./3.*pow(2.,0.6667)*G + 
      8.*sqrt(2)/9./sqrt(-6.*pow(2.,0.3333)*pow(B,2)/G + 3*pow(2.,0.6667)*G + 2.) + 8./9.) - 2./3.;
  }
  else
    r = 0.3333;
  printf("R = %lf\n",r);
  double d = d3*r;
  double alpha = acos(0.5*(d3-d)/d);
  if (A < 0)
    alpha = -alpha;
  vector<double> px_n(4);
  vector<double> py_n(4);
  px_n[0] = 0.0; px_n[1] = 0.5*(d3-d);   px_n[2] = 0.5*(d3-d)+d; px_n[3] = panelX_[id1];
  py_n[0] = 0.0; py_n[1] = d*sin(alpha); py_n[2] = d*sin(alpha); py_n[3] = panelY_[id4]-panelY_[id1];
  for (int i=0; i<4; i++)
    py_n[i] += panelY_[id1];
  // Coordinate rotation according to p1-p4 line segment
  vector<double> pxT_n(4);
  vector<double> pyT_n(4);
  for (int i=0; i<4; i++) {
    pxT_n[i] = cos(theta)*px_n[i] - sin(theta)*py_n[i];
    pyT_n[i] = sin(theta)*px_n[i] + cos(theta)*py_n[i];
  }
  // Update panels
  panelX_[id1] = pxT_n[0]; panelX_[id2] = pxT_n[1]; panelX_[id3] = px_n[2]; panelX_[id4] = px_n[3];
  panelY_[id1] = pyT_n[0]; panelY_[id2] = pyT_n[1]; panelY_[id3] = py_n[2]; panelY_[id4] = py_n[3];

}

double Airfoil::computeJaggednessCriterion(int id1, int id2, int id3, int id4) {
  // Compute jaggedness criterion

  // Compute line segment vectors connecting points
  double dx12 = panelX_[id2] - panelX_[id1];
  double dy12 = panelY_[id2] - panelY_[id1];
  double norm12 = sqrt(pow(dx12,2) + pow(dy12,2));
  dx12 = dx12/norm12; dy12 = dy12/norm12;
  double dx23 = panelX_[id3] - panelX_[id2];
  double dy23 = panelY_[id3] - panelY_[id2];
  double norm23 = sqrt(pow(dx23,2) + pow(dy23,2));
  dx23 = dx23/norm23; dy23 = dy23/norm23;
  double dx34 = panelX_[id4] - panelX_[id3];
  double dy34 = panelY_[id4] - panelY_[id3];
  double norm34 = sqrt(pow(dx34,2) + pow(dy34,2));
  dx34 = dx34/norm34; dy34 = dy34/norm34;
  // Compute inner product
  double ip1 = dx12*dx23 + dy12*dy23;
  double ip2 = dx23*dx34 + dy23*dy34;
  // Compute cross product (l12 x l23)
  double cp1 = dx12*dy23 - dy12*dx23;
  double cp2 = dx23*dy34 - dy23*dx34;
  // Compute signed turning angle
  double alpha23 = atan2(cp1,ip1);
  double alpha34 = atan2(cp2,ip2);
  // Compute jaggedness criterion
  double JAG = std::abs(alpha23 - alpha34)*180.0/M_PI;

  return JAG;

}


void Airfoil::growIce(vector<double>& sTHERMO, vector<double>& mice, double DT, double chord, const char* strSurf) {
  // Function to update XY grid coordinates based on ice growth rate for DT time interval

  vector<double> indAIRFOIL;
  vector<double> indTHERMO;
  double minimum;
  int ind;
  double s_min,s_max;

  // Match s-coordinates of airfoil to those from finely-resolved thermo calculation
  if (strcmp(strSurf,"UPPER")==0) {
    s_min = 0.0;
    s_max = 0.4;
  }
  else if (strcmp(strSurf,"LOWER")==0) {
    s_min = -0.4;
    s_max = 0.0;
  }
  else if (strcmp(strSurf,"ENTIRE")==0) {
    s_min = -0.4;
    s_max = 0.4;
  }
  vector<double> tmp1(sTHERMO.size());
  vector<double> tmp2(sTHERMO.size());
  double sCoord;
  for (int i=0; i<panelS_.size(); i++) {
    sCoord = panelS_(i) - stagPt_;
    if ( ((sCoord >= s_min) && (sCoord <= s_max)) ) {
      for (int j=0; j<sTHERMO.size(); j++)
	tmp1[j] = sTHERMO[j] - sCoord;
      tmp2 = abs(tmp1);
      minimum = min(tmp2,ind);
      indTHERMO.push_back(ind);
      indAIRFOIL.push_back(i);
    }
  }

  // Calculate ice growth
  double rhoICE = 917.0;
  double xNEW,yNEW,dH,dH_old,dH_tmp,ip,theta;
  vector<double> DH(indAIRFOIL.size());
  double ds = sTHERMO[1]-sTHERMO[0];
  for (int i=0; i<indAIRFOIL.size(); i++) {
    dH = mice[indTHERMO[i]]*DT/rhoICE;
    DH[i] = dH;
  }
  // Moving average smoothing
  // vector<double> DH_new = movingAverage(DH,4.0);
  // DH = DH_new;

  // Laplacian smooth normal vectors
  double eps = 0.0;
  int NL = indAIRFOIL.size();
  vector<double> normalX(NL);
  vector<double> normalY(NL);
  for (int i=0; i<NL; i++) {
    normalX[i] = normal_(indAIRFOIL[i],0);
    normalY[i] = normal_(indAIRFOIL[i],1);
  }
  vector<vector<double> > LAPL_NX(NL,vector<double>(NL));
  vector<vector<double> > LAPL_NY(NL,vector<double>(NL));
  vector<double> NX_smooth(NL);
  vector<double> NY_smooth(NL);
  LAPL_NX = LaplacianMatrix(NL,eps);
  LAPL_NY = LaplacianMatrix(NL,eps);
  NX_smooth = tridiagSolve(LAPL_NX,normalX);
  NY_smooth = tridiagSolve(LAPL_NY,normalY);

  // Correction for area oblation
  vector<double> DH_area(NL);
  DH_area[0] = DH[0]; DH_area[1] = DH[0]; DH_area[2] = DH[0];
  dH_old = 0.0;
  for (int i=1; i<NL; i++) {
    sCoord = panelS_(indAIRFOIL[i]) - stagPt_;
    ip = normal_(indAIRFOIL[i],0)*normal_(indAIRFOIL[i]-1,0) + normal_(indAIRFOIL[i],1)*normal_(indAIRFOIL[i]-1,1);
    //ip = NX_smooth[i]*NX_smooth[i-1] + NY_smooth[i]*NY_smooth[i-1];
    theta = acos(ip);
    //DH_area[i] = DH[i]*2*ds/(2*ds + DH[i-1]*sin(theta));
    //DH_area[i] = DH[i]*2*ds/(2*ds + DH_area[i-1]*sin(theta));
    DH_area[i] = DH[i];
  }

  // Implicit Laplacian smoothing
  eps = 5.0;
  vector<vector<double> > LAPL_DH(NL,vector<double>(NL));
  LAPL_DH = LaplacianMatrix(NL,eps);
  vector<double> DH_smooth(NL);
  DH_smooth = tridiagSolve(LAPL_DH,DH_area);
  // Displace each grid point along its normal vector
  double NX,NY;
  for (int i=0; i<NL; i++) {
    if ( (strcmp(strSurf,"UPPER")==0) && (i < 2) )
      dH = DH_smooth[2];
    else if ( (strcmp(strSurf,"LOWER")==0) && (i > NL-3) )
      dH = DH_smooth[NL-3];
    else
      dH = DH_smooth[i];
    NX = NX_smooth[i];
    NY = NY_smooth[i];
    xNEW = panelX_(indAIRFOIL[i]) + dH*NX;
    yNEW = panelY_(indAIRFOIL[i]) + dH*NY;
    panelX_(indAIRFOIL[i]) = xNEW;
    panelY_(indAIRFOIL[i]) = yNEW;
  }

}


vector<double> Airfoil::getBetaBins() {

  return BetaBins_;

}

vector<double> Airfoil::getBeta() {

  return Beta_;

}

vector<double> Airfoil::getX() {
  vector<double> X(panelX_.rows());
  for (int i=0; i<X.size(); i++)
    X[i] = panelX_(i);
  return X;
}

vector<double> Airfoil::getY() {
  vector<double> Y(panelY_.rows());
  for (int i=0; i<Y.size(); i++)
    Y[i] = panelY_(i);
  return Y;
}

void Airfoil::setStagPt(double sLoc) {
  stagPt_ = sLoc;

}

double Airfoil::getStagPt() {
  return stagPt_;
}



vector<double> Airfoil::tridiagSolve(vector<vector<double>>& A, vector<double>& r) {
  // Tridiagonal matrix solver; solves A*x = r

  int N = r.size();
  // Get band-diagonal components of matrix A
  vector<double> a(N);
  vector<double> b(N);
  vector<double> c(N);
  a[0]   = 0; b[0]   = A[0][0];     c[0]   = 0;
  a[N-1] = 0; b[N-1] = A[N-1][N-1]; c[N-1] = 0;
  for (int i=1; i<N-1; i++) {
    a[i] = A[i][i-1];
    b[i] = A[i][i];
    c[i] = A[i][i+1];
  }
  // Main algorithm
  vector<double> x(N);
  double bet;
  vector<double> gam(N);
  // Decomposition and forward substitution
  bet = b[0];
  x[0] = r[0]/bet;
  for (int i=1; i<N; i++) {
    gam[i] = c[i-1]/bet;
    bet = b[i] - a[i]*gam[i];
    x[i] = (r[i] - a[i]*x[i-1])/bet;
  }
  // Backsubstitution
  for (int i=N-2; i>-1; i--) {
    x[i] -= gam[i+1]*x[i+1];
  }

  return x;

} 


vector<vector<double> > Airfoil::LaplacianMatrix(int N, double eps) {
  // Subroutine to compute a Laplacian matrix

  vector<vector<double> > LAPL(N,vector<double>(N));
  // Compute Laplacian matrix
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      LAPL[i][j] = 0.0;
    }
  }
  LAPL[0][0] = 1 + 2*eps;
  LAPL[0][1] = -eps;
  LAPL[N-1][N-1] = 1 + 2*eps;
  LAPL[N-1][N-2] = -eps;
  for (int i=1; i<N-1; i++) {
    LAPL[i][i-1] = -eps;
    LAPL[i][i]   = 1 + 2*eps;
    LAPL[i][i+1] = -eps;
  }

  return LAPL;

} 


vector<double> Airfoil::movingAverage(vector<double>& X, double smooth) {
  // Subroutine to perform a moving average

  int N = X.size();
  vector<double> X_smooth(N);
  double NX_smooth,NY_smooth;
  int startIND = (int) smooth/2;
  int endIND   = N - startIND;
  for (int i=0; i<N; i++) {
    X_smooth[i] = 0.0;
    if ((i > startIND) && (i < endIND)) {
      for (int j=0; j<smooth+1; j++) {
  	X_smooth[i] += X[i+j-startIND]/(smooth+1);
      }
    }
    else
      X_smooth[i] = X[i];
  }

  return X_smooth;
}
