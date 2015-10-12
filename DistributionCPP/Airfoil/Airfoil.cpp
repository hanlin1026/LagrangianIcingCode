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
    normal_(i,0) = -tangent_(i,1);
    normal_(i,1) = tangent_(i,0);
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
  int i1 = floor(0.30*I);
  int i2 = floor(0.70*I);
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

void Airfoil::growIce(vector<double>& sTHERMO, vector<double>& mice, double DT, const char* strSurf) {
  // Function to update XY grid coordinates based on ice growth rate for DT time interval

  vector<double> indAIRFOIL;
  vector<double> indTHERMO;
  // Interpolate s-coordinates from thermo onto airfoil grid
  double minimum;
  int ind;
  double s_min,s_max;
  if (strcmp(strSurf,"UPPER")==0) {
    s_min = 0.0;
    s_max = 0.4;
  }
  else if (strcmp(strSurf,"LOWER")==0) {
    s_min = -0.4;
    s_max = 0.0;
  }
  vector<double> tmp1(sTHERMO.size());
  vector<double> tmp2(sTHERMO.size());
  double sCoord;
  for (int i=0; i<panelS_.size(); i++) {
    sCoord = panelS_(i) - stagPt_;
    if ((sCoord >= s_min) && (sCoord <= s_max)) {
      for (int j=0; j<sTHERMO.size(); j++)
	tmp1[j] = sTHERMO[j] - sCoord;
      tmp2 = abs(tmp1);
      minimum = min(tmp2,ind);
      indTHERMO.push_back(ind);
      indAIRFOIL.push_back(i);
    }
  }
  // Displace each grid point along its normal vector
  double dt = 1.0;
  double rhoICE = 917.0;
  double time = 0.0;
  double xNEW,yNEW,dH;
  while (time <= DT) {
    for (int i=0; i<indAIRFOIL.size(); i++) {
      dH = mice[indTHERMO[i]]*dt/rhoICE;
      xNEW = panelX_(indAIRFOIL[i]) + dH*normal_(indAIRFOIL[i],0);
      yNEW = panelY_(indAIRFOIL[i]) + dH*normal_(indAIRFOIL[i],1);
      panelX_(indAIRFOIL[i]) = xNEW;
      panelY_(indAIRFOIL[i]) = yNEW;
    }
    time += dt;
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
