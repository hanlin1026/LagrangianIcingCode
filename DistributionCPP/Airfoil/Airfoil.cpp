#include "Airfoil.h"
#include <cmath>
#include <stdlib.h>

using namespace std;

Airfoil::Airfoil(Eigen::VectorXd& X, Eigen::VectorXd& Y) {
  // Constructor: takes grid points X,Y
  
  // Set panel center points, tangent/normal vectors
  int gridPts = X.size();
  panelX_.resize(gridPts-1);
  panelY_.resize(gridPts-1);
  tangent_.resize(gridPts-1,2);
  normal_.resize(gridPts-1,2);
  double ds_x, ds_y;
  for (int i=0; i<gridPts-1; i++) {
    ds_x = X(i+1)-X(i);
    ds_y = Y(i+1)-Y(i);
    panelX_(i) = X(i) + 0.5*ds_x;
    panelY_(i) = Y(i) + 0.5*ds_y;
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
  // Output to file
  FILE* fout = fopen("AirfoilXY.out","w");
  for (int i=0; i<gridPts-1; i++) {
    fprintf(fout,"%f\t%f\n",panelX_[i],panelY_[i]);
  }
  fclose(fout);
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

double Airfoil::calcIncidenceAngle = calcIncidenceAngle(XYq,UVq) {
  // Function to calculate the incidence angle of a droplet impinging
  // on the airfoil surface

  // Find closest panel
  std::vector<double> XYa(2);
  std::vector<double> NxNy(2);
  std::vector<double> TxTy(2);
  this->findPanel(XYq,XYa,NxNy,TxTy);
  // Find angle between airfoil surface normal vector and query velocity
  double velNorm = sqrt(pow(UVq[0],2) + pow(UVq[1],2));
  std::vector<double> vel(2);
  velUnitNorm[0] = UVq[0]/velNorm;
  velUnitNorm[1] = UVq[1]/velNorm;
  double projection = velUnitNorm[0]*NxNy[0] + velUnitNorm[1]*NxNy[1];
  double theta = acos(projection);

  return theta;

}
