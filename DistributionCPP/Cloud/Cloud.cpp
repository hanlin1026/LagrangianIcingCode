#include "Cloud.h"
#include <limits>

using namespace std;
using namespace Eigen;

Cloud::Cloud(State& state, Bucket& gridQT, double rhol) {
  // Set initial state of particles
  state_ = state;
  rhol_ = rhol;
  particles_ = state.size_;
  sigma_ = 75.64e-3;
  // Search grid QT for initial cell indices
  indCell_.reserve(particles_);
  double xq, yq, Xnn, Ynn, indCell;
  for (int i=0; i<particles_; i++) {
    xq = state.x_(i);
    yq = state.y_(i);
    gridQT.knnSearch(&xq,&yq,&Xnn,&Ynn,&indCell);
    indCell_[i] = indCell;
  }

}

Cloud::~Cloud() {

}

void Cloud::addParticle(State& state, Bucket& gridQT) {
  // Function to add new particles to the cloud

  // Append new state elements 
  state_.appendState(state);
  // Search grid for new state cell indices
  double xq, yq, Xnn, Ynn, indCell;
  for (int i=0; i<state.size_; i++) {
    xq = state.x_(i);
    yq = state.y_(i);
    gridQT.knnSearch(&xq,&yq,&Xnn,&Ynn,&indCell);
    indCell_.push_back(indCell);
  }
}

State Cloud::getState() {

  return state_;
}

void Cloud::findInSimulation() {
  // Function which calculates which particles are currently being advected

  // Indices of all particles in the simulation
  vector<int> indT(particles_);
  for (int i=0; i<particles_; i++) {
    indT[i] = i;
  }
  // Discount those particles which have already impinged
  int numImp = impingeTotal_.size();
  vector<int> indAdv;
  vector<int> diff;
  if (!impingeTotal_.empty()) {
    sort(impingeTotal_.begin(),impingeTotal_.end());
    set_difference(indT.begin(),indT.end(),impingeTotal_.begin(),impingeTotal_.end(),inserter(diff,diff.begin()) );
  }
  indAdv = diff;
  diff.clear();
  // Find and discount those particles which have passed the airfoil
  vector<int> indPast;
  for (int i=0; i<indAdv.size(); i++) {
    if (state_.x_(indAdv[i]) > 1) {
      indPast.push_back(indAdv[i]);
    }
  }
  set_difference(indAdv.begin(),indAdv.end(),indPast.begin(),indPast.end(),inserter(diff,diff.begin()));
  indAdv.clear();
  indAdv = diff;
  indAdv_ = indAdv;
  
}

void Cloud::computeNewCellLocations(PLOT3D& grid) {
  // Function to compute new cells occupied by particles

  // Get C = indCell(indAdv)
  vector<int> C(indAdv_.size());
  for (int i=0; i<indAdv_.size(); i++) {
    C[i] = indCell_[indAdv_[i]];
  }
  // Particle positions, cell centers
  vector<double> xP(indAdv_.size());
  vector<double> yP(indAdv_.size());
  MatrixXd xC = grid.getXCENT();
  MatrixXd yC = grid.getYCENT();
  for (int i=0; i<indAdv_.size(); i++) {
    xP[i] = state_.x_[indAdv_[i]];
    yP[i] = state_.y_[indAdv_[i]];
  }
  // Get corner points
  vector<int> N(indAdv_.size());
  vector<int> S(indAdv_.size());
  vector<int> E(indAdv_.size());
  vector<int> W(indAdv_.size());
  vector<int> SW(indAdv_.size());
  vector<int> SE(indAdv_.size());
  vector<int> NW(indAdv_.size());
  vector<int> NE(indAdv_.size());
  int nI = grid.getNX();
  int nJ = grid.getNY();
  for (int i=0; i<C.size(); i++) {
    N[i] =  C[i] + nI;
    S[i] =  C[i] - nI;
    E[i] =  C[i] + 1;
    W[i] =  C[i] - 1;
    SW[i] = S[i] - 1;
    SE[i] = S[i] + 1;
    NW[i] = N[i] - 1;
    NE[i] = N[i] + 1;
  }
  // Transform particle and neighbor positions
  vector<double> PI(indAdv_.size()); vector<double> PJ(indAdv_.size());
  vector<double> CI(indAdv_.size()); vector<double> CJ(indAdv_.size());
  vector<double> NI(indAdv_.size()); vector<double> NJ(indAdv_.size());
  vector<double> EI(indAdv_.size()); vector<double> EJ(indAdv_.size());
  vector<double> WI(indAdv_.size()); vector<double> WJ(indAdv_.size());
  vector<double> NWI(indAdv_.size()); vector<double> NWJ(indAdv_.size());
  vector<double> NEI(indAdv_.size()); vector<double> NEJ(indAdv_.size());
  for (int i=0; i<indAdv_.size(); i++) {
    grid.transformXYtoIJ(C[i],xP[i],yP[i],PI[i],PJ[i]);
    grid.transformXYtoIJ(C[i],xC(C[i]),yC(C[i]),CI[i],CJ[i]);
    grid.transformXYtoIJ(C[i],xC(N[i]),yC(N[i]),NI[i],NJ[i]);
    grid.transformXYtoIJ(C[i],xC(E[i]),yC(E[i]),EI[i],EJ[i]);
    grid.transformXYtoIJ(C[i],xC(W[i]),yC(W[i]),WI[i],WJ[i]);
    grid.transformXYtoIJ(C[i],xC(NW[i]),yC(NW[i]),NWI[i],NWJ[i]);
    grid.transformXYtoIJ(C[i],xC(NE[i]),yC(NE[i]),NEI[i],NEJ[i]);
  }
  // Find particles trying to 'glitch' through the airfoil surface
  vector<int> glitch; vector<int> noGlitch;
  for (int i=0; i<indAdv_.size(); i++) {
    if (S[i] <= 0) {
      glitch.push_back(i);
    }
    else {
      noGlitch.push_back(i);
    }
  }
  // Transform neighbor positions for non-glitching particles
  vector<double> SI(indAdv_.size()); vector<double> SJ(indAdv_.size());
  vector<double> SWI(indAdv_.size()); vector<double> SWJ(indAdv_.size());
  vector<double> SEI(indAdv_.size()); vector<double> SEJ(indAdv_.size());
  int indtmp;
  for (int i=0; i<noGlitch.size(); i++) {
    indtmp = noGlitch[i];
    grid.transformXYtoIJ(C[indtmp],xC(S[indtmp]),yC(S[indtmp]),SI[indtmp],SJ[indtmp]);
    grid.transformXYtoIJ(C[indtmp],xC(SW[indtmp]),yC(SW[indtmp]),SWI[indtmp],SWJ[indtmp]);
    grid.transformXYtoIJ(C[indtmp],xC(SE[indtmp]),yC(SE[indtmp]),SEI[indtmp],SEJ[indtmp]);
  }
  // Calculate distances in transformed plane
  vector<double> dC(indAdv_.size());
  vector<double> dN(indAdv_.size());
  vector<double> dE(indAdv_.size());
  vector<double> dW(indAdv_.size());
  vector<double> dNW(indAdv_.size());
  vector<double> dNE(indAdv_.size());
  vector<double> dS(indAdv_.size());
  vector<double> dSW(indAdv_.size());
  vector<double> dSE(indAdv_.size());
  for (int i=0; i<indAdv_.size(); i++) {
    dC[i] = pow(CI[i]-PI[i],2) + pow(CJ[i]-PJ[i],2);
    dN[i] = pow(NI[i]-PI[i],2) + pow(NJ[i]-PJ[i],2);
    dE[i] = pow(EI[i]-PI[i],2) + pow(EJ[i]-PJ[i],2);
    dW[i] = pow(WI[i]-PI[i],2) + pow(WJ[i]-PJ[i],2);
    dNW[i] = pow(NWI[i]-PI[i],2) + pow(NWJ[i]-PJ[i],2);
    dNE[i] = pow(NEI[i]-PI[i],2) + pow(NEJ[i]-PJ[i],2);
  }
  for (int i=0; i<noGlitch.size(); i++) {
    indtmp = noGlitch[i];
    dS[i] = pow(SI[indtmp]-PI[indtmp],2) + pow(SJ[indtmp]-PJ[indtmp],2);
    dSW[i] = pow(SWI[indtmp]-PI[indtmp],2) + pow(SWJ[indtmp]-PJ[indtmp],2);
    dSE[i] = pow(SEI[indtmp]-PI[indtmp],2) + pow(SEJ[indtmp]-PJ[indtmp],2);
  }
  // Set distances for glitching particles to be infinity for dS,dSW,dSE
  for (int i=0; i<glitch.size(); i++) {
    indtmp = glitch[i];
    dS[indtmp] = numeric_limits<double>::infinity();
    dSW[indtmp] = numeric_limits<double>::infinity();
    dSE[indtmp] = numeric_limits<double>::infinity();
  }
  // Find minimum distance
  vector<int> indNN(9);
  int min_index;
  vector<int> indCellNew(indAdv_.size());
  for (int i=0; i<indAdv_.size(); i++) {
    indNN[0] = C[i];
    indNN[1] = N[i];
    indNN[2] = S[i];
    indNN[3] = E[i];
    indNN[4] = W[i];
    indNN[5] = SW[i];
    indNN[6] = SE[i];
    indNN[7] = NW[i];
    indNN[8] = NE[i];
    min_index = min_element(indNN.begin(), indNN.end()) - indNN.begin();
    indCellNew[i] = indNN[min_index];
  }
  indCell_ = indCellNew;
}

void Cloud::calcDtandImpinge(Airfoil& airfoil, PLOT3D& grid) {
  // Function to set local timesteps based on CFL condition

  this->findInSimulation();
  if (!indAdv_.empty()) {
    impinge_.clear();
    this->computeNewCellLocations(grid);
    // Set timesteps
    double x,y,u,v,velMag,normVel,Lmin;
    int ind,nI;
    bool flag,flag1,flag2;
    vector<double> XYq(2);
    vector<double> XYa(2);
    vector<double> NxNy(2);
    vector<double> TxTy(2);
    nI = grid.getNX();
    dt_.reserve(indAdv_.size());
    for (int i=0; i<indAdv_.size(); i++) {
      x = state_.x_(indAdv_[i]);
      y = state_.y_(indAdv_[i]);
      u = state_.u_(indAdv_[i]);
      v = state_.v_(indAdv_[i]);
      // Calculate normal velocities
      XYq[0] = x; XYq[1] = y;
      airfoil.findPanel(XYq,XYa,NxNy,TxTy);
      normVel = u*NxNy[0] + v*NxNy[1];
      ind = indCell_[indAdv_[i]];
      flag1 = (ind-4*nI <= 0);
      flag2 = (normVel < 0);
      flag = flag1 && flag2;
      if (flag==true) {
        impinge_.push_back(indAdv_[i]);
      }
      // Set timesteps based on CFL condition
      velMag = sqrt(pow(u,2) + pow(v,2));
      Lmin = grid.getLMIN(indAdv_[i]);
      dt_[i] = 0.5*Lmin/velMag;
      
    }
  }

}




void Cloud::setIndAdv(vector<int>& indAdv) {
  
  indAdv = indAdv_; 
}

vector<int> Cloud::getIndAdv() {

  return indAdv_;
}
