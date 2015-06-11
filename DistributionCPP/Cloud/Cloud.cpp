#include "Cloud.h"
#include <math.h>
#include <limits>

using namespace std;
using namespace Eigen;

Cloud::Cloud(State& state, PLOT3D& grid, double rhol) {
  // Set initial state of particles
  state_ = state;
  rhoL_ = rhol;
  particles_ = state.size_;
  sigma_ = 75.64e-3;
  // Search grid QT for initial cell indices
  indCell_.reserve(particles_);
  double xq, yq, Xnn, Ynn;
  int indCell;
  for (int i=0; i<particles_; i++) {
    xq = state.x_(i);
    yq = state.y_(i);
    grid.pointSearch(xq,yq,Xnn,Ynn,indCell);
    indCell_[i] = indCell;
  }

}

Cloud::~Cloud() {

}

void Cloud::addParticle(State& state, PLOT3D& grid) {
  // Function to add new particles to the cloud

  // Append new state elements 
  state_.appendState(state);
  // Search grid for new state cell indices
  double xq, yq, Xnn, Ynn;
  int indCell;
  for (int i=0; i<state.size_; i++) {
    xq = state.x_(i);
    yq = state.y_(i);
    grid.pointSearch(xq,yq,Xnn,Ynn,indCell);
    indCell_.push_back(indCell);
  }
}

State Cloud::getState() {

  return state_;
}

vector<int> Cloud::getIMPINGE() {
  return impinge_;
}
vector<int> Cloud::getIMPINGETOTAL() {
  return impingeTotal_;
}
vector<int> Cloud::getINDCELL() {
  return indCell_;
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
  else {
    diff = indT;
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
  vector<double> indDNN(9);
  int min_index;
  vector<int> indCellNew(particles_);
  for (int i=0; i<particles_; i++) {
    indCellNew[i] = indCell_[i];
  }
  for (int i=0; i<indAdv_.size(); i++) {
    indDNN[0] = dC[i];  indNN[0] = C[i];
    indDNN[1] = dN[i];  indNN[1] = N[i];
    indDNN[2] = dS[i];  indNN[2] = S[i];
    indDNN[3] = dE[i];  indNN[3] = E[i];
    indDNN[4] = dW[i];  indNN[4] = W[i];
    indDNN[5] = dSW[i]; indNN[5] = SW[i];
    indDNN[6] = dSE[i]; indNN[6] = SE[i];
    indDNN[7] = dNW[i]; indNN[7] = NW[i];
    indDNN[8] = dNE[i]; indNN[8] = NE[i];
    min_index = min_element(indDNN.begin(), indDNN.end()) - indDNN.begin();
    indCellNew[indAdv_[i]] = indNN[min_index];
  }
  indCell_.clear();
  indCell_ = indCellNew;
  //printf("CELL = %d\n",indCell_[0]);
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
      Lmin = grid.getLMIN(ind);
      dt_[i] = 0.5*Lmin/velMag;
      
    }
    //printf("V = %f\tLMIN = %f\n",velMag,Lmin);
  }

}

void Cloud::transportSLD(PLOT3D& grid) {
  // Function to advect droplets

  if (!indAdv_.empty()) {
    double Tinf = grid.getTINF();
    double x,y,u,v,r,dt;
    double xnn,ynn;
    int indnn;
    double rhoG,uG,vG,muG;
    double xnp1,ynp1,unp1,vnp1;
    double Re,CD,tau,m;
    double g = -9.81;     // m/s
    double C1 = 1.458e-6; // kg/(ms*sqrt(K))
    double S = 110.4;     // K
    for (int i=0; i<indAdv_.size(); i++) {
      x = state_.x_(indAdv_[i]);
      y = state_.y_(indAdv_[i]);
      u = state_.u_(indAdv_[i]);
      v = state_.v_(indAdv_[i]);
      r = state_.r_(indAdv_[i]);
      dt = dt_[i];
      // Get fluid properties at nearest neighbor cell
      rhoG = grid.getRHOCENT(indCell_[indAdv_[i]]);
      uG = grid.getUCENT(indCell_[indAdv_[i]]);
      vG = grid.getVCENT(indCell_[indAdv_[i]]);
      // Sutherland's law
      muG = C1*pow(Tinf,1.5)/(Tinf+S);
      // Force parameter calculations
      Re = 2*rhoG*r/muG*sqrt( pow(u-uG,2) + pow(v-vG,2) );
      CD = 24/Re*(1 + 0.15*pow(Re,0.687));
      if (Re==0) {
        tau = 0;
      }
      else {
        tau = 24/Re/CD*(2*rhoL_*pow(r,2)/9/muG);
      }
      m = 4/3*rhoL_*M_PI*pow(r,3);
      // Advection equations
      xnp1 = x + uG*dt + (u-uG)*(1-exp(-dt/tau))*tau;
      ynp1 = y + vG*dt + (v-vG)*(1-exp(-dt/tau))*tau + (dt-(1-exp(-dt/tau))*tau)*tau*g;
      unp1 = uG + exp(-dt/tau)*(u-uG);
      vnp1 = vG + exp(-dt/tau)*(v-vG) + (1-exp(-dt/tau))*tau*g;
      // Update particle states
      state_.x_(indAdv_[i]) = xnp1;
      state_.y_(indAdv_[i]) = ynp1;
      state_.u_(indAdv_[i]) = unp1;
      state_.v_(indAdv_[i]) = vnp1;
    }
  }

}

void Cloud::computeImpingementRegimes(Airfoil& airfoil) {
  // Function to divide impinging droplets into 3 classes
  // (ie. bounce, spread, splash)
  
  double x,y,u,v,r,t,temp,muL;
  vector<double> XYq(2);
  vector<double> XYnn(2);
  vector<double> NxNy(2);
  vector<double> TxTy(2);
  vector<int> splashSpread;
  double vNormSq,vTang,We,Oh,K;
  double Ks0,Kb0,wsr,wsf,wbr,wbf,hr,hf,R,delta,R_tilda,fs,fb;
  // Parameters used in the impingement regime calculation
  Ks0 = 3000.0; Ks0_ = Ks0;
  Kb0 = 600.0;  Kb0_ = Kb0;
  wsr = 20.0/3.0;
  wsf = 5.0/6.0;
  wbr = 32.0;
  wbf = 1.0;
  hr = 20e-6;
  hf = 0.0;
  // Clear/size vectors as appropriate
  K_.reserve(impinge_.size());
  fs_.reserve(impinge_.size());
  fb_.reserve(impinge_.size());
  vNormSq_.reserve(impinge_.size());
  vTang_.reserve(impinge_.size());
  bounce_.clear();
  spread_.clear();
  splash_.clear();
  for (int i=0; i<impinge_.size(); i++) {
    x = state_.x_(impinge_[i]);
    y = state_.y_(impinge_[i]);
    u = state_.u_(impinge_[i]);
    v = state_.v_(impinge_[i]);
    r = state_.r_(impinge_[i]);
    t = state_.time_(impinge_[i]);
    temp = state_.temp_(impinge_[i]);
    muL = (2.414e-5)*pow( 10 , 247.8/(temp-140) );
    // Find local points of impingement, normal vectors
    XYq[0] = x; XYq[1] = y;
    airfoil.findPanel(XYq,XYnn,NxNy,TxTy);
    // Compute normal velocities at airfoil surface
    vNormSq = pow( u*NxNy[0] + v*NxNy[1], 2 );
    vNormSq_[i] = vNormSq;
    vTang = u*TxTy[0] + v*TxTy[1];
    vTang_[i] = vTang;
    We = 2*rhoL_*r*vNormSq/sigma_;
    Oh = sqrt( muL/2/rhoL_/sigma_/r );
    K  = We*pow(Oh,-0.4);
    K_[i] = K;
    // Intermediate parameter calculations
    R = hr/2/r; // Dimensionless wall roughness height
    delta = hf/2/r; // Dimensionless film thickness
    R_tilda = pow(R,2)/(R+delta); // Modified wall roughness due to presence of film
    fs = (1 + pow(R_tilda,2))*(1 + pow(delta,2))/(1 + wsr*pow(R_tilda,2))/(1 + wsf*pow(delta,2));
    fb = (1 + pow(R_tilda,2))*(1 + pow(delta,2))/(1 + wsr*pow(R_tilda,2))/(1 + wbr*pow(R_tilda,4))/(1 + wbf*pow(delta,2));
    fs_[i] = fs;
    fb_[i] = fb;
    // Calculate impingement regime
    if (K < Kb0*fb) {
      bounce_.push_back(i);
    }
    else if ((K > Kb0*fb) && (K < Ks0*fs)) {
      spread_.push_back(i);
      splashSpread.push_back(impinge_[i]);
    }
    else {
      splash_.push_back(i);
      splashSpread.push_back(impinge_[i]);
    } 
  }
  // Update impingeTotal
  vector<int> diff;
  sort(impingeTotal_.begin(),impingeTotal_.end());
  sort(splashSpread.begin(),splashSpread.end());
  set_difference(splashSpread.begin(),splashSpread.end(),impingeTotal_.begin(),impingeTotal_.end(),inserter(diff,diff.begin()) );
  for (int i=0; i<diff.size(); i++) {
    impingeTotal_.push_back(diff[i]);
  }

}

void Cloud::bounceDynamics(Airfoil& airfoil) {
  // Function to compute bounce dynamics
  
  if (!bounce_.empty()) {
    double x,y,u,v,r;
    double K,Ks,Kb,vNormSq;
    double vN,vT;
    double vNorm,vTang,uNew,vNew;
    vector<double> XYq(2);
    vector<double> XYa(2);
    vector<double> NxNy(2);
    vector<double> TxTy(2);
    int indBounce;
    for (int i=0; i<bounce_.size(); i++) {
      indBounce = impinge_[bounce_[i]];
      x = state_.x_(indBounce);
      y = state_.y_(indBounce);
      u = state_.u_(indBounce);
      v = state_.v_(indBounce);
      r = state_.r_(indBounce);
      K = K_[bounce_[i]];
      Ks = fs_[bounce_[i]]*Ks0_;
      Kb = fb_[bounce_[i]]*Kb0_;
      vNormSq = vNormSq_[bounce_[i]];
      XYq[0] = x; XYq[1] = y;
      airfoil.findPanel(XYq,XYa,NxNy,TxTy);
      vNorm = sqrt(vNormSq);
      vTang = vTang_[bounce_[i]];
      // Calculate post-impact velocities
      vN = 4*vNorm*( sqrt(K/Kb) - K/Kb );
      vT = 0.8*vTang;
      uNew = vN*NxNy[0] + vT*TxTy[0];
      vNew = vN*NxNy[1] + vT*TxTy[1];
      // Set new bouncing velocity
      state_.u_[indBounce] = uNew;
      state_.v_[indBounce] = vNew;
    }
  }
}

void Cloud::splashDynamics(Airfoil& airfoil) {
  // Function to compute splash dynamics

  if (!splash_.empty()) {
    double x,y,u,v,r;
    double K,Ks,Kb,vNormSq;
    double vN,vT;
    double vNorm,vTang,uNew,vNew;
    double theta,sCoord;
    vector<double> XYq(2);
    vector<double> UVq(2);
    vector<double> XYa(2);
    vector<double> NxNy(2);
    vector<double> TxTy(2);
    int indSplash;
    // Index over each splashing parcel
    for (int i=0; i<splash_.size(); i++) {
      // Get splashing parcel properties
      indSplash = impinge_[splash_[i]];
      x = state_.x_(indSplash);
      y = state_.y_(indSplash);
      u = state_.u_(indSplash);
      v = state_.v_(indSplash);
      r = state_.r_(indSplash);
      K = K_[splash_[i]];
      Ks = fs_[splash_[i]]*Ks0_;
      vNormSq = vNormSq_[splash_[i]];
      vTang = vTang_[splash_[i]];
      XYq[0] = x; XYq[1] = y;
      UVq[0] = u; UVq[1] = v;
      // Calculate s-coords of impinging parcel
      sCoord = airfoil.interpXYtoS(XYq);
      // Calculate impinging incidence angle
      airfoil.findPanel(XYq,XYa,NxNy,TxTy);
      theta = airfoil.calcIncidenceAngle(XYq,UVq);

    }
  }
}

void Cloud::setIndAdv(vector<int>& indAdv) {
  
  indAdv = indAdv_; 
}

vector<int> Cloud::getIndAdv() {

  return indAdv_;
}
