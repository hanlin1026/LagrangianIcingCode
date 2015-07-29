#include "Cloud.h"
#include <math.h>
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
using namespace Eigen;

Cloud::Cloud(State& state, PLOT3D& grid, double rhol, ParcelScalars& PARCEL) {
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
  // String specifying whether or not to use splashing at all
  if (PARCEL.SplashFlag_ == 1) {
    SplashFlag_ = true;
  }
  else {
    SplashFlag_ = false;
  }
  // String specifying whether or not to track splash particles
  if (PARCEL.TrackSplashFlag_ == 1) {
    TrackSplashParticles_ = true;
  }
  else {
    TrackSplashParticles_ = false;
  }

}

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

void Cloud::addParticles(State& state, int indCell) {
  // Function to add new particles to the cloud

  // Append new state elements 
  state_.appendState(state);
  // Initialize particles as being in same cell as parent
  for (int i=0; i<state.size_; i++) {
    indCell_.push_back(indCell);
  }
  particles_ += state.size_;
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
  hr = 20.0e-6;
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
    muL = (2.414e-5)*pow( 10.0 , 247.8/(temp-140.0) );
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
    R = (hr/2.0)/r; // Dimensionless wall roughness height
    delta = (hf/2.0)/r; // Dimensionless film thickness
    R_tilda = pow(R,2)/(R+delta); // Modified wall roughness due to presence of film
    fs = (1 + pow(R_tilda,2))*(1 + pow(delta,2))/(1 + wsr*pow(R_tilda,2))/(1 + wsf*pow(delta,2));
    fb = (1 + pow(R_tilda,2))*(1 + pow(delta,2))/(1 + wsr*pow(R_tilda,2))/(1 + wbr*pow(R_tilda,4))/(1 + wbf*pow(delta,2));
    fs_[i] = fs;
    fb_[i] = fb;
    // Calculate impingement regime
    if (K < Kb0*fb) {
      if (SplashFlag_ == true) {
	bounce_.push_back(i);
      }
      else {
	// If we are not using splashing, turn off bouncing by treating it as a spread
	spread_.push_back(i);
	splashSpread.push_back(impinge_[i]);
      }
    }
    else if ((K > Kb0*fb) && (K < Ks0*fs)) {
      spread_.push_back(i);
      splashSpread.push_back(impinge_[i]);
    }
    else {
      if (SplashFlag_ == true) {
	splash_.push_back(i);
      }
      else {
	// If we are not using splashing, turn it off by treating it as a spread
	spread_.push_back(i);
      }
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
  // TEMPORARY: output Cossali number to file
  //plot(s(ind3)-airfoil.stagPt,K(ind3)./(Ks0*fs(ind3)),'b.');

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
      vN = 4.0*vNorm*( sqrt(K/Kb) - K/Kb );
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
    // Declare lots of parameters
    double x,y,u,v,r,temp,Time,numDrop;
    double K,Ks,Kb,vNormSq;
    double vN,vT;
    double vNorm,vTang,uNew,vNew;
    double theta,sCoord;
    double a,b,ms_m0,m0,ms,mStick,rStick;
    vector<double> XYq(2);
    vector<double> UVq(2);
    vector<double> XYa(2);
    vector<double> NxNy(2);
    vector<double> TxTy(2);
    int indSplash;
    double var = 0.2; double A0 = 0.09; double A1 = 0.51; double delK = 1500.0;
    double rm_rd,rm,mu;
    int dropRes = 1000;
    vector<double> dropsize(dropRes);
    vector<double> dropsizeCDF(dropRes);
    vector<double> rnew;
    double diffDropSize,mCHILD,mPARENT,dsamp;
    int numnew,RandIndex,indNN,indCellParent;
    int vRes = 1000;
    vector<double> vratio(vRes);
    vector<double> vratioCDF(vRes);
    vector<double> elev(2);
    vector<double> v1(2);
    vector<double> v2(2);
    double diffVratio,vCDFSamp,vrat,v2mag,e1,e2,elevation,foilAngle;
    int numChildSplash,numDropChild;
    // Initialize uniform random number generators
    default_random_engine generator;
    uniform_real_distribution<double> randCDF(0.05,0.95);
    uniform_real_distribution<double> randVelCDF(0,1);
    uniform_real_distribution<double> randE1(0,25.0*M_PI/180.0);
    uniform_real_distribution<double> randE2(M_PI-25.0*M_PI/180.0,M_PI);
    // Initialize spline interpolation
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, dropRes);
    gsl_spline *splineVel = gsl_spline_alloc(gsl_interp_linear, dropRes);
    // Initialize random seed
    srand ( time(NULL) );
    // Index over each splashing parcel
    for (int i=0; i<splash_.size(); i++) {
      // Get splashing parcel properties
      indSplash = impinge_[splash_[i]];
      x = state_.x_(indSplash);
      y = state_.y_(indSplash);
      u = state_.u_(indSplash);
      v = state_.v_(indSplash);
      r = state_.r_(indSplash);
      temp = state_.temp_(indSplash);
      Time = state_.time_(indSplash);
      numDrop = state_.numDrop_(indSplash);
      indCellParent = indCell_[indSplash];
      K = K_[splash_[i]];
      Ks = fs_[splash_[i]]*Ks0_;
      vNormSq = vNormSq_[splash_[i]];
      vTang = vTang_[splash_[i]];
      XYq[0] = x; XYq[1] = y;
      UVq[0] = u; UVq[1] = v;
      // Calculate s-coords of impinging parcel
      sCoord = airfoil.interpXYtoS(XYq);
      // Check to see whether we are using splashing at all, or not
      if (SplashFlag_ == true) {
	// Calculate impinging incidence angle
	airfoil.findPanel(XYq,XYa,NxNy,TxTy,indNN);
	theta = airfoil.calcIncidenceAngle(XYq,UVq,indNN);
	// Calculate impinging mass loss parameters
	a = 1.0-0.3*sin(theta);
	b = (1.0/8.0)*(1.0+3.0*cos(theta));
	// Calculate splashing ejection mass (ms) and sticking mass
	// (mStick)
	ms_m0 = max(a - pow(Ks/K,b),0.0);
	m0 = (4.0/3.0)*M_PI*pow(r,3);
	ms = ms_m0*m0;
	mStick = rhoL_*(m0-ms);
	// Update parent particle properties (mass/radius)
	rStick = pow(mStick/(rhoL_*(4.0/3.0)*M_PI),1.0/3.0);
	state_.r_(indSplash) = rStick;
      }
      else {
	// Simply treat particles marked as "splash" as spreading
	mStick = rhoL_*(4.0/3.0)*M_PI*pow(r,3);
	ms = 0;
      }
      // Check to see whether we are tracking child splash particles
      // Interpolate analytical expression for the CDF to get child droplet size
      if ((ms != 0) && (TrackSplashParticles_ == true)) {
	numChildSplash = 100;
	// Calculate dropsize CDF
	rm_rd = A0 + A1*exp(-K/delK);
	rm = rm_rd*r;
	mu = log(rm);
	diffDropSize = (r-0.05*rm)/(dropRes-1);
	for (int j=0; j<dropRes; j++) {
	  dropsize[j] = 0.05*rm + j*diffDropSize;
	  dropsizeCDF[j] = 0.5 + 0.5*erf((1/sqrt(2*var))*(log(dropsize[j])-mu));
	}
	// Create spline of CDF for interpolation
	gsl_spline_init(spline, dropsizeCDF.data(), dropsize.data(), dropRes);
	// Draw child particles until mass is conserved
	mCHILD = 0;
	mPARENT = ms*rhoL_;
	numnew = 0;
	while ((mCHILD < mPARENT) && (numnew < numChildSplash)) {
	  dsamp = randCDF(generator);
	  rnew.push_back(gsl_spline_eval(spline, dsamp, acc));
	  mCHILD = mCHILD + (4.0/3.0)*M_PI*pow(rnew[numnew],3)*rhoL_;
	  numnew++;
	}
	// Calculate child parcel number density
	numDropChild = round(mPARENT/mCHILD);
	if (numDropChild == 0) {
	  numDropChild = 1;
	}
	// Calculate post splashing droplet velocities, interpolate
	// analytical expression for the CDF to get magnitude of v2
	// for splash droplets, where v_new = v1 + v2
	diffVratio = 1.0/(vRes-1);
	for (int j=0; j<vRes; j++) {
	  vratio[j] = j*diffVratio;
	  vratioCDF[j] = 1 - exp(-13.7984*pow(vratio[j],2.5));
	}
	// Create spline of velocity CDF for interpolation
	gsl_spline_init(splineVel, vratioCDF.data(), vratio.data(), vRes);
	State stateChildren(numnew);
	for (int j=0; j<numnew; j++) {
	  // Interpolate velocity spline
	  vCDFSamp = randVelCDF(generator);
	  vrat = gsl_spline_eval(splineVel,vCDFSamp,acc);
	  v2mag = vrat*sqrt(vNormSq);
	  // Calculate elevation (relative to surface tangent) of splash droplet rebounds
	  e1 = randE1(generator);
	  e2 = randE2(generator);
	  elev[0] = e1; elev[1] = e2;
	  RandIndex = rand() % 2;
	  elevation = elev[RandIndex];
	  foilAngle = atan2(TxTy[1],TxTy[0]);
	  // Calculate v2
	  v2[0] = v2mag*cos(foilAngle+elevation);
	  v2[1] = v2mag*sin(foilAngle+elevation);
	  // Calculate v1
	  v1[0] = 0.8*vTang*cos(foilAngle);
	  v1[1] = 0.8*vTang*sin(foilAngle);
	  // Total splashed velocity = v1 + v2
	  stateChildren.u_(j) = v1[0] + v2[0];
	  stateChildren.v_(j) = v1[1] + v2[1];
	  // Set other child properties (inherited from parent)
	  stateChildren.x_(j) = x;
	  stateChildren.y_(j) = y;
	  stateChildren.r_(j) = rnew[j];
	  stateChildren.temp_(j) = temp;
	  stateChildren.time_(j) = Time;
	  stateChildren.numDrop_(j) = numDrop*numDropChild;
	}
	// Add child particles to the cloud
	this->addParticles(stateChildren,indCellParent);
      }
      // Add mass which has "stuck" to the airfoil
      airfoil.appendFilm(sCoord,numDrop*mStick);

    }

  }

}

void Cloud::spreadDynamics(Airfoil& airfoil) {
  // Function to compute spreading dynamics (pure stick)

  double x,y,r;
  vector<double> XYq(2);
  double numDrop,mSpread,sCoord,indSpread;
  // Index over each spreading parcel
  for (int i=0; i<spread_.size(); i++) {
    // Get splashing parcel properties
    indSpread = impinge_[spread_[i]];
    x = state_.x_(indSpread);
    y = state_.y_(indSpread);
    r = state_.r_(indSpread);
    numDrop = state_.numDrop_(indSpread);
    mSpread = (4.0/3.0)*M_PI*pow(r,3)*rhoL_;
    XYq[0] = x; XYq[1] = y;
    sCoord = airfoil.interpXYtoS(XYq);
    airfoil.appendFilm(sCoord,numDrop*mSpread);
    
  }

}

void Cloud::setIndAdv(vector<int>& indAdv) {
  
  indAdv = indAdv_; 
}

vector<int> Cloud::getIndAdv() {

  return indAdv_;
}

vector<int> Cloud::getIndSplash() {

  return splash_;
}

void Cloud::setState(State& state, PLOT3D& grid) {
  // Function to reset entire state of cloud

  // Set initial state of particles
  int size = state.size_;
  state_.size_ = size;
  state_.x_.resize(size);
  state_.y_.resize(size);
  state_.u_.resize(size);
  state_.v_.resize(size);
  state_.r_.resize(size);
  state_.temp_.resize(size);
  state_.time_.resize(size);
  state_.numDrop_.resize(size);
  for (int i=0; i<size; i++) {
    state_.x_(i) = state.x_(i);
    state_.y_(i) = state.y_(i); 
    state_.u_(i) = state.u_(i);
    state_.v_(i) = state.v_(i);
    state_.r_(i) = state.r_(i);
    state_.temp_(i) = state.temp_(i);
    state_.time_(i) = 0;
    state_.numDrop_(i) = 1;
  }

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

void Cloud::clearData() {
  // Function to clear any data associated with cloud

  // Clear other elements
  particles_ = 0;
  impingeTotal_.clear();
  indCell_.clear();
  indAdv_.clear();
  dt_.clear();
  impinge_.clear();
  bounce_.clear();
  spread_.clear();
  splash_.clear();
  K_.clear();
  fs_.clear();
  fb_.clear();
  vNormSq_.clear();
  vTang_.clear();

}

double Cloud::calcTotalMass() {
  // Function to calculate and return total mass

  double mass = 0;
  double R,numDrop;
  for (int i=0; i<particles_; i++) {
    R = state_.r_(i);
    numDrop = state_.numDrop_(i);
    mass += (4.0/3.0)*M_PI*rhoL_*pow(R,3)*numDrop;
  }

  return mass;
}
