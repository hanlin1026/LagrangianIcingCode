#include <eigen3/Eigen/Dense>
#include "State.h"
#include <Cloud/ParcelScalars.h>

State::State() {

}
State::State(int size) {
  // Constructor
  
  size_ = size;
  x_.resize(size);
  y_.resize(size);
  u_.resize(size);
  v_.resize(size);
  r_.resize(size);
  temp_.resize(size);
  time_.resize(size);
  numDrop_.resize(size);
}

State::State(const char* distributionType, ParcelScalars& scalars, PLOT3D& p3d) {
  // Constructor

  size_ = scalars.particles_;
  x_.resize(size_);
  y_.resize(size_);
  u_.resize(size_);
  v_.resize(size_);
  r_.resize(size_);
  temp_.resize(size_);
  time_.resize(size_);
  numDrop_.resize(size_);
  double Xnn, Ynn;
  int indnn;
  if (strcmp(distributionType,"MonoDispersed") == 0) {
    if (scalars.Xmin_ != scalars.Xmax_) {
      // Initialize all particles as having the same size
      // Select particle locations randomly and set state
      default_random_engine generator;
      uniform_real_distribution<double> distX(scalars.Xmin_,scalars.Xmax_);
      uniform_real_distribution<double> distY(scalars.Ymin_,scalars.Ymax_);
      for (int i=0; i<size_; i++) {
	x_(i) = distX(generator);
	y_(i) = distY(generator);
	p3d.pointSearch(x_(i),y_(i),Xnn,Ynn,indnn);
	u_(i) = p3d.getUCENT(indnn);
	v_(i) = p3d.getVCENT(indnn);
	r_(i) = scalars.Rmean_;
	temp_(i) = scalars.Tmean_;
	time_(i) = 0;
	numDrop_(i) = 1;
      }
    }
    else {
      // Initialize particles along line of constant X
      double dY = (scalars.Ymax_-scalars.Ymin_)/(size_-1);
      for (int i=0; i<size_; i++) {
	x_(i) = scalars.Xmax_;
	y_(i) = scalars.Ymin_ + i*dY;
	p3d.pointSearch(x_(i),y_(i),Xnn,Ynn,indnn);
	u_(i) = p3d.getUCENT(indnn);
	v_(i) = p3d.getVCENT(indnn);
	r_(i) = scalars.Rmean_;
	temp_(i) = scalars.Tmean_;
	time_(i) = 0;
	numDrop_(i) = 1;
      }
    }

  }

  else if (strcmp(distributionType,"Parcels") == 0) {
    
  }

  else {
    printf("ERROR: Please select either either 'MonoDispersed' or 'Parcels' in the State constructor");
  }

}

State::~State() {

}

void State::appendState(State& addition) {
  // Function to append elements to state
  
  int deltaSize = addition.size_;
  // Resize state element vectors
  x_.conservativeResize(size_ + deltaSize);
  y_.conservativeResize(size_ + deltaSize);
  u_.conservativeResize(size_ + deltaSize);
  v_.conservativeResize(size_ + deltaSize);
  r_.conservativeResize(size_ + deltaSize);
  temp_.conservativeResize(size_ + deltaSize);
  time_.conservativeResize(size_ + deltaSize);
  numDrop_.conservativeResize(size_ + deltaSize);
  // Append elements to state
  for (int i=0; i<deltaSize; i++) {
    x_(size_+i) = addition.x_(i);
    y_(size_+i) = addition.y_(i);
    u_(size_+i) = addition.u_(i);
    v_(size_+i) = addition.v_(i);
    r_(size_+i) = addition.r_(i);
    temp_(size_+i) = addition.temp_(i);
    time_(size_+i) = addition.time_(i);
    numDrop_(size_+i) = addition.numDrop_(i);
  }
  // Update state size
  size_ += deltaSize;
}
