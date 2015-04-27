#include <eigen3/Eigen/Dense>
#include "State.h"

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
}
