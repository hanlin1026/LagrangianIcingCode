#ifndef VECTOROPERATIONS_H_
#define VECTOROPERATIONS_H_

// Simple auxiliary functions and operator overloading for use with std::vector

#include <stdlib.h>

using namespace std;

inline vector<double> operator*(double scal, vector<double> vec) {
  for (int i=0; i<vec.size(); i++) {
    vec[i] *= scal;
  }
  return vec;
}

inline vector<double> operator*(vector<double> vec1, std::vector<double> vec2) {
  for (int i=0; i<vec1.size(); i++) {
    vec1[i] *= vec2[i];
  }
  return vec1;
}

inline vector<double> operator-(vector<double> vec1, vector<double> vec2) {
  vector<double> vec = vec1;
  for (int i=0; i<vec.size(); i++) {
    vec[i] -= vec2[i];
  }
  return vec;
}

inline vector<double> operator-(vector<double>& vec, double scal) {
  for (int i=0; i<vec.size(); i++) {
    vec[i] -= scal;
  }
  return vec;
}

inline double max(vector<double>& vec) {
  double maximum = 0.0;
  for (int i=0; i<vec.size(); i++) {
    if (vec[i]>maximum) {
      maximum = vec[i];
    }
  }

  return maximum;
}

inline double min(vector<double> vec, int& ind) {
  double minimum = vec[0];
  ind = 0;
  for (int i=1; i<vec.size(); i++) {
    if (vec[i]<minimum) {
      minimum = vec[i];
      ind = i;
    }
  }

  return minimum;
}

inline vector<double> abs(vector<double> vec) {
  vector<double> absVec = vec;
  for (int i=0; i<vec.size(); i++) {
    if (vec[i]<0) {
      absVec[i] = -1.0*vec[i];
    }
  }
  
  return absVec;
}

inline vector<double> find(vector<double>& vec, double val, bool& flag) {
  // Function to find indices where vec equals val

  vector<double> indices;
  flag = false;
  for (int i=0; i<vec.size(); i++) {
    if (vec[i] == val) {
      indices.push_back(i);
      flag = true;
    }
  }

  return indices;
}

inline vector<double> find(vector<double>& vec, vector<double>& vec2, bool& flag) {
  // Function to find indices where vec[i] equals vec2[i]

  vector<double> indices;
  flag = false;
  for (int i=0; i<vec.size(); i++) {
    if (vec[i] == vec2[i]) {
      indices.push_back(i);
      flag = true;
    }
  }

  return indices;
}

inline vector<double> flipud(vector<double>& vec) {
  // Function to flip a vector so that vec[0] becomes vec[end] and vice versa

  vector<double> tmp(vec.size());
  for (int i=0; i<vec.size(); i++) {
    tmp[i] = vec[vec.size()-1-i];
  }
  
  return tmp;
}

#endif
