#ifndef CLOUD_H_
#define CLOUD_H_

#include <stdio.h>
#include "State.h"
#include <QuadTree/Bucket.h>
#include <Airfoil/Airfoil.h>
#include <Grid/PLOT3D.h>
#include <findAll.h>
#include <eigen3/Eigen/Dense>
#include "ParcelScalars.h"

class Cloud {
 public:
  Cloud(State& state, PLOT3D& grid, double rhol, ParcelScalars& PARCEL);
  Cloud(State& state, PLOT3D& grid, double rhol);
  ~Cloud();
  void addParticles(State& state, int indCellParent);
  // Methods for SLD dynamics
  void calcDtandImpinge(Airfoil& airfoil, PLOT3D& grid);
  void transportSLD(PLOT3D& grid);
  void computeImpingementRegimes(Airfoil& airfoil);
  void bounceDynamics(Airfoil& airfoil);
  void splashDynamics(Airfoil& airfoil);
  void spreadDynamics(Airfoil& airfoil);
  // Set/get methods
  State getState();
  void setState(State& state, PLOT3D& grid);
  void setIndAdv(std::vector<int>& indAdv);
  std::vector<int> getIndAdv();
  std::vector<int> getIMPINGE();
  std::vector<int> getIMPINGETOTAL();
  std::vector<int> getINDCELL();
  std::vector<int> getIndSplash();
  // Calculate total mass method
  double calcTotalMass();
  // Clear data
  void clearData();

 private:
  State state_;
  double rhoL_;
  int particles_;
  double sigma_;
  std::vector<int> impingeTotal_;
  std::vector<int> indCell_;
  std::vector<int> indAdv_;
  std::vector<double> dt_;
  std::vector<int> impinge_;
  std::vector<int> bounce_;
  std::vector<int> spread_;
  std::vector<int> splash_;
  std::vector<double> K_;
  std::vector<double> fs_;
  std::vector<double> fb_;
  std::vector<double> vNormSq_;
  std::vector<double> vTang_;
  double Ks0_,Kb0_;
  void findInSimulation();
  void computeNewCellLocations(PLOT3D& grid);
  bool TrackSplashParticles_;
  bool SplashFlag_;

};

#endif
