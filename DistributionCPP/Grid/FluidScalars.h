#ifndef __FLUIDSCALARS_H__
#define __FLUIDSCALARS_H__

struct FluidScalars {
  // Scalars needed for the Grid class
  double pinf_;
  double R_;
  double Tinf_;
  double rhoinf_;
  double Ubar_;
  double rhol_;
  // Scalars needed for constructor of Grid class
  int nx_;
  int ny_;
  float mach_;
  float alpha_;
  float reynolds_;
  float time_;
  // Track splashed particles indicator
  int calcImpingementLimits_;

};

#endif
