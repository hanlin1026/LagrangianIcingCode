#ifndef __FLUIDSCALARS_H__
#define __FLUIDSCALARS_H__

struct FluidScalars {
  // Filenames for grid/flow solution
  std::string gridfile_;
  std::string solnfile_;
  std::string heatfile_;
  std::string betafile_;
  std::string outfile_;
  // Scalars needed for the Grid class
  double pinf_;
  double R_;
  double Tinf_;
  double rhoinf_;
  double Ubar_;
  double rhol_;
  double chord_;
  // Scalars needed for constructor of Grid class
  int nx_;
  int ny_;
  float mach_;
  float alpha_;
  float reynolds_;
  float time_;
  // Track splashed particles indicator
  int calcImpingementLimits_;

  // Scalars needed for thermodynamics module
  double Td_;
  double NPts_;
  double Uinf_;
  double LWC_;

};

#endif
