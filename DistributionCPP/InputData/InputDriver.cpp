#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <random>
#include <iterator>
#include "../findAll.h"
#include "../Grid/FluidScalars.h"
#include "../Cloud/ParcelScalars.h"
#include "readInputParams.h"

// Test driver program for reading input 

int main(int argc, const char *argv[]) {
  const char *inFileName = "Input.dat";
  FluidScalars PROPS;
  ParcelScalars PARCEL;
  readInputParams(PROPS, PARCEL, inFileName);
  printf("pinf = %f\n",PROPS.pinf_);
  printf("R = %f\n",PROPS.R_);
  printf("Tinf = %f\n",PROPS.Tinf_);
  printf("rhol = %f\n",PROPS.rhol_);
  printf("particles = %d\n",PARCEL.particles_);
  printf("Rmean = %f\n",PARCEL.Rmean_);
  printf("Tmean = %f\n",PARCEL.Tmean_);
}
