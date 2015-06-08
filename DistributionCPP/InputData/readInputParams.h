#ifndef __READINPUTPARAMS_H__
#define __READINPUTPARAMS_H__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "../Grid/FluidScalars.h"
#include "../Cloud/ParcelScalars.h"

void readInputParams(FluidScalars& PROPS, ParcelScalars& PARCEL, const char *inFileName);


#endif
