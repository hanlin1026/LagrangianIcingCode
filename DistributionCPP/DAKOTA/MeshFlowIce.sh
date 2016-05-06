#!/bin/bash

# *********************************
# SHELL SCRIPT DRIVER
# *********************************

# $1 is params.in.(fn_eval_num) FROM Dakota
# $2 is results.out.(fn_eval_num) returned to Dakota

# *********************************
# INPUT PARSING
# *********************************

# Get run number
num=${1//[^0-9]/}

# Master directory
masterDir="/home/adegenna/LagrangianIcingCode/DistributionCPP"

# Parse params.in for input icing parameters
NPARAMS=$(awk 'FNR==1 { print $1 }' $1)
TINF=$(awk 'FNR==2 {print $1 }' $1)
LWC=$(awk 'FNR==3 {print $1 }' $1)
printf -v TINF "%.*f" 2 "$TINF"
printf -v LWC "%.*E" 3 "$LWC"

# *********************************
# MESH-->FLOW-->ICE
# *********************************

# Initialize T0 directory with appropriate files
BASE=$(pwd)

 # Fill in blanks for input files (FLO.d, INP.inp)
find $BASE/"FLO.d" -type f -print0 | xargs -0 sed -i 's/tinf/'$TINF'/g'
find $BASE/"INP.inp" -type f -print0 | xargs -0 sed -i 's/TINF/'$TINF'/g'
find $BASE/"INP.inp" -type f -print0 | xargs -0 sed -i 's/LWC/'$LWC'/g'

# Qsub mesh/flow/icing script
qsub -v par_name=par_value\[,var1=$masterDir,var2=$BASE,var3=$num,\] $BASE/DRIVE.sh
