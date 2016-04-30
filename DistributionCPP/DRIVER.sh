#!/bin/bash

# *********************************
# SHELL SCRIPT DRIVER
# *********************************

# Set up multiple runs in serial
inFile="/home/adegenna/LagrangianIcingCode/DistributionCPP/405.inp"
BASE="/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/RUN405"
workDir="T_ROE"
motherDir=$(pwd)

for i in {6..7}
do
    # I/O directories
    inDir=$BASE/$workDir$((i-1))
    outDir=$BASE/$workDir$i
    mkdir $outDir
    cp $BASE/horn.d $outDir

    # Run icing binary
    ./IcingDriver $inFile $inDir $outDir
    echo INPUT:  $inFile
    echo INDIR:  $inDir  
    echo OUTDIR: $outDir

    # Run GAIR/HYPERG
    cd $outDir
    ~/Mesh2D/GAIR/gair
    ~/Mesh2D/HYPERG/hyperg
    ~/Flo103/flo103_sa < horn.d

    # Go back to motherDir
    cd $motherDir
done
