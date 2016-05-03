#!/bin/bash

# *********************************
# SHELL SCRIPT DRIVER
# *********************************

# Set up multiple runs in serial
inFile="/home/adegenna/LagrangianIcingCode/DistributionCPP/404.inp"
BASE="/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/RUN404"
workDir="T_SIMUL_ROE"
motherDir="/home/adegenna/LagrangianIcingCode/DistributionCPP"

for i in {1..7}
do
    # I/O directories
    inDir=$BASE/$workDir$((i-1))
    outDir=$BASE/$workDir$i
    mkdir $outDir
    cp $BASE/horn.d $outDir

    # Run icing binary
    echo INPUT:  $inFile
    echo INDIR:  $inDir  
    echo OUTDIR: $outDir
    cp $motherDir/header1 $outDir
    cp $motherDir/header2 $outDir
    $motherDir/IcingDriver $inFile $inDir $outDir

    # Run GAIR/HYPERG
    cd $outDir
    /home/adegenna/Mesh2D/GAIR/gair
    /home/adegenna/Mesh2D/HYPERG/hyperg
    /home/adegenna/Flo103/flo103_sa < horn.d

    # Go back to motherDir
    cd $motherDir
done
