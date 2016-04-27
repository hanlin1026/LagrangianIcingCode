#!/bin/bash

# *********************************
# SHELL SCRIPT DRIVER
# *********************************

# I/O files/directories
inFile="/home/adegenna/LagrangianIcingCode/DistributionCPP/308.inp"
inDir="/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/RUN308/T3area2"
outDir="/home/adegenna/LagrangianIcingCode/DistributionCPP/Grid/RUN308/T4area2"

# Run icing binary
./IcingDriver $inFile $inDir $outDir
echo INPUT:  $inFile
echo INDIR:  $inDir  
echo OUTDIR: $outDir
mv XY_NEW.out $outDir

# Run GAIR/HYPERG
~/Mesh2D/GAIR/gair
~/Mesh2D/HYPERG/hyperg
~/Flo103/flo103_sa < horn.d

# Move files from working directory to output directory
mv fort.30 $outDir
mv fort.31 $outDir
mv fort.32 $outDir
mv fort.77 $outDir
mv fort.71 $outDir
mv fort.66 $outDir
mv MESH.D $outDir
mv MESH.P3D $outDir
mv cnv.plot $outDir
mv cp.plot $outDir
mv xyz.bin $outDir
mv heatflux $outDir
mv q103.0.40E+01.bin $outDir
