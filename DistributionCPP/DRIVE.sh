#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -m abe

# *********************************
# SHELL SCRIPT DRIVER (INDIVIDUAL RUN)
# *********************************

# Files/directories
masterDir=$var1
BASE=$var2
workDir="T"
inFile=$BASE/*.inp

# Create initial directory/solution
inDir="$BASE/T0"
outDir=$BASE/T1
mkdir $inDir
cp $BASE/"FLO.d" $inDir
cp $BASE/"fort.30" $inDir

for i in {1..7}
do
    # I/O directories
    inDir=$BASE/$workDir$((i-1))
    outDir=$BASE/$workDir$i
    mkdir $outDir
    cp $inFile $inDir

    # Run GAIR/HYPERG
    cd $inDir
    /home/adegenna/Mesh2D/GAIR/gair
    /home/adegenna/Mesh2D/HYPERG/hyperg
    /home/adegenna/Flo103/flo103_sa < "FLO.d"

    # Go back to masterDir
    cd $masterDir

    # Run icing binary
    echo INPUT:  $inFile
    echo INDIR:  $inDir  
    echo OUTDIR: $outDir
    cp $masterDir/header1 $outDir
    cp $masterDir/header2 $outDir
    $masterDir/IcingDriver $inFile $inDir $outDir

done
