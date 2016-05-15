#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -m abe

# *********************************
# SHELL SCRIPT DRIVER (INDIVIDUAL RUN)
# *********************************

# Files/directories
masterDir=$var1
BASE=$var2
steps=$var3
workDir="T"
inFile=$BASE/*.inp

# Create initial directory/solution (use initial fort.30 file)
inDir="$BASE/T0"
outDir=$BASE/T1
mkdir $inDir
cp $BASE/"fort.30" $inDir

for (( i=8; i<=$steps; i++ ))
do
    # I/O directories
    inDir=$BASE/$workDir$((i-1))
    outDir=$BASE/$workDir$i
    mkdir $outDir
    cp $inFile $inDir
    cp $BASE/"FLO.d" $inDir

    # Run GAIR/HYPERG
    cd $inDir
    /home/adegenna/Mesh2D/GAIR/gair
    /home/adegenna/Mesh2D/HYPERG/hyperg
    /home/adegenna/Flo103/flo103_sa < "FLO.d" > "FLO.out"

    # Go back to masterDir
    cd $masterDir

    # Run icing binary
    echo INPUT:  $inFile
    echo INDIR:  $inDir  
    echo OUTDIR: $outDir
    cp $masterDir/header1 $outDir
    cp $masterDir/header2 $outDir
    $masterDir/CATFISH $inFile $inDir $outDir > $inDir/"CATFISH.out"

done
