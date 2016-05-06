#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -M adegenna@princeton.edu

export PATH=$PBS_O_PATH

# *********************************
# SHELL SCRIPT DRIVER (INDIVIDUAL RUN)
# *********************************

# *********************************
# INITIALIZATION
# *********************************

# Files/directories
masterDir=$var1
BASE=$var2
num=$var3
inFile=$BASE/INP.inp

# Make initial T0 directory
step="T"
mkdir "$BASE/$step"0
cp "$BASE/fort.30" "$BASE/$step"0

# *********************************
# MESH-->FLOW-->ICE
# *********************************

for i in {1..5} # Indexes over timesteps of ice accretion
do
    # I/O directories
    inDir=$BASE/$step$((i-1))
    outDir=$BASE/$step$i
    mkdir $outDir
    echo "$BASE/FLO.d"
    cp "$BASE/FLO.d" $inDir
    cp "$BASE/INP.inp" $inDir
    
    # Run GAIR/HYPERG/FLO103
    cd $inDir
    /home/adegenna/Mesh2D/GAIR/gair
    /home/adegenna/Mesh2D/HYPERG/hyperg
    /home/adegenna/Flo103/flo103_sa < FLO.d

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

# *********************************
# OUTPUT
# *********************************

# Go to last timestep directory
cd $BASE/$step/5

# Parse out results from fort.66 (CL,CD)
CL=$(awk 'FNR==1 { print $2 }' fort.66)
CD=$(awk 'FNR==1 { print $3 }' fort.66)

# Write to results file
echo -e "$CL\tCL\t$CD\tCD" >> $BASE/results.out.$num
cp $BASE/results.out.$num $BASE/results.out
