#!/bin/bash

# *********************************
# SHELL SCRIPT DRIVER
# *********************************

# Tab-delimited file describing multiple job submissions
masterDir="/home/adegenna/CATFISH/Validations"
masterList=$masterDir/"RUNS.list"
masterInp=$masterDir/"INP.inp"
masterFLO=$masterDir/"FLO.d"
masterDrive=$masterDir/"DRIVE.sh"
runsDir=$masterDir/"RUNS"
fort=$masterDir/"fort.30"
header1=$masterDir/"header1"
header2=$masterDir/"header2"

for k in {1..26} # Indexes over individual runs
do
    # Read in run number
    RUN=$(awk 'FNR == '$((k+1))' { print $1 }' $masterList)
    
    # Make base directory for current run, copy template icing input files and template driver into it
    BASE=$runsDir/$RUN
    mkdir $BASE
    inFile=$BASE/$RUN".inp"
    cp $masterInp $inFile
    cp $masterDrive $BASE
    cp $masterFLO $BASE
    cp $fort $BASE
    cp $header1 $BASE
    cp $header2 $BASE

    # Get flow/icing input parameters for current run
    PINF=$(awk 'FNR == '$((k+1))' { print $2 }' $masterList)
    TINF=$(awk 'FNR == '$((k+1))' { print $3 }' $masterList)
    RDROP=$(awk 'FNR == '$((k+1))' { print $4 }' $masterList)
    UINF=$(awk 'FNR == '$((k+1))' { print $5 }' $masterList)
    LWC=$(awk 'FNR == '$((k+1))' { print $6 }' $masterList)
    CHORD=$(awk 'FNR == '$((k+1))' { print $7 }' $masterList)
    MACH=$(awk 'FNR == '$((k+1))' { print $8 }' $masterList)
    RE=$(awk 'FNR == '$((k+1))' { print $9 }' $masterList)
    DT=$(awk 'FNR == '$((k+1))' { print $10 }' $masterList)
    ROUGH=$(awk 'FNR == '$((k+1))' { print $11 }' $masterList)
    TIME=$(awk 'FNR == '$((k+1))' { print $12 }' $masterList)

    # Fill in blanks for template icing input file
    find $inFile -type f -print0 | xargs -0 sed -i 's/TINF/'$TINF'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/PINF/'$PINF'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/RDROP/'$RDROP'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/UINF/'$UINF'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/LWC/'$LWC'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/CHORD/'$CHORD'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/MACH/'$MACH'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/RE/'$RE'/g'
    find $inFile -type f -print0 | xargs -0 sed -i 's/DT/'$DT'/g'
    
    # Fill in blanks for GAIR/HYPERG/FLO103 input file
    find $BASE/"FLO.d" -type f -print0 | xargs -0 sed -i 's/mach/'$MACH'/g'
    find $BASE/"FLO.d" -type f -print0 | xargs -0 sed -i 's/re/'$RE'/g'
    find $BASE/"FLO.d" -type f -print0 | xargs -0 sed -i 's/tinf/'$TINF'/g'
    find $BASE/"FLO.d" -type f -print0 | xargs -0 sed -i 's/rough/'$ROUGH'/g'

    # Qsub job driver scripts
    qsub -v par_name=par_value\[,var1=$masterDir,var2=$BASE,var3=$TIME,\] $BASE/DRIVE.sh

done

