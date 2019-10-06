#!/bin/bash

# directories with all cell code
# NOTE: IF IN CELLS DIR, DO NOT CHANGE
cellsdir=~/cells
srcdir=$cellsdir/src
maindir=$cellsdir/main

# directory for all output for cell simulations
outputdir=~/project/cells

# directory for simulations specific to hopper flow
simtypedir=$outputdir/hopper

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NCELLS=$1
NT=$2
sizeDisp=$3
g=$4
w0=$5
w=$6
partition=$7
time=$8

let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=hopper2D_SP_N"$NCELLS"_NT"$NT"_sig"$sizeDisp"_g"$g"_w0"$w0"_w"$w"
runstr="$basestr"_si"$startSeed"_sf"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/hopper2D_SP.cpp
echo Running $numSeeds sims of $NCELLS SPM particles, size dispersion = $sizeDisp,
echo flowing through hopper with width $w0, orifice $w, with forcing strength g = $g

# run compiler
rm -f $binf
g++ --std=c++11 -I $srcdir $mainf $srcdir/*.cpp -o $binf 
echo compiling with : g++ --std=c++11 -I $srcdir $mainf $srcdir/*.cpp -o $binf 

# check compilation
if [[ ! -f $binf ]]
then
    echo -- binary file does not exist, compilation failed.
    exit 1
fi

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# loop over files
let fcount=0

# LOOP OVER FILES. 
for seed in `seq $startSeed $numSeedsPerRun $endSeed`; do
    # count files
    let fcount=$fcount+1

    # echo to console
    echo On base seed $seed

    # echo string of numSeedPerRun commands to task file
    runString="cd `pwd`"

    # loop over seeds to go into runString
    let ssMax=$numSeedsPerRun-1

    for ss in `seq 0 $ssMax`; do
        # get seed for actual run
        let runseed=$seed+ss

        # get file str
        filestr="$basestr"_seed"$seed"

        # make output directory
        specificdir=$simdatadir/"$filestr"
        mkdir -p $specificdir

        # create output files
        posf=$specificdir/$filestr.pos
        enf=$specificdir/$filestr.en

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NT $sizeDisp $g $w0 $w $runseed $posf"
    done

    # finish off run string
    runString="$runString ;"

    # echo to task file
    echo "$runString" >> $taskf
done

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# get number of jobs to submit to each array
let arraynum=$fcount
echo -- total number of array runs = $arraynum

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$arraynum >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. NT
# 3. size dispersion of particles
# 4. forcing strength g
# 5. hopper reservoir width w0
# 6. hopper orifice width w
# 7. partition
# 8. time
# 9. num seeds per run (for each entry in array)
# 10. number of runs (number of array entries, i.e. arraynum)
# 11. start seed (end seed determined by number of runs)




