#!/bin/bash
# directories with code

#example call: bash bash/seq/seqJamFractureSubmit.sh 24 24 1.08 0.2 1e-7 1.0 0 0.0 0.05 4e6 pi_ohern,day,scavenge 0-12:00:00 1 1

cellsdir=~/cells
srcdir=$cellsdir/src
maindir=$cellsdir/sequential

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/at965/cells

# directory for simulations specific to fracture
simtypedir=$outputdir/fracture

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NCELLS=$1
NV=$2
calA0=$3
phiMin=$4
Ptol=$5
kl=$6
kb=$7
att=$8
B=$9
NT="${10}"
partition="${11}"
time="${12}"
numRuns="${13}"
startSeed="${14}"

numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=fracture_N"$NCELLS"_NV"$NV"_calA"$calA0"_kl"$kl"_kb"$kb"_att"$att"_B"$B"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/fracture/jamFracture.cpp
echo Running $numSeeds fracture sims of $NCELLS cells with $NV verts, bidisperse , calA0 = $calA0 and attraction parameter att = $att, damping coefficient = $B with $NT timesteps

# run compiler
rm -f $binf
g++ --std=c++11 -O3 $mainf -o $binf
echo compiling with : g++ --std=c++11 -O3 $mainf -o $binf

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

        # create output files
        posf=$simdatadir/$filestr.pos
        shapef=$simdatadir/$filestr.shape
        energyf=$simdatadir/$filestr.energy

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $calA0 $phiMin $Ptol $kl $kb $att $B $runseed $NT $posf $shapef $energyf"
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
echo \#SBATCH --mail-type=END,FAIL >> $slurmf
echo \#SBATCH --mail-user=andrewtondata@gmail.com >> $slurmf
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
# 2. NV
# 3. calA0
# 4. phiMin
# 5. Ptol
# 6. kl
# 7. kb
# 8. att
# 9. number of timesteps
# 10. partition
# 11. time
# 12. number of runs (number of array entries, i.e. arraynum)
# 13. start seed (end seed determined by number of runs)
