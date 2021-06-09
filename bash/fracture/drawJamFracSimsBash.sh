#!/bin/bash
# directories with code
#should print movies to the directory matlab_out/
#example call: bash bash/seq/drawJamFracSimsBash.sh 24 24 1.08 1.0 0.0 0.5 0.05 1 1 0 1 1 pi_ohern,day,scavenge 0-12:00:00

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
kl=$4
kb=$5
att=$6
B=$7
startSeed=$8
numRuns=$9
makePlots="${10}"
makeAMovie="${11}"
onCluster="${12}"
partition="${13}"
time="${14}"

numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=matlabfracture_N"$NCELLS"_NV"$NV"_calA"$calA0"_kl"$kl"_kb"$kb"_att"$att"_B"$B"
runstr=matlab"$basestr"_startseed"$startSeed"_endseed"$endSeed"

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
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH --mail-type=END,FAIL >> $slurmf
echo \#SBATCH --mail-user=andrewtondata@gmail.com >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo module load MATLAB/2020b >> $slurmf
echo matlab -batch "drawJamFracSims($NCELLS, $NV, $calA0, $kl, $kb, $att, $B, $startSeed, $numRuns, $makePlots, $makeAMovie, $onCluster)" >> $slurmf
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
# 4. kl
# 5. kb
# 6. attraction strength
# 7. damping coeff
# 8. start seed
# 9. number of runs
# 10. boolean plots or no
# 11. boolean movie or no
# 12. boolean cluster or no
# 13. partition to run on
# 14. how much time to allocate
