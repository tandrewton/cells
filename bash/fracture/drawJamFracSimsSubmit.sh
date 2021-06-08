#!/bin/bash
#SBATCH -J matlabDraw
#SBATCH --output=out/matlabDraw.out
#SBATCH --error=out/matlabDraw.err
#SBATCH --chdir=/home/at965/cells
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=andrewtondata@gmail.com
#SBATCH -c 4
#SBATCH -t 6:00:00
#SBATCH -p pi_ohern,day,scavenge

module load MATLAB/2019b
#matlab -nodisplay -nosplash -r YourFunction < /dev/null
matlab -nodisplay -nosplash -r drawJamFracSims < /dev/null
