#!/bin/sh
#SBATCH --job-name=vq
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nd2560@columbia.edu
#SBATCH --account=stats
# SBATCH --qos=statistics
#SBATCH -c 1                     # The number of cpu cores to use.
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=9:00:00
#SBATCH --output=Out/main_run_%A-%a.out
#SBATCH --array=1-1000
pwd; hostname; date

RUN="${SLURM_ARRAY_TASK_ID}"

echo "This is task ${RUN}"

module load R/4.3.1
R

# Out dir
# Rscript "--args seed=${RUN}" Tab1.R # > Out/test_%j-%a.out
R  --vanilla --no-save --args seed=${RUN} < Testingtwosamexec.R > /moto/home/nd2560/codeOutput/main_run${RUN}.log

# in R script

date
                                                                                                             
