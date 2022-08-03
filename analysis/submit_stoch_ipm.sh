#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J Draco_ipm

#Resources
# max running time
#SBATCH -t 48:00:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=25

# Load modules
module load foss/2019b R/4.0.0-2

setwd="/gpfs0/home/evers/Dracocephalum_pop_analysis/"

export MC_CORES=${SLURM_CPUS_PER_TASK:-1}

Rscript --vanilla "$HOME"/Dracocephalum_pop_analysis/analysis/ipm_stoch_analysis.R \
"$setwd"
