#!/bin/bash

#SBATCH -D /work/evers/

#SBATCH -o /work/%u/%x-%A.log

#Specify job name
#SBATCH -J Draco_ipm

#Resources
# max running time
#SBATCH -t 1:30:00

# memory per core (hard limit)
#SBATCH --mem-per-cpu=2G

# create output direcotry per job
OUTPUT_DIR="/work/$USER/$SLURM_JOB_NAME-$SLURM_ARRAY_JOB_ID"
mkdir -p $OUTPUT_DIR

# Load modules
module load foss/2019b R/4.0.0-2

setwd="/gpfs0/home/evers/Dracocephalum_pop_analysis/"

Rscript --vanilla "$HOME"/Dracocephalum_pop_analysis/analysis/ipm_stoch_analysis.R \
"$OUTPUT_DIR" \
"$setwd"
