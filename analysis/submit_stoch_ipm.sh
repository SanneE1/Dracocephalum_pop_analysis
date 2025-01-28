#!/bin/bash

#SBATCH -D /datos_severs/Dracocephalum_pop_analysis/

#SBATCH -o /datos_severs/Dracocephalum_pop_analysis/results/%x-%A_%a.log

#Specify job name
#SBATCH -J Draco_ipm

# N of array jobs
#SBATCH --array=1-1530

#Resources
# max running time

# memory per core (hard limit)
#SBATCH --mem-per-cpu=2G

# create output direcotry per job
OUTPUT_DIR="/datos_severs/Dracocephalum_pop_analysis/results/$SLURM_JOB_NAME-$SLURM_ARRAY_JOB_ID"
mkdir -p $OUTPUT_DIR

# Load modules
# module load foss/2019b R/4.0.0-2

setwd="/datos_severs/Dracocephalum_pop_analysis/"

Rscript --vanilla "$HOME"/Dracocephalum_pop_analysis/analysis/ipm_stoch_analysis.R \
"$OUTPUT_DIR" \
"$setwd"
