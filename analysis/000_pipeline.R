### Pipeline for the Dracocephalum population analysis

library(tidyverse)
library(RColorBrewer)
library(lme4)
library(patchwork)
library(rethinking)
library(mgcv)
library(parallel)
library(forecast)


# ----------------------------------------------------------------------
# Demographic data
# ----------------------------------------------------------------------

# Format raw data into long format + add vitalrates
source("analysis/data_formatting.R")
# Files produced: 
# data/Dracocephalum_long_format.csv 
# data/Dracocephalum_with_vital_rates.csv

# ----------------------------------------------------------------------
# Climate data
# ----------------------------------------------------------------------

# Run "R/CHELSA_data_get_using_Python.R" manually. You'll need to switch to python
# and copy paste the lines produced in the 1st half of the script. After downloading through 
# python you continue on with data processing in R (more details in script)
# File produced: 
# "data/CHELSA_data.csv"

# Download .nc files from here first: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V1%2Fchelsa_cmip5_ts
# then run: 
# source("R/CHELSA_future_ts_extract_and_format.R")
# File produced: 
# "data/CHELSA_future_ts_formatted.csv"


# ----------------------------------------------------------------------
# Vital rates
# ----------------------------------------------------------------------

# Size dependent vital rates
source("analysis/vital_rates_models_elastnet.R")
# Files produced:
# results/rds/seasons_surv.rds
# results/rds/seasons_growth.rds
# results/rds/seasons_growth_sd.rds
# results/rds/seasons_flowp.rds
# results/rds/seasons_seedp.rds
# results/rds/seasons_seedn.rds
# results/rds/seasons_seedn_sd.rds


# State independent variables
source("analysis/state_independent_variables.R")
# File produced:
# results/rds/state_independent_VR.rds

# ----------------------------------------------------------------------
# Population dynamics model
# ----------------------------------------------------------------------

# Short-term IBM Calculate extinction probability from 2022-2100
# Takes a long time. with submit_extinct_prob.sh it can be run on the UFZ HPC
source("analysis/extinction_prop.R")
# File produced:
# results/rds/extinction_probability.rds

# ----------------------------------------------------------------------
# IBM diagnostics
# ----------------------------------------------------------------------

# Run model during observational data period for in sample model validation
# (out sample plotted in Plot_results.R)
source("analysis/IBM_validation.R")

# Run same future scenarios for small selection of shading values to evaluate 
# effect of different starting values
source("analysis/IBM_seedbank_sensitivity.R")

# ----------------------------------------------------------------------
# Process results
# ----------------------------------------------------------------------

# Plot population model results
source("analysis/Plot_results.R")





