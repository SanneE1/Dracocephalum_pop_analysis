### Pipeline for the Dracocephalum population analysis

library(tidyverse)
library(RColorBrewer)
library(lme4)
library(patchwork)
library(rethinking)
library(ipmr)
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
# python you continueu on with data processing in R (more details in script)
# File produced: 
# "data/CHELSA_data.csv"

# Download .nc files from here first: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V1%2Fchelsa_cmip5_ts
# then run: 
# source("R/CHELSA_future_ts_extract_and_format.R")
# File produced: 
# "data/CHELSA_future_ts_formatted.csv"


rmarkdown::render('analysis/climate_stations.Rmd', 
                  output_dir = "results/exploratory plots/")
rmarkdown::render('analysis/climate_exploration.Rmd', 
                  output_dir = "results/exploratory plots/")
rmarkdown::render('analysis/check_states_and_state_variables.Rmd', 
                  output_dir = "results/exploratory plots/")

# Create Arima model that can be used to simulate climate for IPMs
source("analysis/Arima_climate_models.R")
# File produced:
# results/rds/ARIMA_clim_mods.rds


# ----------------------------------------------------------------------
# Vital rates
# ----------------------------------------------------------------------

# Functional linear models for state dependent variables
source("analysis/FLM_climate_models.R")
# Files produced:
# results/rds/VR_FLM.rds
# results/rds/VR_mod_infos.rds

# State independent variables
source("analysis/state_independent_variables.R")
# File produced:
# results/rds/state_independent_VR.rds

# ----------------------------------------------------------------------
# Population model
# ----------------------------------------------------------------------

# deterministic ipm (locality/year specific. no other covariates)
source("analysis/ipm_det_analysis.R")
# File produced:
# result/rds/deterministic_ipm.rds

# Long-term IPM with climate & environmental variables
# takes a long time. With submit_stoch_ipm.sh it can be run on the UFZ HPC
source("analysis/ipm_stoch_analysis.R")
# File produced:
# result/overview_lambda_env_levels.csv

# Calculate extinction probability from 2022-2100
# Takes a long time. with submit_extinct_prob.sh it can be run on the UFZ HPC
source("analysis/extinction_prop.R")
# File produced:
# results/rds/extinction_probability.rds

# Model possible management options -> planting adults or adding seeds
source("analysis/management_strategies.R")
# Files produced:
# results/rds/management_projections_transplants.rds
# results/rds/management_projections_seedaddition.rds



# ----------------------------------------------------------------------
# Process results
# ----------------------------------------------------------------------

# Plot population model results
source("analysis/Plot_results.R")


# Create document sumarizing results
rmarkdown::render('results/Summarise results.Rmd', 
                  output_format = "pdf_document",
                  output_dir = "results/")




