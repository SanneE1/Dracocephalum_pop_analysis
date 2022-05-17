library(targets)
library(tarchetypes)

  lapply(
    as.list(
      grep(
        append(list.files("R/", pattern = "functions", full.names = T), 
               list.files("analysis/", full.names = T)), 
        pattern = ".R$", value = T, ignore.case = T)),
    source)
 

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "lme4", "patchwork", 
                            "rethinking", "ipmr", "mgcv",
                            "parallel"))

list(
  tar_target(
    raw_data_file_location,
    "data/Dracocephalum_03_21.xlsx",
    format = "file"
  ),
  tar_target(
    raw_data_long_format,
    format_raw_excel_data_to_long_format(raw_data_file_location)
  ),
  tar_target(
    explore_plot_raw_data,
    explore_plot(raw_data_long_format)
    ),
  tar_target(
    data_for_modeling,
    transform_raw_data_and_save(raw_data_long_format),
    format = "file"
  ),
  tar_render(
    state_dependency,
    "analysis/check_states_and_state_variables.Rmd",
    params = list(raw_data_location = data_for_modeling),
    output_dir = "results/exploratory plots/"
  ),
  tar_target(
    raw_climate_data,
    "data/climatic data Dobrichovice.xlsx",
    format = "file"
  ),
  tar_render(
    climate_stations,
    "analysis/climate_stations.Rmd",
    params = list(gps_locations = "data/Draco clima GPS.xlsx",
                  data_loggers_folder = "data/clima/",
                  climate_dobrichovice = "data/climatic data Dobrichovice.xlsx",
                  chelsa_data = "data/CHELSA_data.csv"),
    output_dir = "results/exploratory plots/"
  ),
  tar_target(
    transformed_climate,
    transform_and_save_climate_data(raw_climate_data),
    format = "file"
  ), 
  tar_render(
    climate_exploration,
    "analysis/climate_exploration.Rmd",
    params = list(climate_data_location = transformed_climate),
    output_dir = "results/exploratory plots/"
  ),
  tar_target(
   VR_FLM,
   VR_climate_FLM(data_for_modeling, transformed_climate, lag = 48)
  ),
  tar_target(
    state_independent_variables,
    estimate_non_state_dependent_variables(data_for_modeling)
  )
)
