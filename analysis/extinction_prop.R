### script to run extinction probability analysis
library(tidyverse)
library(glmnet)
library(lme4)

source("R/functions_ibm.R")
# source("R/functions_ipmr.R")

# Projections will be from 2022 to 2100
n_it = 79

## -------------------------------------------------------------------------------------------
## Format CHELSA's time series for IBM
## -------------------------------------------------------------------------------------------

fut_clim <- read.csv("data/CHELSA_future_ts_formatted.csv") %>% 
  filter(complete.cases(.))

clim_spring <- fut_clim %>% 
  filter(month >= 3 & month <= 5) %>%
  mutate(year_t0 = year - 1, .keep = "unused") %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_spring = mean(pet_scaled, na.rm = T),
            pr_spring = mean(pr_scaled, na.rm = T),
            tas_spring = mean(tas_scaled, na.rm = T)) %>%
  ungroup()
clim_summer <- fut_clim %>% 
  filter(month >= 6 & month <= 8) %>%
  rename(year_t0 = year) %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_summer = mean(pet_scaled, na.rm = T),
            pr_summer = mean(pr_scaled, na.rm = T),
            tas_summer = mean(tas_scaled, na.rm = T)) %>%
  ungroup()
clim_dormant <- fut_clim %>% 
  filter(month <= 2 | month >= 9) %>%
  mutate(year_t0 = ifelse(month < 3, year - 1, year)) %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_dormant = mean(pet_scaled, na.rm = T),
            pr_dormant = mean(pr_scaled, na.rm = T),
            tas_dormant = mean(tas_scaled, na.rm = T)) %>%
  ungroup() 


fut_clim <- left_join(clim_spring, clim_summer)
fut_clim <- left_join(fut_clim, clim_dormant)


## -------------------------------------------------------------------------------------------
## Other variables for IBM
## -------------------------------------------------------------------------------------------

VR_mods <- list(surv = readRDS("results/rds/seasons_surv.rds"),
                growth = readRDS("results/rds/seasons_growth.rds"),
                flower_p = readRDS("results/rds/seasons_flowp.rds"),
                seedp = readRDS("results/rds/seasons_seedp.rds"),
                seedn = readRDS("results/rds/seasons_seedn.rds")
)

state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")

demo_data <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
  mutate(herb_shading_t0 = scale(herb_shading_t0),
         shrub_shading_t0 = scale(shrub_shading_t0))

params <- list(
  surv_mod = VR_mods$surv,
  
  grow_mod = VR_mods$growth,
  grow_sd = readRDS(file = "results/rds/seasons_growth_sd.rds"),
  
  pflower_mod = VR_mods$flower_p,
  
  pseed_mod = VR_mods$seedp,
  
  nseed_mod = VR_mods$seedn,
  nseed_sd = readRDS(file = "results/rds/seasons_seedn_sd.rds"),
  
  seed_surv1 = 0.45,  ## Probability of seed being viable at the next census 
  seed_surv2 = 0.089,  ## Probability of viable seed surviving first year in seed bank. 
  ## (i.e. produced in t, surv to census t+1 (with % above), not germinated in t+1 and 
  seed_surv3 = 0.663,  ## Probability of viable seed in seedbank (yr 1) surviving 2nd year
  
  germ_mean = mean(state_independent_variables$est_germination_rate$germ),     
  
  sdl_surv_mean = boot::inv.logit(fixef(state_independent_variables$sdl_surv)),
  
  sdl_d_int = lme4::fixef(state_independent_variables$sdl_size_d)[1],
  sdl_d_sd = sd(resid(state_independent_variables$sdl_size_d)),
  
  soil_depth = calc_stats(demo_data, "soil_depth"),
  slope = calc_stats(demo_data, "slope"),
  rock = calc_stats(demo_data, "rock"),
  
  herb_center = attr(demo_data$herb_shading_t0, "scaled:center"),
  herb_scale = attr(demo_data$herb_shading_t0, "scaled:scale"),
  shrub_center = attr(demo_data$shrub_shading_t0, "scaled:center"),
  shrub_scale = attr(demo_data$shrub_shading_t0, "scaled:scale")
)


## Starting population vector
data <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
  filter(year_t0 == "2021")
data <- split(data, data$population)
pop_n <- list(CR = 45,
              HK = 202,
              KS = 25,
              RU = 43
)

pop_vec_start <- lapply(data, function(x) (x %>% filter(stage_t0 %in% c("veg","flow")))$ln_stems_t0) %>%
  purrr::map2(., pop_n, ~ rep(.x, round(.y/length(.x),0)))

sdl_n <- lapply(data, function(x) nrow(x %>% filter(stage_t0 == "sdl")))


## -------------------------------------------------------------------------------------------
## Set up environmental values
## -------------------------------------------------------------------------------------------


### Loop through different populations and env_param levels 
localities <- c("Cr", "Hk", "Ks", "Ru")
shrub_shading <- seq(0,6, length.out = 4)
herb_shading <- seq(0,6, length.out = 4)

model <- c("ACCESS1-3", "CESM1-BGC", "CMCC-CM", "MIROC5")
scenario <- c("rcp45", "rcp85")


df_env <- expand.grid(localities = localities, 
                      shrub_shading = shrub_shading, 
                      herb_shading = herb_shading,
                      scenario = scenario,
                      model = model
) 

df_env <- df_env %>% rowid_to_column()

## -------------------------------------------------------------------------------------------
## IBM
## -------------------------------------------------------------------------------------------
gc()
print("start ibm")

# args <- commandArgs(trailingOnly = TRUE)
# taskID = as.integer(args[1])

for(taskID in c(1:512)) {
results <- lapply(as.list(c(1:10)), function(x) {
  a <- ibm_ext_p(i = taskID, df_env = df_env,
                 params = params,
                 fut_clim = fut_clim,
                 clim_hist = hist_clim,
                 pop_vec = pop_vec_start,
                 sdl_n = sdl_n,
                 n_it = n_it)
  return(a)}) %>% bind_rows()

saveRDS(results, file = paste0("results/ibm/df_", taskID, ".rds"))
}







