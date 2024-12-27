### script to run extinction probability analysis
rm(list = ls())

# getwd()

library(tidyverse)
library(lme4)
library(patchwork)
library(ipmr)
library(mgcv)
library(parallel)
library(forecast)

# args = commandArgs(trailingOnly = T)

# setwd(args[1])

source("R/functions_ibm.R")
source("R/functions_ipmr.R")

# Projections will be from 2022 to 2100
n_it = 79
lag = 24

## -------------------------------------------------------------------------------------------
## Format CHELSA's time series for IBM
## -------------------------------------------------------------------------------------------

fut_clim <- read.csv("data/CHELSA_future_ts_formatted.csv") %>% 
  filter(complete.cases(.))

clim_ts <- fut_clim %>% 
  filter(scenario != "historical") %>%
  dplyr::select(c(locality, model, scenario, month, year, pr_scaled, tas_scaled, pet_scaled)) %>%
  pivot_longer(., contains("scaled"), values_to = "value", names_to = "variable") %>%
  split(., list(.$locality, .$model, .$scenario, .$variable)) %>%
  lapply(., function(x) {
    x %>% mutate(value = ts(value, frequency = 12 , start = c(2006,1)))
  })

hist_clim <- readRDS("results/rds/ARIMA_clim_mods.rds")$clim_hist_model %>%
  lapply(., function(x)
    simulate(x, nsim = ((n_it * 12) + (3*lag))) %>%
      ts(., start= c(2019,1), frequency = 12))

## -------------------------------------------------------------------------------------------
## Other variables for IBM
## -------------------------------------------------------------------------------------------

VR_FLM <- readRDS("results/rds/VR_FLM.rds")
state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")

demo_data <- read.csv("data/Dracocephalum_with_vital_rates.csv") 

params <- list(
  surv_mod = VR_FLM$surv,
  
  grow_mod = VR_FLM$growth,
  grow_sd = sd(resid(VR_FLM$growth)),
  
  pflower_mod = VR_FLM$flower_p,
  
  pseed_mod = VR_FLM$seedp,
  
  nseed_mod = VR_FLM$seedn,
  nseed_sd = sd(resid(VR_FLM$seedn)),
  
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
  rock = calc_stats(demo_data, "rock")
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

pop_vec <- lapply(data, function(x) (x %>% filter(stage_t0 %in% c("veg","flow")))$ln_stems_t0) %>%
  purrr::map2(., pop_n, ~ rep(.x, round(.y/length(.x),0)))

sdl_n <- lapply(data, function(x) nrow(x %>% filter(stage_t0 == "sdl")))


## -------------------------------------------------------------------------------------------
## Set up environmental values
## -------------------------------------------------------------------------------------------


### Loop through different populations and env_param levels 
localities <- c("Cr", "Hk", "Ks", "Ru")
shading <- seq(0,6, length.out = 4)

model <- c("ACCESS1", "CESM1", "CMCC", "MIROC5")
scenario <- c("rcp45", "rcp85")


df_env <- expand.grid(localities = localities, 
                      shading = shading, 
                      scenario = scenario,
                      model = model
) %>% 
  rbind(.,
        expand.grid(
          localities = localities, 
          shading = shading,  
          scenario = NA,
          model = "No change"
        )) %>% 
  mutate(localities = as.character(localities)) 

df_env <- df_env[rep(1:nrow(df_env), 100),] %>% rowid_to_column()

## -------------------------------------------------------------------------------------------
## IBM
## -------------------------------------------------------------------------------------------
gc()
print("start ibm")
### Set up parallel    --------------------- For some reason not working. Easier right now to just let lapply run overnight
# cl <- makeCluster(2, outfile = "")
# # cl <- makeForkCluster(outfile = "")
# 
# clusterExport(cl=cl, c("df_env", "ibm_ext_p", "yearly_loop", "mod_pred",
#                        "params", "clim_ts", "hist_clim",
#                        "pop_vec", "sdl_n",
#                        "n_it", "lag", "sampling_env"))
# 
# clusterEvalQ(cl, {
#   library(tidyverse)
#   library(lme4)
#   library(patchwork)
#   library(ipmr)
#   library(mgcv)
#   library(parallel)
#   library(forecast) }
#   )

for(x in c(1:nrow(df_env))) {
  a <- ibm_ext_p(i = x, df_env = df_env,
                 params = params,
                 clim_ts = clim_ts,
                 clim_hist = hist_clim,
                 pop_vec = pop_vec,
                 sdl_n = sdl_n,
                 n_it = n_it)
  
  write.csv(a, paste0("results/imb/df_", x, ".csv"))
  
}
# 
# df <- lapply(as.list(c(1:nrow(df_env))),
#                   function(x) ibm_ext_p(i = x, df_env = df_env,
#                                                  params = params,
#                                                  clim_ts = clim_ts,
#                                                  clim_hist = hist_clim,
#                                                  pop_vec = pop_vec,
#                                                  sdl_n = sdl_n,
#                                                  n_it = n_it)) %>%
#   bind_rows()

# stopCluster(cl)

saveRDS(df, file = "results/rds/extinction_probability.rds")

