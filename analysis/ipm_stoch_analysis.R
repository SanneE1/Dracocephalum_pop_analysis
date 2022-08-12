## This script calculates lambda for different environmental variable levels
## under asymptotic behaviour (i.e. long run simulations)

getwd()

library(tidyverse)
library(lme4)
library(patchwork)
library(ipmr)
library(mgcv)
library(parallel)
library(forecast)

args = commandArgs(trailingOnly = T)

setwd(args[1])

source("R/functions_ipmr.R")

VR_FLM <- readRDS("results/rds/VR_FLM.rds")
state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")
climate_models <- readRDS("results/rds/ARIMA_clim_mods.rds")

lag = 24
n_it = 10000
# param/model list 
params <- list(
  surv_mod = VR_FLM$surv,
  s_int = coef(VR_FLM$surv)[1],
  s_stems = coef(VR_FLM$surv)[2],
  s_site_CR = 0,
  s_site_HK = coef(VR_FLM$surv)[3],
  s_site_KS = coef(VR_FLM$surv)[4],
  s_site_RU = coef(VR_FLM$surv)[5],
  s_shading = coef(VR_FLM$surv)[6],
  
  grow_mod = VR_FLM$growth,
  g_int = coef(VR_FLM$growth)[1],
  g_stems = coef(VR_FLM$growth)[2],
  g_site_CR = 0,
  g_site_HK = coef(VR_FLM$growth)[3],
  g_site_KS = coef(VR_FLM$growth)[4],
  g_site_RU = coef(VR_FLM$growth)[5],
  grow_sd = sd(resid(VR_FLM$growth)),
  
  pflower_mod = VR_FLM$flower_p,
  fp_int = coef(VR_FLM$flower_p)[1],
  fp_stems = coef(VR_FLM$flower_p)[2],
  fp_site_CR = 0,
  fp_site_HK = coef(VR_FLM$flower_p)[3],
  fp_site_KS = coef(VR_FLM$flower_p)[4],
  fp_site_RU = coef(VR_FLM$flower_p)[5],
  fp_slope = coef(VR_FLM$flower_p)[6],
  
  pabort_mod = VR_FLM$abort_p,
  ab_int = coef(VR_FLM$abort_p)[1],
  ab_stems = coef(VR_FLM$abort_p)[2],
  ab_site_CR = 0,
  ab_site_HK = coef(VR_FLM$abort_p)[3],
  ab_site_KS = coef(VR_FLM$abort_p)[4],
  ab_site_RU = coef(VR_FLM$abort_p)[5],

  nseed_mod = VR_FLM$n_seeds,
  ns_int = coef(VR_FLM$n_seeds)[1],
  ns_stems = coef(VR_FLM$n_seeds)[2],
  ns_site_CR = 0,
  ns_site_HK = coef(VR_FLM$n_seeds)[3],
  ns_site_KS = coef(VR_FLM$n_seeds)[4],
  ns_site_RU = coef(VR_FLM$n_seeds)[5],
  
  seed_surv1 = 0.57826,  ## Probability of seed being viable at the next census 
  seed_surv2 = 0.15374,  ## Probability of viable seed surviving first year in seed bank. 
  ## (i.e. produced in t, surv to census t+1 (with % above), not germinated in t+1 and 
  ## then remain viable in seedbank till census t+2 with % of 0.15374) 
  seed_surv3 = 0.66658,  ## Probability of viable seed in seedbank (yr 1) surviving 2nd year
  
  germ_mean = mean(state_independent_variables$est_germination_rate$germ),     
  
  sdl_surv_mod = state_independent_variables$sdl_surv,
  sdl_s_int = lme4::fixef(state_independent_variables$sdl_surv)[1],
  
  sdl_d_int = lme4::fixef(state_independent_variables$sdl_size_d)[1],
  sdl_size_d_sd = sd(resid(state_independent_variables$sdl_size_d))
)



## Set integration params
L <- min(VR_FLM$growth$model$ln_stems_t0, na.rm = T)
U <- max(VR_FLM$growth$model$ln_stems_t0, na.rm = T) * 1.1
n = 100


##------------------------------------------------------------------------------
## Calculate lambda at different co-variate levels
##------------------------------------------------------------------------------

### Loop through different populations and env_param levels 
localities <- c("Cr", "Hk", "Ks", "Ru")
shading <- seq(0,8, length.out = 9)
slope <- seq(0,80, length.out = 4)

model <- c("ACCESS1", "CESM1", "CMCC", "MIROC5")
scenario <- c("rcp45", "rcp85")

df_env_hist <- expand.grid(localities = localities, 
                           shading = shading, 
                           slope = slope,
                           time = "hist",
                           scenario = NA,
                           model = NA
) 
df_env_fut <- expand.grid(localities = localities, 
                          shading = shading, 
                          slope = slope,
                          time = "future",
                          scenario = scenario,
                          model = model
) 

df_env <- rbind(df_env_hist, df_env_fut) %>% 
  mutate(localities = as.character(localities))

rep <- rep(c(1:nrow(df_env)), each = 10)


### Set up parallel
cl <- makeForkCluster(outfile = "")  
clusterExport(cl=cl, c("df_env", "ipm_loop", "run_ipm", 
                       "params", "climate_models", 
                       "U", "L", "n", "n_it", "lag",
                       "sampling_env", "FLM_clim_predict"))

clusterEvalQ(cl, c(library("ipmr"), library("dplyr"), library("forecast")))

df <- parLapplyLB(cl,
                  as.list(rep),
                  function(x) tryCatch(ipm_loop(i = x, df_env = df_env,
                                                params = params, 
                                                climate_models = climate_models,
                                                n_it = n_it, 
                                                U = U, L = L, n = n),
                                       error = function(e) NULL)) %>% 
  bind_rows()

stopCluster(cl)

write.csv(df, file = "results/overview_lambda_env_levels.csv", 
          row.names = F)

