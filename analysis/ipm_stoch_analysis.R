## This script calculates lambda for different environmental variable levels
## under asymptotic behaviour (i.e. long run simulations)

library(dplyr)
library(lme4)
library(ipmr)
library(mgcv)
library(parallel)
library(forecast)

source("R/functions_ipmr.R")

VR_FLM <- readRDS("results/rds/VR_FLM.rds")
state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")
climate_models <- readRDS("results/rds/ARIMA_clim_mods.rds")

lag = 24
n_it = 5000

# param/model list 
params <- list(
  surv_mod = VR_FLM$surv,
  s_int = coef(VR_FLM$surv)[1],
  s_stems = coef(VR_FLM$surv)[2],
  s_site_CR = 0,
  s_site_HK = coef(VR_FLM$surv)[3],
  s_site_KS = coef(VR_FLM$surv)[4],
  s_site_RU = coef(VR_FLM$surv)[5],
  s_slope = coef(VR_FLM$surv)[6],
  s_rock = coef(VR_FLM$surv)[7],
  s_soil_depth = coef(VR_FLM$surv)[8],
  
  grow_mod = VR_FLM$growth,
  grow_sd = sd(resid(VR_FLM$growth)),
  g_int = coef(VR_FLM$growth)[1],
  g_stems = coef(VR_FLM$growth)[2],
  g_site_CR = 0,
  g_site_HK = coef(VR_FLM$growth)[3],
  g_site_KS = coef(VR_FLM$growth)[4],
  g_site_RU = coef(VR_FLM$growth)[5],
  g_slope = coef(VR_FLM$growth)[6],
  g_rock = coef(VR_FLM$growth)[7],
  g_soil_depth = coef(VR_FLM$growth)[8],
  
  pflower_mod = VR_FLM$flower_p,
  fp_int = coef(VR_FLM$flower_p)[1],
  fp_stems = coef(VR_FLM$flower_p)[2],
  fp_site_CR = 0,
  fp_site_HK = coef(VR_FLM$flower_p)[3],
  fp_site_KS = coef(VR_FLM$flower_p)[4],
  fp_site_RU = coef(VR_FLM$flower_p)[5],
  fp_slope = coef(VR_FLM$flower_p)[6],
  fp_rock = coef(VR_FLM$flower_p)[7],
  fp_soil_depth = coef(VR_FLM$flower_p)[8],
  
  seedp_mod = VR_FLM$seedp,
  sp_int = coef(VR_FLM$seedp)[1],
  sp_stems = coef(VR_FLM$seedp)[2],
  sp_site_CR = 0,
  sp_site_HK = coef(VR_FLM$seedp)[3],
  sp_site_KS = coef(VR_FLM$seedp)[4],
  sp_site_RU = coef(VR_FLM$seedp)[5],
  sp_slope = coef(VR_FLM$seedp)[6],
  sp_rock = coef(VR_FLM$seedp)[7],
  sp_soil_depth = coef(VR_FLM$seedp)[8],
  
  seedn_mod = VR_FLM$seedn,
  sn_int = coef(VR_FLM$seedn)[1],
  sn_stems = coef(VR_FLM$seedn)[2],
  sn_site_CR = 0,
  sn_site_HK = coef(VR_FLM$seedn)[3],
  sn_site_KS = coef(VR_FLM$seedn)[4],
  sn_site_RU = coef(VR_FLM$seedn)[5],
  sn_slope = coef(VR_FLM$seedn)[6],
  sn_rock = coef(VR_FLM$seedn)[7],
  sn_soil_depth = coef(VR_FLM$seedn)[8],
  
  seed_surv1 = 0.45,  ## Probability of seed being viable at the next census 
  seed_surv2 = 0.089,  ## Probability of viable seed surviving first year in seed bank. 
  ## (i.e. produced in t, surv to census t+1 (with % above), not germinated in t+1 and 
  seed_surv3 = 0.663,  ## Probability of viable seed in seedbank (yr 1) surviving 2nd year
  
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
shading <- seq(0,6, length.out = 4)
slope <- seq(0, 50, length.out = 6)
rock <- seq(0, 80, length.out = 6)
soil_depth <- seq(0,10, length.out = 6)

model <- c("ACCESS1", "CESM1", "CMCC", "MIROC5")
scenario <- c("rcp45", "rcp85")

hist <- expand.grid(localities = localities, 
                           shading = shading,
                    slope = 14,
                    rock = 30,
                    soil_depth = 5,
                           time = "hist",
                           scenario = NA,
                           model = NA
) 
fut <- expand.grid(localities = localities, 
                          shading = shading, 
                          slope = 14,
                          rock = 30,
                          soil_depth = 5,
                          time = "future",
                          scenario = scenario,
                          model = model
) 

slope <- expand.grid(localities = "CR", 
                   shading = 3, 
                   slope = slope,
                   rock = 30,
                   soil_depth = 5,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

rock <- expand.grid(localities = localities, 
                   shading = 3, 
                   slope = 14,
                   rock = rock,
                   soil_depth = 5,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

soil <- expand.grid(localities = localities, 
                   shading = 3, 
                   slope = 14,
                   rock = 30,
                   soil_depth = soil_depth,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

df_env <- rbind(hist, fut, rock, slope, soil) %>% 
  mutate(localities = as.character(localities))

rep <- rep(c(1:nrow(df_env)))


# ### Set up parallel
cl <- makeCluster(detectCores() - 1)
clusterExport(cl=cl, c("df_env", "ipm_loop", "run_ipm",
                       "params", "climate_models",
                       "U", "L", "n", "n_it", "lag",
                       "sampling_env", "FLM_clim_predict"))

clusterEvalQ(cl, c(library("ipmr"), library("dplyr"), library("forecast")))

df <- parLapply(cl,
                  as.list(rep),
                  function(x) ipm_loop(i = x, df_env = df_env,
                                                params = params,
                                                climate_models = climate_models,
                                                n_it = n_it,
                                                U = U, L = L, n = n,
                                                save = T)) %>%
  bind_rows()

stopCluster(cl)

write.csv(df, file = "results/overview_lambda_env_levels.csv",
          row.names = F)
