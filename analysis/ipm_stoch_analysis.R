## This script calculates lambda for different environmental variable levels
## under asymptotic behaviour (i.e. long run simulations)

args = commandArgs(trailingOnly = T)
setwd(args[2])

library(dplyr)
library(lme4)
library(ipmr)
library(glmnet)
library(parallel)
library(forecast)

source("R/functions_ipmr.R")

state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")
climate_models <- readRDS("results/rds/ARIMA_clim_mods.rds")

demo_data <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
  mutate(herb_shading_t0 = scale(herb_shading_t0),
         shrub_shading_t0 = scale(shrub_shading_t0))

lag = 24
n_it = 100 #5000

# param/model list 
params <- list(
  surv_mod = readRDS("results/rds/seasons_surv.rds"),
  
  grow_mod = readRDS('results/rds/seasons_growth.rds'),
  grow_sd = readRDS("results/rds/seasons_growth_sd.rds"),
  
  pflower_mod = readRDS("results/rds/seasons_flowp.rds"),
  
  seedp_mod = readRDS("results/rds/seasons_seedp.rds"),
  
  seedn_mod = readRDS('results/rds/seasons_seedn.rds'),
  
  seed_surv1 = 0.45,  ## Probability of seed being viable at the next census 
  seed_surv2 = 0.089,  ## Probability of viable seed surviving first year in seed bank. 
  ## (i.e. produced in t, surv to census t+1 (with % above), not germinated in t+1 and 
  seed_surv3 = 0.663,  ## Probability of viable seed in seedbank (yr 1) surviving 2nd year
  
  germ_mean = mean(state_independent_variables$est_germination_rate$germ),     
  
  sdl_surv_mod = state_independent_variables$sdl_surv,
  sdl_s_int = lme4::fixef(state_independent_variables$sdl_surv)[1],
  
  sdl_d_int = lme4::fixef(state_independent_variables$sdl_size_d)[1],
  sdl_size_d_sd = sd(resid(state_independent_variables$sdl_size_d)),
  
  herb_center = attr(demo_data$herb_shading_t0, "scaled:center"),
  herb_scale = attr(demo_data$herb_shading_t0, "scaled:scale"),
  shrub_center = attr(demo_data$shrub_shading_t0, "scaled:center"),
  shrub_scale = attr(demo_data$shrub_shading_t0, "scaled:scale"))


## Set integration params
L <- 0
U <- 4.672829 * 1.1
n = 100


##------------------------------------------------------------------------------
## Calculate lambda at different co-variate levels
##------------------------------------------------------------------------------

### Loop through different populations and env_param levels 
localities <- c("Cr", "Hk", "Ks", "Ru")

shrub_shading <- seq(0,6, length.out = 4)
herb_shading <- seq(0,6, length.out = 4)

slope <- seq(0, 50, length.out = 4)
rock <- seq(0, 80, length.out = 4)
soil_depth <- seq(0,10, length.out = 4)

model <- c("ACCESS1-3", "CESM1-BGC", "CMCC-CM", "MIROC5")
scenario <- c("rcp45", "rcp85")

hist_shrub <- expand.grid(localities = localities,
                    shrub_shading = shrub_shading,
                    herb_shading = 2,
                    slope = 14,
                    rock = 30,
                    soil_depth = 5,
                    time = "hist",
                    scenario = NA,
                    model = NA
) 
hist_herb <- expand.grid(localities = localities,
                         shrub_shading = 2,
                         herb_shading = herb_shading,
                         slope = 14,
                         rock = 30,
                         soil_depth = 5,
                         time = "hist",
                         scenario = NA,
                         model = NA
)
fut_herb <- expand.grid(localities = localities, 
                   shrub_shading = 2,
                   herb_shading = herb_shading,
                   slope = 14,
                   rock = 30,
                   soil_depth = 5,
                   time = "future",
                   scenario = scenario,
                   model = model
) 
fut_shrub <- expand.grid(localities = localities, 
                        shrub_shading = shrub_shading,
                        herb_shading = 2,
                        slope = 14,
                        rock = 30,
                        soil_depth = 5,
                        time = "future",
                        scenario = scenario,
                        model = model
) 
slope <- expand.grid(localities = "Cr", 
                     shrub_shading = 2,
                     herb_shading = 2,
                     slope = slope,
                   rock = 30,
                   soil_depth = 5,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

rock <- expand.grid(localities = "Cr", 
                    shrub_shading = 2,
                    herb_shading = 2,
                    slope = 14,
                   rock = rock,
                   soil_depth = 5,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

soil <- expand.grid(localities = "Cr", 
                    shrub_shading = 2,
                    herb_shading = 2,
                    slope = 14,
                   rock = 30,
                   soil_depth = soil_depth,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

df_env <- rbind(hist_herb, hist_shrub, fut_herb, fut_shrub, rock, slope, soil) %>% 
  mutate(localities = as.character(localities))

run <- c(1:nrow(df_env))


### Set up parallel

# # Local machine
cl <- makeCluster(detectCores()-4)

clusterExport(cl=cl, c("df_env", "ipm_loop", "run_ipm",
                       "params", "climate_models", "run",
                       "U", "L", "n", "n_it",
                       "sampling_env", "mod_pred"))

clusterEvalQ(cl, c(library("ipmr"), library("dplyr"), library("forecast"), library("glmnet") ))

df <- clusterApply(cl,
                  as.list(c(289:length(run))),
                  function(x) tryCatch(
                    ipm_loop(i = run[x], df_env = df_env,
                                                params = params,
                                                climate_models = climate_models,
                                                n_it = n_it,
                                                U = U, L = L, n = n,
                                                save = T),
                                       error = function(e) NULL)
                  ) %>%
  bind_rows()

stopCluster(cl)

write.csv(df, file = "results/overview_lambda_env_levels.csv",
          row.names = F)


## Run as an array job
# taskID <- as.integer(sys.getenv("SLURM_ARRAY_TASK_ID"))
# 
# df <- ipm_loop(i = run[taskID], df_env = df_env,
#                params = params,
#                climate_models = climate_models,
#                n_it = n_it,
#                U = U, L = L, n = n,
#                save = T)






