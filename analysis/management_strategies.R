### script to investigate effort/transplant needed to maintain population

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

source("R/functions_ipmr.R")

## -------------------------------------------------------------------------------------------
## Format CHELSA's time series for IPM
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


# Projections will be from 2022 to 2100

n_it = 79


## -------------------------------------------------------------------------------------------
## Other variables for IPM
## -------------------------------------------------------------------------------------------

VR_FLM <- readRDS("results/VR_FLM.rds")
state_independent_variables <- readRDS("results/state_independent_VR.rds")
lag = 24

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


## Starting population vector

pop_vec <- table(cut(VR_FLM$surv$model$ln_stems_t0[which(VR_FLM$surv$model$year_t0 == "2017")], 
                     breaks=seq(L,U,length.out=(n+1)))) %>%
  as.numeric(as.matrix(.)) * 5

sdl_n <- length(state_independent_variables$sdl_surv@frame$survival_t1[which(state_independent_variables$sdl_surv@frame$year_t0 == "2017")])

## -------------------------------------------------------------------------------------------
## Set up environmental values
## -------------------------------------------------------------------------------------------


### Loop through different populations and env_param levels 
localities <- c("Cr", "Hk", "Ks", "Ru")
shading <- 2
slope <- 20
yrs_between_man <- seq(2, 10, by = 2)
effort_trans <- seq(5, 50, by = 5)
effort_seed <- seq(1000, 10000, by = 1000)

model <- c("ACCESS1", "CESM1", "CMCC", "MIROC5")
scenario <- c("rcp45", "rcp85")


df_env_trans <- expand.grid(localities = localities, 
                      shading = shading, 
                      slope = slope,
                      scenario = scenario,
                      model = model,
                      effort = effort_trans,
                      yrs_between_man = yrs_between_man
) %>% 
  mutate(localities = as.character(localities))

df_env_seed <- expand.grid(localities = localities, 
                            shading = shading, 
                            slope = slope,
                            scenario = scenario,
                            model = model,
                            effort = effort_seed,
                            yrs_between_man = yrs_between_man
) %>% 
  mutate(localities = as.character(localities))

## -------------------------------------------------------------------------------------------
## IPM
## -------------------------------------------------------------------------------------------


### Set up parallel
cl <- makeCluster(detectCores() - 2 ) 
# cl <- makeForkCluster(outfile = "")

clusterExport(cl=cl, c("df_env_trans", "df_env_seed",
                       "man_trans", "man_seedadd", "proto_ipm", 
                       "params", "clim_ts", "pop_vec", "sdl_n", 
                       "U", "L", "n", "n_it", "lag",
                       "sampling_env", "FLM_clim_predict"))

clusterEvalQ(cl, c(library("ipmr"), library("dplyr"), library("forecast")))

df_trans <- parLapplyLB(cl,
                  as.list(c(1:nrow(df_env_trans))),
                  function(x) tryCatch(man_trans(i = x, df_env = df_env_trans,
                                                 params = params, 
                                                 clim_ts = clim_ts,
                                                 n_it = n_it, 
                                                 U = U, L = L, n = n), 
                                       error = function(e) NULL)) %>% 
  bind_rows()

saveRDS(df_trans, file = "results/rds/managment_projections_transplants.rds")

df_seedadd <- parLapplyLB(cl,
                        as.list(c(1:nrow(df_env))),
                        function(x) tryCatch(man_seedadd(i = x, df_env = df_env_seed,
                                                       params = params, 
                                                       clim_ts = clim_ts,
                                                       n_it = n_it, 
                                                       U = U, L = L, n = n), 
                                             error = function(e) NULL)) %>% 
  bind_rows()

saveRDS(df_seedadd, file = "results/rds/managment_projections_seedaddition.rds")

