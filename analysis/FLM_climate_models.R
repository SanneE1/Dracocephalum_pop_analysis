library(tidyverse)
library(lme4)
library(mgcv)
library(patchwork)
library(MuMIn)
library(parallel)

### FLM climate model for Growth

data_for_modeling = "data/Dracocephalum_with_vital_rates.csv"
CHELSA_data = "data/CHELSA_data.csv"
lag = 24  # number of months prior to census to include

source("R/functions_GAM.R")
source("R/functions_cross_validation.R")

## Format dataframe:
### response & covariates and add climate variable as lagg since census in columns 
### i.e. temp.00 = month of census, temp.01 = temperature one month before census etc.
### Looks like for the GAM you need to use a variable name ("temp") and then "." and then add the numbers

demo_data <- read.csv(data_for_modeling) %>% rowwise() %>%
  mutate(population = population,
         tot_shading_tm2 = ifelse(any(!is.na(c(herb_shading_tm2, shrub_shading_tm2))), 
                                  sum(herb_shading_tm2, shrub_shading_tm2, na.rm = T),
                                  NA),
         tot_shading_tm1 = ifelse(any(!is.na(c(herb_shading_tm1, shrub_shading_tm1))), 
                                  sum(herb_shading_tm1, shrub_shading_tm1, na.rm = T),
                                  NA),
         tot_shading_t0 = ifelse(any(!is.na(c(herb_shading_t0, shrub_shading_t0))), 
                                 sum(herb_shading_t0, shrub_shading_t0, na.rm = T),
                                 NA),
         tot_shading_t1 = ifelse(any(!is.na(c(herb_shading_t1, shrub_shading_t1))), 
                                 sum(herb_shading_t1, shrub_shading_t1, na.rm = T),
                                 NA)
  ) %>% ungroup() %>%
  filter(is.finite(ln_stems_t0) & !(stage_t0 == "sdl"))


# climate data frame for VRs with responses measured in t1
clim_data_t1 <- read.csv(CHELSA_data) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("tas_scaled", "pr_scaled", "pet_scaled"), 
                        demo_data = demo_data, 
                        response_t1 = T,
                        lag = lag)
# Climate data frame for VRs with responses measured in t0 
clim_data_t0 <- read.csv(CHELSA_data) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("tas_scaled", "pr_scaled", "pet_scaled"), 
                        demo_data = demo_data, 
                        response_t1 = F,
                        lag = lag)

### Merge the two dataframes and keep only entries for which we have full climate data available

data_t1 <- left_join(demo_data, clim_data_t1) %>% 
  mutate(year_t0 = factor(year_t0))

data_t0 <- left_join(demo_data, clim_data_t0) %>% 
  mutate(year_t0 = factor(year_t0))

# Add shading matrices for te() functions below
data_t1$tot_shading_m <- apply(data_t1, 1, function(x) as.integer(c(rep(x['tot_shading_tm1'], 6),
                                                                    rep(x['tot_shading_t0'], 12),
                                                                    rep(x['tot_shading_t1'], 7)))) %>% t
data_t1$slope_m <- apply(data_t1, 1, function(x) as.integer(c(rep(x['slope'], 25)))) %>% t
data_t1$rock_m <- apply(data_t1, 1, function(x) as.integer(c(rep(x['rock'], 25)))) %>% t
data_t1$soil_m <- apply(data_t1, 1, function(x) as.integer(c(rep(x['soil_depth'], 25)))) %>% t


data_t0$tot_shading_m <- apply(data_t0, 1, function(x) as.integer(c(rep(x['tot_shading_tm2'], 6),
                                                                    rep(x['tot_shading_tm1'], 12),
                                                                    rep(x['tot_shading_t0'], 7)))) %>% t
data_t0$slope_m <- apply(data_t0, 1, function(x) as.integer(c(rep(x['slope'], 25)))) %>% t
data_t0$rock_m <- apply(data_t0, 1, function(x) as.integer(c(rep(x['rock'], 25)))) %>% t
data_t0$soil_m <- apply(data_t0, 1, function(x) as.integer(c(rep(x['soil_depth'], 25)))) %>% t

formula_parts <- c(
  tas = "s(lags, by= tas_scaledcovar, bs = 'cs')", 
  tas_int = "ti(lags, tot_shading_t0, by= tas_scaledcovar, bs = 'cs')",
  pr = "s(lags, by= pr_scaledcovar, bs = 'cs')",
  pr_int = "ti(lags, tot_shading_t0, by= pr_scaledcovar, bs = 'cs')",
  pet = "s(lags, by= pet_scaledcovar, bs = 'cs')", 
  pet_int = "ti(lags, tot_shading_t0, by= pet_scaledcovar, bs = 'cs')",
  shade = "s(tot_shading_t0, k = 8, bs = 'cs')", 
  slope = "s(slope, k = 8, bs = 'cs')",
  rock = "s(rock, k = 8, bs = 'cs')", 
  soil = "s(soil_depth, k = 8, bs = 'cs')"
)


## -------------------------------------------------------
## Start parallel
## -------------------------------------------------------

cores <- detectCores() - 4

cl <- makeCluster(cores)

clusterEvalQ(cl, c(library(mgcv), library(MuMIn), library(dplyr)))
clusterExport(cl=cl, c("gam.crossvalidation"))

## -------------------------------------------------------
## Survival
## -------------------------------------------------------


surv_data <- data_t1 %>%
  select(survival_t1, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth, tot_shading_m) %>%
  filter(complete.cases(.))

# Basemodel
base_surv <- "survival_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
formula_surv_full <- paste(c(base_surv, formula_parts), collapse = "+")  

#FLM models
# surv_full <- gam(as.formula(formula_surv_full),
#                  data=surv_data ,
#                  family = binomial(link = "logit"),
#                  method="GCV.Cp",
#                  gamma=1.4,
#                  select = T,
#                  na.action = "na.fail")
# 
# summary(surv_full)

# plot(mgcViz::getViz(surv_full))

surv_combinations <- generate_combinations(formula_parts) %>% 
  paste(base_surv, ., sep = "+")

clusterExport(cl=cl, c("surv_data"))

surv_mods <- parLapply(cl, 
                         surv_combinations, 
                         function(x) gam(as.formula(x),
                                         data=surv_data ,
                                         family = binomial(link = "logit"),
                                         method="GCV.Cp",
                                         gamma=1.4,
                                         select = T,
                                         na.action = "na.fail")
)
saveRDS(surv_mods, file = "results/rds/surv_flm_mods.rds")


surv_cross <- parLapply(cl,
                          surv_mods,
                          function(x) gam.crossvalidation(mod = x, 
                                                          response_column =  "survival_t1", 
                                                          subset_length = 10)) %>%
  bind_rows()

saveRDS(surv_cross, file = "results/rds/surv_flm_cross.rds")

 
# View(surv_cross)
# 
# surv_combinations[c(232, 176,84)]
# 
# plot(mgcViz::getViz(surv_mods[[232]]))
# plot(mgcViz::getViz(surv_mods[[250]]))
# plot(mgcViz::getViz(surv_mods[[176]]))
# plot(mgcViz::getViz(surv_mods[[84]]))
# 
# 
# saveRDS(surv_mods[[176]], file = "results/rds/surv_flm.rds")


## -------------------------------------------------------
## Growth
## -------------------------------------------------------

growth_data <- data_t1 %>%
  filter(survival_t1 == 1) %>%
  select(ln_stems_t1, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth, tot_shading_m) %>%
  filter(complete.cases(.))

# Basemodel
base_growth <- "ln_stems_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"

formula_growth_full <- paste(c(base_growth, formula_parts), collapse = "+")  

#FLM models
# 
# growth_full <- gam(as.formula(formula_growth_full),
#                    data=growth_data,
#                    method="GCV.Cp",
#                    gamma=1.4,
#                    select = T,
#                    na.action = "na.fail")
# 
# summary(growth_full)
# 
# plot(mgcViz::getViz(growth_full))

growth_combinations <- generate_combinations(formula_parts) %>% 
  paste(base_growth, ., sep = "+")

clusterExport(cl=cl, c("growth_data"))

growth_mods <- parLapply(cl, 
                         growth_combinations, 
                      function(x) gam(as.formula(x),
                                      data=growth_data,
                                      method="GCV.Cp",
                                      gamma=1.4,
                                      select = T,
                                      na.action = "na.fail")
  )
saveRDS(growth_mods, file = "results/rds/growth_flm_mods.rds")


growth_cross <- parLapply(cl,
                          growth_mods,
                          function(x) gam.crossvalidation(mod = x, 
                                                          response_column =  "ln_stems_t1", 
                                                          subset_length = 3))

saveRDS(growth_cross, file = "results/rds/growth_flm_cross.rds")


# View(growth_cross)
# 
# growth_combinations[c(189,202,221,234)]
# 
# plot(mgcViz::getViz(growth_mods[[189]]))
# plot(mgcViz::getViz(growth_mods[[202]]))
# 
# saveRDS(growth_mods[[189]], file = "results/rds/growth_flm.rds")


## -------------------------------------------------------
## Flower probability
## -------------------------------------------------------

flowp_data <- data_t0 %>% 
  select(flower_p_t0, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth, tot_shading_m) %>%
  filter(complete.cases(.))

# Basemodel
base_flowp <- "flower_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"

formula_flowp_full <- paste(c(base_flowp, formula_parts), collapse = "+")  

#FLM models

# flowp_full <- gam(as.formula(formula_flowp_full),
#                   data=flowp_data ,
#                   family = binomial(link = "logit"),
#                   method="GCV.Cp",
#                   gamma=1.4,
#                   select = T,
#                   na.action = "na.fail")
# 
# 
# summary(flowp_full)
# 
# plot(mgcViz::getViz(flowp_full))


flowp_combinations <- generate_combinations(formula_parts[-c(7,10)]) %>% 
  paste(base_flowp, ., sep = "+")

clusterExport(cl=cl, c("flowp_data"))

flowp_mods <- parLapply(cl, 
                         flowp_combinations, 
                         function(x) gam(as.formula(x),
                                         data=flowp_data ,
                                         family = binomial(link = "logit"),
                                         method="GCV.Cp",
                                         gamma=1.4,
                                         select = T,
                                         na.action = "na.fail")
)
saveRDS(flowp_mods, file = "results/rds/flowp_flm_mods.rds")


flowp_cross <- parLapply(cl,
                         flowp_mods,
                          function(x) gam.crossvalidation(mod = x, 
                                                          response_column =  "flower_p_t0", 
                                                          subset_length = 3)) %>%
  bind_rows()


saveRDS(flowp_cross, file = "results/rds/flowp_flm_cross.rds")


# View(flowp_cross)
# 
# flowp_combinations[c(249,227,191,168,189,113)]
# 
# summary(flowp_mods[[249]])
# plot(mgcViz::getViz(flowp_mods[[249]]))

summary(flowp_mods[[227]])
plot(mgcViz::getViz(flowp_mods[[227]]))

# summary(flowp_mods[[191]])
# plot(mgcViz::getViz(flowp_mods[[191]]))
# 
# summary(flowp_mods[[168]])
# plot(mgcViz::getViz(flowp_mods[[168]]))
# 
# summary(flowp_mods[[189]])
# plot(mgcViz::getViz(flowp_mods[[189]]))
# 
# summary(flowp_mods[[113]])
# plot(mgcViz::getViz(flowp_mods[[113]]))



saveRDS(flowp_mods[[227]], file = "results/rds/flowp_flm.rds")


## -------------------------------------------------------
## Seed probability
## -------------------------------------------------------

seedp_data <- data_t0 %>%
  filter(flower_p_t0 == 1) %>%
  select(seed_p_t0, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth) %>%
  filter(complete.cases(.))
  
# Basemodel
base_seedp <- "seed_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
formula_seedp_full <- paste(c(base_seedp, formula_parts), collapse = "+")  

#FLM models
# 
# seedp_full <- gam(as.formula(formula_seedp_full),
#                   data=seedp_data,
#                   family = binomial(link = "logit"),
#                   method="GCV.Cp",
#                   gamma=1.4,
#                   select = T,
#                   na.action = "na.fail")
# 
# summary(seedp_full)
# 
# plot(mgcViz::getViz(seedp_full))


seedp_combinations <- generate_combinations(formula_parts[-c(8,10)]) %>% 
  paste(base_seedp, ., sep = "+")

clusterExport(cl=cl, c("seedp_data"))

seedp_mods <- parLapply(cl, 
                        seedp_combinations, 
                        function(x) gam(as.formula(x),
                                        data=seedp_data,
                                        family = binomial(link = "logit"),
                                        method="GCV.Cp",
                                        gamma=1.4,
                                        select = T,
                                        na.action = "na.fail")
)
saveRDS(seedp_mods, file = "results/rds/seedp_flm_mods.rds")


seedp_cross <- parLapply(cl,
                         seedp_mods,
                         function(x) gam.crossvalidation(mod = x, 
                                                         response_column =  "seed_p_t0", 
                                                         subset_length = 3)) %>%
  bind_rows()


saveRDS(seedp_cross, file = "results/rds/seedp_flm_cross.rds")


# View(seedp_cross)
# 
# seedp_combinations[c(253,238,251,255,249,190)]
# 
# summary(seedp_mods[[253]])
# plot(mgcViz::getViz(seedp_mods[[253]]))
# 
# summary(seedp_mods[[238]])
# plot(mgcViz::getViz(seedp_mods[[238]]))
# 
# summary(seedp_mods[[251]])
# plot(mgcViz::getViz(seedp_mods[[251]]))
# 
# summary(seedp_mods[[255]])
# plot(mgcViz::getViz(seedp_mods[[255]]))
# 
# summary(seedp_mods[[249]])
# plot(mgcViz::getViz(seedp_mods[[249]]))

summary(seedp_mods[[190]])
plot(mgcViz::getViz(seedp_mods[[190]]))


saveRDS(seedp_mods[[190]], file = "results/rds/seedp_flm.rds")

## -------------------------------------------------------
## Seed numbers
## -------------------------------------------------------

seedn_data <- data_t0 %>%
  filter(flower_p_t0 == 1 & seed_p_t0 == 1) %>%
  select(est_seed_n_t0, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth) %>%
  filter(complete.cases(.))

# Basemodel
base_seedn <- "est_seed_n_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
formula_seedn_full <- paste(c(base_seedn, formula_parts), collapse = "+")  

#FLM models
# 
# seedn_full <- gam(as.formula(formula_seedn_full),
#                   data=seedn_data,
#                   family = Gamma(link = "log"),
#                   method="GCV.Cp",
#                   gamma=1.4,
#                   select = T,
#                   na.action = "na.fail")
# 
# summary(seedn_full)
# plot(mgcViz::getViz(seedn_full))


seedn_combinations <- generate_combinations(formula_parts[-c(3,9)]) %>% 
  paste(base_seedn, ., sep = "+")

clusterExport(cl=cl, c("seedn_data"))

seedn_mods <- parLapply(cl, 
                        seedn_combinations, 
                        function(x) gam(as.formula(x),
                                        data=seedn_data,
                                        family = Gamma(link = "log"),
                                        method="GCV.Cp",
                                        gamma=1.4,
                                        select = T,
                                        na.action = "na.fail")
)
saveRDS(seedn_mods, file = "results/rds/seedn_flm_mods.rds")


seedn_cross <- parLapply(cl,
                         seedn_mods,
                         function(x) gam.crossvalidation(mod = x, 
                                                         response_column =  "seed_p_t0", 
                                                         subset_length = 3)) %>%
  bind_rows()


saveRDS(seedn_cross, file = "results/rds/seedn_flm_cross.rds")


# View(seedn_cross)

# seedn_combinations[c(56, 64, 15)]

# summary(seedn_mods[[56]])
# plot(mgcViz::getViz(seedn_mods[[56]]))

# summary(seedn_mods[[64]])
# plot(mgcViz::getViz(seedn_mods[[64]]))

summary(seedn_mods[[15]])
plot(mgcViz::getViz(seedn_mods[[15]]))


saveRDS(seedn_mods[[15]], file = "results/rds/seedn_flm.rds")

## -------------------------------------------------------
## Save results
## -------------------------------------------------------

saveRDS(list(surv = readRDS(file = "results/rds/surv_flm.rds"),
             growth = readRDS(file = "results/rds/growth_flm.rds"),
             flower_p = readRDS(file = "results/rds/flowp_flm.rds"),
             seedp = readRDS(file = "results/rds/seedp_flm.rds"),
             seedn = readRDS(file = "results/rds/seedn_flm.rds")),
        file = "results/rds/VR_FLM.rds")


