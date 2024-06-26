library(tidyverse)
library(lme4)
library(mgcv)

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



## -------------------------------------------------------
## Survival
## -------------------------------------------------------


# Basemodel
base_surv <- "survival_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're') + 
s(lags, k=lag, by= tas_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pr_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pet_scaledcovar, bs = 'cs') +
s(tot_shading_t0, k = 8, bs = 'cs') + slope + rock + soil_depth" 

#FLM models

surv_mod <- gam(as.formula(base_surv),
      data=data_t1 ,
      family = binomial(link = "logit"),
      method="GCV.Cp",
      gamma=1.4,
      select = T)


plot_spline_coeff(best_model = surv_mod, tas = T, vital_rate = "Survival")
plot_spline_coeff(best_model = surv_mod, pr = T, vital_rate = "Survival")
plot_spline_coeff(best_model = surv_mod, pet = T, vital_rate = "Survival")
plot_spline_coeff(best_model = surv_mod, shade = T, vital_rate = "Survival")


## -------------------------------------------------------
## Growth
## -------------------------------------------------------

# Basemodel
base_growth <- "ln_stems_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're') + 
s(lags, k=lag, by= tas_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pr_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pet_scaledcovar, bs = 'cs') +
s(tot_shading_t0, k = 8, bs = 'cs') + slope + rock + soil_depth"

#FLM models

growth_mod <- gam(as.formula(base_growth),
      data=data_t1 %>%
        filter(survival_t1 == 1),
      method="GCV.Cp",
      gamma=1.4,
      select = T)


plot_spline_coeff(best_model = growth_mod, tas = T, vital_rate = "Survival")
plot_spline_coeff(best_model = growth_mod, pr = T, vital_rate = "Survival")
plot_spline_coeff(best_model = growth_mod, pet = T, vital_rate = "Survival")
plot_spline_coeff(best_model = growth_mod, shade = T, vital_rate = "Survival")




## -------------------------------------------------------
## Flower probability
## -------------------------------------------------------

# Basemodel
base_flowp <- "flower_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're') + 
s(lags, k=lag, by= tas_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pr_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pet_scaledcovar, bs = 'cs') +
s(tot_shading_t0, k = 8, bs = 'cs') + slope + rock + soil_depth"

#FLM models

flowp_mod <- gam(as.formula(base_flowp),
      data=data_t0 ,
      family = binomial(link = "logit"),
      method="GCV.Cp",
      gamma=1.4,
      select = T)


plot_spline_coeff(best_model = flowp_mod, tas = T, vital_rate = "Survival")
plot_spline_coeff(best_model = flowp_mod, pr = T, vital_rate = "Survival")
plot_spline_coeff(best_model = flowp_mod, pet = T, vital_rate = "Survival")
plot_spline_coeff(best_model = flowp_mod, shade = T, vital_rate = "Survival")


## -------------------------------------------------------
## Abortion probability
## -------------------------------------------------------

# Basemodel
base_seedp <- "seed_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're') + 
s(lags, k=lag, by= tas_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pr_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pet_scaledcovar, bs = 'cs') +
s(tot_shading_t0, k = 8, bs = 'cs') + slope + rock + soil_depth"

#FLM models

seedp_mod <- gam(as.formula(base_seedp),
      data=data_t0 %>%
        filter(flower_p_t0 == 1),
      family = binomial(link = "logit"),
      method="GCV.Cp",
      gamma=1.4,
      select = T)


plot_spline_coeff(best_model = seedp_mod, tas = T, vital_rate = "Survival")
plot_spline_coeff(best_model = seedp_mod, pr = T, vital_rate = "Survival")
plot_spline_coeff(best_model = seedp_mod, pet = T, vital_rate = "Survival")
plot_spline_coeff(best_model = seedp_mod, shade = T, vital_rate = "Survival")



## -------------------------------------------------------
## Seed numbers
## -------------------------------------------------------

# Basemodel
base_seedn <- "est_seed_n_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're') + 
s(lags, k=lag, by= tas_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pr_scaledcovar, bs = 'cs') + s(lags, k=lag, by= pet_scaledcovar, bs = 'cs') +
s(tot_shading_t0, k = 8, bs = 'cs') + slope + rock + soil_depth"

#FLM models

seedn_mod <- gam(as.formula(base_seedn),
      data=data_t0 %>%
        filter(flower_p_t0 == 1 & seed_p_t0 == 1),
      family = Gamma(link = "log"),
      method="GCV.Cp",
      gamma=1.4,
      select = T)


plot_spline_coeff(best_model = seedn_mod, tas = T, vital_rate = "Survival")
plot_spline_coeff(best_model = seedn_mod, pr = T, vital_rate = "Survival")
plot_spline_coeff(best_model = seedn_mod, pet = T, vital_rate = "Survival")
plot_spline_coeff(best_model = seedn_mod, shade = T, vital_rate = "Survival")



## -------------------------------------------------------
## Save results
## -------------------------------------------------------


saveRDS(list(surv = surv_mods[[as.symbol(surv_aic.cv$model[1])]],
             growth = growth_mods[[as.symbol(growth_aic.cv$model[1])]],
             flower_p = flowp_mods[[as.symbol(flowp_aic.cv$model[1])]],
             abort_p = abp_mods[[as.symbol(abp_aic.cv$model[1])]],
             n_seeds = seed_mods[[as.symbol(seed_aic.cv$model[1])]]),
        file = "results/rds/VR_FLM.rds")

saveRDS(list(surv = surv_mods,
             growth = growth_mods,
             flower_p = flowp_mods,
             abort_p = abp_mods,
             n_seeds = seed_mods),
        file = "results/rds/VR_FLM_all.rds")

saveRDS(
  list(
    surv = surv_aic.cv,
    growth = growth_aic.cv,
    flowp = flowp_aic.cv,
    abp = abp_aic.cv,
    seeds = seed_aic.cv
  ),
  file = "results/rds/VR_mod_AIC_CV.rds"
)

## Check if crossvalidation would select a different model from AIC
surv_aic.cv$model[which.min(surv_aic.cv$AIC)]
surv_aic.cv$model[which.min(surv_aic.cv$RMSE)]

growth_aic.cv$model[which.min(growth_aic.cv$AIC)]
growth_aic.cv$model[which.min(growth_aic.cv$RMSE)]

flowp_aic.cv$model[which.min(flowp_aic.cv$AIC)]
flowp_aic.cv$model[which.min(flowp_aic.cv$RMSE)]

abp_aic.cv$model[which.min(abp_aic.cv$AIC)]
abp_aic.cv$model[which.min(abp_aic.cv$RMSE)]

seed_aic.cv$model[which.min(seed_aic.cv$AIC)]
seed_aic.cv$model[which.min(seed_aic.cv$RMSE)]




