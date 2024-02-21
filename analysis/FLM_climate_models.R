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



smooth_forms <- list(
  null = "",
  tot = "+ tot_shading_t0",
  slope = "+ slope",
  rock = "+ rock",
  soil = "+ soil_depth",
  tas24 = '+ s(lags, k=lag, by= tas_scaledcovar)',
  tas8 = '+ s(lags, k=lag/3, by= tas_scaledcovar)',
  tas_tot = '+ te(lags, tot_shading_m, k=lag/3, by= tas_scaledcovar)',
  tas_slope = '+ te(lags, slope_m, k=lag/3, by= tas_scaledcovar)',
  tas_rock = '+ te(lags, rock_m, k=lag/3, by= tas_scaledcovar)',
  tas_soil = '+ te(lags, soil_m, k=lag/3, by= tas_scaledcovar)',
  pr24 = '+ s(lags, k=lag, by= pr_scaledcovar)',
  pr8 = '+ s(lags, k=lag/3, by= pr_scaledcovar)',
  pr_tot = '+ te(lags, tot_shading_m, k=lag/3, by= pr_scaledcovar)',
  pr_slope = '+ te(lags, slope_m, k=lag/3, by= pr_scaledcovar)',
  pr_rock = '+ te(lags, rock_m, k=lag/3, by= pr_scaledcovar)',
  pr_soil = '+ te(lags, soil_m, k=lag/3, by= pr_scaledcovar)',
  pet24 = '+ s(lags, k=lag, by= pet_scaledcovar)',
  pet8 = '+ s(lags, k=lag/3, by= pet_scaledcovar)',
  pet_tot = '+ te(lags, tot_shading_m, k=lag/3, by= pet_scaledcovar)',
  pet_slope = '+ te(lags, slope_m, k=lag/3, by= pet_scaledcovar)',
  pet_rock = '+ te(lags, rock_m, k=lag/3, by= pr_scaledcovar)',
  pet_soil = '+ te(lags, soil_m, k=lag/3, by= pr_scaledcovar)'
)

smooth_forms0 <- list(
  null = "",
  tot = "+ tot_shading_t0",
  slope = "+ slope",
  rock = "+ rock",
  soil = "+ soil_depth",
  tas24 = '+ s(lags, k=lag, by= tas_scaledcovar)',
  tas5 = '+ s(lags, k=lag/5, by= tas_scaledcovar)',
  tas_tot = '+ te(lags, tot_shading_m, k=lag/5, by= tas_scaledcovar)',
  tas_slope = '+ te(lags, slope_m, k=lag/5, by= tas_scaledcovar)',
  tas_rock = '+ te(lags, rock_m, k=lag/5, by= tas_scaledcovar)',
  tas_soil = '+ te(lags, soil_m, k=lag/5, by= tas_scaledcovar)',
  pr24 = '+ s(lags, k=lag, by= pr_scaledcovar)',
  pr8 = '+ s(lags, k=lag/5, by= pr_scaledcovar)',
  pr_tot = '+ te(lags, tot_shading_m, k=lag/5, by= pr_scaledcovar)',
  pr_slope = '+ te(lags, slope_m, k=lag/5, by= pr_scaledcovar)',
  pr_rock = '+ te(lags, rock_m, k=lag/5, by= pr_scaledcovar)',
  pr_soil = '+ te(lags, soil_m, k=lag/5, by= pr_scaledcovar)',
  pet24 = '+ s(lags, k=lag, by= pet_scaledcovar)',
  pet8 = '+ s(lags, k=lag/5, by= pet_scaledcovar)',
  pet_tot = '+ te(lags, tot_shading_m, k=lag/5, by= pet_scaledcovar)',
  pet_slope = '+ te(lags, slope_m, k=lag/5, by= pet_scaledcovar)',
  pet_rock = '+ te(lags, rock_m, k=lag/5, by= pr_scaledcovar)',
  pet_soil = '+ te(lags, soil_m, k=lag/5, by= pr_scaledcovar)'
)

VR_FLM_all <- readRDS("results/rds/VR_FLM_all.rds")

surv_mods <- VR_FLM_all$surv
growth_mods <- VR_FLM_all$growth
flowp_mods <- VR_FLM_all$flower_p
abp_mods <- VR_FLM_all$abort_p
seed_mods <- VR_FLM_all$n_seeds

## -------------------------------------------------------
## Survival
## -------------------------------------------------------


# # Basemodel
# base_surv <- "survival_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
# 
# #FLM models
# 
# surv_mods <- lapply(smooth_forms, function(x)
#   
#   gam(formula(paste(base_surv, x)),
#       data=data_t1 ,
#       family = binomial(link = "logit"),
#       method="GCV.Cp", 
#       gamma=1.2)
# )


surv_aic <- bbmle::AICtab(surv_mods, weights = T, base = T) %>%
  as.data.frame %>%
  rownames_to_column

surv_cv <- sapply(surv_mods, gam.crossvalidation)  %>% 
  t %>%
  as.data.frame %>%
  rownames_to_column

surv_aic.cv <- left_join(surv_aic, surv_cv) %>% 
  arrange(dAIC) %>%
  rename(model = rowname)


## -------------------------------------------------------
## Growth
## -------------------------------------------------------

# # Basemodel
# base_growth <- "ln_stems_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
# 
# #FLM models
# 
# growth_mods <- lapply(smooth_forms, function(x)
#   
#   gam(formula(paste(base_growth, x)),
#       data=data_t1 %>%
#         filter(survival_t1 == 1),
#       method="GCV.Cp",gamma=1.2)
# )


growth_aic <- bbmle::AICtab(growth_mods,
                            weights = T, base = T) %>%
  as.data.frame %>% 
  rownames_to_column

growth_cv <- sapply(growth_mods, gam.crossvalidation) %>% 
  t %>%
  as.data.frame %>% 
  rownames_to_column

growth_aic.cv <- left_join(growth_aic, growth_cv) %>% 
  arrange(dAIC) %>%
  rename(model = rowname)



## -------------------------------------------------------
## Flower probability
## -------------------------------------------------------

# # Basemodel
# base_flowp <- "flower_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
# 
# #FLM models
# 
# flowp_mods <- lapply(smooth_forms0, function(x)
#   
#   gam(formula(paste(base_flowp, x)),
#       data=data_t0 ,
#       family = binomial(link = "logit"),
#       method="GCV.Cp", 
#       gamma=1.2)
# )



flowp_aic <- bbmle::AICtab(flowp_mods,
                           weights = T, base = T) %>%
  as.data.frame %>% 
  rownames_to_column

flowp_cv <- sapply(flowp_mods, gam.crossvalidation) %>% 
  t %>%
  as.data.frame %>% 
  rownames_to_column

flowp_aic.cv <- left_join(flowp_aic, flowp_cv) %>% 
  arrange(dAIC) %>%
  rename(model = rowname)



## -------------------------------------------------------
## Abortion probability
## -------------------------------------------------------

# # Basemodel
# base_abp <- "seed_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
# 
# #FLM models
# 
# abp_mods <- lapply(smooth_forms0, function(x)
#   
#   gam(formula(paste(base_abp, x)),
#       data=data_t0 %>% 
#         filter(flower_p_t0 == 1),
#       family = binomial(link = "logit"),
#       method="GCV.Cp",
#       gamma=1.2)
# )

abp_aic <- bbmle::AICtab(abp_mods,
                         weights = T, base = T) %>%
  as.data.frame %>% 
  rownames_to_column

abp_cv <- sapply(abp_mods, gam.crossvalidation) %>% 
  t %>%
  as.data.frame %>% 
  rownames_to_column

abp_aic.cv <- left_join(abp_aic, abp_cv) %>% 
  arrange(dAIC) %>%
  rename(model = rowname)



## -------------------------------------------------------
## Seed numbers
## -------------------------------------------------------

# # Basemodel
# base_seed <- "est_seed_n_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
# 
# #FLM models
# 
# seed_mods <- lapply(smooth_forms0, function(x)
#   
#   gam(formula(paste(base_seed, x)),
#       data=data_t0 %>% 
#         filter(flower_p_t0 == 1 & seed_p_t0 == 1),
#       family = Gamma(link = "log"),
#       method="GCV.Cp",
#       gamma=1.2)
# )

seed_aic <- bbmle::AICtab(seed_mods,
                          weights = T, base = T) %>%
  as.data.frame %>% 
  rownames_to_column

seed_cv <- sapply(seed_mods,gam.crossvalidation) %>% 
  t %>%
  as.data.frame %>% 
  rownames_to_column

seed_aic.cv <- left_join(seed_aic, seed_cv) %>% 
  arrange(dAIC) %>%
  rename(model = rowname)



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




