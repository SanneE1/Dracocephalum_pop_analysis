### FLM climate model for Growth

data_for_modeling = "data/Dracocephalum_with_vital_rates.csv"
CHELSA_data = "data/CHELSA_data.csv"
lag = 24  # number of months prior to census to include

source("R/functions_GAM.R")

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
  tas24 = '+ s(lags, k=lag, by= tas_scaledcovar)',
  tas8 = '+ s(lags, k=lag/3, by= tas_scaledcovar)',
  tas_tot = '+ te(lags, tot_shading_m, k=lag/5, by= tas_scaledcovar)',
  tas_slope = '+ te(lags, slope_m, k=lag/3, by= tas_scaledcovar)',
  tas_rock = '+ te(lags, rock_m, k=lag/3, by= tas_scaledcovar)',
  tas_soil = '+ te(lags, soil_m, k=lag/3, by= tas_scaledcovar)',
  pr24 = '+ s(lags, k=lag, by= pr_scaledcovar)',
  pr8 = '+ s(lags, k=lag/3, by= pr_scaledcovar)',
  pr_tot = '+ te(lags, tot_shading_m, k=lag/5, by= pr_scaledcovar)',
  pr_slope = '+ te(lags, slope_m, k=lag/3, by= pr_scaledcovar)',
  pr_rock = '+ te(lags, rock_m, k=lag/3, by= pr_scaledcovar)',
  pr_soil = '+ te(lags, soil_m, k=lag/3, by= pr_scaledcovar)',
  pet24 = '+ s(lags, k=lag, by= pet_scaledcovar)',
  pet8 = '+ s(lags, k=lag/3, by= pet_scaledcovar)',
  pet_tot = '+ te(lags, tot_shading_m, k=lag/5, by= pet_scaledcovar)',
  pet_slope = '+ te(lags, slope_m, k=lag/3, by= pet_scaledcovar)',
  pet_rock = '+ te(lags, rock_m, k=lag/3, by= pr_scaledcovar)',
  pet_soil = '+ te(lags, soil_m, k=lag/3, by= pr_scaledcovar)'
)

## -------------------------------------------------------
## Survival
## -------------------------------------------------------


# Basemodel
surv <- glmer(survival_t1 ~ ln_stems_t0 + population + soil_depth + rock + slope + 
                tot_shading_t0 + (1|year_t0),
              data = na.omit(data_t1[c("survival_t1", "ln_stems_t0", "year_t0", "population", 
                                       "soil_depth", "rock", "slope", 
                                       "tot_shading_t0")]), family = "binomial")

surv_drop <- drop1(surv, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.05, "yes", "no")
  ) %>% as_tibble()

base_surv <- surv_drop %>% filter(include == "yes" | dropped_covariate == "ln_stems_t0") %>%
  pull(dropped_covariate) %>% paste(., collapse = "+") %>%
  paste0("survival_t1 ~", ., '+ s(year_t0, bs="re")')

#FLM models

surv_mods <- lapply(smooth_forms, function(x)
  
  gam(formula(paste(base_surv, x)),
      data=data_t1 ,
      family = binomial(link = "logit"),
      method="GCV.Cp", 
      gamma=1.2)
)

surv_aic <- bbmle::AICtab(surv_mods,
                          weights = T, base = T, sort = T)


summary(surv_mods[[as.symbol(attributes(surv_aic)$row.names[1])]])   ## This eval(as.symbol) is maybe a bit much but will let this analysis run in case input data or something changes the resulting best model

if(attributes(surv_aic)$row.names[1] == "pet_tot") {
  plot_spline_coeff(best_model = surv_mods[[as.symbol(attributes(surv_aic)$row.names[1])]],
                    lag = lag,
                    pet = T, shade = T,
                    vital_rate = "surv",
                    save_plot = T
  )
} else {
  warning("Different survival model than expected with the lowest AIC - Plotting of spline skipped")
}

## -------------------------------------------------------
## Growth
## -------------------------------------------------------

# Basemodel
growth <- lmer(ln_stems_t1 ~ ln_stems_t0 + population + soil_depth + rock + slope + 
                 tot_shading_t0 + (1|year_t0),
               data = na.omit(data_t1[c("ln_stems_t1", "ln_stems_t0", "year_t0", "population", 
                                        "soil_depth", "rock", "slope", 
                                        "tot_shading_t0")]) %>% filter(is.finite(ln_stems_t0) & is.finite(ln_stems_t1))
)

growth_drop <- drop1(growth, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.05, "yes", "no")
  ) %>% as_tibble()

base_growth <- growth_drop %>% filter(include == "yes" | dropped_covariate == "ln_stems_t0") %>%
  pull(dropped_covariate) %>% paste(., collapse = "+") %>%
  paste0("ln_stems_t1 ~", ., '+ s(year_t0, bs="re")')


#FLM models

growth_mods <- lapply(smooth_forms, function(x)
  
  gam(formula(paste(base_growth, x)),
      data=data_t1 %>%
        filter(survival_t1 == 1),
      method="GCV.Cp",gamma=1.2)
)


growth_aic <- bbmle::AICtab(growth_mods,
                            weights = T, base = T, sort = T)


summary(growth_mods[[as.symbol(attributes(growth_aic)$row.names[1])]])   ## This eval(as.symbol) is maybe a bit much but will let this analysis run in case input data or something changes the resulting best model

if(attributes(growth_aic)$row.names[1] == "pet_tot") {
  ## Taking the model with the lowest AIC that doesn't use PET  
  plot_spline_coeff(best_model = growth_mods[[as.symbol(attributes(growth_aic)$row.names[1])]],
                    lag = lag,
                    pet = T, shade = T,
                    vital_rate = "growth",
                    save_plot = T
  )
} else {
  warning("Different growth model than expected with the lowest AIC - Plotting of spline skipped")
}


## -------------------------------------------------------
## Flower probability
## -------------------------------------------------------

# Basemodel
flowp <- glmer(flower_p_t0 ~ ln_stems_t0 + population + soil_depth + rock + slope + 
                 tot_shading_t0 + (1|year_t0),
               data = na.omit(data_t0[c("flower_p_t0", "ln_stems_t0", "year_t0", "population", 
                                        "soil_depth", "rock", "slope", 
                                        "tot_shading_t0")]), 
               family = "binomial")

flowp_drop <- drop1(flowp, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.05, "yes", "no")
  ) %>% as_tibble()

base_flowp <- flowp_drop %>% filter(include == "yes" | dropped_covariate == "ln_stems_t0") %>%
  pull(dropped_covariate) %>% paste(., collapse = "+") %>%
  paste0("flower_p_t0 ~", ., '+ s(year_t0, bs="re")')

#FLM models

flowp_mods <- lapply(smooth_forms0, function(x)
  
  gam(formula(paste(base_flowp, x)),
      data=data_t0 ,
      family = binomial(link = "logit"),
      method="GCV.Cp", 
      gamma=1.2)
)



flowp_aic <- bbmle::AICtab(flowp_mods,
                           weights = T, base = T, sort = T)


summary(flowp_mods[[as.symbol(attributes(flowp_aic)$row.names[1])]])   ## This eval(as.symbol) is maybe a bit much but will let this analysis run in case input data or something changes the resulting best model


if(attributes(flowp_aic)$row.names[1] == "pet_tot") {
  plot_spline_coeff(best_model = flowp_mods[[as.symbol(attributes(flowp_aic)$row.names[1])]],
                    lag = lag,
                    pet = T, shade = T,
                    vital_rate = "flowp",
                    save_plot = T
  )
} else {
  warning("Different flowp model than expected with the lowest AIC - Plotting of spline skipped")
}

## -------------------------------------------------------
## Abortion probability
## -------------------------------------------------------

# Basemodel
abp <- glmer(seed_p_t0 ~ ln_stems_t0 + population + soil_depth + rock + slope + 
               tot_shading_t0 + (1|year_t0),
             data = na.omit(data_t0[c("seed_p_t0", "ln_stems_t0", "year_t0", "population", 
                                      "soil_depth", "rock", "slope", 
                                      "tot_shading_t0")]), 
             family = "binomial")

abp_drop <- drop1(abp, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.05, "yes", "no")
  ) %>% as_tibble()

base_abp <- abp_drop %>% filter(include == "yes" | dropped_covariate == "ln_stems_t0" ) %>%
  pull(dropped_covariate) %>% paste(., collapse = "+") %>%
  paste0("seed_p_t0 ~", ., '+ s(year_t0, bs="re")')

#FLM models

abp_mods <- lapply(smooth_forms0, function(x)
  
  gam(formula(paste(base_abp, x)),
      data=data_t0 %>% 
        filter(flower_p_t0 == 1),
      family = binomial(link = "logit"),
      method="GCV.Cp",
      gamma=1.2)
)

abp_aic <- bbmle::AICtab(abp_mods,
                         weights = T, base = T, sort = T)


summary(abp_mods[[as.symbol(attributes(abp_aic)$row.names[1])]])   


if(attributes(abp_aic)$row.names[1] == "pr_tot") {
  plot_spline_coeff(best_model = abp_mods[[as.symbol(attributes(abp_aic)$row.names[1])]],
                    lag = lag,
                    pr = T, shade = T,
                    vital_rate = "abortion prob.",
                    save_plot = T
  )
} else {
  warning("Different abortion probability model than expected with the lowest AIC - Plotting of spline skipped")
}

## -------------------------------------------------------
## Seed numbers
## -------------------------------------------------------

# Basemodel
seed <- glmer(est_seed_n_t0 ~ ln_stems_t0 + population + soil_depth + rock + slope + 
                tot_shading_t0 + (1|year_t0),
              data = data_t0 %>%
                filter(flower_p_t0 == 1 & seed_p_t0 == 1), 
              family = Gamma(link = "log"))

seed_drop <- drop1(seed, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.05, "yes", "no")
  ) %>% as_tibble()

base_seed <- seed_drop %>% filter(include == "yes" | dropped_covariate == "ln_stems_t0" ) %>%
  pull(dropped_covariate) %>% paste(., collapse = "+") %>%
  paste0("est_seed_n_t0 ~", ., '+ s(year_t0, bs="re")')

#FLM models

seed_mods <- lapply(smooth_forms0, function(x)
  
  gam(formula(paste(base_seed, x)),
      data=data_t0 %>% 
        filter(flower_p_t0 == 1 & seed_p_t0 == 1),
      family = Gamma(link = "log"),
      method="GCV.Cp",
      gamma=1.2)
)

seed_aic <- bbmle::AICtab(seed_mods,
                          weights = T, base = T, sort = T)


summary(seed_mods[[as.symbol(attributes(seed_aic)$row.names[1])]])   


if(attributes(seed_aic)$row.names[1] == "pet_tot") {
  plot_spline_coeff(best_model = seed_mods[[as.symbol(attributes(seed_aic)$row.names[1])]],
                    lag = lag,
                    pet = T, shade = T,
                    vital_rate = "seed production",
                    save_plot = T
  )
} else {
  warning("Different seed production model than expected with the lowest AIC - Plotting of spline skipped")
}



saveRDS(list(surv = surv_mods[[as.symbol(attributes(surv_aic)$row.names[1])]],
             growth = growth_mods[[as.symbol(attributes(growth_aic)$row.names[1])]],
             flower_p = flowp_mods[[as.symbol(attributes(flowp_aic)$row.names[1])]],
             abort_p = abp_mods[[as.symbol(attributes(abp_aic)$row.names[1])]],
             n_seeds = seed_mods[[as.symbol(attributes(seed_aic)$row.names[1])]]),
        file = "results/VR_FLM.rds")


rm(list = ls())
