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
  tas = "s(lags, by= tas_scaledcovar, bs = 'cs', k = 24)", 
  pr = "s(lags, by= pr_scaledcovar, bs = 'cs', k = 24)",
  pet = "s(lags, by= pet_scaledcovar, bs = 'cs', k = 24)",
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
clusterExport(cl=cl, c("gam.selection.criteria"))

## -------------------------------------------------------
## Survival
## -------------------------------------------------------

surv_data <- data_t1 %>%
  select(survival_t1, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth, tot_shading_m) %>%
  filter(complete.cases(.))

clusterExport(cl=cl, c("surv_data"))

# Basemodel
base_surv <- "survival_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
surv_combinations <- generate_combinations(formula_parts) %>%
  paste(base_surv, ., sep = "+")

#FLM models
# surv_mods <- parLapply(cl,
#                          surv_combinations,
#                          function(x) gam(as.formula(x),
#                                          data=surv_data ,
#                                          family = binomial(link = "logit"),
#                                          method="GCV.Cp",
#                                          gamma=1.4,
#                                          select = T,
#                                          na.action = "na.fail")
# )
# saveRDS(surv_mods, file = "results/rds/surv_flm_mods.rds")
surv_mods <- readRDS("results/rds/surv_flm_mods.rds")

surv_select <- parLapply(cl,
                          surv_mods,
                          function(x) gam.selection.criteria(mod = x, 
                                                          response_column =  "survival_t1")) %>%
  bind_rows()


View(surv_select)

surv_combinations[c(125,110,74,120,127)]

# plot(mgcViz::getViz(surv_mods[[125]]))
# plot(mgcViz::getViz(surv_mods[[110]]))
plot(mgcViz::getViz(surv_mods[[74]]))
plot(mgcViz::getViz(surv_mods[[120]]))
# plot(mgcViz::getViz(surv_mods[[127]]))


surv_cross <- parLapply(cl, surv_mods[c(74,120)], 
                        function(x) gam.selection.criteria(mod = x, 
                                        response_column =  "survival_t1", 
                                        do.cross.validation = T, subset_length = 10)) %>%
  bind_rows()
saveRDS(surv_cross, file = "results/rds/surv_flm_cross.rds")

View(surv_cross)
summary(surv_mods[[74]])

saveRDS(surv_mods[[74]], file = "results/rds/surv_flm.rds")


## -------------------------------------------------------
## Growth
## -------------------------------------------------------

growth_data <- data_t1 %>%
  filter(survival_t1 == 1) %>%
  select(ln_stems_t1, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth, tot_shading_m) %>%
  filter(complete.cases(.))

clusterExport(cl=cl, c("growth_data"))

# Basemodel
base_growth <- "ln_stems_t1 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
growth_combinations <- generate_combinations(formula_parts) %>%
  paste(base_growth, ., sep = "+")

# #FLM models
# growth_mods <- parLapply(cl,
#                          growth_combinations,
#                       function(x) gam(as.formula(x),
#                                       data=growth_data,
#                                       method="GCV.Cp",
#                                       gamma=1.4,
#                                       select = T,
#                                       na.action = "na.fail")
#   )
# saveRDS(growth_mods, file = "results/rds/growth_flm_mods.rds")
growth_mods <- readRDS("results/rds/growth_flm_mods.rds")

growth_select <- parLapply(cl,
                          growth_mods,
                          function(x) gam.selection.criteria(mod = x, 
                                                          response_column =  "ln_stems_t1")
                          ) %>% bind_rows()

View(growth_select)

growth_combinations[c(122,104,127,123,67,101)]

# plot(mgcViz::getViz(growth_mods[[122]]))
plot(mgcViz::getViz(growth_mods[[104]]))
# plot(mgcViz::getViz(growth_mods[[127]]))
# plot(mgcViz::getViz(growth_mods[[123]]))
plot(mgcViz::getViz(growth_mods[[67]]))
# plot(mgcViz::getViz(growth_mods[[101]]))

growth_cross <- parLapply(cl,
                           growth_mods[c(104, 67)],
                           function(x) gam.selection.criteria(mod = x,
                                                              response_column =  "ln_stems_t1",
                                                              do.cross.validation = T, subset_length = 10)
) %>% bind_rows()
saveRDS(growth_cross, file = "results/rds/growth_flm_cross.rds")

View(growth_cross)

summary(growth_mods[[104]])
summary(growth_mods[[67]])

saveRDS(growth_mods[[67]], file = "results/rds/growth_flm.rds")


## -------------------------------------------------------
## Flower probability
## -------------------------------------------------------

flowp_data <- data_t0 %>% 
  select(flower_p_t0, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth, tot_shading_m) %>%
  filter(complete.cases(.))

clusterExport(cl=cl, c("flowp_data"))

# Basemodel
base_flowp <- "flower_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
flowp_combinations <- generate_combinations(formula_parts) %>% 
  paste(base_flowp, ., sep = "+")


#FLM models
# flowp_mods <- parLapply(cl,
#                          flowp_combinations,
#                          function(x) gam(as.formula(x),
#                                          data=flowp_data ,
#                                          family = binomial(link = "logit"),
#                                          method="GCV.Cp",
#                                          gamma=1.4,
#                                          select = T,
#                                          na.action = "na.fail")
# )
# saveRDS(flowp_mods, file = "results/rds/flowp_flm_mods.rds")
flowp_mods <- readRDS("results/rds/flowp_flm_mods.rds")


flowp_select <- parLapply(cl,
                         flowp_mods,
                          function(x) gam.selection.criteria(mod = x, 
                                                          response_column =  "flower_p_t0")) %>%
  bind_rows()


View(flowp_select)

flowp_combinations[c(123,127,117,126,124,125,112,102,105,120,115,87,77)]

plot(mgcViz::getViz(flowp_mods[[123]]))
# plot(mgcViz::getViz(flowp_mods[[127]]))
plot(mgcViz::getViz(flowp_mods[[117]]))
plot(mgcViz::getViz(flowp_mods[[126]]))
# plot(mgcViz::getViz(flowp_mods[[124]]))
# plot(mgcViz::getViz(flowp_mods[[125]]))
plot(mgcViz::getViz(flowp_mods[[112]]))
# plot(mgcViz::getViz(flowp_mods[[105]]))
# plot(mgcViz::getViz(flowp_mods[[120]]))
# plot(mgcViz::getViz(flowp_mods[[115]]))
plot(mgcViz::getViz(flowp_mods[[87]]))
plot(mgcViz::getViz(flowp_mods[[77]]))


flowp_cross <- parLapply(cl,
                          flowp_mods[c(123,117,126,112,87,77)],
                          function(x) gam.selection.criteria(mod = x, 
                                                             response_column =  "flower_p_t0", 
                                                             do.cross.validation = T, subset_length = 10)) %>%
  bind_rows()
saveRDS(flowp_cross, file = "results/rds/flowp_flm_cross.rds")

View(flowp_cross)
summary(flowp_mods[[123]])
summary(flowp_mods[[112]])
summary(flowp_mods[[77]])

saveRDS(flowp_mods[[123]], file = "results/rds/flowp_flm.rds")


## -------------------------------------------------------
## Seed probability
## -------------------------------------------------------

seedp_data <- data_t0 %>%
  filter(flower_p_t0 == 1) %>%
  select(seed_p_t0, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth) %>%
  filter(complete.cases(.))

clusterExport(cl=cl, c("seedp_data"))

# Basemodel
base_seedp <- "seed_p_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"
seedp_combinations <- generate_combinations(formula_parts) %>% 
  paste(base_seedp, ., sep = "+")

#FLM models
# seedp_mods <- parLapply(cl, 
#                         seedp_combinations, 
#                         function(x) gam(as.formula(x),
#                                         data=seedp_data,
#                                         family = binomial(link = "logit"),
#                                         method="GCV.Cp",
#                                         gamma=1.4,
#                                         select = T,
#                                         na.action = "na.fail")
# )
# saveRDS(seedp_mods, file = "results/rds/seedp_flm_mods.rds")
seedp_mods <- readRDS("results/rds/seedp_flm_mods.rds")

seedp_select <- parLapply(cl,
                         seedp_mods,
                         function(x) gam.selection.criteria(mod = x, 
                                                         response_column =  "seed_p_t0")) %>%
  bind_rows()

View(seedp_select)


seedp_combinations[c(126,75,100,120,79,77,66,36)]
# plot(mgcViz::getViz(seedp_mods[[126]]))
# plot(mgcViz::getViz(seedp_mods[[75]]))
plot(mgcViz::getViz(seedp_mods[[100]]))
# plot(mgcViz::getViz(seedp_mods[[120]]))
# plot(mgcViz::getViz(seedp_mods[[79]]))
# plot(mgcViz::getViz(seedp_mods[[77]]))
plot(mgcViz::getViz(seedp_mods[[66]]))
# plot(mgcViz::getViz(seedp_mods[[36]]))
plot(mgcViz::getViz(seedp_mods[[49]]))

seedp_cross <- parLapply(cl,
                          seedp_mods[c(100,66,49)],
                          function(x) gam.selection.criteria(mod = x, 
                                                             response_column =  "seed_p_t0", 
                                                             do.cross.validation = T, subset_length = 10)) %>%
  bind_rows()

saveRDS(seedp_cross, file = "results/rds/seedp_flm_cross.rds")

View(seedp_cross)

summary(seedp_mods[[49]])

saveRDS(seedp_mods[[49]], file = "results/rds/seedp_flm.rds")

## -------------------------------------------------------
## Seed numbers
## -------------------------------------------------------

seedn_data <- data_t0 %>%
  filter(flower_p_t0 == 1 & seed_p_t0 == 1) %>%
  select(est_seed_n_t0, ln_stems_t0, population, year_t0, lags, contains("scaledcovar"), 
         tot_shading_t0, slope, rock, soil_depth) %>%
  filter(complete.cases(.))

clusterExport(cl=cl, c("seedn_data"))

# Basemodel
base_seedn <- "est_seed_n_t0 ~ ln_stems_t0 + population + s(year_t0, bs = 're')"

seedn_combinations <- generate_combinations(formula_parts) %>% 
  paste(base_seedn, ., sep = "+")

#FLM models
# seedn_mods <- parLapply(cl, 
#                         seedn_combinations, 
#                         function(x) gam(as.formula(x),
#                                         data=seedn_data,
#                                         family = Gamma(link = "log"),
#                                         method="GCV.Cp",
#                                         gamma=1.4,
#                                         select = T,
#                                         na.action = "na.fail")
# )
# saveRDS(seedn_mods, file = "results/rds/seedn_flm_mods.rds")
seedn_mods <- readRDS("results/rds/seedn_flm_mods.rds")

seedn_select <- parLapply(cl,
                         seedn_mods,
                         function(x) gam.selection.criteria(mod = x, 
                                                         response_column =  "est_seed_n_t0")) %>%
  bind_rows()

View(seedn_select)

seedn_combinations[c(118,91,92,50,106,124,119,95,96,56)]
# plot(mgcViz::getViz(seedn_mods[[118]]))
plot(mgcViz::getViz(seedn_mods[[91]]))
# plot(mgcViz::getViz(seedn_mods[[92]]))
plot(mgcViz::getViz(seedn_mods[[50]]))
plot(mgcViz::getViz(seedn_mods[[106]]))
# plot(mgcViz::getViz(seedn_mods[[124]]))
# plot(mgcViz::getViz(seedn_mods[[119]]))
plot(mgcViz::getViz(seedn_mods[[95]]))
# plot(mgcViz::getViz(seedn_mods[[96]]))
plot(mgcViz::getViz(seedn_mods[[56]]))


seedn_cross <- parLapply(cl,
                          seedn_mods[c(91,50,106,95,56)],
                          function(x) gam.selection.criteria(mod = x, 
                                                             response_column =  "est_seed_n_t0", 
                                                             do.cross.validation = T)) %>%
  bind_rows()
saveRDS(seedn_cross, file = "results/rds/seedn_flm_cross.rds")

View(seedn_cross)

summary(seedn_mods[[91]])

saveRDS(seedn_mods[[91]], file = "results/rds/seedn_flm.rds")

## -------------------------------------------------------
## Save results
## -------------------------------------------------------

stopCluster(cl)

saveRDS(list(surv = readRDS(file = "results/rds/surv_flm.rds"),
             growth = readRDS(file = "results/rds/growth_flm.rds"),
             flower_p = readRDS(file = "results/rds/flowp_flm.rds"),
             seedp = readRDS(file = "results/rds/seedp_flm.rds"),
             seedn = readRDS(file = "results/rds/seedn_flm.rds")),
        file = "results/rds/VR_FLM.rds")


