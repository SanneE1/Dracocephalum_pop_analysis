library(lme4)
library(dplyr)
library(ggplot2)
library(brms)
library(MuMIn)
library(glmnet)

data_for_modeling = "data/Dracocephalum_with_vital_rates.csv"
CHELSA_data = "data/CHELSA_data.csv"

demo_data <- read.csv(data_for_modeling) %>%
  mutate(population = population,
         rock = scale(rock),
         slope = scale(slope),
         soil_depth = scale(soil_depth),
         herb_shading_t0 = scale(herb_shading_t0),
         shrub_shading_t0 = scale(shrub_shading_t0)) %>% ungroup() %>%
  filter(is.finite(ln_stems_t0) & !(stage_t0 == "sdl")) %>%
  select(-c(contains("seeds_n_"), contains("tm1"), contains("tm2")))



clim_data <- read.csv(CHELSA_data)

clim_spring <- clim_data %>% 
  filter(month > 2 & month < 6) %>%
  mutate(year_t0 = year - 1, .keep = "unused") %>%
  group_by(locality, year_t0) %>%
  summarise(pet_spring = mean(pet_scaled, na.rm = T),
            pr_spring = mean(pr_scaled, na.rm = T),
            tas_spring = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")
clim_summer <- clim_data %>% 
  filter(month > 2 & month < 6) %>%
  rename(year_t0 = year) %>%
  group_by(locality, year_t0) %>%
  summarise(pet_summer = mean(pet_scaled, na.rm = T),
            pr_summer = mean(pr_scaled, na.rm = T),
            tas_summer = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")
clim_dormant <- clim_data %>% 
  filter(month < 3 | month > 8) %>%
  mutate(year_t0 = ifelse(month < 3, year - 1, year)) %>%
  group_by(locality, year_t0) %>%
  summarise(pet_dormant = mean(pet_scaled, na.rm = T),
            pr_dormant = mean(pr_scaled, na.rm = T),
            tas_dormant = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")


df <- left_join(demo_data, clim_dormant)
df <- left_join(df, clim_summer)
df <- left_join(df, clim_spring)


surv_df <- df %>%
  filter(!is.na(survival_t1) & !is.na(ln_stems_t0)) %>%
  filter(complete.cases(population, rock, slope, herb_shading_t0, shrub_shading_t0, 
                        pr_dormant, tas_dormant, pr_summer, tas_summer,
                        pr_spring, tas_spring))

#-------------------------------------------------------------------------------------------
# SURVIVAL
#-------------------------------------------------------------------------------------------

surv_gam <- gam(survival_t1 ~ ln_stems_t0 + population  +
                  s(rock) + s(slope) + s(soil_depth) +
                  s(herb_shading_t0) + s(shrub_shading_t0) +
                  s(pr_dormant) + s(tas_dormant) +
                  s(pr_summer)+ s(tas_summer) +
                  s(pr_spring) + s(tas_spring) +
                  ti(pr_dormant, herb_shading_t0) + ti(pr_dormant, shrub_shading_t0) +
                  ti(pr_summer, herb_shading_t0) + ti(pr_summer, shrub_shading_t0) +
                  ti(pr_spring, herb_shading_t0) + ti(pr_spring, shrub_shading_t0) +
                  ti(tas_dormant, herb_shading_t0) + ti(tas_dormant, shrub_shading_t0) +
                  ti(tas_summer, herb_shading_t0) + ti(tas_summer, shrub_shading_t0) +
                  ti(tas_spring, herb_shading_t0) + ti(tas_spring, shrub_shading_t0),
                data=surv_df ,
                family = binomial(link = "logit"),
                method="GCV.Cp",
                gamma=5,
                select = T,
                na.action = "na.fail")

edf_surv <- summary(surv_gam)$s.table[, "edf"]
p_surv <- summary(surv_gam)$s.table[, "p-value"]

surv_selected <- names(edf_surv)[edf_surv > 1 & p_surv < 0.05]

surv_mod <- gam(as.formula(paste0("survival_t1 ~ ln_stems_t0 + population  +",
                                  paste(surv_selected, collapse = "+"))),
                data=surv_df ,
                family = binomial(link = "logit"),
                method="GCV.Cp",
                gamma=5,
                select = T,
                na.action = "na.fail")


saveRDS(surv_final, file = "results/rds/seasons_gam_surv.rds")
rm(surv_gam, edf_surv, p_surv, surv_selected, surv_mod)
gc()

#-------------------------------------------------------------------------------------------
# GROWTH
#-------------------------------------------------------------------------------------------

growth_gam <- gam(ln_stems_t1 ~ ln_stems_t0 + population  + s(rock) + s(slope) + s(soil_depth) +
                    s(herb_shading_t0) + s(shrub_shading_t0) +
                    s(pr_dormant) + s(tas_dormant) +
                    s(pr_summer)+ s(tas_summer) +
                    s(pr_spring) + s(tas_spring) +
                    ti(pr_dormant, herb_shading_t0) + ti(pr_dormant, shrub_shading_t0) +
                    ti(pr_summer, herb_shading_t0) + ti(pr_summer, shrub_shading_t0) +
                    ti(pr_spring, herb_shading_t0) + ti(pr_spring, shrub_shading_t0) +
                    ti(tas_dormant, herb_shading_t0) + ti(tas_dormant, shrub_shading_t0) +
                    ti(tas_summer, herb_shading_t0) + ti(tas_summer, shrub_shading_t0) +
                    ti(tas_spring, herb_shading_t0) + ti(tas_spring, shrub_shading_t0),
                  data = surv_df %>% filter(survival_t1 == 1),
                  family = gaussian(),
                  method="GCV.Cp",
                  gamma=2,
                  select = T,
                  na.action = "na.fail")

edf_growth <- summary(growth_gam)$s.table[, "edf"]
p_growth <- summary(growth_gam)$s.table[, "p-value"]

growth_selected <- names(edf_growth)[edf_growth > 1 & p_growth < 0.05]

growth_mod <- gam(as.formula(paste0("ln_stems_t1 ~ ln_stems_t0 + population  +",
                                    paste(growth_selected, collapse = "+"))),
                  data = surv_df %>% filter(survival_t1 == 1),
                  family = gaussian(),
                  method="GCV.Cp",
                  gamma=2,
                  select = T,
                  na.action = "na.fail")

saveRDS(growth_mod, file = "results/rds/seasons_gam_growth.rds")
rm(growth_gam, edf_growth, p_growth, growth_selected, growth_mod)
gc()

#-------------------------------------------------------------------------------------------
# FLOWER PROBABILITY
#-------------------------------------------------------------------------------------------

flowp_gam <- gam(flower_p_t0 ~ ln_stems_t0 + population  + s(rock) + s(slope) + s(soil_depth) + 
                   s(herb_shading_t0) + s(shrub_shading_t0) + 
                   s(pr_dormant) + s(tas_dormant) +
                   s(pr_summer)+ s(tas_summer) + 
                   s(pr_spring) + s(tas_spring) +
                   ti(pr_dormant, herb_shading_t0) + ti(pr_dormant, shrub_shading_t0) + 
                   ti(pr_summer, herb_shading_t0) + ti(pr_summer, shrub_shading_t0) + 
                   ti(pr_spring, herb_shading_t0) + ti(pr_spring, shrub_shading_t0) + 
                   ti(tas_dormant, herb_shading_t0) + ti(tas_dormant, shrub_shading_t0) + 
                   ti(tas_summer, herb_shading_t0) + ti(tas_summer, shrub_shading_t0) + 
                   ti(tas_spring, herb_shading_t0) + ti(tas_spring, shrub_shading_t0),
                 data = surv_df %>% filter(!is.na(flower_p_t0)),
                 family = binomial(link = "logit"),
                 method="GCV.Cp",
                 gamma=2,
                 select = T,
                 na.action = "na.fail")

edf_flowp <- summary(flowp_gam)$s.table[, "edf"]
p_flowp <- summary(flowp_gam)$s.table[, "p-value"]

flowp_selected <- names(edf_flowp)[edf_flowp > 1 & p_flowp < 0.05]

flowp_mod <- gam(as.formula(paste0("flower_p_t0 ~ ln_stems_t0 + population  +", 
                                   paste(flowp_selected, collapse = "+"))),
                 data = surv_df %>% filter(!is.na(flower_p_t0)),
                 family = binomial(link = "logit"),
                 method="GCV.Cp",
                 gamma=2,
                 select = T,
                 na.action = "na.fail")

saveRDS(flowp_mod, file = "results/rds/seasons_gam_flowp.rds")
rm(flowp_gam, edf_flowp, p_flowp, flowp_selected, flowp_mod)
gc()

#-------------------------------------------------------------------------------------------
# SEED PROBABILITY
#-------------------------------------------------------------------------------------------

seedp_gam <- gam(seed_p_t0 ~ ln_stems_t0 + population  + s(rock) + s(slope) + s(soil_depth) + 
                   s(herb_shading_t0, k = 4) + s(shrub_shading_t0, k = 4) +
                   s(pr_dormant) + s(tas_dormant) +
                   s(pr_summer)+ s(tas_summer) +
                   s(pr_spring) + s(tas_spring) +
                   ti(pr_dormant, herb_shading_t0) + ti(pr_dormant, shrub_shading_t0) +
                   ti(pr_summer, herb_shading_t0) + ti(pr_summer, shrub_shading_t0) +
                   ti(pr_spring, herb_shading_t0) + ti(pr_spring, shrub_shading_t0) +
                   ti(tas_dormant, herb_shading_t0) + ti(tas_dormant, shrub_shading_t0) +
                   ti(tas_summer, herb_shading_t0) + ti(tas_summer, shrub_shading_t0) +
                   ti(tas_spring, herb_shading_t0) + ti(tas_spring, shrub_shading_t0),
                 data = surv_df %>% filter(flower_p_t0 == 1),
                 family = binomial(link = "logit"),
                 method="GCV.Cp",
                 gamma=2,
                 select = T,
                 na.action = "na.fail")

edf_seedp <- summary(seedp_gam)$s.table[, "edf"]
p_seedp <- summary(seedp_gam)$s.table[, "p-value"]

seedp_selected <- names(edf_seedp)[edf_seedp > 1 & p_seedp < 0.05]

seedp_mod <- gam(as.formula(paste0("seed_p_t0 ~ ln_stems_t0 + population  +", 
                                   paste(seedp_selected, collapse = "+"))),
                 data = surv_df %>% filter(flower_p_t0 == 1),
                 family = binomial(link = "logit"),
                 method="GCV.Cp",
                 gamma=2,
                 select = T,
                 na.action = "na.fail")

saveRDS(seedp_mod, file = "results/rds/seasons_gam_seedp.rds")
rm(seedp_gam, edf_seedp, p_seedp, seedp_selected, seedp_mod)
gc()

#-------------------------------------------------------------------------------------------
# SEED NUMBERS
#-------------------------------------------------------------------------------------------

seedn_gam <- gam(est_seed_n_t0 ~ ln_stems_t0 + population  + s(rock) + s(slope) + s(soil_depth) + 
                   s(herb_shading_t0, k = 4) + s(shrub_shading_t0, k = 4) +
                   s(pr_dormant) + s(tas_dormant) +
                   s(pr_summer)+ s(tas_summer) + 
                   s(pr_spring) + s(tas_spring) +
                   ti(pr_dormant, herb_shading_t0) + ti(pr_dormant, shrub_shading_t0) + 
                   ti(pr_summer, herb_shading_t0) + ti(pr_summer, shrub_shading_t0) + 
                   ti(pr_spring, herb_shading_t0) + ti(pr_spring, shrub_shading_t0) + 
                   ti(tas_dormant, herb_shading_t0) + ti(tas_dormant, shrub_shading_t0) + 
                   ti(tas_summer, herb_shading_t0) + ti(tas_summer, shrub_shading_t0) + 
                   ti(tas_spring, herb_shading_t0) + ti(tas_spring, shrub_shading_t0),
                 data = surv_df %>% filter(seed_p_t0 == 1),
                 family = Gamma(link = 'log'),
                 method="GCV.Cp",
                 gamma=2,
                 select = T,
                 na.action = "na.fail")


edf_seedn <- summary(seedn_gam)$s.table[, "edf"]
p_seedn <- summary(seedn_gam)$s.table[, "p-value"]

seedn_selected <- names(p_seedn)[edf_seedn > 1 & p_seedn < 0.05]

seedn_mod <- gam(as.formula(paste0("est_seed_n_t0 ~ ln_stems_t0 + population  +", 
                                   paste(seedn_selected, collapse = "+"))),
                 data = surv_df %>% filter(seed_p_t0 == 1),
                 family = Gamma(link = 'log'),
                 method="GCV.Cp",
                 gamma=2,
                 select = T,
                 na.action = "na.fail")


saveRDS(seedn_mod, file = "results/rds/seasons_gam_seedn.rds")
rm(seedn_gam, edf_seedn, p_seedn, seedn_selected, seedn_mod)
gc()
