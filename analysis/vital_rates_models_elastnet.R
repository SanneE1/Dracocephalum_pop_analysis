library(lme4)
library(dplyr)
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
  filter(month >= 3 & month <= 5) %>%
  mutate(year_t0 = year - 1, .keep = "unused") %>%
  group_by(locality, year_t0) %>%
  summarise(pet_spring = mean(pet_scaled, na.rm = T),
            pr_spring = mean(pr_scaled, na.rm = T),
            tas_spring = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")
clim_summer <- clim_data %>% 
  filter(month >= 6 & month <= 8) %>%
  rename(year_t0 = year) %>%
  group_by(locality, year_t0) %>%
  summarise(pet_summer = mean(pet_scaled, na.rm = T),
            pr_summer = mean(pr_scaled, na.rm = T),
            tas_summer = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")
clim_dormant <- clim_data %>% 
  filter(month <= 2 | month >= 9) %>%
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

#-------------------------------------------------------------------------------------------
# SURVIVAL
#-------------------------------------------------------------------------------------------

surv_df <- df %>%
  filter(!is.na(survival_t1) & !is.na(ln_stems_t0)) %>%
  filter(complete.cases(population, rock, slope, herb_shading_t0, shrub_shading_t0 + pet_dormant, 
                        pr_dormant, tas_dormant, pet_summer, pr_summer, tas_summer,
                        pet_spring, pr_spring, tas_spring))


y_surv <- surv_df$survival_t1
X_surv <- model.matrix(survival_t1 ~ ln_stems_t0 + population + rock + slope + soil_depth + 
                         herb_shading_t0 + shrub_shading_t0 +
                         pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                         herb_shading_t0:pr_dormant + shrub_shading_t0:pr_dormant +
                         herb_shading_t0:pr_summer + shrub_shading_t0:pr_summer +
                         herb_shading_t0:pr_spring + shrub_shading_t0:pr_spring + 
                         herb_shading_t0:tas_dormant + shrub_shading_t0:tas_dormant +
                         herb_shading_t0:tas_summer + shrub_shading_t0:tas_summer +
                         herb_shading_t0:tas_spring + shrub_shading_t0:tas_spring +
                         herb_shading_t0:pet_dormant + shrub_shading_t0:pet_dormant +
                         herb_shading_t0:pet_summer + shrub_shading_t0:pet_summer +
                         herb_shading_t0:pet_spring + shrub_shading_t0:pet_spring, 
                       data = surv_df)[, -1]  
cv_surv <- cv.glmnet(X_surv, y_surv, family = "binomial", alpha = 0.5)  # alpha = 0.5 for elastic net

# Plot cross-validation results
plot(cv_surv)

# Get the optimal lambda values
surv_l_min <- cv_surv$lambda.min      # Lambda that minimizes CV error
surv_l_1se <- cv_surv$lambda.1se      # Lambda within 1 SE of minimum


surv_final <- glmnet(X_surv, y_surv, family = "binomial", alpha = 0.5, lambda = surv_l_1se)

coef_surv <- coef(surv_final)
non_zero_surv <- coef_surv[coef_surv[,1] != 0, , drop = FALSE]


saveRDS(surv_final, file = "results/rds/seasons_surv.rds")
gc()

#-------------------------------------------------------------------------------------------
# GROWTH
#-------------------------------------------------------------------------------------------
growth_df <- surv_df %>% filter(survival_t1 == 1)

y_growth <- growth_df$ln_stems_t1
X_growth <- model.matrix(ln_stems_t1 ~ ln_stems_t0 + population + rock + slope + soil_depth + 
                           herb_shading_t0 + shrub_shading_t0 +
                           pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                           herb_shading_t0:pr_dormant + shrub_shading_t0:pr_dormant +
                           herb_shading_t0:pr_summer + shrub_shading_t0:pr_summer +
                           herb_shading_t0:pr_spring + shrub_shading_t0:pr_spring + 
                           herb_shading_t0:tas_dormant + shrub_shading_t0:tas_dormant +
                           herb_shading_t0:tas_summer + shrub_shading_t0:tas_summer +
                           herb_shading_t0:tas_spring + shrub_shading_t0:tas_spring +
                           herb_shading_t0:pet_dormant + shrub_shading_t0:pet_dormant +
                           herb_shading_t0:pet_summer + shrub_shading_t0:pet_summer +
                           herb_shading_t0:pet_spring + shrub_shading_t0:pet_spring, 
                         data = growth_df)[, -1]  
cv_growth <- cv.glmnet(X_growth, y_growth, family = "gaussian", alpha = 0.5)  # alpha = 0.5 for elastic net

# Plot cross-validation results
plot(cv_growth)

# Get the optimal lambda values
growth_l_min <- cv_growth$lambda.min      # Lambda that minimizes CV error
growth_l_1se <- cv_growth$lambda.1se      # Lambda within 1 SE of minimum


growth_final <- glmnet(X_growth, y_growth, family = "gaussian", alpha = 0.5, lambda = growth_l_1se)

coef_growth <- coef(growth_final)
non_zero_growth <- coef_growth[coef_growth[,1] != 0, , drop = FALSE]


saveRDS(growth_final, file = "results/rds/seasons_growth.rds")
gc()

## Get residuals of model for growth prediction in IBM
obs <- growth_df$ln_stems_t1
pred <- predict(growth_final, newx = X_growth)

sd_growth <- sd(obs - pred)
saveRDS(sd_growth, file = "results/rds/seasons_growth_sd.rds")

#-------------------------------------------------------------------------------------------
# FLOWER PROBABILITY
#-------------------------------------------------------------------------------------------
flowp_df <- surv_df %>% filter(!is.na(flower_p_t0))

y_flowp <- flowp_df$flower_p_t0
X_flowp <- model.matrix(flower_p_t0 ~ ln_stems_t0 + population + rock + slope + soil_depth + 
                          herb_shading_t0 + shrub_shading_t0 +
                          pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                          herb_shading_t0:pr_dormant + shrub_shading_t0:pr_dormant +
                          herb_shading_t0:pr_summer + shrub_shading_t0:pr_summer +
                          herb_shading_t0:pr_spring + shrub_shading_t0:pr_spring + 
                          herb_shading_t0:tas_dormant + shrub_shading_t0:tas_dormant +
                          herb_shading_t0:tas_summer + shrub_shading_t0:tas_summer +
                          herb_shading_t0:tas_spring + shrub_shading_t0:tas_spring +
                          herb_shading_t0:pet_dormant + shrub_shading_t0:pet_dormant +
                          herb_shading_t0:pet_summer + shrub_shading_t0:pet_summer +
                          herb_shading_t0:pet_spring + shrub_shading_t0:pet_spring, 
                        data = flowp_df)[, -1]  
cv_flowp <- cv.glmnet(X_flowp, y_flowp, family = "binomial", alpha = 0.5)  # alpha = 0.5 for elastic net

# Plot cross-validation results
plot(cv_flowp)

# Get the optimal lambda values
flowp_l_min <- cv_flowp$lambda.min      # Lambda that minimizes CV error
flowp_l_1se <- cv_flowp$lambda.1se      # Lambda within 1 SE of minimum


flowp_final <- glmnet(X_flowp, y_flowp, family = "binomial", alpha = 0.5, lambda = flowp_l_1se)

coef_flowp <- coef(flowp_final)
non_zero_flowp <- coef_flowp[coef_flowp[,1] != 0, , drop = FALSE]


saveRDS(flowp_final, file = "results/rds/seasons_flowp.rds")
gc()

#-------------------------------------------------------------------------------------------
# SEED PROBABILITY
#-------------------------------------------------------------------------------------------
seedp_df <- surv_df %>% filter(flower_p_t0 == 1)

y_seedp <- seedp_df$seed_p_t0
X_seedp <- model.matrix(seed_p_t0 ~ ln_stems_t0 + population + rock + slope + soil_depth + 
                          herb_shading_t0 + shrub_shading_t0 +
                          pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                          herb_shading_t0:pr_dormant + shrub_shading_t0:pr_dormant +
                          herb_shading_t0:pr_summer + shrub_shading_t0:pr_summer +
                          herb_shading_t0:pr_spring + shrub_shading_t0:pr_spring + 
                          herb_shading_t0:tas_dormant + shrub_shading_t0:tas_dormant +
                          herb_shading_t0:tas_summer + shrub_shading_t0:tas_summer +
                          herb_shading_t0:tas_spring + shrub_shading_t0:tas_spring +
                          herb_shading_t0:pet_dormant + shrub_shading_t0:pet_dormant +
                          herb_shading_t0:pet_summer + shrub_shading_t0:pet_summer +
                          herb_shading_t0:pet_spring + shrub_shading_t0:pet_spring, 
                        data = seedp_df)[, -1]  
cv_seedp <- cv.glmnet(X_seedp, y_seedp, family = "binomial", alpha = 0.5)  # alpha = 0.5 for elastic net

# Plot cross-validation results
plot(cv_seedp)

# Get the optimal lambda values
seedp_l_min <- cv_seedp$lambda.min      # Lambda that minimizes CV error
seedp_l_1se <- cv_seedp$lambda.1se      # Lambda within 1 SE of minimum


seedp_final <- glmnet(X_seedp, y_seedp, family = "binomial", alpha = 0.5, lambda = seedp_l_1se)

coef_seedp <- coef(seedp_final)
non_zero_seedp <- coef_seedp[coef_seedp[,1] != 0, , drop = FALSE]


saveRDS(seedp_final, file = "results/rds/seasons_seedp.rds")
gc()

#-------------------------------------------------------------------------------------------
# SEED NUMBERS
#-------------------------------------------------------------------------------------------

seedn_df <- surv_df %>% filter(seed_p_t0 == 1)

y_seedn <- seedn_df$est_seed_n_t0
X_seedn <- model.matrix(est_seed_n_t0 ~ ln_stems_t0 + population + rock + slope + soil_depth + 
                          herb_shading_t0 + shrub_shading_t0 +
                          pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                          herb_shading_t0:pr_dormant + shrub_shading_t0:pr_dormant +
                          herb_shading_t0:pr_summer + shrub_shading_t0:pr_summer +
                          herb_shading_t0:pr_spring + shrub_shading_t0:pr_spring + 
                          herb_shading_t0:tas_dormant + shrub_shading_t0:tas_dormant +
                          herb_shading_t0:tas_summer + shrub_shading_t0:tas_summer +
                          herb_shading_t0:tas_spring + shrub_shading_t0:tas_spring +
                          herb_shading_t0:pet_dormant + shrub_shading_t0:pet_dormant +
                          herb_shading_t0:pet_summer + shrub_shading_t0:pet_summer +
                          herb_shading_t0:pet_spring + shrub_shading_t0:pet_spring, 
                        data = seedn_df)[, -1]  
cv_seedn <- cv.glmnet(X_seedn, y_seedn, family = Gamma(link = 'log'), alpha = 0.5)  # alpha = 0.5 for elastic net

# Plot cross-validation results
plot(cv_seedn)

# Get the optimal lambda values
seedn_l_min <- cv_seedn$lambda.min      # Lambda that minimizes CV error
seedn_l_1se <- cv_seedn$lambda.1se      # Lambda within 1 SE of minimum


seedn_final <- glmnet(X_seedn, y_seedn, family = Gamma(link = 'log'), alpha = 0.5, lambda = seedn_l_1se)

coef_seedn <- coef(seedn_final)
non_zero_seedn <- coef_seedn[coef_seedn[,1] != 0, , drop = FALSE]


saveRDS(seedn_final, file = "results/rds/seasons_seedn.rds")
gc()


## Get residuals of model for growth prediction in IBM
obs <- seedn_df$est_seed_n_t0
pred <- predict(seedn_final, newx = X_seedn)

sd_seedn <- sd(obs - pred)
saveRDS(sd_seedn, file = "results/rds/seasons_seedn_sd.rds")




