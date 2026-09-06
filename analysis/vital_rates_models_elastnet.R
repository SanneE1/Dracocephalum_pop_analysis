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
  select(-c(contains("seeds_n_"), contains("tm1"), contains("tm2"))) %>%
  arrange(plant_ID, year_t0) %>%
  group_by(plant_ID) %>%
  mutate(flower_p_t1 = lead(flower_p_t0),
         seed_p_t1 = lead(seed_p_t0),
         est_seed_n_t1 = lead(est_seed_n_t0))



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
  filter(complete.cases(population, rock, slope, soil_depth, herb_shading_t0, shrub_shading_t0 + pet_dormant, 
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
                         herb_shading_t0:pet_spring + shrub_shading_t0:pet_spring +
                         ln_stems_t0:slope + 
                         pr_dormant:slope + tas_dormant:slope + pet_dormant:slope + pr_summer:slope + tas_summer:slope + pet_summer:slope + pr_spring:slope + tas_spring:slope + pet_spring:slope +
                         pr_dormant:ln_stems_t0 + tas_dormant:ln_stems_t0 + pet_dormant:ln_stems_t0 + pr_summer:ln_stems_t0 + tas_summer:ln_stems_t0 + pet_summer:ln_stems_t0 + pr_spring:ln_stems_t0 + tas_spring:ln_stems_t0 + pet_spring:ln_stems_t0, 
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
growth_df <- df %>%
  filter(survival_t1 == 1 & !is.na(ln_stems_t0)) %>%
  filter(complete.cases(population, rock, slope, soil_depth, herb_shading_t0, shrub_shading_t0 + pet_dormant, 
                        pr_dormant, tas_dormant, pet_summer, pr_summer, tas_summer,
                        pet_spring, pr_spring, tas_spring))

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
                           herb_shading_t0:pet_spring + shrub_shading_t0:pet_spring +
                           ln_stems_t0:slope + 
                           pr_dormant:slope + tas_dormant:slope + pet_dormant:slope + pr_summer:slope + tas_summer:slope + pet_summer:slope + pr_spring:slope + tas_spring:slope + pet_spring:slope +
                           pr_dormant:ln_stems_t0 + tas_dormant:ln_stems_t0 + pet_dormant:ln_stems_t0 + pr_summer:ln_stems_t0 + tas_summer:ln_stems_t0 + pet_summer:ln_stems_t0 + pr_spring:ln_stems_t0 + tas_spring:ln_stems_t0 + pet_spring:ln_stems_t0, 
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
flowp_df <- df %>%
  filter(!is.na(flower_p_t1) & is.finite(ln_stems_t1)) %>%
  filter(complete.cases(flower_p_t1, ln_stems_t1, population, rock, slope, soil_depth, herb_shading_t1, shrub_shading_t1, pet_dormant, 
                        pr_dormant, tas_dormant, pet_summer, pr_summer, tas_summer,
                        pet_spring, pr_spring, tas_spring))


y_flowp <- flowp_df$flower_p_t1
X_flowp <- model.matrix(flower_p_t1 ~ ln_stems_t1 + population + rock + slope + soil_depth + 
                          herb_shading_t1 + shrub_shading_t1 +
                          pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                          herb_shading_t1:pr_dormant + shrub_shading_t1:pr_dormant +
                          herb_shading_t1:pr_summer + shrub_shading_t1:pr_summer +
                          herb_shading_t1:pr_spring + shrub_shading_t1:pr_spring + 
                          herb_shading_t1:tas_dormant + shrub_shading_t1:tas_dormant +
                          herb_shading_t1:tas_summer + shrub_shading_t1:tas_summer +
                          herb_shading_t1:tas_spring + shrub_shading_t1:tas_spring +
                          herb_shading_t1:pet_dormant + shrub_shading_t1:pet_dormant +
                          herb_shading_t1:pet_summer + shrub_shading_t1:pet_summer +
                          herb_shading_t1:pet_spring + shrub_shading_t1:pet_spring +
                          ln_stems_t1:slope + 
                          pr_dormant:slope + tas_dormant:slope + pet_dormant:slope + pr_summer:slope + tas_summer:slope + pet_summer:slope + pr_spring:slope + tas_spring:slope + pet_spring:slope +
                          pr_dormant:ln_stems_t1 + tas_dormant:ln_stems_t1 + pet_dormant:ln_stems_t1 + pr_summer:ln_stems_t1 + tas_summer:ln_stems_t1 + pet_summer:ln_stems_t1 + pr_spring:ln_stems_t1 + tas_spring:ln_stems_t1 + pet_spring:ln_stems_t1, 
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
seedp_df <- df %>%
  filter(flower_p_t1 == 1 & is.finite(ln_stems_t1)) %>%
  filter(complete.cases(seed_p_t1, ln_stems_t1, population, rock, slope, soil_depth, herb_shading_t1, shrub_shading_t1, pet_dormant, 
                        pr_dormant, tas_dormant, pet_summer, pr_summer, tas_summer,
                        pet_spring, pr_spring, tas_spring))
  
y_seedp <- seedp_df$seed_p_t1
X_seedp <- model.matrix(seed_p_t1 ~ ln_stems_t1 + population + rock + slope + soil_depth + 
                          herb_shading_t1 + shrub_shading_t1 +
                          pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                          herb_shading_t1:pr_dormant + shrub_shading_t1:pr_dormant +
                          herb_shading_t1:pr_summer + shrub_shading_t1:pr_summer +
                          herb_shading_t1:pr_spring + shrub_shading_t1:pr_spring + 
                          herb_shading_t1:tas_dormant + shrub_shading_t1:tas_dormant +
                          herb_shading_t1:tas_summer + shrub_shading_t1:tas_summer +
                          herb_shading_t1:tas_spring + shrub_shading_t1:tas_spring +
                          herb_shading_t1:pet_dormant + shrub_shading_t1:pet_dormant +
                          herb_shading_t1:pet_summer + shrub_shading_t1:pet_summer +
                          herb_shading_t1:pet_spring + shrub_shading_t1:pet_spring +
                          ln_stems_t1:slope + 
                          pr_dormant:slope + tas_dormant:slope + pet_dormant:slope + pr_summer:slope + tas_summer:slope + pet_summer:slope + pr_spring:slope + tas_spring:slope + pet_spring:slope +
                          pr_dormant:ln_stems_t1 + tas_dormant:ln_stems_t1 + pet_dormant:ln_stems_t1 + pr_summer:ln_stems_t1 + tas_summer:ln_stems_t1 + pet_summer:ln_stems_t1 + pr_spring:ln_stems_t1 + tas_spring:ln_stems_t1 + pet_spring:ln_stems_t1, 
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

seedn_df <- df %>%
  filter(seed_p_t1 == 1 & is.finite(ln_stems_t1)) %>%
  filter(complete.cases(est_seed_n_t1, ln_stems_t1, population, rock, slope, soil_depth, herb_shading_t1, shrub_shading_t1, pet_dormant, 
                        pr_dormant, tas_dormant, pet_summer, pr_summer, tas_summer,
                        pet_spring, pr_spring, tas_spring))

y_seedn <- seedn_df$est_seed_n_t1
X_seedn <- model.matrix(est_seed_n_t1 ~ ln_stems_t1 + population + rock + slope + soil_depth + 
                          herb_shading_t1 + shrub_shading_t1 +
                          pr_dormant + tas_dormant + pet_dormant + pr_summer + tas_summer + pet_summer + pr_spring + tas_spring + pet_spring +
                          herb_shading_t1:pr_dormant + shrub_shading_t1:pr_dormant +
                          herb_shading_t1:pr_summer + shrub_shading_t1:pr_summer +
                          herb_shading_t1:pr_spring + shrub_shading_t1:pr_spring + 
                          herb_shading_t1:tas_dormant + shrub_shading_t1:tas_dormant +
                          herb_shading_t1:tas_summer + shrub_shading_t1:tas_summer +
                          herb_shading_t1:tas_spring + shrub_shading_t1:tas_spring +
                          herb_shading_t1:pet_dormant + shrub_shading_t1:pet_dormant +
                          herb_shading_t1:pet_summer + shrub_shading_t1:pet_summer +
                          herb_shading_t1:pet_spring + shrub_shading_t1:pet_spring +
                          ln_stems_t1:slope + 
                          pr_dormant:slope + tas_dormant:slope + pet_dormant:slope + pr_summer:slope + tas_summer:slope + pet_summer:slope + pr_spring:slope + tas_spring:slope + pet_spring:slope +
                          pr_dormant:ln_stems_t1 + tas_dormant:ln_stems_t1 + pet_dormant:ln_stems_t1 + pr_summer:ln_stems_t1 + tas_summer:ln_stems_t1 + pet_summer:ln_stems_t1 + pr_spring:ln_stems_t1 + tas_spring:ln_stems_t1 + pet_spring:ln_stems_t1, 
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




#-------------------------------------------------------------------------------------------
# Plot covariate estimates (requires the X and y matrices and final models to be loaded)
#-------------------------------------------------------------------------------------------
source("R/functions_plotting.R")

boot_surv <- bootstrap_glmnet_coefs(x = X_surv, y = y_surv,
                                     alpha = 0.5,          # match your original fit
                                     lambda = surv_final$lambda,
                                     n_boot = 1000,
                                     family = "binomial")

plot_coef_clim_shading(model = surv_final, boot_coefs = boot_surv,
                       save = T, name = "Survival_clim_shading_coef", plot_title = 'Survival')


boot_growth <- bootstrap_glmnet_coefs(x = X_growth, y = y_growth,
                                    alpha = 0.5,          # match your original fit
                                    lambda = growth_final$lambda,
                                    n_boot = 1000,
                                    family = "gaussian")

plot_coef_clim_shading(model = growth_final, boot_coefs = boot_growth,
                       save = T, name = "Growth_clim_shading_coef", plot_title = 'Growth')


boot_flowp <- bootstrap_glmnet_coefs(x = X_flowp, y = y_flowp,
                                    alpha = 0.5,          # match your original fit
                                    lambda = flowp_final$lambda,
                                    n_boot = 1000,
                                    family = "binomial")

plot_coef_clim_shading(model = flowp_final, boot_coefs = boot_flowp,
                       save = T, name = "Flower_P_clim_shading_coef", plot_title = 'Flower probability')

boot_seedp <- bootstrap_glmnet_coefs(x = X_seedp, y = y_seedp,
                                     alpha = 0.5,          # match your original fit
                                     lambda = seedp_final$lambda,
                                     n_boot = 1000,
                                     family = "binomial")

plot_coef_clim_shading(model = seedp_final, boot_coefs = boot_seedp,
                       save = T, name = "Seed_P_clim_shading_coef", plot_title = 'Seed probability')

boot_seedn <- bootstrap_glmnet_coefs(x = X_seedn, y = y_seedn,
                                     alpha = 0.5,          # match your original fit
                                     lambda = seedn_final$lambda,
                                     n_boot = 1000,
                                     family = Gamma(link = 'log'))

plot_coef_clim_shading(model = seedn_final, boot_coefs = boot_seedn,
                       save = T, name = "Seed_N_clim_shading_coef", plot_title = 'Number of seeds')


View(build_effect_size_table(model = surv_final, boot_coefs = boot_surv) %>%
  arrange(shading_type, variable, period) %>%
  filter(sign_reversal == T & ci_excludes_baseline == T & ci_excludes_zero == T))

View(build_effect_size_table(model = surv_final, boot_coefs = boot_surv) %>%
       arrange(shading_type, variable, period))


View(build_effect_size_table(model = growth_final, boot_coefs = boot_growth) %>%
       arrange(shading_type, variable, period) %>%
       filter(sign_reversal == T & ci_excludes_baseline == T & ci_excludes_zero == T))

View(build_effect_size_table(model = growth_final, boot_coefs = boot_growth) %>%
       arrange(shading_type, variable, period))


View(build_effect_size_table(model = flowp_final, boot_coefs = boot_flowp) %>%
       arrange(shading_type, variable, period) %>%
       filter(sign_reversal == T & ci_excludes_baseline == T & ci_excludes_zero == T))

View(build_effect_size_table(model = flowp_final, boot_coefs = boot_flowp) %>%
       arrange(shading_type, variable, period))


