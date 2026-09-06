set.seed(42)

### script to run extinction probability analysis
library(tidyverse)
library(glmnet)
library(lme4)

source("R/functions_ibm.R")

if(!dir.exists("results/ibm_validation")) dir.create("results/ibm_validation", recursive = T)

# Projections will be from 2022 to 2100
n_it = 15

## -------------------------------------------------------------------------------------------
## Format CHELSA's time series for IBM
## -------------------------------------------------------------------------------------------

hist_clim <- read.csv("data/CHELSA_data.csv") %>% 
  filter(complete.cases(.))

clim_spring <- hist_clim %>% 
  filter(month >= 3 & month <= 5) %>%
  mutate(year_t0 = year - 1, .keep = "unused") %>%
  group_by(locality, year_t0) %>%
  summarise(pet_spring = mean(pet_scaled, na.rm = T),
            pr_spring = mean(pr_scaled, na.rm = T),
            tas_spring = mean(tas_scaled, na.rm = T)) %>%
  ungroup()
clim_summer <- hist_clim %>% 
  filter(month >= 6 & month <= 8) %>%
  rename(year_t0 = year) %>%
  group_by(locality, year_t0) %>%
  summarise(pet_summer = mean(pet_scaled, na.rm = T),
            pr_summer = mean(pr_scaled, na.rm = T),
            tas_summer = mean(tas_scaled, na.rm = T)) %>%
  ungroup()
clim_dormant <- hist_clim %>% 
  filter(month <= 2 | month >= 9) %>%
  mutate(year_t0 = ifelse(month < 3, year - 1, year)) %>%
  group_by(locality, year_t0) %>%
  summarise(pet_dormant = mean(pet_scaled, na.rm = T),
            pr_dormant = mean(pr_scaled, na.rm = T),
            tas_dormant = mean(tas_scaled, na.rm = T)) %>%
  ungroup() 


hist_clim <- left_join(clim_spring, clim_summer)
hist_clim <- left_join(hist_clim, clim_dormant)


## -------------------------------------------------------------------------------------------
## Other variables for IBM
## -------------------------------------------------------------------------------------------

VR_mods <- list(surv = readRDS("results/rds/seasons_surv.rds"),
                growth = readRDS("results/rds/seasons_growth.rds"),
                flower_p = readRDS("results/rds/seasons_flowp.rds"),
                seedp = readRDS("results/rds/seasons_seedp.rds"),
                seedn = readRDS("results/rds/seasons_seedn.rds")
)

state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")

demo_data <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
  mutate(herb_shading_t0 = scale(herb_shading_t0),
         shrub_shading_t0 = scale(shrub_shading_t0))

params <- list(
  surv_mod = VR_mods$surv,
  
  grow_mod = VR_mods$growth,
  grow_sd = readRDS(file = "results/rds/seasons_growth_sd.rds"),
  
  pflower_mod = VR_mods$flower_p,
  
  pseed_mod = VR_mods$seedp,
  
  nseed_mod = VR_mods$seedn,
  nseed_sd = readRDS(file = "results/rds/seasons_seedn_sd.rds"),
  
  seed_surv1 = 0.45,  ## Probability of seed being viable at the next census 
  seed_surv2 = 0.089,  ## Probability of viable seed surviving first year in seed bank. 
  ## (i.e. produced in t, surv to census t+1 (with % above), not germinated in t+1 and 
  seed_surv3 = 0.663,  ## Probability of viable seed in seedbank (yr 1) surviving 2nd year
  
  germ_mean = mean(state_independent_variables$est_germination_rate$germ),     
  
  sdl_surv_mean = boot::inv.logit(fixef(state_independent_variables$sdl_surv)),
  
  sdl_d_int = lme4::fixef(state_independent_variables$sdl_size_d)[1],
  sdl_d_sd = sd(resid(state_independent_variables$sdl_size_d)),
  
  soil_depth = calc_stats(demo_data, "soil_depth"),
  slope = calc_stats(demo_data, "slope"),
  rock = calc_stats(demo_data, "rock"),
  
  herb_center = attr(demo_data$herb_shading_t0, "scaled:center"),
  herb_scale = attr(demo_data$herb_shading_t0, "scaled:scale"),
  shrub_center = attr(demo_data$shrub_shading_t0, "scaled:center"),
  shrub_scale = attr(demo_data$shrub_shading_t0, "scaled:scale")
)


## Starting population vector
data <- read.csv("data/Dracocephalum_with_vital_rates.csv") 
data <- split(data, data$population)
data <- lapply(data, function(x){x %>% filter(year_t0 == min(year_t0))})

obs_n <- readxl::read_xlsx("data/Population_sizes_validation.xlsx") %>%
  filter(year == 2003 & stage == "total") %>%
  select(-c(stage, `Haknovec 2`)) %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "obs_size") %>%
  mutate(population = population |> 
           recode_values(
             "Císařská rokle" ~ "Cr",
             "Haknovec 1" ~ "Hk",
             "Kodská stěna" ~ "Ks",
             "Radotínské údolí" ~ "Ru"
           )) %>%
  select(-year)

pop_n <- as.list(obs_n$obs_size)
names(pop_n) <- obs_n$population


## -------------------------------------------------------------------------------------------
## Set up environmental values
## -------------------------------------------------------------------------------------------


### Loop through different populations and env_param levels 
localities <- c("Cr", "Hk", "Ks", "Ru")
shrub_shading <- c(2, 4)
herb_shading <- c(2, 4)
sb1_vec <- seq(5000, 20000, by = 2500)
sb2_vec <- seq(2500, 10000, by = 2500)

df_env <- expand.grid(localities = localities, 
                       shrub_shading = shrub_shading, 
                       herb_shading = herb_shading,
                      sb1 = sb1_vec,
                      sb2 = sb2_vec
) 


## -------------------------------------------------------------------------------------------
## IBM function
## -------------------------------------------------------------------------------------------


ibm_validation <- function(i, df_env, params,
                           clim_hist,
                           start_data, start_n, 
                           sdl_n,
                           start_yr = 2003, end_yr = 2018) {
  # gc()
  message(i)
  
  n_it = end_yr - start_yr
  
  # Simulation settings
  locality_i <- df_env$localities[i]
  
  Hshading_i <- df_env$herb_shading[i]
  Sshading_i <- df_env$shrub_shading[i]
  
  clim <- clim_hist %>% filter(locality == locality_i)
  
  # Select pop_vectors
  pop_vec <-  sample((start_data[[toupper(locality_i)]] %>% filter(stage_t0 %in% c("veg","flow")))$ln_stems_t0,
                     start_n[[locality_i]], replace = TRUE)
  
  sdl <- ceiling(nrow(start_data[[locality_i]] %>% filter(stage_t0 == "sdl")) *
                   (start_n[[locality_i]] / nrow(start_data[[locality_i]] %>% filter(stage_t0 %in% c("veg","flow")))))
  
  sb1 = df_env$sb1[i] 
  sb2 = df_env$sb2[i] 
  
  # set up starting shading & slope
  pop_Hshading <- rbinom(n = length(pop_vec), 20, prob = (Hshading_i/20))
  pop_Sshading <- rbinom(n = length(pop_vec), 20, prob = (Sshading_i/20))
  pop_slope <- cadence_rbgamma(n = length(pop_vec), prob = params$slope$prob_non_zero[which(params$slope$population == toupper(locality_i))],
                               scale = params$slope$gamma_scale[which(params$slope$population == toupper(locality_i))], 
                               shape = params$slope$gamma_shape[which(params$slope$population == toupper(locality_i))])
  
  pop_rock <- cadence_rbgamma(n = length(pop_vec), prob = params$rock$prob_non_zero[which(params$rock$population == toupper(locality_i))],
                              scale = params$rock$gamma_scale[which(params$rock$population == toupper(locality_i))], 
                              shape = params$rock$gamma_shape[which(params$rock$population == toupper(locality_i))])
  
  pop_sd <- cadence_rbgamma(n = length(pop_vec), prob = params$soil_depth$prob_non_zero[which(params$soil_depth$population == toupper(locality_i))],
                            scale = params$soil_depth$gamma_scale[which(params$soil_depth$population == toupper(locality_i))], 
                            shape = params$soil_depth$gamma_shape[which(params$soil_depth$population == toupper(locality_i))])
  
  # Set up start population
  indi <- data.frame(
    size = pop_vec,
    Hshading = pop_Hshading,
    Sshading = pop_Sshading,
    slope = pop_slope,
    rock = pop_rock,
    soil_depth = pop_sd
  )
  
  pop_size <- data.frame(
    yr = c(start_yr:end_yr),
    n_plants = c(nrow(indi), rep(NA, n_it)),
    n_sdl = c(sdl, rep(NA, n_it)),
    n_sb1 = c(sb1, rep(NA, n_it)),
    n_sb2 = c(sb2, rep(NA, n_it))
  )
  
  ## IBM
  yr <- start_yr
  
  while (yr < end_yr & nrow(indi) > 0){
    
    ibm_yr <- yearly_loop(yr, clim, locality = locality_i,
                          indi, sdl, sb1, sb2,
                          Hshading = Hshading_i,
                          Sshading = Sshading_i)
    
    yr <- ibm_yr$yr
    indi <- ibm_yr$indi
    sdl <- ibm_yr$sdl
    sb1 <- ibm_yr$sb1
    sb2 <- ibm_yr$sb2
    
    pop_size$n_plants[yr - start_yr] <- nrow(indi)
    pop_size$n_sdl[yr - start_yr] <- sdl
    pop_size$n_sb1[yr - start_yr] <- sb1
    pop_size$n_sb2[yr - start_yr] <- sb2
    
  }
  
  b <- data.frame(
    locality = locality_i,
    sb1_start = df_env$sb1[i],
    sb2_start = df_env$sb2[i],
    Hshading = Hshading_i,
    Sshading = Sshading_i,
    below_ext_size = any(pop_size$n_plants < 10), 
    yr_of_ext = ifelse(any(pop_size$n_plants < 10), min(which(pop_size$n_plants < 10)), NA)
  )
  b$pop_size = list(list(plants = pop_size$n_plants,
                         sdl = pop_size$n_sdl,
                         sb1 = pop_size$n_sb1,
                         sb2 = pop_size$n_sb2))
  
  return(b)
  
}

## -------------------------------------------------------------------------------------------
## Run IBM
## -------------------------------------------------------------------------------------------


gc()
print("start ibm")

# args <- commandArgs(trailingOnly = TRUE)
# taskID = as.integer(args[1])

for(taskID in c(1:nrow(df_env))) {
  results <- lapply(as.list(c(1:50)), function(x) {
    a <- ibm_validation(i = taskID, df_env = df_env,
                   params = params,
                   clim_hist = hist_clim,
                   start_data = data,
                   start_n = pop_n)
    return(a)}) %>% bind_rows()
  
  saveRDS(results, file = paste0("results/ibm_validation/df_", taskID, ".rds"))
}


## -------------------------------------------------------------------------------------------
## Plotting of results
## -------------------------------------------------------------------------------------------
rm(list = ls())
gc()

df <- lapply(list.files("results/ibm_validation/", full.names = T), readRDS) %>%
  bind_rows() %>%
  mutate(locality = as.factor(locality),
         yr_of_ext = yr_of_ext + 2003,
         yr_of_ext = replace_na(yr_of_ext, 2019))

names(df$pop_size) <- paste(df$locality, df$sb1_start, df$sb2_start, c(1:nrow(df)), sep = "_")


ext_pl <- lapply(df$pop_size, function(x) x$plants) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>% tibble::rowid_to_column() %>%
  pivot_longer(num_range(prefix = "", range = c(1:15)), names_to = "time_step", values_to = "pop_size_plant") %>%
  mutate(time_step = as.numeric(time_step),
         year = time_step + 2003)

ext_sdl <- lapply(df$pop_size, function(x) x$sdl) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:15)), names_to = "time_step", values_to = "pop_size_sdl") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2003)

ext_sb1 <- lapply(df$pop_size, function(x) x$sb1) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:15)), names_to = "time_step", values_to = "pop_size_sb1") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2003)

ext_sb2 <- lapply(df$pop_size, function(x) x$sb2) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:(15))), names_to = "time_step", values_to = "pop_size_sb2") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2003)

ext <- ext_pl %>% cbind(.,
                        pop_size_sdl = ext_sdl$pop_size_sdl,
                        pop_size_sb1 = ext_sb1$pop_size_sb1,
                        pop_size_sb2 = ext_sb2$pop_size_sb2) %>%
  mutate(below_ext_size = ifelse(pop_size_plant < 10, T, F))

rm(ext_pl, ext_sdl, ext_sb1, ext_sb2)


ext_summary <- ext %>%
  group_by(locality, year) %>% #, sb1_start, sb2_start) %>%
  summarise(
    mean_pop = mean(pop_size_plant, na.rm = TRUE),
    q1 = quantile(pop_size_plant, 0.25, na.rm = TRUE),
    q3 = quantile(pop_size_plant, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

val_df <- readxl::read_xlsx("data/Population_sizes_validation.xlsx") %>%
  filter(stage == "total") %>%
  dplyr::select(-c(`Haknovec 2`, stage)) %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "obs_size") %>%
  mutate(locality = population |>
           recode_values(
             "Císařská rokle" ~ "Cr",
             "Haknovec 1" ~ "Hk",
             "Kodská stěna" ~ "Ks",
             "Radotínské údolí" ~ "Ru"
           ))

ggplot(ext_summary, aes(x = year)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.3, colour = NA) +
  geom_line(aes(y = mean_pop),
            linewidth = 0.8) +
  geom_line(data = val_df, aes(x = year, y = obs_size, shape = "Observed data"), colour = "red") +
  coord_cartesian(xlim = c(2003, 2018)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  scale_shape_manual(name = NULL, values = c("Observed data" = 18)) +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size")


ggplot() +
  geom_line(data = ext, aes(x = year, y = pop_size_plant, group = rowid)) +
  geom_line(data = val_df, aes(x = year, y = obs_size, shape = "Observed data"), colour = "red") +
  facet_grid(rows = vars(locality))
