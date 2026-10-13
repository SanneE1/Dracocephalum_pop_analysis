set.seed(42)

### script to run the effect of different seedbank survival on the IBM projections
library(tidyverse)
library(glmnet)
library(lme4)

source("R/functions_ibm.R")

# Projections will be from 2019 to 2100
start_yr = 2018
end_yr = 2100
n_it = end_yr - start_yr

if(!dir.exists(file.path("results", "ibm_seedbank_survival"))) dir.create(file.path("results", "ibm_seedbank_survival"), recursive = T)

## -------------------------------------------------------------------------------------------
## Format CHELSA's time series for IBM
## -------------------------------------------------------------------------------------------

fut_clim <- read.csv("data/CHELSA_future_ts_formatted.csv") %>%
  filter(complete.cases(.))

clim_spring <- fut_clim %>%
  filter(month >= 3 & month <= 5) %>%
  mutate(year_t0 = year - 1, .keep = "unused") %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_spring = mean(pet_scaled, na.rm = T),
            pr_spring = mean(pr_scaled, na.rm = T),
            tas_spring = mean(tas_scaled, na.rm = T)) %>%
  ungroup()
clim_summer <- fut_clim %>%
  filter(month >= 6 & month <= 8) %>%
  rename(year_t0 = year) %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_summer = mean(pet_scaled, na.rm = T),
            pr_summer = mean(pr_scaled, na.rm = T),
            tas_summer = mean(tas_scaled, na.rm = T)) %>%
  ungroup()
clim_dormant <- fut_clim %>%
  filter(month <= 2 | month >= 9) %>%
  mutate(year_t0 = ifelse(month < 3, year - 1, year)) %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_dormant = mean(pet_scaled, na.rm = T),
            pr_dormant = mean(pr_scaled, na.rm = T),
            tas_dormant = mean(tas_scaled, na.rm = T)) %>%
  ungroup()


fut_clim <- left_join(clim_spring, clim_summer)
fut_clim <- left_join(fut_clim, clim_dormant)


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

  seed_viability = 0.45,  ## Probability of seed produced by the plants, being viable
  seed_surv1 = 0.089,     ## Probability of viable seed surviving first year in seed bank.
  ## (i.e. produced in t, surv to census t+1 (with % above), not germinated in t+1 and
  seed_surv2 = 0.663,  ## Probability of viable seed in seedbank (yr 1) surviving 2nd year

  sdl_est_mean = mean(state_independent_variables$est_germination_rate$germ),

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
data <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
  filter(year_t0 == start_yr)
data <- split(data, data$population)


obs_n <- readxl::read_xlsx("data/Population_sizes_validation.xlsx") %>%
  filter(year == start_yr & stage == "total") %>%
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


### Loop through different populations and seedbank survival levels
localities <- c("Cr", "Ks", "Ru")
shrub_shading <- c(2)
herb_shading <- c(2)

## Multiplier applied to both seedbank survival probabilities (seed_surv1 & seed_surv2),
## 1 = estimated survival (as in IBM_projections.R). Survival is capped at 1.
sb_mult <- c(0.5, 0.75, 1, 1.25, 1.5)

## Starting seedbank sizes (as in ibm_ext_p)
sb1_start <- 10000
sb2_start <- 5000

n_rep <- 20

model <- c("ACCESS1-3", "CESM1-BGC", "CMCC-CM", "MIROC5")
scenario <- c("rcp45", "rcp85")


df_env <- expand.grid(localities = localities,
                      shrub_shading = shrub_shading,
                      herb_shading = herb_shading,
                      scenario = scenario,
                      model = model,
                      sb_mult = sb_mult
) %>% rowid_to_column()

## -------------------------------------------------------------------------------------------
## IBM functions
## -------------------------------------------------------------------------------------------

## Same as yearly_loop() in functions_ibm.R, but takes params as an argument (instead of
## using the global params) so seedbank survival can differ between runs
yearly_loop_sbs <- function(yr, clim, locality,
                            indi, sdl, sb1, sb2,
                            Hshading, Sshading, params){

  # retrieve this year's climate
  clim_yr <- clim %>% filter(year_t0 == yr)
  clim_fec <- clim %>% filter(year_t0 == (yr - 1))

  # calculate population size
  Npop <- nrow(indi)

  # Plants ------------------------------------------------------

  # generate binomial random number for survival
  surv <- rbinom(n = Npop, prob = mod_pred(indi = indi, model = params$surv_mod,
                                           params = params, locality = locality,
                                           clim = clim_yr),
                 size = 1)


  # estimate size for surviving individuals
  z1 <- rep(NA, Npop)
  if (sum(surv) > 0) {
    z1[which(surv == 1)] <- rnorm(n = length(which(surv == 1)),
                                  mean = mod_pred(indi = indi[which(surv == 1),], model = params$grow_mod,
                                                  params = params, locality = locality,
                                                  clim = clim_yr),
                                  sd = params$grow_sd)
  }
  indi_z1 <- indi %>%
    mutate(surv = surv,
           z1 = z1)

  # calculate total number of seeds produced
  flowp <- rbinom(n = Npop, prob = mod_pred(indi = indi, model = params$pflower_mod,
                                            params = params, locality = locality,
                                            clim = clim_fec),
                  size = 1)

  seedp <- rep(NA, Npop)
  if (length(which(flowp == 1)) > 0) {
    seedp[which(flowp == 1)] <- rbinom(n = length(which(flowp == 1)),
                                       prob = mod_pred(indi = indi[which(flowp == 1),],
                                                       model = params$pseed_mod,
                                                       params = params, locality = locality,
                                                       clim = clim_fec),
                                       size = 1)
  }

  nseeds <- rep(NA, Npop)
  if(length(which(seedp == 1)) > 0 ) {
    seed_mean <- mod_pred(indi = indi[which(seedp == 1),], model = params$nseed_mod,
                          params = params, locality = locality,
                          clim = clim_fec)

    seed_gamma_shape1 <- (seed_mean^2)/(params$nseed_sd^2)
    seed_gamma_shape2 <- (seed_mean)/(params$nseed_sd^2)
    nseeds[which(seedp == 1)] <- round(
      rgamma(n = length(which(seedp == 1)), shape = seed_gamma_shape1, scale = seed_gamma_shape2),
      0)

    ## cap individual seed production to 515
    nseeds[which(nseeds > 515)] <- 515
  }

  tot_seeds <- sum(nseeds, na.rm = T)


  # Discrete stages  ------------------------------------------------------

  # survival of discrete stages
  seeds <- rbinom(1, prob = params$seed_viability, size = tot_seeds)
  sb1 <- rbinom(1, prob = params$seed_surv1, size = sb1)
  sb2 <- rbinom(1, prob = params$seed_surv2, size = sb2)
  sdl_to_plant <- rbinom(1, prob = params$sdl_surv_mean, size = sdl)

  # discrete stage transitions

  seed_to_sdl <- rbinom(1, prob = params$sdl_est_mean, size = seeds)
  seed_to_sb1 <- seeds - seed_to_sdl

  sb1_to_sdl <- rbinom(1, prob = params$sdl_est_mean, size = sb1)
  sb1_to_sb2 <- sb1 - sb1_to_sdl

  sb2_to_sdl <- rbinom(1, prob = params$sdl_est_mean, size = sb2)

  sb1 = seed_to_sb1
  sb2 = sb1_to_sb2
  sdl = seed_to_sdl + sb1_to_sdl + sb2_to_sdl

  # recruitment into continuous state
  new_indi <- data.frame(
    size = rnorm(sdl_to_plant, mean = params$sdl_d_int, sd = params$sdl_d_sd),
    Hshading = rpois(sdl_to_plant, Hshading),
    Sshading = rpois(sdl_to_plant, Sshading),
    slope = cadence_rbgamma(n = sdl_to_plant, prob = params$slope$prob_non_zero[which(params$slope$population == toupper(locality))],
                            scale = params$slope$gamma_scale[which(params$slope$population == toupper(locality))],
                            shape = params$slope$gamma_shape[which(params$slope$population == toupper(locality))]),
    rock = cadence_rbgamma(n = sdl_to_plant, prob = params$rock$prob_non_zero[which(params$rock$population == toupper(locality))],
                           scale = params$rock$gamma_scale[which(params$rock$population == toupper(locality))],
                           shape = params$rock$gamma_shape[which(params$rock$population == toupper(locality))]),
    soil_depth = cadence_rbgamma(n = sdl_to_plant, prob = params$soil_depth$prob_non_zero[which(params$soil_depth$population == toupper(locality))],
                                 scale = params$soil_depth$gamma_scale[which(params$soil_depth$population == toupper(locality))],
                                 shape = params$soil_depth$gamma_shape[which(params$soil_depth$population == toupper(locality))])
  )

  # Set next steps individuals dataframe  ------------------------------------------------------
  indi <- bind_rows(indi_z1 %>% dplyr::select(-size) %>%
                      filter(!is.na(z1)) %>% rename(size = z1),
                    new_indi)


  return(list(indi = indi,
              sdl = sdl,
              sb1 = sb1,
              sb2 = sb2,
              yr = yr + 1))

}


## Same as ibm_ext_p() in functions_ibm.R, but scales seedbank survival with df_env$sb_mult
ibm_sbs <- function(i, df_env, params,
                    fut_clim,
                    start_data, start_n,
                    start_yr, end_yr,
                    sb1_start, sb2_start) {
  # gc()
  message(i)

  n_it = end_yr - start_yr

  # Simulation settings
  scenario_i <- df_env$scenario[i]
  model_i <- df_env$model[i]
  locality_i <- df_env$localities[i]

  Hshading_i <- df_env$herb_shading[i]
  Sshading_i <- df_env$shrub_shading[i]

  sb_mult_i <- df_env$sb_mult[i]

  # scale seedbank survival (probabilities can't exceed 1)
  params$seed_surv1 <- min(params$seed_surv1 * sb_mult_i, 1)
  params$seed_surv2 <- min(params$seed_surv2 * sb_mult_i, 1)

  clim <- fut_clim %>% filter(scenario == scenario_i & model == model_i & locality == locality_i)


  # Select pop_vectors
  pop_vec <-  sample((start_data[[toupper(locality_i)]] %>% filter(stage_t0 %in% c("veg","flow")))$ln_stems_t0,
                     start_n[[locality_i]], replace = TRUE)

  sdl <- ceiling(nrow(start_data[[locality_i]] %>% filter(stage_t0 == "sdl")) *
                   (start_n[[locality_i]] / nrow(start_data[[locality_i]] %>% filter(stage_t0 %in% c("veg","flow")))))

  sb1 = sb1_start
  sb2 = sb2_start

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

    ibm_yr <- yearly_loop_sbs(yr, clim, locality = locality_i,
                              indi, sdl, sb1, sb2,
                              Hshading = Hshading_i,
                              Sshading = Sshading_i,
                              params = params)

    yr <- ibm_yr$yr
    indi <- ibm_yr$indi
    sdl <- ibm_yr$sdl
    sb1 <- ibm_yr$sb1
    sb2 <- ibm_yr$sb2

    pop_size$n_plants[yr - start_yr + 1] <- nrow(indi)
    pop_size$n_sdl[yr - start_yr + 1] <- sdl
    pop_size$n_sb1[yr - start_yr + 1] <- sb1
    pop_size$n_sb2[yr - start_yr + 1] <- sb2

  }

  b <- data.frame(
    locality = locality_i,
    model = model_i,
    scenario = scenario_i,
    sb_mult = sb_mult_i,
    Hshading = Hshading_i,
    Sshading = Sshading_i,
    below_ext_size = any(pop_size$n_plants < 10, na.rm = T),
    yr_of_ext = ifelse(any(pop_size$n_plants < 10, na.rm = T), min(which(pop_size$n_plants < 10)), NA)
  )
  b$pop_size = list(list(plants = pop_size$n_plants,
                         sdl = pop_size$n_sdl,
                         sb1 = pop_size$n_sb1,
                         sb2 = pop_size$n_sb2))

  return(b)

}

## -------------------------------------------------------------------------------------------
## IBM run
## -------------------------------------------------------------------------------------------

gc()
print("start ibm")

for(taskID in c(1:nrow(df_env))) {
  results <- lapply(as.list(c(1:n_rep)), function(x) {
    a <- ibm_sbs(i = taskID, df_env = df_env,
                 params = params,
                 fut_clim = fut_clim,
                 start_data = data, start_n = pop_n,
                 start_yr = start_yr, end_yr = end_yr,
                 sb1_start = sb1_start, sb2_start = sb2_start)
    return(a)}) %>% bind_rows()

  saveRDS(results, file = paste0("results/ibm_seedbank_survival/df_", taskID, ".rds"))
}

## -------------------------------------------------------------------------------------------
## Format IBM results
## -------------------------------------------------------------------------------------------

gc()

df <- lapply(list.files("results/ibm_seedbank_survival/", full.names = T), readRDS) %>%
  bind_rows() %>%
  mutate(locality = as.factor(locality),
         sb_mult = factor(sb_mult, levels = .env$sb_mult, labels = paste0(.env$sb_mult * 100, "%"))) %>%
  rowid_to_column("run")

## One row per run and year. Runs that went extinct stop being simulated (NA), these are
## set to 0 plants. A run is counted as extinct from the first year with < 10 plants
ext <- df %>%
  mutate(year = list(start_yr:end_yr),
         plants = lapply(pop_size, function(x) x$plants),
         sb1 = lapply(pop_size, function(x) x$sb1),
         sb2 = lapply(pop_size, function(x) x$sb2)) %>%
  select(-pop_size) %>%
  unnest(c(year, plants, sb1, sb2)) %>%
  mutate(across(c(plants, sb1, sb2), ~ replace_na(.x, 0))) %>%
  group_by(run) %>%
  mutate(below_ext_size = cumsum(plants < 10) > 0) %>%
  ungroup()

ext_summary <- ext %>%
  group_by(locality, scenario, year, sb_mult) %>%
  summarise(
    median_pop = median(plants),
    q1 = quantile(plants, 0.25),
    q3 = quantile(plants, 0.75),
    median_sb = median(sb1 + sb2),
    ext_prob = mean(below_ext_size),
    .groups = "drop"
  )

## -------------------------------------------------------------------------------------------
## Plots
## -------------------------------------------------------------------------------------------

## Each non-estimated seedbank survival level gets its own column, with the estimated
## survival (100%) repeated in every column as the reference
ref_sb <- "100%"
alt_levels <- setdiff(levels(ext_summary$sb_mult), ref_sb)

ext_alt <- ext_summary %>% filter(sb_mult != ref_sb) %>%
  mutate(sb_mult = droplevels(sb_mult))
ext_ref <- ext_summary %>% filter(sb_mult == ref_sb) %>%
  select(-sb_mult) %>%
  cross_join(data.frame(sb_mult = factor(alt_levels, levels = alt_levels)))

## below the estimated survival = orange, above = blue
sb_cols <- setNames(c("#D55E00", "#E69F00", "#56B4E9", "#0072B2"), alt_levels)

## Projected population sizes
ggplot() +
  geom_ribbon(data = ext_ref, aes(x = year, ymin = q1, ymax = q3), fill = "grey40", alpha = 0.2) +
  geom_line(data = ext_ref, aes(x = year, y = median_pop), colour = "grey20") +
  geom_ribbon(data = ext_alt, aes(x = year, ymin = q1, ymax = q3, fill = sb_mult), alpha = 0.25) +
  geom_line(data = ext_alt, aes(x = year, y = median_pop, colour = sb_mult)) +
  facet_grid(rows = vars(locality, scenario), cols = vars(sb_mult), scales = "free_y") +
  scale_colour_manual(values = sb_cols) + scale_fill_manual(values = sb_cols) +
  theme(legend.position = "none") +
  ylab("Population size") + ylim(c(0, NA)) +
  ggtitle("Projected population sizes", 
          subtitle = "Seedbank survival as % of estimated survival (grey = 100%)")

ggsave("results/population_trend_seedbank_survival.png", width = 10, height = 12)


## Projected extinction probability
ggplot() +
  geom_line(data = ext_ref, aes(x = year, y = ext_prob), colour = "grey20") +
  geom_line(data = ext_alt, aes(x = year, y = ext_prob, colour = sb_mult)) +
  facet_grid(rows = vars(locality, scenario), cols = vars(sb_mult)) +
  scale_colour_manual(values = sb_cols) +
  theme(legend.position = "none") +
  ylab("Extinction probability (< 10 plants)") + ylim(c(0, 1)) +
  ggtitle("Projected extinction probability", 
          subtitle = "Seedbank survival as % of estimated survival (grey = 100%)")

ggsave("results/extinction_probability_seedbank_survival.png", width = 10, height = 12)


## Total seedbank size
ggplot() +
  geom_line(data = ext_ref, aes(x = year, y = median_sb + 1), colour = "grey20") +
  geom_line(data = ext_alt, aes(x = year, y = median_sb + 1, colour = sb_mult)) +
  scale_y_log10() +
  facet_grid(rows = vars(locality, scenario), cols = vars(sb_mult)) +
  scale_colour_manual(values = sb_cols) +
  theme(legend.position = "none") +
  ylab("Seedbank size (+1, log scale)") +
  ggtitle("Projected seedbank sizes", 
          subtitle = "Seedbank survival as % of estimated survival (grey = 100%)")

ggsave("results/seedbank_size_seedbank_survival.png", width = 10, height = 12)
