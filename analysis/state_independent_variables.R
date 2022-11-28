## Script to estimate germination rate

data <- read.csv('data/Dracocephalum_with_vital_rates.csv') %>%
  rowwise() %>%
  mutate(tot_shading_t0 = herb_shading_t0 + shrub_shading_t0)

n_sdl_t3 <- data %>% filter(stage_t0 == "sdl") %>% 
  group_by(population, year_t0) %>% summarise(n_sdl_t3 = n()) %>%
  mutate(year_t0 = year_t0 - 3)
n_seeds_t0 <- data %>% group_by(population, year_t0) %>% 
  summarise(n_seeds_t0 = sum(est_seed_n_t0, na.rm = T))
n_seeds_t1 <- data %>% group_by(population, year_t0) %>% 
  summarise(n_seeds_t1 = sum(est_seed_n_t0, na.rm = T)) %>% mutate(year_t0 = year_t0 -1)
n_seeds_t2 <- data %>% group_by(population, year_t0) %>% 
  summarise(n_seeds_t2 = sum(est_seed_n_t0, na.rm = T)) %>% mutate(year_t0 = year_t0 - 2)


df <- plyr::join_all(list(n_seeds_t0, n_seeds_t1, n_seeds_t2, n_sdl_t3)) %>% ungroup() %>% 
  filter(complete.cases(.)) %>% 
  mutate(ln_n_sdl_t3 = log(n_sdl_t3),
         ln_seeds_t0 = log(n_seeds_t0),
         ln_seeds_t1 = log(n_seeds_t1),
         ln_seeds_t2 = log(n_seeds_t2)
  ) 


### sdl_t2 normal distribution - tried a poisson distribution, but this had a higher n_eff
germ_model <- rethinking::ulam(
  alist(
    ln_n_sdl_t3 ~ dnorm(sdl_mean, sdl_sd),
    sdl_mean <- (n_seeds_t2 * 0.45 * germ) + (n_seeds_t1 * (0.45 * (1-germ)) * (0.089 * germ)) + ((n_seeds_t1 * (0.45 * (1-germ)) * (0.089 * germ)) * (0.663 * germ)), 
    germ ~ dbeta(alpha, kappa),
    alpha ~ dnorm(2, 0.25),
    kappa ~ dnorm(10, 1),
    sdl_sd ~ dunif(0,100)
  ), data = df ,
  chains = 8, cores = 4,
  log_lik = T)


sdl_data <- data %>% filter(stage_t0 == "sdl")

## herb & shrub shading not significant for seedling survival and CI of estimates includes 0
# sdl_surv <- glmer(survival_t1 ~ population + herb_shading_t0 + shrub_shading_t0 + (1|year_t0), family = "binomial", data = sdl_data)
# confint.merMod(sdl_surv)
sdl_surv <- glmer(survival_t1 ~ population + soil_depth + rock + slope + 
                    tot_shading_t0 + (1|year_t0), 
                  family = "binomial", 
                  data = na.omit(sdl_data[c('survival_t1', 'population', 'soil_depth',
                                          'rock', 'slope', 'tot_shading_t0', 'year_t0')])
)

surv_drop <- drop1(sdl_surv, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.5, "yes", "no")
  ) %>% as_tibble()


base_surv <- surv_drop %>% filter(include == "yes" | dropped_covariate == "ln_stems_t0") %>%
  pull(dropped_covariate) %>% paste(., collapse = "+") %>%
  paste0("survival_t1 ~", ., '(1|year_t0)')

sdl_surv <- glmer(formula(base_surv), 
                  family = "binomial", data = sdl_data)


## herb & shrub shading not significant for seedling size distribution and CI of estimates includes 0
# sdl_size_d <- lmer(ln_stems_t1 ~ population + herb_shading_t0 + shrub_shading_t0 + (1|year_t0), data = sdl_data %>% filter(survival_t1 == 1))
# confint.merMod(sdl_size_d)
sdl_size_d <- lmer(ln_stems_t1 ~ population + soil_depth + rock + slope + 
                     tot_shading_t0 + (1|year_t0), 
                   data = na.omit(sdl_data[c('survival_t1', 'ln_stems_t1', 'population', 
                                             'soil_depth', 'rock', 'slope', 
                                             'tot_shading_t0', 'year_t0')]) %>% 
                                    filter(survival_t1 == 1)
                   )

sdld_drop <- drop1(sdl_size_d, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.5, "yes", "no")
  ) %>% as_tibble()


base_sdld <- sdld_drop %>% filter(include == "yes" | dropped_covariate == "ln_stems_t0") %>%
  pull(dropped_covariate) %>% paste(., collapse = "+") %>%
  paste0("ln_stems_t1 ~", ., '(1|year_t0)')

sdl_size_d <- lmer(formula(base_sdld), 
                   data = sdl_data %>% filter(survival_t1 == 1))


SI_mods <- list(est_germination_rate = rethinking::extract.samples(germ_model, n = 1e5),
                sdl_surv = sdl_surv,
                sdl_size_d = sdl_size_d)

saveRDS(SI_mods, file = "results/rds/state_independent_VR.rds")

rm(list = ls())
