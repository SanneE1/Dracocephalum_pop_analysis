## Script to estimate germination rate
estimate_non_state_dependent_variables <- function(data_for_modeling) {

data <- read.csv(data_for_modeling)
  
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
    sdl_mean <- (n_seeds_t2 * 0.578 * germ) + (n_seeds_t1 * (0.578 * (1-germ)) * (0.154 * germ)) + ((n_seeds_t1 * (0.578 * (1-germ)) * (0.154 * germ)) * (0.667 * germ)), 
    germ ~ dbeta(alpha, kappa),
    alpha ~ dnorm(2, 0.25),
    kappa ~ dnorm(10, 1),
    sdl_sd ~ dunif(0,100)
  ), data = df ,
  chains = 8, cores = 4,
  log_lik = T)


sdl_data <- data %>% filter(stage_t0 == "sdl") %>%
  mutate(n_stems_t1 = n_fl_stems_t1 + n_veg_stems_t1)

## herb & shrub shading not significant for seedling survival, but herb_shading's CI does not span 0, so could be included
sdl_surv <- glmer(survival_t1 ~ population + herb_shading_t0 + (1|year_t0), family = "binomial", data = sdl_data)

## herb shading's CI does not spann 0, shrub shading's CI does include 0
sdl_size_d <- lmer(n_stems_t1 ~ population + herb_shading_t0 + (1|year_t0), data = sdl_data %>% filter(survival_t1 == 1))



return(list(est_germination_rate = rethinking::extract.samples(germ_model, n = 1e5),
            sdl_surv = sdl_surv,
            sdl_size_d = sdl_size_d))

}

