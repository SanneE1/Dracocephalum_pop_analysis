# Script to run deterministic, locality/year specific IPM

data_for_modeling = "data/Dracocephalum_with_vital_rates.csv"
state_independent_variables <- readRDS("results/state_independent_VR.rds")

demo_data <- read.csv(data_for_modeling) %>%
  mutate(population = factor(population)) %>%
  filter(is.finite(ln_stems_t0) & stage_t0 != "sdl") 

sdl_data <- read.csv(data_for_modeling) %>% 
  filter(stage_t0 == "sdl") 

#----------------------------------------------------------------
# Vital rate models
#----------------------------------------------------------------

surv <- glmer(survival_t1 ~ ln_stems_t0 + (population|year_t0), 
              family = "binomial", data = demo_data)

growth <- lmer(ln_stems_t1 ~ ln_stems_t0 + (population|year_t0), 
               data = demo_data %>% 
                 filter(is.finite(ln_stems_t1)))

pflow <- glmer(flower_p_t0 ~ ln_stems_t0 + (population|year_t0), 
               family = "binomial", data = demo_data)

pab <- glmer(seed_p_t0 ~ ln_stems_t0 + (population|year_t0), 
             family = "binomial", data = demo_data)

nseed <- glmer(est_seed_n_t0 ~ ln_stems_t0 + (population|year_t0), 
               family = Gamma(link = "log"), 
               data = demo_data %>% 
                 filter(est_seed_n_t0 > 0))

sdl_surv <- glmer(survival_t1 ~ (population|year_t0), 
                  family = "binomial", data = sdl_data)

sdl_size_d <- lmer(ln_stems_t1 ~  (population|year_t0), 
                   data = sdl_data %>% filter(survival_t1 == 1))


#----------------------------------------------------------------
# Set up variables for IPM
#----------------------------------------------------------------
yrs <- paste0("r_", 2003:2020)

surv_params <- as.list(c(
  fixef(surv), 
  ranef(surv)$year_t0$`(Intercept)`,
  rep(0,18),
  ranef(surv)$year_t0$populationHK,
  ranef(surv)$year_t0$populationKS,
  ranef(surv)$year_t0$populationRU
) )
names(surv_params) <- c("s_int", "s_stems", 
                        paste0("s_int_", yrs),
                        paste0("s_int_", yrs, "_CR"),
                        paste0("s_int_", yrs, "_HK"),
                        paste0("s_int_", yrs, "_KS"),
                        paste0("s_int_", yrs, "_RU"))

growth_params <- as.list(c(
  fixef(growth), 
  ranef(growth)$year_t0$`(Intercept)`,
  rep(0,18),
  ranef(growth)$year_t0$populationHK,
  ranef(growth)$year_t0$populationKS,
  ranef(growth)$year_t0$populationRU,
  sd(resid(growth))
) )
names(growth_params) <- c("g_int", "g_stems", 
                          paste0("g_int_", yrs),
                          paste0("g_int_", yrs, "_CR"),
                          paste0("g_int_", yrs, "_HK"),
                          paste0("g_int_", yrs, "_KS"),
                          paste0("g_int_", yrs, "_RU"),
                          "grow_sd")

pflow_params <- as.list(c(     ## this model also contains 2021
  fixef(pflow), 
  ranef(pflow)$year_t0$`(Intercept)`[1:18],
  rep(0,18),
  ranef(pflow)$year_t0$populationHK[1:18],
  ranef(pflow)$year_t0$populationKS[1:18],
  ranef(pflow)$year_t0$populationRU[1:18]
) )
names(pflow_params) <- c("fp_int", "fp_stems", 
                         paste0("fp_int_", yrs),
                         paste0("fp_int_", yrs, "_CR"),
                         paste0("fp_int_", yrs, "_HK"),
                         paste0("fp_int_", yrs, "_KS"),
                         paste0("fp_int_", yrs, "_RU"))


pab_params <- as.list(c(     ## this model only contains data starting from 2010 as well as 2021
  fixef(pab), 
  rep(0,7), ranef(pab)$year_t0$`(Intercept)`[1:11],
  rep(0,18),
  rep(0,7), ranef(pab)$year_t0$populationHK[1:11],
  rep(0,7), ranef(pab)$year_t0$populationKS[1:11],
  rep(0,7), ranef(pab)$year_t0$populationRU[1:11]
) )
names(pab_params) <- c("pab_int", "pab_stems", 
                       paste0("pab_int_", yrs),
                       paste0("pab_int_", yrs, "_CR"),
                       paste0("pab_int_", yrs, "_HK"),
                       paste0("pab_int_", yrs, "_KS"),
                       paste0("pab_int_", yrs, "_RU"))

nseed_params <- as.list(c(     ## this model only contains data starting from 2010 as well as 2021
  fixef(nseed), 
  rep(0,7), ranef(nseed)$year_t0$`(Intercept)`[1:11],
  rep(0,18),
  rep(0,7), ranef(nseed)$year_t0$populationHK[1:11],
  rep(0,7), ranef(nseed)$year_t0$populationKS[1:11],
  rep(0,7), ranef(nseed)$year_t0$populationRU[1:11]
) )
names(nseed_params) <- c("nseed_int", "nseed_stems", 
                         paste0("nseed_int_", yrs),
                         paste0("nseed_int_", yrs, "_CR"),
                         paste0("nseed_int_", yrs, "_HK"),
                         paste0("nseed_int_", yrs, "_KS"),
                         paste0("nseed_int_", yrs, "_RU"))

sdl_s_params <- as.list(c(     
  fixef(sdl_surv), 
  ranef(sdl_surv)$year_t0$`(Intercept)`,
  rep(0,18),
  ranef(sdl_surv)$year_t0$populationHK,
  ranef(sdl_surv)$year_t0$populationKS,
  ranef(sdl_surv)$year_t0$populationRU
) )
names(sdl_s_params) <- c("sdl_s_int", 
                         paste0("sdl_s_int_", yrs),
                         paste0("sdl_s_int_", yrs, "_CR"),
                         paste0("sdl_s_int_", yrs, "_HK"),
                         paste0("sdl_s_int_", yrs, "_KS"),
                         paste0("sdl_s_int_", yrs, "_RU"))

sdl_d_params <- as.list(c(     
  fixef(sdl_size_d), 
  ranef(sdl_size_d)$year_t0$`(Intercept)`[1:15], 0 , ranef(sdl_size_d)$year_t0$`(Intercept)`[16:17],
  rep(0,18),
  ranef(sdl_size_d)$year_t0$populationHK[1:15], 0 , ranef(sdl_size_d)$year_t0$populationHK[16:17],
  ranef(sdl_size_d)$year_t0$populationKS[1:15], 0 , ranef(sdl_size_d)$year_t0$populationKS[16:17],
  ranef(sdl_size_d)$year_t0$populationRU[1:15], 0 , ranef(sdl_size_d)$year_t0$populationRU[16:17],
  sd(resid(sdl_size_d))
) )
names(sdl_d_params) <- c("sdl_d_int", 
                         paste0("sdl_d_int_", yrs),
                         paste0("sdl_d_int_", yrs, "_CR"),
                         paste0("sdl_d_int_", yrs, "_HK"),
                         paste0("sdl_d_int_", yrs, "_KS"),
                         paste0("sdl_d_int_", yrs, "_RU"),
                         "sdl_d_sd")


other_params <- list(
  
  seed_surv1 = 0.57826,  
  seed_surv2 = 0.15374,  
  seed_surv3 = 0.66658,  
  
  germ_mean = mean(state_independent_variables$est_germination_rate$germ)    
  
)

params <- c(surv_params, growth_params, pflow_params, pab_params, nseed_params,
            sdl_s_params, sdl_d_params, other_params
)



L <- min(demo_data$ln_stems_t0, na.rm = T)
U <- max(demo_data$ln_stems_t0, na.rm = T) * 1.2
n = 100

#----------------------------------------------------------------
# Year/Locality specific IPM
#----------------------------------------------------------------

det_ipm <- init_ipm("general", "di", "det") %>%
  define_kernel(
    name = "P_yr_loc",
    family = "CC",
    formula = s_yr_loc * g_yr_loc * d_stems,
    
    s_yr_loc = plogis(s_lin_yr_loc),
    s_lin_yr_loc = s_int + s_stems * stems_1 + s_int_r_yr + s_int_r_yr_loc,
    g_yr_loc = dnorm(stems_2, g_mu_yr_loc, grow_sd),
    g_mu_yr_loc = g_int + g_stems * stems_1 + g_int_r_yr + g_int_r_yr_loc,
    
    data_list = params,
    states = list(c("stems")),
    
    uses_par_sets = TRUE,
    par_set_indices = list(yr = c(2003:2020),
                           loc = c("CR", "HK", "KS", "RU")),
    
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "g_yr_loc")
  ) %>%
  define_kernel(
    name = "F_to_SDL_yr_loc",
    family = "CD",
    formula = fp_yr_loc * pabort_yr_loc * n_seeds_yr_loc * surv1_seed * germ * d_stems,
    
    fp_yr_loc = plogis(fp_lin_yr_loc),
    fp_lin_yr_loc = fp_int + fp_stems * stems_1 + fp_int_r_yr + fp_int_r_yr_loc,
    pabort_yr_loc = plogis(ab_lin_yr_loc),
    ab_lin_yr_loc = pab_int + pab_stems * stems_1 + pab_int_r_yr + pab_int_r_yr_loc,
    n_seeds_yr_loc = exp(n_seeds_lin_yr_loc),
    n_seeds_lin_yr_loc = nseed_int + nseed_stems * stems_1 + nseed_int_r_yr + nseed_int_r_yr_loc,
    germ = germ_mean,
    surv1_seed = seed_surv1,
    
    data_list = params,
    states = list(c("stems", "sdl")),
    
    uses_par_sets = TRUE,
    par_set_indices = list(yr = c(2003:2020),
                           loc = c("CR", "HK", "KS", "RU")),
    
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name = "F_to_SB1_yr_loc",
    family = "CD",
    formula = fp_yr_loc * pabort_yr_loc * n_seeds_yr_loc * surv1_seed * (1-germ) * d_stems,
    
    fp_yr_loc = plogis(fp_lin_yr_loc),
    fp_lin_yr_loc = fp_int + fp_stems * stems_1 + fp_int_r_yr + fp_int_r_yr_loc,
    pabort_yr_loc = plogis(ab_lin_yr_loc),
    ab_lin_yr_loc = pab_int + pab_stems * stems_1 + pab_int_r_yr + pab_int_r_yr_loc,
    n_seeds_yr_loc = exp(n_seeds_lin_yr_loc),
    n_seeds_lin_yr_loc = nseed_int + nseed_stems * stems_1 + nseed_int_r_yr + nseed_int_r_yr_loc,
    germ = germ_mean,
    surv1_seed = seed_surv1,
    
    data_list = params,
    states = list(c("stems", "sb1")),
    
    uses_par_sets = TRUE,
    par_set_indices = list(yr = c(2003:2020),
                           loc = c("CR", "HK", "KS", "RU")),
    
    evict_cor = FALSE
  )%>%
  define_kernel(
    name = "SB1_to_SB2",
    family = "DD",
    formula = (1-germ) * surv2_seed,
    
    germ = germ_mean,
    surv2_seed = seed_surv2,
    
    data_list = params,
    states = list(c("sb1", "sb2")),
    
    uses_par_sets = FALSE,
    
    evict_cor = FALSE
    
  ) %>%
  define_kernel(
    name = "SB1_to_SDL",
    family = "DD",
    formula = germ * surv2_seed,
    
    germ = germ_mean,
    surv2_seed = seed_surv2,
    
    data_list = params,
    states = list(c("sb1", "sdl")),
    
    uses_par_sets = FALSE,
    
    evict_cor = FALSE
    
  ) %>%
  define_kernel(
    name = "SB2_to_SDL",
    family = "DD",
    formula = germ * surv3_seed,
    
    germ       = germ_mean,
    surv3_seed = seed_surv3,
    
    data_list = params,
    states = list(c("sb2", "sdl")),
    
    uses_par_sets = FALSE,
    
    evict_cor = FALSE
    
  ) %>%
  define_kernel(
    name = "SDL_to_Plant_yr_loc",
    family = "DC",
    formula = sdl_s_yr_loc * d_size_yr_loc * d_stems,
    
    sdl_s_yr_loc = plogis(sdl_s_lin_yr_loc),
    sdl_s_lin_yr_loc = sdl_s_int + sdl_s_int_r_yr + sdl_s_int_r_yr_loc,
    d_size_yr_loc = dnorm(stems_2, sdl_d_lin_yr_loc, sdl_d_sd),
    sdl_d_lin_yr_loc = sdl_d_int + sdl_d_int_r_yr + sdl_d_int_r_yr_loc,
    
    data_list = params,
    states = list(c("sdl", "stems")),
    
    uses_par_sets = TRUE,
    par_set_indices = list(yr = c(2003:2020),
                           loc = c("CR", "HK", "KS", "RU")),
    
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "d_size_yr_loc")
  ) %>%
  define_impl(
    list(
      P_yr_loc =        list(int_rule = 'midpoint',
                             state_start = "stems",
                             state_end = "stems"),
      F_to_SDL_yr_loc =     list(int_rule = 'midpoint',
                                 state_start = "stems",
                                 state_end = "sdl"),
      F_to_SB1_yr_loc =     list(int_rule = 'midpoint',
                                 state_start = "stems",
                                 state_end = "sb1"),
      SB1_to_SB2 =   list(int_rule = "midpoint",
                          state_start = "sb1",
                          state_end = "sb2"),
      SB1_to_SDL =   list(int_rule = "midpoint",
                          state_start = "sb1",
                          state_end = "sdl"),
      SB2_to_SDL =   list(int_rule = "midpoint",
                          state_start = "sb2",
                          state_end = "sdl"),
      SDL_to_Plant_yr_loc = list(int_rule = "midpoint",
                                 state_start = "sdl",
                                 state_end = "stems")
    )
  ) %>% define_domains(
    stems = c(L, U, n)
  ) %>% 
  define_pop_state(
    pop_vectors = list(
      n_stems_yr_loc = runif(n),
      n_sb1_yr_loc = 20,
      n_sb2_yr_loc = 20,
      n_sdl_yr_loc = 20
    )
  ) %>%
  make_ipm(
    iterate = T,
    iterations = 1000
  )


lam <- stack(lambda(det_ipm) ) %>% 
  tidyr::separate(col = ind, into = c(NA, "year_t0", "locality"), sep = "_") %>%
  rename(lambda = values)


ind_lambda <- ggplot(lam, aes(x = year_t0, y = lambda, colour = locality)) + geom_point()

ggsave(ind_lambda, filename = here::here("results", "individual_lambdas.png"),
       width = 6.75, height = 4, units = "in", dpi = 300, type = "cairo")

deterministic_ipm = list(det_ipm = det_ipm,
                         individual_lambdas = lam)

saveRDS(deterministic_ipm, file = "results/deterministic_ipm.rds")

rm(list = ls())
