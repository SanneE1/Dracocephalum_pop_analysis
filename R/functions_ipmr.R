## simple environmental sampler function. 
## Now this is re-calculated for each year (so last years temp doesn't line up with 
## the lagged temp in current year)
sampling_env <- function(iteration, env_params, start_year = 2023) {
  suppressWarnings({
    
    ### get values for window function
    mon <- seq(from = 0, to = 0.95, length.out = 12)
    year_t1 <- (start_year + iteration)
    year_t0 <- year_t1 - 1
    
    pr_spring <- mean(as.numeric(window(env_params[[grep("pr", names(env_params), value = T)]], start = (year_t1+mon[2]), end = (year_t1+mon[5]))))
    tas_spring <- mean(as.numeric(window(env_params[[grep("tas", names(env_params), value = T)]], start = (year_t1+mon[2]), end = (year_t1+mon[5]))))
    pet_spring <- mean(as.numeric(window(env_params[[grep("pet", names(env_params), value = T)]], start = (year_t1+mon[2]), end = (year_t1+mon[5]))))
    
    pr_dormant <- mean(as.numeric(window(env_params[[grep("pr", names(env_params), value = T)]], start = (year_t0+mon[8]), end = (year_t1+mon[2]))))
    tas_dormant <- mean(as.numeric(window(env_params[[grep("tas", names(env_params), value = T)]], start = (year_t0+mon[8]), end = (year_t1+mon[2]))))
    pet_dormant <- mean(as.numeric(window(env_params[[grep("pet", names(env_params), value = T)]], start = (year_t0+mon[8]), end = (year_t1+mon[2]))))
    
    pr_summer <- mean(as.numeric(window(env_params[[grep("pr", names(env_params), value = T)]], start = (year_t0+mon[5]), end = (year_t0+mon[8]))))
    tas_summer <- mean(as.numeric(window(env_params[[grep("tas", names(env_params), value = T)]], start = (year_t0+mon[5]), end = (year_t0+mon[8]))))
    pet_summer <- mean(as.numeric(window(env_params[[grep("pet", names(env_params), value = T)]], start = (year_t0+mon[5]), end = (year_t0+mon[8]))))
    
  })
  
  return(list(pr_spring = pr_spring,
              tas_spring = tas_spring,
              pet_spring = pet_spring,
              pr_dormant = pr_dormant,
              tas_dormant = tas_dormant,
              pet_dormant = pet_dormant,
              pr_summer = pr_summer,
              tas_summer = tas_summer,
              pet_summer = pet_summer
  ))
}

mod_pred <- function(model, size,
                     env_params,  
                     pr_spring, tas_spring, pet_spring,
                     pr_dormant, tas_dormant, pet_dormant,
                     pr_summer, tas_summer, pet_summer) {
  ## New datalist, 
  new_data <- data.frame(
    survival_t1 = 1,
    ln_stems_t0 = size,  
    population = factor(rep(toupper(env_params$locality), n), levels = c("CR", "HK", "KS", "RU")),
    slope = env_params$slope,
    rock = env_params$rock,
    soil_depth = env_params$soil_depth,
    ## scaling these values given the center and scale attributes
    herb_shading_t0 = (env_params$herb_shading - env_params$herb_center) / env_params$herb_scale,
    shrub_shading_t0 = (env_params$shrub_shading - env_params$shrub_center) / env_params$shrub_scale,
    
    pr_dormant = pr_dormant,
    tas_dormant = tas_dormant,
    pet_dormant = pet_dormant,
    pr_spring = pr_spring,
    tas_spring = tas_spring,
    pet_spring = pet_spring,
    pr_summer = pr_summer,
    tas_summer = tas_summer,
    pet_summer = pet_summer
  )
  
  new_data = model.matrix( ~ ln_stems_t0 + population + rock + slope + soil_depth + 
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
                           data = new_data)[, -1]
  
  pt <- predict(model, newx = new_data, type = "response")
  return(pt)
}


run_ipm <- function(params, env_params, locality, 
                    n_it = n_it, U, L, n){
  
  init_ipm("general", "di", "stoch", "param") %>%
    define_kernel(
      name      = "P_loc",
      family    = "CC",
      formula   = s_loc * g_loc * d_stems,
      
      s_loc         = mod_pred(model = surv_mod, size = stems_1, 
                               env_params = env_params, 
                               pr_spring = pr_spring, tas_spring = tas_spring, pet_spring = pet_spring,
                               pr_dormant = pr_dormant, tas_dormant = tas_dormant, pet_dormant = pet_dormant,
                               pr_summer = pr_summer, tas_summer = tas_summer, pet_summer = pet_summer),
      g_loc         =  dnorm(stems_2, g_mu_loc, grow_sd),
      g_mu_loc      =  mod_pred(model = grow_mod, size = stems_1, 
                                env_params, 
                                pr_spring, tas_spring, pet_spring,
                                pr_dormant, tas_dormant, pet_dormant,
                                pr_summer, tas_summer, pet_summer),
      
      data_list     = params,
      states        = list(c("stems")),
      
      uses_par_sets = TRUE,
      par_set_indices = list(loc = c("CR", "HK", "KS", "RU")),
      
      evict_cor     = TRUE,
      evict_fun     = discrete_extrema("g_loc", "stems", 100, 100)
    ) %>%
    define_kernel(
      name = "F_to_SDL_loc",
      family = "CD",
      formula = fp_loc * seedp_loc * seedn_loc * surv1_seed * germ * d_stems,
      
      fp_loc          = mod_pred(model = pflower_mod, size = stems_1, 
                                 env_params, 
                                 pr_spring, tas_spring, pet_spring,
                                 pr_dormant, tas_dormant, pet_dormant,
                                 pr_summer, tas_summer, pet_summer),
      
      seedp_loc     = mod_pred(model = seedp_mod, size = stems_1, 
                               env_params, 
                               pr_spring, tas_spring, pet_spring,
                               pr_dormant, tas_dormant, pet_dormant,
                               pr_summer, tas_summer, pet_summer),
      
      seedn_loc = ifelse(seedn_mod_loc > 515, 515, seedn_mod_loc),
      seedn_mod_loc     = mod_pred(model = seedn_mod, size = stems_1, 
                                   env_params, 
                                   pr_spring, tas_spring, pet_spring,
                                   pr_dormant, tas_dormant, pet_dormant,
                                   pr_summer, tas_summer, pet_summer),
      
      germ        = germ_mean, 
      surv1_seed  = seed_surv1,
      
      data_list = params,
      states = list(c("stems", "sdl")),
      
      uses_par_sets = TRUE,
      par_set_indices = list(loc = c("CR", "HK", "KS", "RU")),
      
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name = "F_to_SB1_loc",
      family = "CD",
      formula = fp_loc * seedp_loc * seedn_loc * surv1_seed * (1-germ) * d_stems,
      
      fp_loc          = mod_pred(model = pflower_mod, size = stems_1, 
                                 env_params, 
                                 pr_spring, tas_spring, pet_spring,
                                 pr_dormant, tas_dormant, pet_dormant,
                                 pr_summer, tas_summer, pet_summer),
      
      seedp_loc     = mod_pred(model = seedp_mod, size = stems_1, 
                               env_params, 
                               pr_spring, tas_spring, pet_spring,
                               pr_dormant, tas_dormant, pet_dormant,
                               pr_summer, tas_summer, pet_summer),
      
      seedn_loc = ifelse(seedn_mod_loc > 515, 515, seedn_mod_loc),
      seedn_mod_loc     = mod_pred(model = seedn_mod, size = stems_1, 
                                   env_params, 
                                   pr_spring, tas_spring, pet_spring,
                                   pr_dormant, tas_dormant, pet_dormant,
                                   pr_summer, tas_summer, pet_summer),
      
      germ        = germ_mean, 
      surv1_seed  = seed_surv1,
      
      data_list = params,
      states = list(c("stems", "sb1")),
      
      uses_par_sets = TRUE,
      par_set_indices = list(loc = c("CR", "HK", "KS", "RU")),
      
      evict_cor = FALSE
      
    ) %>%
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
      name = "SDL_to_Plant_loc",
      family = "DC",
      formula = sdl_surv_loc * d_size_loc * d_stems,
      
      sdl_surv_loc     = plogis(sdl_s_linear),
      sdl_s_linear = sdl_s_int,
      d_size_loc       = dnorm(stems_2, sdl_d_linear, sdl_size_d_sd),
      sdl_d_linear = sdl_d_int,
      
      data_list = params,
      states = list(c("sdl", "stems")),
      
      uses_par_sets = TRUE,
      par_set_indices = list(loc = c("CR", "HK", "KS", "RU")),
      
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "d_size_loc")
      
    ) %>%
    define_impl(
      list(
        P_loc =        list(int_rule = 'midpoint',
                            state_start = "stems",
                            state_end = "stems"),
        F_to_SDL_loc =     list(int_rule = 'midpoint',
                                state_start = "stems",
                                state_end = "sdl"),
        F_to_SB1_loc =     list(int_rule = 'midpoint',
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
        SDL_to_Plant_loc = list(int_rule = "midpoint",
                                state_start = "sdl",
                                state_end = "stems")
      )
    ) %>%
    define_domains(
      stems = c(L, U, n)
    ) %>%
    define_env_state(
      env_covs = sample_env(iteration = t, env_params = env_params),
      data_list = list(env_params = env_params,
                       sample_env = sampling_env)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_stems = runif(n),
        n_sb1 = 20,
        n_sb2 = 20,
        n_sdl = 20
      )
    ) %>%
    make_ipm(
      iterate = T,
      iterations = n_it,
      kernel_seq = rep(locality, n_it),
      usr_funs = list(mod_pred = mod_pred),
      return_sub_kernels = T
    )
  
}

ipm_loop <- function(i, df_env, params,
                     climate_models,
                     n_it, U = U, L = L, n = n, save = F) {
  print(i)
  loc <- df_env$localities[i]
  
  clim_mod <- climate_models[[grep(as.character(df_env$time[i]), 
                                   names(climate_models), value = T)]]
  clim_mod <- clim_mod[grep(as.character(loc), names(clim_mod), value = T)]
  
  if(as.character(df_env$time[i]) == "future") {
    
    txt <- paste0(df_env$model[i], '.*', df_env$scenario[i])
    clim_mod <- clim_mod[grep(txt, names(clim_mod), value = T)]
    
  }
  
  clim_sim <- lapply(clim_mod, function(x)
    simulate(x, nsim = (((n_it + 1) * 12))) %>%
      ts(., start= c(2023,1), frequency = 12))
  
  # environmental params
  env_params <- append(
    clim_sim,
    list(
      herb_shading = df_env$herb_shading[i],
      herb_center = params$herb_center,
      herb_scale = params$herb_scale,
      shrub_shading = df_env$shrub_shading[i],
      shrub_center = params$shrub_center,
      shrub_scale = params$shrub_scale,
      slope = df_env$slope[i],
      rock = df_env$rock[i],
      soil_depth = df_env$soil_depth[i],
      locality = df_env$localities[i]
    ))
  
  ipm <- run_ipm(params = params, env_params = env_params, 
                 locality = toupper(loc), 
                 n_it = n_it, U = U, L = L, n = n)
  
  df1 <- tibble(time = df_env$time[i],
                locality = df_env$localities[i],
                model = df_env$model[i],
                scenario = df_env$scenario[i],
                herb_shading = df_env$herb_shading[i],
                shrub_shading = df_env$shrub_shading[i],
                slope = env_params$slope,
                rock = env_params$rock,
                soil_depth = env_params$soil_depth,
                lambda = ipmr::lambda(ipm)
  )
  
  if(save){
    write.csv(df1, file = paste0("results/stoch_ipms/lambda_env_levels_", i, ".csv"), row.names = F)
  }
  
  return(df1)
}


