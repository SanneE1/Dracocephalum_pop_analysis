## simple environmental sampler function. 
## Now this is re-calculated for each year (so last years temp doesn't line up with 
## the lagged temp in current year)
sampling_env <- function(iteration, env_params, start_year = 2023) {
  
  lags_sim <- matrix(c(1:(env_params$lags + 1)), nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  
  ### get values for window function
  mon <- seq(from = 0, to = 0.95, length.out = 12)
  end <- (start_year + iteration) + mon[7]
  start <- end - (1/12 * (env_params$lags + 1))
  
  pet_sim <- rev(as.numeric(window(env_params[[grep("pet", names(env_params), value = T)]], start, end))) %>%
    matrix(., nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  pr_sim <- rev(as.numeric(window(env_params[[grep("pr", names(env_params), value = T)]], start, end))) %>%
    matrix(., nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  tas_sim <- rev(as.numeric(window(env_params[[grep("tas", names(env_params), value = T)]], start, end))) %>%
    matrix(., nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  
  
  return(list(lags = lags_sim,
              temp = tas_sim,
              precip = pr_sim,
              pet = pet_sim,
              shading = env_params$shading,
              rock = env_params$rock,
              slope = env_params$slope,
              soil_depth = env_params$soil_depth
  ))
}

FLM_clim_predict <- function(model, lag_vec, 
                             temp_vec, precip_vec, pet_vec, 
                             shading, rock, slope, soil) {
  
  ## dummy data list, won't be using the predictions for n_stems : year_t0 but these need to be provided for predict function
  new_data <- list(
    ln_stems_t0 = 1,  
    population = "CR",
    tot_shading_t0 = shading,
    slope = slope,
    rock = rock,
    soil_depth = soil,
    lags = matrix(lag_vec, nrow = 1),
    tas_scaledcovar = matrix(temp_vec, nrow = 1),
    pr_scaledcovar = matrix(precip_vec, nrow = 1),
    pet_scaledcovar = matrix(pet_vec, nrow = 1),
    year_t0 = 2016
  )
  
  pt <- mgcv::predict.gam(model, new_data, type = "terms", exclude = "s(year_t0)")
  
  a <-  sum(pt[, grep("s\\(", attributes(pt)[[2]][[2]], value = T)])
  
  return(a)
}

run_ipm <- function(params, env_params, locality, 
                    n_it = n_it, U, L, n){
  
  init_ipm("general", "di", "stoch", "param") %>%
    define_kernel(
      name      = "P_loc",
      family    = "CC",
      formula   = s_loc * g_loc * d_stems,
      
      s_loc         =  plogis(s_linear_loc),
      s_linear_loc  =  s_int + s_stems * stems_1 + s_site_loc + 
        FLM_clim_predict(model = surv_mod, lag_vec = lags, shading = shading,
                         temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock),
      g_loc         =  dnorm(stems_2, g_mu_loc, grow_sd),
      g_mu_loc      =  g_int + g_stems * stems_1 + g_site_loc +
        FLM_clim_predict(model = grow_mod, shading = shading, lag_vec = lags, 
                         temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock), 
      
      data_list     = params,
      states        = list(c("stems")),
      
      uses_par_sets = TRUE,
      par_set_indices = list(loc = c("CR", "HK", "KS", "RU")),
      
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", "g_loc")
    ) %>%
    define_kernel(
      name = "F_to_SDL_loc",
      family = "CD",
      formula = fp_loc * seedp_loc * seedn_loc * surv1_seed * germ * d_stems,
      
      fp_loc          = plogis(fp_linear_loc),
      fp_linear_loc   = fp_int + fp_stems * stems_1 + fp_site_loc +
        FLM_clim_predict(model = pflower_mod, lag_vec = lags, shading = shading,
                                         temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock),
      
      seedp_loc     = plogis(sp_linear_loc),
      sp_linear_loc   = sp_int + sp_stems * stems_1 + sp_site_loc +
        FLM_clim_predict(model = seedp_mod, lag_vec = lags, shading = shading,
                                       temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock),
      
      seedn_loc = ifelse(seedn_mod_loc > 515, 515, seedn_mod_loc),
      seedn_mod_loc     = exp(seedn_linear_loc),
      seedn_linear_loc = sn_int + sn_stems * stems_1 + sn_site_loc +
        FLM_clim_predict(model = seedn_mod, lag_vec = lags, shading = shading,
                         temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock),
      
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
      
      fp_loc          = plogis(fp_linear_loc),
      fp_linear_loc   = fp_int + fp_stems * stems_1 + fp_site_loc +
        FLM_clim_predict(model = pflower_mod, lag_vec = lags, shading = shading,
                         temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock),
      
      seedp_loc     = plogis(sp_linear_loc),
      sp_linear_loc   = sp_int + sp_stems * stems_1 + sp_site_loc +
        FLM_clim_predict(model = seedp_mod, lag_vec = lags, shading = shading,
                         temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock),
      
      seedn_loc = ifelse(seedn_mod_loc > 515, 515, seedn_mod_loc),
      seedn_mod_loc     = exp(seedn_linear_loc),
      seedn_linear_loc = sn_int + sn_stems * stems_1 + sn_site_loc +
        FLM_clim_predict(model = seedn_mod, lag_vec = lags, shading = shading,
                         temp_vec = temp, precip_vec = precip, pet_vec = pet, 
                         slope = slope, soil = soil_depth, rock = rock),
      
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
      usr_funs = list(FLM_clim_predict = FLM_clim_predict),
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
    simulate(x, nsim = ((n_it * 12) + (3*lag))) %>%
      ts(., start= c(2023,1), frequency = 12))
  
  # environmental params
  env_params <- append(
    clim_sim,
    list(
      lags = lag,
      shading = df_env$shading[i],
      slope = df_env$slope[i],
      rock = df_env$rock[i],
      soil_depth = df_env$soil_depth[i]
    ))
  
  ipm <- run_ipm(params = params, env_params = env_params, 
                 locality = toupper(loc), 
                 n_it = n_it, U = U, L = L, n = n)
  
  df1 <- tibble(time = df_env$time[i],
                    locality = loc,
                    model = df_env$model[i],
                    scenario = df_env$scenario[i],
                    shading = df_env$shading[i],
                slope = df_env$slope[i],
                rock = df_env$rock[i],
                soil_depth = df_env$soil_depth[i],
                    lambda = ipmr::lambda(ipm)#,
                    # all_lambda = paste(ipmr::lambda(ipm, type_lambda = "all"), collapse = ",")
                    )
  
  if(save){
    write.csv(df1, file = paste0("results/stoch_ipms/lambda_env_levels_", i, ".csv"), row.names = F)
  }
  
  return(df1)
}


