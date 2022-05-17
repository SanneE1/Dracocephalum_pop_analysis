## simple environmental sampler function. 
## Now this is re-calculated for each year (so last years temp doesn't line up with 
## the lagged temp in current year)
sampling_env <- function(time, env_params, init_temp = 0, init_prcp = 0) {
  
  lags_est <- matrix(c(1:(env_params$lags + 1)), nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  
  t_est <-  matrix(c(init_temp + env_params$yr_increase_temp * time + rnorm(13,0,1),
                     init_temp + env_params$yr_increase_temp * (time - 1) + rnorm(12,0,1),
                     init_temp + env_params$yr_increase_temp * (time - 2) + rnorm(12,0,1),
                     init_temp + env_params$yr_increase_temp * (time - 3) + rnorm(12,0,1)), 
                   nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  
  
  p_est <-  matrix(c(init_prcp + env_params$yr_increase_prcp * time + rnorm(13,0,1),
                     init_prcp + env_params$yr_increase_prcp * (time - 1) + rnorm(12,0,1),
                     init_prcp + env_params$yr_increase_prcp * (time - 2) + rnorm(12,0,1),
                     init_prcp + env_params$yr_increase_prcp * (time - 3) + rnorm(12,0,1)), 
                   nrow = 1, ncol = (env_params$lags + 1), byrow = T)
  
  return(list(lags = lags_est,
              temp = t_est,
              precip = p_est,
              herb_shading = env_params$herb_shading,
              shrub_shading = env_params$shrub_shading,
              slope = env_params$slope,
              soil = env_params$soil))
}

FLM_clim_predict <- function(model, lag_vec, temp_vec, precip_vec) {
  
  ## dummy datalist, won't be using the predictions for n_stems : year_t0 but these need to be provided for predict function
  new_data <- list(
    ln_stems_t0 = 2,  
    population = "CR",
    herb_shading_t0 = 2,
    shrub_shading_t0 = 2,
    slope = 20,
    year_t0 = 2010,
    soil_d3 = 4,
    lags = matrix(lag_vec, nrow = 1),
    Tavg_scaledcovar = matrix(temp_vec, nrow = 1),
    Prec_scaledcovar = matrix(precip_vec, nrow = 1)
  )
  
  pt <- mgcv::predict.gam(model, new_data, type = "terms")
  
  return(
    sum(pt[, c('s(lags):Tavg_scaledcovar', 's(lags):Prec_scaledcovar')])
  )
}


run_ipm <- function(params, env_params, locality, 
                    n_it = 200, U, L, n){
  
  init_ipm("general", "di", "stoch", "param") %>%
    define_kernel(
      name      = "P_loc",
      family    = "CC",
      formula   = s_loc * g * d_stems,
      
      s_loc         =  plogis(s_linear_loc),
      s_linear_loc  =  s_int + s_stems * stems_1 + s_site_loc + s_herb * herb_shading +
        FLM_clim_predict(model = surv_mod, ### spline model prediction
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip), 
      g         =  dnorm(stems_2, g_mu, grow_sd),
      g_mu      =  g_int + g_stems * stems_1 + g_herb * herb_shading +
        FLM_clim_predict(model = grow_mod, 
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip), ### model prediction
      
      data_list     = params,
      states        = list(c("stems")),
      
      uses_par_sets = TRUE,
      par_set_indices = list(loc = c("CR", "HK", "KS", "RU")),
      
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", "g")
    ) %>%
    define_kernel(
      name = "F_to_SDL_loc",
      family = "CD",
      formula = fp_loc * pabort_loc * n_seeds_loc * surv1_seed * germ * d_stems,
      
      fp_loc          = plogis(fp_linear_loc),
      fp_linear_loc   = fp_int + fp_stems * stems_1 + 
        fp_shrub * shrub_shading + fp_slope * slope + fp_soil * soil +
        FLM_clim_predict(model = pflower_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip),
      pabort_loc     = plogis(ab_linear_loc),
      ab_linear_loc   = ab_int + ab_stems * stems_1 + ab_shrub * shrub_shading + ab_soil * soil +
        FLM_clim_predict(model = pabort_mod,
                         lag_vec = lags,
                         temp_vec = temp,
                         precip_vec = precip),
      n_seeds_loc     = exp(n_seeds_linear_loc),
      n_seeds_linear_loc = ns_int + ns_stems * stems_1 + ns_slope * slope +
        FLM_clim_predict(model = nseed_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip),
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
      formula = fp_loc * pabort_loc * n_seeds_loc *  surv1_seed * (1-germ) * d_stems,
      
      fp_loc          = plogis(fp_linear_loc),
      fp_linear_loc   = fp_int + fp_stems * stems_1 + 
        fp_shrub * shrub_shading + fp_slope * slope + fp_soil * soil +
        FLM_clim_predict(model = pflower_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip),
      pabort_loc     = plogis(ab_linear_loc),
      ab_linear_loc   = ab_int + ab_stems * stems_1 + ab_shrub * shrub_shading + ab_soil * soil +
        FLM_clim_predict(model = pabort_mod,
                         lag_vec = lags,
                         temp_vec = temp,
                         precip_vec = precip),
      n_seeds_loc     = exp(n_seeds_linear_loc),
      n_seeds_linear_loc = ns_int + ns_stems * stems_1 + ns_slope * slope +
        FLM_clim_predict(model = nseed_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip),
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
      sdl_s_linear = sdl_s_int + sdl_s_site_loc + sdl_s_herb * herb_shading,
      d_size_loc       = dnorm(stems_2, sdl_d_linear, sdl_size_d_sd),
      sdl_d_linear = sdl_d_int + sdl_d_site_loc + sdl_d_herb * herb_shading,
      
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
      env_covs = sample_env(time = t, env_params = env_params),
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

ipm_loop <- function(i, df_env, params = params,
                     n_it, U = U, L = L, n = n) {
  
  # environmental params
  env_params <- list(
    lags = lag,
    yr_increase_temp = 0,
    yr_increase_prcp = 0,
    herb_shading = df_env$herb_shading[i],
    shrub_shading = df_env$shrub_shading[i],
    slope = df_env$slope[i],
    soil = df_env$soil[i]
  )
  
  ipm <- run_ipm(params = params, env_params = env_params, 
                 locality = as.character(df_env$localities[i]), 
                 n_it = n_it, U = U, L = L, n = n)
  
  
  df1 <- as.data.frame(env_params)
  df1$lambda <- lambda(ipm)
  
  return(df1)
}

ipm_elast <- function(params = params,
                      mean_env_params, mean_lambda, 
                      pertub_size = 0.01, parameter_perturb,
                      n_it, U = U, L = L, n = n) {
  
  ## perturb parameter in params list
  params_high <- params
  params_high[[parameter_perturb]] <- params[[parameter_perturb]] + pertub_size
  
  params_low <- params
  params_low[[parameter_perturb]] <- params[[parameter_perturb]] - pertub_size
  
  
  ipm_high <- run_ipm(params = params_high, env_params = mean_env_params, 
                      locality = "CR", 
                      n_it = n_it, U = U, L = L, n = n)
  
  ipm_low <- run_ipm(params = params_low, env_params = mean_env_params, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  lambda_high <- lambda(ipm_high)
  lambda_low <- lambda(ipm_low)
  sensitivity <- (lambda_high - lambda_low)/(2*pertub_size)
  elasticity <- sensitivity * params[[parameter_perturb]] / mean_lambda
  
  
  return(
    data.frame(parameter = parameter_perturb,
               sensitivity = sensitivity,
               elasticity = elasticity)
  )
  
}

