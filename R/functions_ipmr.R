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
              slope = env_params$slope
  ))
}

FLM_clim_predict <- function(model, lag_vec, 
                             temp_vec, precip_vec, pet_vec, 
                             shading) {
  
  ## dummy datalist, won't be using the predictions for n_stems : year_t0 but these need to be provided for predict function
  new_data <- list(
    ln_stems_t0 = 2,  
    population = "CR",
    tot_shading_t0 = 0,
    slope = 0,
    soil_depth = 0,
    year_t0 = 2016,
    tot_shading_m = matrix(rep(shading, length(lag_vec)), nrow = 1),
    lags = matrix(lag_vec, nrow = 1),
    tas_scaledcovar = matrix(temp_vec, nrow = 1),
    pr_scaledcovar = matrix(precip_vec, nrow = 1),
    pet_scaledcovar = matrix(pet_vec, nrow = 1)
  )
  
  pt <- mgcv::predict.gam(model, new_data, type = "terms")
  
  
  
  return(
    sum(pt[, grep("scaledcovar", attributes(pt)[[2]][[2]], value = T)])
  )
}

proto_ipm <- function(params, env_params, locality, 
                      n_it = n_it, U, L, n) {
  init_ipm("general", "di", "stoch", "param") %>%
    define_kernel(
      name      = "P_loc",
      family    = "CC",
      formula   = s_loc * g_loc * d_stems,
      
      s_loc         =  plogis(s_linear_loc),
      s_linear_loc  =  s_int + s_stems * stems_1 + s_site_loc + 
        s_shading * shading + 
        FLM_clim_predict(model = surv_mod, ### spline model prediction
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading
        ), 
      g_loc         =  dnorm(stems_2, g_mu_loc, grow_sd),
      g_mu_loc      =  g_int + g_stems * stems_1 + g_site_loc +
        FLM_clim_predict(model = grow_mod, 
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading), ### model prediction
      
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
      formula = fp_loc * pabort_loc * n_seeds_loc * surv1_seed * germ * d_stems,
      
      fp_loc          = plogis(fp_linear_loc),
      fp_linear_loc   = fp_int + fp_stems * stems_1 + fp_site_loc +
        fp_slope * slope +
        FLM_clim_predict(model = pflower_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      pabort_loc     = plogis(ab_linear_loc),
      ab_linear_loc   = ab_int + ab_stems * stems_1 + ab_site_loc +
        FLM_clim_predict(model = pabort_mod,
                         lag_vec = lags,
                         temp_vec = temp,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      n_seeds_loc     = exp(n_seeds_linear_loc),
      n_seeds_linear_loc = ns_int + ns_stems * stems_1 + ns_site_loc +
        FLM_clim_predict(model = nseed_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
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
      formula = fp_loc * pabort_loc * n_seeds_loc * surv1_seed * (1-germ) * d_stems,
      
      fp_loc          = plogis(fp_linear_loc),
      fp_linear_loc   = fp_int + fp_stems * stems_1 + fp_site_loc +
        fp_slope * slope + 
        FLM_clim_predict(model = pflower_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      pabort_loc     = plogis(ab_linear_loc),
      ab_linear_loc   = ab_int + ab_stems * stems_1 + ab_site_loc +
        FLM_clim_predict(model = pabort_mod,
                         lag_vec = lags,
                         temp_vec = temp,
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
      n_seeds_loc     = exp(n_seeds_linear_loc),
      n_seeds_linear_loc = ns_int + ns_stems * stems_1 + ns_site_loc +
        FLM_clim_predict(model = nseed_mod,
                         lag_vec = lags, 
                         temp_vec = temp, 
                         precip_vec = precip,
                         pet_vec = pet,
                         shading = shading),
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
    ) 
}


run_ipm <- function(params, env_params, locality, 
                    n_it = n_it, U, L, n){
  
   proto_ipm(params, env_params, locality, 
             n_it, U, L, n) %>%
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
                     n_it, U = U, L = L, n = n) {
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
      slope = df_env$slope[i]
    ))
  ipm <- run_ipm(params = params, env_params = env_params, 
                 locality = toupper(loc), 
                 n_it = n_it, U = U, L = L, n = n)
  
  df1 <- data.frame(time = df_env$time[i],
                    locality = loc,
                    model = df_env$model[i],
                    scenario = df_env$scenario[i],
                    shading = df_env$shading[i],
                    slope = df_env$slope[i],
                    lambda = lambda(ipm))
  
  return(df1)
}

ipm_ext_p <- function(i, df_env, params,
                      clim_ts,
                      n_it, U = U, L = L, n = n) {
  
  
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  
  txt <- paste(locality, model, scenario, sep = '.*')
  clim_mod <- clim_ts[grep(txt, names(clim_ts), value = T)] %>%
    lapply(., function(x) x$value)
  
  env_params <- append(
    clim_mod,
    list(lags = lag,
         shading = df_env$shading[i],
         slope = df_env$slope[i]
    ))
  
  
  a <- proto_ipm(params = params, env_params = env_params, locality = toupper(locality),
                 n_it = n_it, U = U, L = L, n = n) %>%
    define_env_state(
      env_covs = sample_env(iteration = t, env_params = env_params, start_year = 2022),
      data_list = list(env_params = env_params,
                       sample_env = sampling_env)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_stems = pop_vec,
        n_sb1 = 2000,
        n_sb2 = 2000,
        n_sdl = sdl_n
      )
    ) %>%
    make_ipm(
      iterate = TRUE,
      iterations = n_it,
      kernel_seq = rep(toupper(locality), n_it),
      usr_funs = list(FLM_clim_predict = FLM_clim_predict),
      return_sub_kernels = TRUE,
      normalize_pop_size = FALSE
    )
  
  pop_time <- apply(a$pop_state$n_stems, 2, sum)
  
  b <- data.frame(
    locality = locality,
    model = model,
    scenario = scenario,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    below_ext_size = any(pop_time < 10), 
    yr_of_ext = ifelse(any(pop_time < 10), min(which(pop_time < 10)), NA)
  )
  b$pop_size = list(apply(a$pop_state$n_stems, 2, sum))
  
  return(b)
  
}


man_trans <- function(i, df_env, params,
                        clim_ts, n_it,
                        U, L, n) {

  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  yrs_between_man <- df_env$yrs_between_man[i]
  
  
  txt <- paste(locality, model, scenario, sep = '.*')
  clim_mod <- clim_ts[grep(txt, names(clim_ts), value = T)] %>%
    lapply(., function(x) x$value)
  
  env_params <- append(
    clim_mod,
    list(lags = lag,
         shading = df_env$shading[i],
         slope = df_env$slope[i]
    ))
  
  
pop_vec1 <- pop_vec
sdl_n1 <- sdl_n
sb1_n1 <- 2000
sb2_n1 <- 2000

pop_size <- c(sum(pop_vec))
lambdas <- c()

for(j in c(0:round(n_it/yrs_between_man, 0))) {

  
  start_yr <- (2022 + (j * yrs_between_man))
  
  if(start_yr >= 2100) next
  if(start_yr + yrs_between_man >= 2100) {
    yrs_between_man = 2100 - start_yr
  }
  
  
  a <- proto_ipm(params = params, env_params = env_params, locality = toupper(locality),
                 n_it = yrs_between_man, U = U, L = L, n = n) %>%
    define_env_state(
      env_covs = sample_env(iteration = t, env_params = env_params, 
                            start_year = start_yr),
      data_list = list(env_params = env_params,
                       sample_env = sampling_env,
                       start_yr = start_yr)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_stems = pop_vec1,
        n_sb1 = sb1_n1,
        n_sb2 = sb2_n1,
        n_sdl = sdl_n1
      )
    ) %>%
    make_ipm(
      iterate = TRUE,
      iterations = yrs_between_man,
      kernel_seq = rep(toupper(locality), yrs_between_man),
      usr_funs = list(FLM_clim_predict = FLM_clim_predict),
      return_sub_kernels = TRUE,
      normalize_pop_size = FALSE
    )
  
  pop_size = c(pop_size, 
               apply(a$pop_state$n_stems, 2, sum)[2:(yrs_between_man+1)])
  lambdas = c(lambdas,
              a$pop_state$lambda)
  
  
  pop_vec1 <- a$pop_state$n_stems[,(yrs_between_man+1)]
  sdl_n1 <- a$pop_state$n_sdl[,(yrs_between_man+1)]
  sb1_n1 <- a$pop_state$n_sb1[,(yrs_between_man+1)]
  sb2_n1 <- a$pop_state$n_sb2[,(yrs_between_man+1)]
  
  pop_vec1[19] <- pop_vec1[19] + (df_env$effort[i]*0.7)   ## 19 is meshpoint for individuals sized ~ log(2)
  pop_vec1[30] <- pop_vec1[30] + (df_env$effort[i]*0.3)   ## 30 is meshpoint for individuals sized ~ log(3)
  
}
  
  b <- data.frame(
    locality = locality,
    model = model,
    scenario = scenario,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    effort = df_env$effort[i],
    yrs_between_man = yrs_between_man,
    below_ext_size = any(pop_size < 10), 
    yr_of_ext = ifelse(any(pop_size < 10), min(which(pop_size < 10)), NA)
  )
  b$pop_size = list(pop_size)
  b$lambdas = list(lambdas) 
  
  
  return(b)
  
  
}

man_seedadd <- function(i, df_env, params,
                      clim_ts, n_it,
                      U, L, n) {
  
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  yrs_between_man <- df_env$yrs_between_man[i]
  
  
  txt <- paste(locality, model, scenario, sep = '.*')
  clim_mod <- clim_ts[grep(txt, names(clim_ts), value = T)] %>%
    lapply(., function(x) x$value)
  
  env_params <- append(
    clim_mod,
    list(lags = lag,
         shading = df_env$shading[i],
         slope = df_env$slope[i]
    ))
  
  
  pop_vec1 <- pop_vec
  sdl_n1 <- sdl_n
  sb1_n1 <- 2000
  sb2_n1 <- 2000
  
  pop_size <- c(sum(pop_vec))
  lambdas <- c()
  
  for(j in c(0:round(n_it/yrs_between_man, 0))) {
    
    
    start_yr <- (2022 + (j * yrs_between_man))
    
    if(start_yr >= 2100) next
    if(start_yr + yrs_between_man >= 2100) {
      yrs_between_man = 2100 - start_yr
    }
    
    
    a <- proto_ipm(params = params, env_params = env_params, locality = toupper(locality),
                   n_it = yrs_between_man, U = U, L = L, n = n) %>%
      define_env_state(
        env_covs = sample_env(iteration = t, env_params = env_params, 
                              start_year = start_yr),
        data_list = list(env_params = env_params,
                         sample_env = sampling_env,
                         start_yr = start_yr)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_stems = pop_vec1,
          n_sb1 = sb1_n1,
          n_sb2 = sb2_n1,
          n_sdl = sdl_n1
        )
      ) %>%
      make_ipm(
        iterate = TRUE,
        iterations = yrs_between_man,
        kernel_seq = rep(toupper(locality), yrs_between_man),
        usr_funs = list(FLM_clim_predict = FLM_clim_predict),
        return_sub_kernels = TRUE,
        normalize_pop_size = FALSE
      )
    
    pop_size = c(pop_size, 
                 apply(a$pop_state$n_stems, 2, sum)[2:(yrs_between_man+1)])
    lambdas = c(lambdas,
                a$pop_state$lambda)
    
    
    pop_vec1 <- a$pop_state$n_stems[,(yrs_between_man+1)]
    sdl_n1 <- a$pop_state$n_sdl[,(yrs_between_man+1)]
    sb1_n1 <- a$pop_state$n_sb1[,(yrs_between_man+1)]
    sb2_n1 <- a$pop_state$n_sb2[,(yrs_between_man+1)]
    
    sb1_n1 <- sb1_n1 + df_env$effort[i]
    
  }
  
  b <- data.frame(
    locality = locality,
    model = model,
    scenario = scenario,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    effort = df_env$effort[i],
    yrs_between_man = yrs_between_man,
    below_ext_size = any(pop_size < 10), 
    yr_of_ext = ifelse(any(pop_size < 10), min(which(pop_size < 10)), NA)
  )
  b$pop_size = list(pop_size)
  b$lambdas = list(lambdas) 
  
  
  return(b)
  
  
}
