






ipm_ext_p <- function(i, df_env, params,
                      clim_ts, clim_hist,
                      pop_vec, sdl_n,
                      n_it, U = U, L = L, n = n) {
  
  
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  
  # Select pop_vectors
  pop_vec <- pop_vec[[grep(toupper(locality), names(pop_vec), value = T)]] %>% as.vector()
  sdl_n <- sdl_n[[grep(toupper(locality), names(sdl_n), value = T)]] %>% as.vector()
  
  if(model == "No change") {
    
    clim_mod <- hist_clim[grep(locality, names(hist_clim), value = T, ignore.case = T)]
    
  } else{
    
    txt <- paste(locality, model, scenario, sep = '.*')
    clim_mod <- clim_ts[grep(txt, names(clim_ts), value = T)] %>%
      lapply(., function(x) x$value)
    
  }
  
  env_params <- append(
    clim_mod,
    list(lags = lag,
         shading = df_env$shading[i],
         slope = df_env$slope[i]
    ))
  
  
  ## IPM
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
                      clim_ts, 
                      pop_vec, sdl_n,
                      n_it,
                      U, L, n) {
  
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  yrs_between_man <- df_env$yrs_between_man[i]
  
  # Select pop_vectors
  pop_vec <- pop_vec[[grep(toupper(locality), names(pop_vec), value = T)]] %>% as.vector()
  sdl_n <- sdl_n[[grep(toupper(locality), names(sdl_n), value = T)]] %>% as.vector()
  
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
                        pop_vec, sdl_n,
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
  
  
  # Select pop_vectors
  pop_vec <- pop_vec[[grep(toupper(locality), names(pop_vec), value = T)]] %>% as.vector()
  sdl_n <- sdl_n[[grep(toupper(locality), names(sdl_n), value = T)]] %>% as.vector()
  
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