## Vital rate functions

mod_pred <- function(indi, model,
                     params = params, env_params = env_params, 
                     locality = locality, clim = clim,
                     lag = 24) {
  n = nrow(indi)
  ## New datalist, 
  new_data <- list(
    ln_stems_t0 = indi$size,  
    population = rep(toupper(locality), n),
    tot_shading_t0 = indi$shading,
    # slope = indi$slope,
    year_t0 = rep(2016, n),   ## just a dummy, should be ignored when predicting below
    tot_shading_m = matrix(indi$shading, ncol = (lag+1), nrow = n),
    lags = matrix(clim$lags, ncol = (lag+1), nrow = n, byrow = T),
    tas_scaledcovar = matrix(clim$temp, ncol = (lag+1), nrow = n, byrow = T),
    pr_scaledcovar = matrix(clim$precip, ncol = (lag+1), nrow = n, byrow = T),
    pet_scaledcovar = matrix(clim$pet, ncol = (lag+1), nrow = n, byrow = T)
  )
  
  pt <- mgcv::predict.gam(model, new_data, 
                          exclude = 'year_t0',
                          type = "response")
  return(pt)
}

yearly_loop <- function(yr, env_params, locality,
                        indi, sdl, sb1, sb2){ #,
                        # slope_mean, slope_sd) {
  
    # retrieve this year's climate
    clim <- sampling_env(iteration = yr, env_params = env_params, start_year = 2021)
    
    # calculate population size 
    pop_size <- nrow(indi)
    
    # Plants ------------------------------------------------------
    
    # generate binomial random number for survival
    surv <- rbinom(n = pop_size, prob = mod_pred(indi = indi, model = params$surv_mod,
                                                 params = params, env_params = env_params, 
                                                 locality = locality, clim = clim), 
                   size = 1)
    
    # estimate size for surviving individuals
    z1 <- rep(NA, pop_size)
    z1[which(surv == 1)] <- rnorm(n = length(which(surv == 1)),
                                  mean = mod_pred(indi = indi[which(surv == 1),], model = params$grow_mod,
                                                  params = params, env_params = env_params, 
                                                  locality = locality, clim = clim),
                                  sd = sd(resid(params$grow_mod)))
    
    indi_z1 <- cbind(indi, surv, z1)
    
    # calculate total number of seeds produced
    flowp <- rbinom(n = pop_size, prob = mod_pred(indi = indi, model = params$pflower_mod,
                                                  params = params, env_params = env_params, 
                                                  locality = locality, clim = clim), 
                    size = 1)
    
    abp <- rep(NA, pop_size)
    abp[which(flowp == 1)] <- rbinom(n = length(which(flowp == 1)), 
                                     prob = mod_pred(indi = indi[which(flowp == 1),], 
                                                     model = params$pabort_mod,
                                                     params = params, env_params = env_params, 
                                                     locality = locality, clim = clim), 
                                     size = 1)
    
    nseeds <- rep(NA, pop_size)
    seed_mean <- mod_pred(indi = indi[which(abp == 1),], model = params$nseed_mod,
                          params = params, env_params = env_params, 
                          locality = locality, clim = clim)
    
    seed_gamma_shape1 <- (seed_mean^2)/(params$nseed_sd^2)
    seed_gamma_shape2 <- (seed_mean)/(params$nseed_sd^2)
    nseeds[which(abp == 1)] <- round(
      rgamma(n = length(which(abp == 1)), shape = seed_gamma_shape1, scale = seed_gamma_shape2),
      0)
    
    ## cap individual seed production to 515
    nseeds[which(nseeds > 515)] <- 515
    
    
    tot_seeds <- sum(nseeds, na.rm = T)
    
    # Divide seeds into seedling and seedbank
    
    to_sb1 <- rbinom(1, prob = (1-params$germ_mean), size = tot_seeds)
    to_sdl <- tot_seeds - to_sb1
    
    # Discrete stages  ------------------------------------------------------
    
    ## do seeds in seedbanks germinate or stay sb
    sb1_to_sdl <- rbinom(1, prob = params$germ_mean, size = sb1)
    sb1_to_sb2 <- sb1 - sb1_to_sdl
    
    sb2_to_sdl <- rbinom(1, prob = params$germ_mean, size = sb2) 
    
    sdl_to_plant <- rbinom(1, prob = params$sdl_surv_mean, size = sdl)
    
    ## how many transitioning seeds survive to next census in each stage  
    sb1 <- rbinom(1, prob = params$seed_surv1, size = to_sb1)
    sb2 <- rbinom(1, prob = params$seed_surv2, size = sb1_to_sb2)
    
    sdl <- rbinom(1, prob = params$seed_surv2, size = sb1_to_sdl) +
      rbinom(1, prob = params$seed_surv3, size = sb2_to_sdl)
    
    # New recruits into plant stage
    new_indi <- data.frame(
      size = rnorm(sdl_to_plant, mean = params$sdl_d_int, sd = params$sdl_d_sd),
      shading = rpois(sdl_to_plant, env_params$shading) #,
      # slope = rgamma(sdl_to_plant, (slope_mean^2)/(slope_sd^2), (slope_mean)/(slope_sd^2))
    )
    
    # Set next steps individuals dataframe  ------------------------------------------------------
    indi <- rbind(indi_z1 %>% select(z1, shading) %>% #, slope) %>% 
                    filter(!is.na(z1)) %>% rename(size = z1),
                  new_indi)
    
    
    return(list(indi = indi,
                sdl = sdl,
                sb1 = sb1,
                sb2 = sb2,
                yr = yr + 1))
    
  
}

ibm_ext_p <- function(i, df_env, params,
                      clim_ts, clim_hist,
                      pop_vec, sdl_n,
                      n_it) {
  message(i)
  
  # Simulation settings
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  
  # slope_mean <- df_env$mean_slope[i]
  # slope_sd <- df_env$sd_slope[i]
  
  if(model == "No change") {
    
    clim_mod <- hist_clim[grep(locality, names(hist_clim), value = T, ignore.case = T)]
    
  } else{
    
    txt <- paste(locality, model, scenario, sep = '.*')
    clim_mod <- clim_ts[grep(txt, names(clim_ts), value = T)] %>%
      lapply(., function(x) x$value)
    
  }
  
  # Environment parameters
  env_params <- append(
    clim_mod,
    list(lags = lag,
         shading = df_env$shading[i] #,
         # slope = df_env$slope[i]
    ))
  
  
  # Select pop_vectors
  pop_vec <- pop_vec[[grep(toupper(locality), names(pop_vec), value = T)]] 
  sdl <- sdl_n[[grep(toupper(locality), names(sdl_n), value = T)]]
  
  sb1 = 2000
  sb2 = 1000
  
  # set up starting shading & slope
  pop_shading <- rbinom(n = length(pop_vec), 20, prob = (env_params$shading/20))

  # Set up start population
  indi <- data.frame(
    size = pop_vec,
    shading = pop_shading
  )
  
  pop_size <- data.frame(
    yr = c(2021:2100),
    n_plants = c(nrow(indi), rep(NA, 79)),
    n_sdl = c(sdl, rep(NA, 79)),
    n_sb1 = c(sb1, rep(NA, 79)),
    n_sb2 = c(sb2, rep(NA, 79))
  )

  ## IBM
  yr <- 1
  
  while (yr < 80 & nrow(indi) > 0){
  
    ibm_yr <- yearly_loop(yr, env_params, locality,
                          indi, sdl, sb1, sb2)
    
    yr <- ibm_yr$yr
    indi <- ibm_yr$indi
    sdl <- ibm_yr$sdl
    sb1 <- ibm_yr$sb1
    sb2 <- ibm_yr$sb2
    
    pop_size$n_plants[yr] <- nrow(indi)
    pop_size$n_sdl[yr] <- sdl
    pop_size$n_sb1[yr] <- sb1
    pop_size$n_sb2[yr] <- sb2
    
  }
  
  b <- data.frame(
    locality = locality,
    model = model,
    scenario = scenario,
    shading = df_env$shading[i],
    below_ext_size = any(pop_size$n_plants < 10), 
    yr_of_ext = ifelse(any(pop_size$n_plants < 10), min(which(pop_size$n_plants < 10)), NA)
  )
  b$pop_size = list(list(plants = pop_size$n_plants,
                    sdl = pop_size$n_sdl,
                    sb1 = pop_size$n_sb1,
                    sb2 = pop_size$n_sb2))
  
  return(b)
  
}


man_trans <- function(i, df_env, params,
                      clim_ts,
                      pop_vec, sdl_n,
                      n_it,
                      U, L, n) {

  print(i)
  
  # Simulation settings
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  
  yrs_between_man <- df_env$yrs_between_man[i]
  effort <- df_env$effort[i]
  
  # slope_mean <- df_env$mean_slope[i]
  # slope_sd <- df_env$sd_slope[i]
  
  txt <- paste(locality, model, scenario, sep = '.*')
  clim_mod <- clim_ts[grep(txt, names(clim_ts), value = T)] %>%
    lapply(., function(x) x$value)
  
  
  # Environment parameters
  env_params <- append(
    clim_mod,
    list(lags = lag,
         shading = df_env$shading[i] #,
         # slope = df_env$slope[i]
    ))
  
  
  # Select pop_vectors
  pop_vec <- pop_vec[[grep(toupper(locality), names(pop_vec), value = T)]] 
  sdl <- sdl_n[[grep(toupper(locality), names(sdl_n), value = T)]]
  
  sb1 = 2000
  sb2 = 1000
  
  # set up starting shading & slope
  pop_shading <- rpois(length(pop_vec), env_params$shading)
  # pop_slope <- rgamma(length(pop_vec), (slope_mean^2)/(slope_sd^2), (slope_mean)/(slope_sd^2))
  
  # Set up management
  management_yrs <- seq(from = 1, to = 80, by = yrs_between_man)
  
  # Set up start population
  indi <- data.frame(
    size = pop_vec,
    shading = pop_shading#,
    # slope = pop_slope
  )
  
  pop_size <- data.frame(
    yr = c(2021:2100),
    n_plants = c(nrow(indi), rep(NA, 79)),
    n_sdl = c(sdl, rep(NA, 79)),
    n_sb1 = c(sb1, rep(NA, 79)),
    n_sb2 = c(sb2, rep(NA, 79))
  )
  
  ## IBM
  yr <- 1
  
  while (yr < 80 & nrow(indi) > 0){
    
    if(yr %in% management_yrs) {
      transplant_indi <- data.frame(
        size = rnorm(effort, mean = params$sdl_d_int, sd = params$sdl_d_sd),
        shading =  rpois(effort, env_params$shading)#,
        # slope = rgamma(effort, (slope_mean^2)/(slope_sd^2), (slope_mean)/(slope_sd^2))
      )
      indi <- rbind(indi, transplant_indi)
    }
    
    ibm_yr <- yearly_loop(yr, env_params, locality,
                          indi, sdl, sb1, sb2) #,
                          # slope_mean, slope_sd)
    
    
    yr <- ibm_yr$yr
    indi <- ibm_yr$indi
    sdl <- ibm_yr$sdl
    sb1 <- ibm_yr$sb1
    sb2 <- ibm_yr$sb2
    
    pop_size$n_plants[yr] <- nrow(indi)
    pop_size$n_sdl[yr] <- sdl
    pop_size$n_sb1[yr] <- sb1
    pop_size$n_sb2[yr] <- sb2
    
  }

  b <- data.frame(
    locality = locality,
    model = model,
    scenario = scenario,
    shading = df_env$shading[i],
    # slope = df_env$slope[i],
    effort = df_env$effort[i],
    yrs_between_man = yrs_between_man,
    below_ext_size = any(pop_size$n_plants < 10), 
    yr_of_ext = ifelse(any(pop_size$n_plants < 10), min(which(pop_size$n_plants < 10)), NA)
  )
  b$pop_size = list(list(plants = pop_size$n_plants,
                         sdl = pop_size$n_sdl,
                         sb1 = pop_size$n_sb1,
                         sb2 = pop_size$n_sb2))


  return(b)


}


man_seedadd <- function(i, df_env, params,
                      clim_ts,
                      pop_vec, sdl_n,
                      n_it,
                      U, L, n) {
  
  print(i)
  
  # Simulation settings
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  
  yrs_between_man <- df_env$yrs_between_man[i]
  effort <- df_env$effort[i]
  
  # slope_mean <- df_env$mean_slope[i]
  # slope_sd <- df_env$sd_slope[i]
  
  txt <- paste(locality, model, scenario, sep = '.*')
  clim_mod <- clim_ts[grep(txt, names(clim_ts), value = T)] %>%
    lapply(., function(x) x$value)
  
  
  # Environment parameters
  env_params <- append(
    clim_mod,
    list(lags = lag,
         shading = df_env$shading[i],
         slope = df_env$slope[i]
    ))
  
  
  # Select pop_vectors
  pop_vec <- pop_vec[[grep(toupper(locality), names(pop_vec), value = T)]] 
  sdl <- sdl_n[[grep(toupper(locality), names(sdl_n), value = T)]]
  
  sb1 = 2000
  sb2 = 1000
  
  # set up starting shading & slope
  pop_shading <- rpois(length(pop_vec), env_params$shading)
  # pop_slope <- rgamma(length(pop_vec), (slope_mean^2)/(slope_sd^2), (slope_mean)/(slope_sd^2))
  
  # Set up management
  management_yrs <- seq(from = 1, to = 80, by = yrs_between_man)
  
  # Set up start population
  indi <- data.frame(
    size = pop_vec,
    shading = pop_shading#,
    # slope = pop_slope
  )
  
  pop_size <- data.frame(
    yr = c(2021:2100),
    n_plants = c(nrow(indi), rep(NA, 79)),
    n_sdl = c(sdl, rep(NA, 79)),
    n_sb1 = c(sb1, rep(NA, 79)),
    n_sb2 = c(sb2, rep(NA, 79))
  )
  
  ## IBM
  yr <- 1
  
  while (yr < 80 & nrow(indi) > 0){
    
    if(yr %in% management_yrs) {
      sb1 <- sb1 + df_env$effort
    }
    
    ibm_yr <- yearly_loop(yr, env_params, locality,
                          indi, sdl, sb1, sb2) #,
                          # slope_mean, slope_sd)
    
    
    yr <- ibm_yr$yr
    indi <- ibm_yr$indi
    sdl <- ibm_yr$sdl
    sb1 <- ibm_yr$sb1
    sb2 <- ibm_yr$sb2
    
    pop_size$n_plants[yr] <- nrow(indi)
    pop_size$n_sdl[yr] <- sdl
    pop_size$n_sb1[yr] <- sb1
    pop_size$n_sb2[yr] <- sb2
    
  }
  
  b <- data.frame(
    locality = locality,
    model = model,
    scenario = scenario,
    shading = df_env$shading[i],
    # slope = df_env$slope[i],
    effort = df_env$effort[i],
    yrs_between_man = yrs_between_man,
    below_ext_size = any(pop_size$n_plants < 10), 
    yr_of_ext = ifelse(any(pop_size$n_plants < 10), min(which(pop_size$n_plants < 10)), NA)
  )
  b$pop_size = list(list(plants = pop_size$n_plants,
                         sdl = pop_size$n_sdl,
                         sb1 = pop_size$n_sb1,
                         sb2 = pop_size$n_sb2))
  
  
  return(b)
  
  
}