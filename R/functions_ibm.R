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
    slope = indi$slope,
    rock = indi$rock,
    soil_depth = indi$soil_depth,
    year_t0 = rep(2016, n),   ## just a dummy, will be ignored when predicting below
    lags = matrix(clim$lags, ncol = (lag+1), nrow = n, byrow = T),
    tas_scaledcovar = matrix(clim$temp, ncol = (lag+1), nrow = n, byrow = T),
    pr_scaledcovar = matrix(clim$precip, ncol = (lag+1), nrow = n, byrow = T),
    pet_scaledcovar = matrix(clim$pet, ncol = (lag+1), nrow = n, byrow = T)
  )
  
  pt <- mgcv::predict.gam(model,newdata = new_data, 
                          exclude = 's(year_t0)',
                          type = "response")
  return(pt)
}

yearly_loop <- function(yr, env_params, locality,
                        indi, sdl, sb1, sb2){ 
  
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
    
    indi_z1 <- indi %>%
      mutate(surv = surv, 
             z1 = z1)
    
    # calculate total number of seeds produced
    flowp <- rbinom(n = pop_size, prob = mod_pred(indi = indi, model = params$pflower_mod,
                                                  params = params, env_params = env_params, 
                                                  locality = locality, clim = clim), 
                    size = 1)
    
    seedp <- rep(NA, pop_size)
    seedp[which(flowp == 1)] <- rbinom(n = length(which(flowp == 1)), 
                                     prob = mod_pred(indi = indi[which(flowp == 1),], 
                                                     model = params$pseed_mod,
                                                     params = params, env_params = env_params, 
                                                     locality = locality, clim = clim), 
                                     size = 1)
    
    nseeds <- rep(NA, pop_size)
    seed_mean <- mod_pred(indi = indi[which(seedp == 1),], model = params$nseed_mod,
                          params = params, env_params = env_params, 
                          locality = locality, clim = clim)
    
    seed_gamma_shape1 <- (seed_mean^2)/(params$nseed_sd^2)
    seed_gamma_shape2 <- (seed_mean)/(params$nseed_sd^2)
    nseeds[which(seedp == 1)] <- round(
      rgamma(n = length(which(seedp == 1)), shape = seed_gamma_shape1, scale = seed_gamma_shape2),
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
      shading = rpois(sdl_to_plant, env_params$shading),
      slope = CaDENCE::rbgamma(n = sdl_to_plant, prob = params$slope$prob_non_zero[which(params$slope$population == toupper(locality))],
                               scale = params$slope$gamma_scale[which(params$slope$population == toupper(locality))], 
                               shape = params$slope$gamma_shape[which(params$slope$population == toupper(locality))]),
      rock = CaDENCE::rbgamma(n = sdl_to_plant, prob = params$rock$prob_non_zero[which(params$rock$population == toupper(locality))],
                              scale = params$rock$gamma_scale[which(params$rock$population == toupper(locality))], 
                              shape = params$rock$gamma_shape[which(params$rock$population == toupper(locality))]),
      soil_depth = CaDENCE::rbgamma(n = sdl_to_plant, prob = params$soil_depth$prob_non_zero[which(params$soil_depth$population == toupper(locality))],
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

ibm_ext_p <- function(i, df_env, params,
                      clim_ts, clim_hist,
                      pop_vec, sdl_n,
                      n_it) {
  # gc()
  message(i)
  
  # Simulation settings
  scenario <- df_env$scenario[i]
  model <- df_env$model[i]
  locality <- df_env$localities[i]
  
  
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
         shading = df_env$shading[i]
         ))
  
  
  # Select pop_vectors
  pop_vec <- pop_vec[[grep(toupper(locality), names(pop_vec), value = T)]] 
  sdl <- sdl_n[[grep(toupper(locality), names(sdl_n), value = T)]]
  
  sb1 = 2000
  sb2 = 1000
  
  # set up starting shading & slope
  pop_shading <- rbinom(n = length(pop_vec), 20, prob = (env_params$shading/20))
  pop_slope <- CaDENCE::rbgamma(n = length(pop_vec), prob = params$slope$prob_non_zero[which(params$slope$population == toupper(locality))],
                                scale = params$slope$gamma_scale[which(params$slope$population == toupper(locality))], 
                                shape = params$slope$gamma_shape[which(params$slope$population == toupper(locality))])
  
  pop_rock <- CaDENCE::rbgamma(n = length(pop_vec), prob = params$rock$prob_non_zero[which(params$rock$population == toupper(locality))],
                               scale = params$rock$gamma_scale[which(params$rock$population == toupper(locality))], 
                               shape = params$rock$gamma_shape[which(params$rock$population == toupper(locality))])
  
  pop_sd <- CaDENCE::rbgamma(n = length(pop_vec), prob = params$soil_depth$prob_non_zero[which(params$soil_depth$population == toupper(locality))],
                             scale = params$soil_depth$gamma_scale[which(params$soil_depth$population == toupper(locality))], 
                             shape = params$soil_depth$gamma_shape[which(params$soil_depth$population == toupper(locality))])

  # Set up start population
  indi <- data.frame(
    size = pop_vec,
    shading = pop_shading,
    slope = pop_slope,
    rock = pop_rock,
    soil_depth = pop_sd
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



calc_stats <- function(data, variable) {
  data %>%
    group_by(population) %>%
    summarise(
      prob_non_zero = mean(get(variable) != 0, na.rm = T),
      gamma_shape = ifelse(sum(!is.na(get(variable)) & get(variable) > 0) > 1, 
                           MASS::fitdistr(get(variable)[!is.na(get(variable)) & get(variable) > 0], "gamma")$estimate["shape"], NA),
      gamma_scale = ifelse(sum(!is.na(get(variable)) & get(variable) > 0) > 1, 
                           1 / MASS::fitdistr(get(variable)[!is.na(get(variable)) & get(variable) > 0], "gamma")$estimate["rate"], NA)
    )
}
