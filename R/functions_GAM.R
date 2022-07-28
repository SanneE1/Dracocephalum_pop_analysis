

climate_wider_for_gam <- function(clim_data, demo_data, variables, response_t1 = T,
                                  lag,
                                  censusMonth = 7) {
  
  ## For which years do we need the climate
  year_t0 <- unique(demo_data$year_t0)
  ## We want to get climate starting from the month of census where response is measured
  ## So if response happens in t1, we need to add a year to the year_t0
  ## However,to later link this back to year_t0 which is in the demographic data we use year_t0 for the dataframe column
  if(response_t1) {
    responseYear <- year_t0 + 1
  } else {
    responseYear <- year_t0
  }
  
  ## Arrange climate data chronologically "backwards"/ descending
  desc_clim <- clim_data %>% dplyr::arrange(desc(year), desc(month))
  
  lagged_climate <- data.frame(year_t0 = year_t0,
                               responseYear = responseYear)
  
  for(j in variables) {
    ## Set up lag dataframe
    ## ncol: lag starts in census month (0) and need extra column for year_t0 and responseYear
    clim_wide <- data.frame(matrix(nrow = length(year_t0), ncol = (as.integer(lag) + 3)))  
    names(clim_wide) <- c("year_t0", "responseYear", paste(j, c(0:as.integer(lag)), sep = "."))
    clim_wide$year_t0 <- year_t0
    clim_wide$responseYear <- responseYear
    
    # fill in lag values for each of the observation years
    for(i in responseYear) {
      n <- which(desc_clim$year == i & desc_clim$month == censusMonth)
      
      if(length(n)==0) next
      
      clim_wide[which(clim_wide$responseYear == i),-c(1,2)] <- desc_clim[c(n:(n+lag)), j]
    }
    
    clim_wide <- clim_wide %>% filter(complete.cases(.[,-c(1,2)]))
    
    ## define and enter covariates the way the fitting function wants 
    ## -- See Teller et al. 2016 Rcode/Climate/AuthenticAnalysis/Utilities/fetchDemoData.R L223-233
    vars <- which(names(clim_wide) %in% names(clim_wide %>% dplyr::select(contains(j))))
    clim_wide$a <- as.matrix(clim_wide[,vars])
    names(clim_wide)[which(names(clim_wide)=="a")] <- paste0(j,"covar")
    
    lagged_climate <- left_join(lagged_climate, clim_wide) %>% dplyr::select(-responseYear)
  }
  
  ## define and enter covariates the way the fitting function wants 
  ## -- See Teller et al. 2016 Rcode/Climate/AuthenticAnalysis/Utilities/fetchDemoData.R L223-233
  lags <- matrix(0,nrow(lagged_climate),length(vars)) 
  for(k in 1:ncol(lags)) lags[,k]=k 
  lagged_climate$lags=as.matrix(lags)
  
  
  if(any(is.na(lagged_climate))) warning("NA's have been filter out of the lagged climate dataframe. Most likely, climate data provided doesn't cover full demography periods + lag period")
  
  
  return(lagged_climate)
}



plot_spline_coeff <- function(best_model, 
                              lag = lag,
                              tas = F, pr = F, pet = F,
                              shade = F, slope = F, rock = F, soil = F,
                              vital_rate, save_plot = T
                              ) {
  
   lags= c(0:lag); ln_stems_t0=1; population = factor("CR")
   year_t0 = 2015

   tas_scaledcovar= 0; pr_scaledcovar = 0; pet_scaledcovar = 0
   tot_shading_m = 0; slope_m = 0; rock_m = 0; soil_m= 0
   interact = ""; covar=""
   
   if(tas) {
     tas_scaledcovar = 1
     yaxis_title = "temperature anomaly coefficient"
     covar = "tas"}
   if(pr) {
     pr_scaledcovar = 1
     yaxis_title = "precipitation anomaly coefficient"
     covar = "pr"}
   if(pet) {
     pet_scaledcovar = 1
     yaxis_title = "potential evapotranspiration anomaly coefficient"
     covar = "pet"}
   
   if(shade) {
     tot_shading_m = c(0,3,6,9)
     legendtitle <- "shading level"
     interact = "shading_m"
   }
   if(slope) {
     slope_m = c(0,25,50,75)
     legendtitle <- "slope (degree)"
     interact = "slope_m"
   }
   if(rock) {
     rock_m = c(0, 33, 66, 99)
     legendtitle <- "rock cover (%)"
     interact = "rock_m"
   }
   if(soil) {
     soil_m = c(0,3,6,9)
     legendtitle <- "soil depth (cm)"
     interact = "soil_m"
   }
   
   getBeta.data <- expand.grid(lags= c(0:lag), 
                               ln_stems_t0=1,
                               population = factor("CR"),
                               year_t0 = 2015,
                               soil_depth = 1,
                               rock = 20,
                               slope = 20,
                               tot_shading_t0 = 3,
                               tas_scaledcovar = tas_scaledcovar, 
                               pr_scaledcovar = pr_scaledcovar, 
                               pet_scaledcovar = pet_scaledcovar,
                               tot_shading_m = tot_shading_m,
                               slope_m = slope_m,
                               rock_m = rock_m,
                               soil_m = soil_m)
   terms.data <- mgcv::predict.gam(best_model, newdata=getBeta.data, type="terms", se=TRUE)
   betas_surv_temp <- data.frame(lag = getBeta.data$lags,
                                interact = getBeta.data[,grep(interact, colnames(getBeta.data))],
                                beta= terms.data[[1]][,grep(paste0(interact, ".*", covar), colnames(terms.data[[1]]), value = T)],
                                se= terms.data[[2]][,grep(paste0(interact, ".*", covar), colnames(terms.data[[2]]), value = TRUE)])

  plot <- ggplot() +
    geom_ribbon(data = betas_surv_temp, aes(x = lag * -1, ymin = beta - se, ymax = beta + se, fill = as.factor(interact)), alpha = 0.2) +
    geom_line(data = betas_surv_temp, aes(x = lag * -1, y = beta, colour = as.factor(interact))) +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = seq(from = -1 * (lag), to = 0, by = 6),
                       limits = c((lag * -1), 0),
                       expand = c(0,0)) +
    scale_fill_discrete(name = legendtitle) +
    scale_colour_discrete(name = legendtitle) +
    xlab("months") + ylab(yaxis_title) +
    theme(legend.position = "bottom")
  

  if(save_plot == T){
  ggsave(plot = plot,
    filename = here::here("results", paste0(vital_rate, "_spline.png")),
         width = 7, height = 5, units = "in", dpi = 300, type = "cairo")
  }
  
  return(plot)
}
