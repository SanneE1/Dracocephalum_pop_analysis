

climate_wider_for_gam <- function(clim_data, demo_data, variables, response_t1 = T,
                                  lag,
                                  censusMonth = 7) {
  
  localities <- unique(clim_data$locality)
  
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
  desc_clim_all <- clim_data %>% dplyr::arrange(desc(year), desc(month))
  
  climate <- data.frame(year_t0 = year_t0,
                        responseYear = responseYear)
  lagged_climate <- c()
  
  for(k in localities) {
    desc_clim <- desc_clim_all %>% filter(locality == k)
    
    for(j in variables) {
      ## Set up lag dataframe
      ## ncol: lag starts in census month (0) and need extra column for year_t0 and responseYear
      clim_wide <- data.frame(matrix(nrow = length(year_t0), ncol = (as.integer(lag) + 4)))  
      names(clim_wide) <- c("year_t0", "responseYear", "population", paste(j, c(0:as.integer(lag)), sep = "."))
      clim_wide$year_t0 <- year_t0
      clim_wide$responseYear <- responseYear
      clim_wide$population <- toupper(k)
      
      # fill in lag values for each of the observation years
      for(i in responseYear) {
        n <- which(desc_clim$year == i & desc_clim$month == censusMonth)
        
        if(length(n)==0) next
        
        clim_wide[which(clim_wide$responseYear == i),-c(1:3)] <- desc_clim[c(n:(n+lag)), j]
      }
      
      ## define and enter covariates the way the fitting function wants 
      ## -- See Teller et al. 2016 Rcode/Climate/AuthenticAnalysis/Utilities/fetchDemoData.R L223-233
      vars <- which(names(clim_wide) %in% names(clim_wide %>% dplyr::select(contains(j))))
      clim_wide$a <- as.matrix(clim_wide[,vars])
      names(clim_wide)[which(names(clim_wide)=="a")] <- paste0(j,"covar")
      
      if(j == variables[1]){
        y <- left_join(climate, clim_wide) %>% dplyr::select(-responseYear)
      } else {
        y <- cbind(y, clim_wide[,-c(1:3)])
      }
      
    }
    
    lagged_climate <- rbind(lagged_climate, y)
  }
  ## define and enter covariates the way the fitting function wants 
  ## -- See Teller et al. 2016 Rcode/Climate/AuthenticAnalysis/Utilities/fetchDemoData.R L223-233
  lags <- matrix(0,nrow(lagged_climate),length(vars)) 
  for(k in 1:ncol(lags)) lags[,k]=k 
  lagged_climate$lags=as.matrix(lags)
  
  return(lagged_climate)
}



plot_spline_coeff <- function(best_model, 
                              lag = 24,
                              tas = F, pr = F, pet = F,
                              shade = F, ro = F, sl = F, sd = F,
                              vital_rate, 
                              save_plot = F
) {
  
  lags= c(0:lag); ln_stems_t0=1; population = factor("CR")

  tas_scaledcovar= 0; pr_scaledcovar = 0; pet_scaledcovar = 0
  tot_shading_t0 = 0; rock = 0; slope = 0; soil_depth = 0
  covar=""; legendtitle = ""
  
  xaxis_title = "months"
  
  if(shade) {
    tot_shading_t0 = c(0:10)
    yaxis_title = "Coefficient estimates"
    xaxis_title = "Shading level"
    covar = "tot_shading_t0"
  }
  if(ro) {
    rock = c(0:70)
    yaxis_title = "Coefficient estimates"
    xaxis_title = "Rock %"
    covar = "rock"
  }
  if(sl) {
    slope = c(0:70)
    yaxis_title = "Coefficient estimates"
    xaxis_title = "Local slope"
    covar = "slope"
  }
  if(sd) {
    soil_depth = c(0:10)
    yaxis_title = "Coefficient estimates"
    xaxis_title = "Soil Depth"
    covar = "soil_depth"
  }
  
  if(tas) {
    tas_scaledcovar = 1
    yaxis_title = "temperature \n coefficient estimates"
    covar = "tas"}
  if(pr) {
    pr_scaledcovar = 1
    yaxis_title = "precipitation \n coefficient estimates"
    covar = "pr"}
  if(pet) {
    pet_scaledcovar = 1
    yaxis_title = "PET \n coefficient estimates"
    covar = "pet"}
  
  
  getBeta.data <- expand.grid(lags= c(0:lag), 
                              ln_stems_t0=1,
                              population = factor("CR"),
                              year_t0 = 2018,
                              soil_depth = soil_depth,
                              rock = rock,
                              slope = slope,
                              tot_shading_t0 = tot_shading_t0,
                              tas_scaledcovar = tas_scaledcovar, 
                              pr_scaledcovar = pr_scaledcovar, 
                              pet_scaledcovar = pet_scaledcovar)
  
  terms.data <- mgcv::predict.gam(best_model, newdata=getBeta.data, type="terms", se=TRUE)
  
  if(tas|pet|pr){
    xaxis = (-1 * getBeta.data$lags)
  } else if(shade){
    xaxis = getBeta.data$tot_shading_t0
  } else if(ro){
    xaxis = getBeta.data$rock
  } else if(sl){
    xaxis = getBeta.data$slope
  } else if(sd){
    xaxis = getBeta.data$soil_depth
  }
  
  
  betas_surv_temp <- data.frame(x = xaxis,
                                beta= terms.data[[1]][,grep(covar, colnames(terms.data[[1]]), value = T)],
                                se= terms.data[[2]][,grep(covar, colnames(terms.data[[2]]), value = TRUE)])
  
  plot <- ggplot() +
    geom_ribbon(data = betas_surv_temp, aes(x = x, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_surv_temp, aes(x = x, y = beta)) +
    geom_hline(aes(yintercept = 0)) +
    viridis::scale_fill_viridis(name = legendtitle, option = "D", discrete = T, direction = -1) +  
    viridis::scale_color_viridis(name = legendtitle, option = "D", discrete = T, direction = -1) +  
    xlab(xaxis_title) + ylab(yaxis_title) + ggtitle(vital_rate) +
    theme(legend.position = "bottom")
 
  if(tas|pr|pet) {
    plot <- plot +
      scale_x_continuous(breaks = seq(from = -1 * (lag), to = 0, by = 6),
                         limits = c((lag * -1), 0),
                         expand = c(0,0)) 
  } 
  
  if(save_plot == T){
    ggsave(plot = plot,
           filename = here::here("results", paste0(vital_rate, "_spline.png")),
           width = 7, height = 5, units = "in", dpi = 300, type = "cairo")
  }
  
  return(plot)
}

generate_combinations <- function(terms) {
  # Create an empty list to store the combinations
  combinations <- list()
  
  # Loop over the range of possible lengths (1 to length of terms)
  for (i in 1:length(terms)) {
    # Generate all combinations of length i
    comb <- combn(terms, i, simplify = FALSE)
    comb <- lapply(comb, function(x) paste(x, collapse = "+"))
    # Append these combinations to the list
    combinations <- c(combinations, comb)
  }
  
  return(combinations)
}
