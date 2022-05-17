

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
  desc_clim <- clim_data %>% dplyr::arrange(desc(Year), desc(Month))
  
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
      n <- which(desc_clim$Year == i & desc_clim$Month == censusMonth)
      
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



