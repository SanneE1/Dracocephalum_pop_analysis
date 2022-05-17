
transform_and_save_climate_data <- function(location_clim_xlsx_file) {
  
  save_location <- "data/climate_data_Dobrichovice.csv"
  
  clim <- readxl::read_xlsx(location_clim_xlsx_file) 
  #for some reason Mesíc was causing trouble when using dplyr::rename, but minimal example worked fine, not sure where it came from but this works too 
  names(clim)[1:3] <- c("Year", "Month", "Day")
  
  clim %>%
    mutate(across(c("Tavg", "Tmin", "Tmax", "Prec_Dobrichovice"), ~gsub(pattern = ",", replacement = ".",.) %>% as.numeric())) %>%
    mutate(GDD_5_contribution = ((Tmin + Tmax)/2) - 5,
           GDD_10_contribution = ((Tmin + Tmax)/2) - 10) %>%
    group_by(Year, Month) %>%
    summarise(Tavg = mean(Tavg),
              Tmin = mean(Tmin),
              Tmax = mean(Tmax),
              Prec = sum(Prec_Dobrichovice),
              GDD_5 = sum(GDD_5_contribution),
              GDD_10 = sum(GDD_10_contribution)) %>%
    ## The next line makes the scaled variables across months
    ## In other words: Is it a cold or warm Januari compared to other Januaries, 
    ## rather than is this Januari cold across all the time series
    ungroup() %>% group_by(Month) %>%  
    mutate(Tavg_scaled = scale(Tavg),
           Tmin_scaled = scale(Tmin),
           Tmax_scaled = scale(Tmax),
           Prec_scaled = scale(Prec)) %>%
    write.csv(., file = save_location, row.names = F)
  
  return(save_location)
  
}
