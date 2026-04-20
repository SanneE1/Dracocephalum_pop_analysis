# Extract values from future climate models for populations

library(raster)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(patchwork)
# Function to process a single file
format_future_ts_climate <- function(file) {
  # message(paste("Processing file:", basename(file)))
  
  tryCatch({
    
    opt <- stringr::str_split(file, "/")[[1]]
    opt <- stringr::str_split(opt[3], "_")[[1]]
    
    nc <- raster::brick(file)
    
    df_extr <- raster::extract(nc, populations[,c("Longitude", "Latitude")], method="bilinear", df = TRUE) 
    
    df <- pivot_longer(df_extr, cols = starts_with("X"), 
                       names_to = c("year", "month"), names_pattern = "X(\\d*)[[:punct:]](\\d*)[[:punct:]].*") %>%
      mutate(variable = opt[1],
             model = opt[3],
             scenario = opt[4])
    
    suppressMessages(df <- full_join(populations, df))
    
    return(df)
    
  },
  error = function(cond) {
    message(paste("there's a problem with file:", file))
    message(conditionMessage(cond))
    NULL})
  
}

files <- list.files(path = "data/future_climat_nc/", pattern = ".nc", recursive = T, full.names = T)

populations <- read.csv("data/Draco_locations.csv")[c(5:8),] %>%
  dplyr::select(-c("X", "Name")) %>%
  rename(locality = Station.ID.Draco.locality.ID) %>%
  tibble::rowid_to_column(var = "ID")

hist_clim <- read.csv("data/CHELSA_data.csv") %>%
  mutate(date = make_date(year, month, 1))

hist_scale <- hist_clim  %>%
  group_by(locality, month) %>%
  summarise(pr_mean = mean(pr, na.rm = T),
            pr_sd = sd(pr, na.rm = T),
            tas_mean = mean(tas, na.rm = T),
            tas_sd = sd(tas, na.rm = T),
            pet_mean = mean(pet, na.rm = T),
            pet_sd = sd(pet, na.rm = T)) 

climA <- lapply(files, format_future_ts_climate) %>%
  bind_rows %>%
  mutate(value = case_when(variable == "pr" ~ value * 86400 * lubridate::days_in_month(as.Date(paste("01", month, year, sep = "/"), "%d/%m/%Y")),               ## pr unit is kg/m2/s and the pet calculations needs kg/m2/month
                           variable == "tas" ~ value - 273.15),
         year = as.integer(year),
         month = as.integer(month)) %>%
  pivot_wider(values_from = value, names_from = variable) %>% 
  arrange(year, month) %>%
  group_by(locality, model, scenario) %>%
  group_split()

clim <- lapply(climA, function(x) {x %>%  mutate(pet = as.numeric(SPEI::thornthwaite(Tave = tas, lat = unique(Latitude), na.rm = T, verbose = F)))}) %>%
  bind_rows() %>%
  mutate(date = make_date(year, month, 1)) %>%
  left_join(., hist_scale) %>%
  rowwise() %>%
  mutate(pr_scaled = (pr-pr_mean)/pr_sd,
         tas_scaled = (tas - tas_mean)/tas_sd,
         pet_scaled = (pet - pet_mean)/pet_sd)


chelsa_clim <- read.csv("data/CHELSA_future_ts_formatted.csv") %>%
  mutate(date = make_date(year, month, 1)) %>%
  filter(scenario != "historical")

station_clim <- read.csv("data/climate_data_Dobrichovice.csv") %>%
  mutate(date = make_date(Year, Month, 1))


chelsahist_vs_cmip5 <- ggplot() +
  geom_line(data = clim, aes(x = date, y = pr, colour = interaction(model, scenario)), alpha = 0.2) +
  geom_line(data = hist_clim, aes(x = date, y = pr), colour = "black") +
  geom_smooth(data = clim, aes(x = date, y = pr, colour = interaction(model, scenario))) +
  geom_smooth(data = hist_clim, aes(x = date, y = pr), colour = "black") + 
  coord_cartesian(xlim = as.Date(c("1970-01-01", "2050-12-31"))) +
  theme(legend.position = "none") + ggtitle("CHELSA historical data (black) vs \nCMIP5 models")

chelsahist_vs_chelsafut <- ggplot() +  # Yes, this is also based on the cmip5 but had to call it something
  geom_line(data = chelsa_clim, aes(x = date, y = precip, colour = interaction(model, scenario)), alpha = 0.2) +
  geom_line(data = hist_clim, aes(x = date, y = pr), colour = "black") +
  geom_smooth(data = chelsa_clim, aes(x = date, y = precip, colour = interaction(model, scenario))) +
  geom_smooth(data = hist_clim, aes(x = date, y = pr), colour = "black") +
  coord_cartesian(xlim = as.Date(c("1970-01-01", "2050-12-31"))) +
  theme(legend.position = "none") + ggtitle("CHELSA historical data (black) vs \nCHELSA's downscaled CMIP5 forecasts")

chelsahist_vs_localstation <- ggplot() +
  geom_line(data = hist_clim, aes(x = date, y = pr), colour = "black") +
  geom_line(data = station_clim, aes(x = date, y = Prec, colour = "red")) +
  geom_smooth(data = station_clim, aes(x = date, y = Prec, colour = "red")) +
  geom_smooth(data = hist_clim, aes(x = date, y = pr), colour = "black") +
  coord_cartesian(xlim = as.Date(c("1970-01-01", "2050-12-31"))) +
  theme(legend.position = "none") + ggtitle("CHELSA historical data (black) vs \nDobrichovice climate station")

chelsahist_vs_cmip5 + chelsahist_vs_chelsafut + chelsahist_vs_localstation

chelsahist_vs_cmip52 <- ggplot() +
  geom_line(data = clim, aes(x = date, y = tas, colour = interaction(model, scenario)), alpha = 0.2) +
  geom_line(data = hist_clim, aes(x = date, y = tas), colour = "black") +
  geom_smooth(data = clim, aes(x = date, y = tas, colour = interaction(model, scenario))) +
  geom_smooth(data = hist_clim, aes(x = date, y = tas), colour = "black") + 
  coord_cartesian(xlim = as.Date(c("1970-01-01", "2050-12-31"))) +
  theme(legend.position = "none") + ggtitle("CHELSA historical data (black) vs \nCMIP5 models")

chelsahist_vs_chelsafut2 <- ggplot() +  # Yes, this is also based on the cmip5 but had to call it something
  geom_line(data = chelsa_clim, aes(x = date, y = tas, colour = interaction(model, scenario)), alpha = 0.2) +
  geom_line(data = hist_clim, aes(x = date, y = tas), colour = "black") +
  geom_smooth(data = chelsa_clim, aes(x = date, y = tas, colour = interaction(model, scenario))) +
  geom_smooth(data = hist_clim, aes(x = date, y = tas), colour = "black") +
  coord_cartesian(xlim = as.Date(c("1970-01-01", "2050-12-31"))) +
  theme(legend.position = "none") + ggtitle("CHELSA historical data (black) vs \nCHELSA's downscaled CMIP5 forecasts")

chelsahist_vs_localstation2 <- ggplot() +
  geom_line(data = hist_clim, aes(x = date, y = tas), colour = "black") +
  geom_line(data = station_clim, aes(x = date, y = Tavg, colour = "red")) +
  geom_smooth(data = station_clim, aes(x = date, y = Tavg, colour = "red")) +
  geom_smooth(data = hist_clim, aes(x = date, y = tas), colour = "black") +
  coord_cartesian(xlim = as.Date(c("1970-01-01", "2050-12-31"))) +
  theme(legend.position = "none") + ggtitle("CHELSA historical data (black) vs \nDobrichovice climate station")

chelsahist_vs_cmip52 + chelsahist_vs_chelsafut2 + chelsahist_vs_localstation2




# write.csv(clim, "cds_future_ts_formatted.csv")

