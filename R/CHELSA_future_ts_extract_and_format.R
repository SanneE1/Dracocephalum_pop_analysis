# This script extracts the climate data from the CHELSA cmip5 time series data from local source.
# As I haven't been able to figure out how to extract from a .nc ts format without fully 
# downloading the files first, download the .nc files from this link first 
# (https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V1%2Fchelsa_cmip5_ts)
# then change the location in L18 to your local folder

# module load foss/2019b R/4.0.0-2

##------------------------------------------------------
## Future Climate
##------------------------------------------------------


### On EVE (The UFZ HPC)
library(ncdf4)
library(raster)
library(tidyverse)
library(lubridate)

gps <- read.csv("/data/lagged/Dracocephalum/Draco_locations.csv")

files <- list.files("/data/lagged/CHELSA/", pattern = "CHELSAcmip5ts", full.names = T, recursive = T)


format_future_ts_climate <- function(file) {
  
  opt <- stringr::str_split(file, "_")
  points <- gps
  coordinates(points) <- ~ Longitude + Latitude
  nc <- raster::brick(file)
  df <- raster::extract(nc,points , method="bilinear", df = T) %>% pivot_longer(-ID, names_to = "date", values_to = "value") %>%
    mutate(
      date = case_when(
        grepl("\\.", date) ~ as.Date(gsub("X", "", date), format = "%Y.%m.%d"), 
        TRUE ~ as.Date(paste0("15/01/",regmatches(opt[[1]][8], regexpr('^\\d{4}', opt[[1]][8]))), "%d/%m/%Y") %m+% 
          months(as.integer(gsub("X", "", date)))),
      variable = opt[[1]][5],
      model = opt[[1]][6],
      scenario = opt[[1]][7]
    )
  
  return(df)
}

df <- lapply(files, format_future_ts_climate) %>% bind_rows() %>%
  left_join(., cbind(ID = c(1:8), gps %>% dplyr::select(`Station.ID.Draco.locality.ID`))) %>%
  mutate(ID = `Station.ID.Draco.locality.ID`) %>%
  filter(ID %in% c("Hk", "Cr", "Ks", "Ru")) 

write.csv(df, "/data/lagged/Dracocephalum/CHELSA_future_data.csv", row.names = F)

## On Local machine
gps <- read.csv("data/Draco_locations.csv") %>%
  rename(locality = Station.ID.Draco.locality.ID) %>%
  dplyr::select(locality, Latitude)

df <- read.csv("data/CHELSA_future_data.csv") %>%
  rename(locality = ID) %>%
  left_join(., gps) %>%
  mutate(
    month = lubridate::month(as.Date(date, "%Y-%m-%d")),
    year = lubridate::year(as.Date(date, "%Y-%m-%d")),
    days_per_month = lubridate::days_in_month(as.Date(paste("01", month, year, sep = "/"), "%d/%m/%Y")),
    value = case_when(variable == "pr" ~ value * 86400 * days_per_month,               ## pr unit is kg/m2/s and the pet calculations needs kg/m2/month
                      variable %in% c("tasmax", "tasmin") ~ value - 273.15)) %>%       ## tasmax/min is in K, converted to C for comparison with monthly sequence (historical) 
  dplyr::select(-c(days_per_month, date, contains("Station"))) %>%
  pivot_wider(values_from = value, names_from = variable) %>%
  rowwise() %>%
  mutate(tas = (tasmin + tasmax) / 2,
         pet = as.numeric(SPEI::hargreaves(Tmin = tasmin, Tmax = tasmax, lat = Latitude, Pre = pr, na.rm = T))
         ) %>% 
  pivot_longer(cols = c("pr", "tasmin", "tasmax", "tas", "pet"), names_to = "variable", values_to = "value")

## Get mean and standard deviation from historical data so those can be used to scale future values
hist <- df %>% filter(scenario == "historical") %>%
  group_by(locality, month, variable) %>%
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T))

df1 <- left_join(df, hist) %>%
  mutate(scaled = (value - mean) / sd) %>% ## scale values, using the historical mean and standard deviation
  dplyr::select(-c(mean, sd)) %>%
  pivot_wider(names_from = variable, values_from = c(value, scaled)) %>%
  rename(
    pr = value_pr,
    tas = value_tas,
    tmin = value_tasmin,
    tmax = value_tasmax,
    pr_scaled = scaled_pr,
    tas_scaled = scaled_tas,
    tmin_scaled = scaled_tasmin,
    tmax_scaled = scaled_tasmax,
    pet_scaled = scaled_pet
  ) 

write.csv(df1, "data/CHELSA_future_ts_formatted.csv", row.names = F)
