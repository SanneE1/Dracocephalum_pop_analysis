# Adapted from a script written by Aldo Compagnoni

## Run this part in R

library(tidyverse)

##------------------------------------------------------
## Current Climate
##------------------------------------------------------

# Root directory
root_dir    <- 'C:/Users/se44heqo/Documents/CHELSA/'

# file with coordinates. Text because I don't know how to read CSV in Anaconda powershell prompt
## file with two columns ("Longitude", "Latitude")
coord_file  <- 'coordinates.txt'
# path to file
read_cmd    <- paste0('cat ', root_dir, coord_file)

# URLS: I stored the WGET paths of chelsa files I downloaded from "https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F"
# Visit that URL, select a folder, then click the icon "Download only the selected files via Wget" on the right!
url_roots   <- read.table( paste0( root_dir, 'CHELSA_monthly_paths.txt') )[,1]
# Add vsicurl - enables to read files online
chelsa_cmd  <- paste0('/vsicurl/', url_roots) 

# Pipe, path to gdallocationinfo.exe, and arguments (look this up online)
gdal_cmd    <- ' | C:/Users/se44heqo/Anaconda3/Library/bin/gdallocationinfo.exe -wgs84 -valonly '


# Produce the names of output files starting from url_roots
file_names  <- gsub( 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/',
                     '',url_roots ) %>% 
  gsub( '.tif', '_coord.txt', .)

# store values: one file at a time
output_cmd  <- paste0(' > ', root_dir, 'results/', file_names )
# these are the lines of code
store       <- paste0(read_cmd, gdal_cmd, chelsa_cmd, output_cmd )
store <- gsub('C:/', '"C:/', store) %>%
  gsub('.txt', '.txt"', .) %>%
  gsub('.exe', '.exe"', .)

# write a file
write.table( store,
             paste0(root_dir,'CHELSA_get.txt'),
             col.names = F, row.names = F, quote = F )


## Copy past the resulting file into the Anaconda Prompt. Once the downloading is complete continue on with the script
## I know this really isn't good practice and not really reproducible, but I haven't figured out how to do this from here
download_location = 'C:/Users/se44heqo/Documents/CHELSA/results/'

coord <- read.table(here::here(root_dir, coord_file), header = F, sep = " ",
                    col.names = c("Longitude", "Latitude"))

files <- list.files(download_location, full.names = T, recursive = T)

df <- lapply(as.list(files),
             function(x) cbind(coord,
                               month = stringr::str_extract(x, '_\\d{2}_') %>% gsub("_", "",.),
                               year = stringr::str_extract(x, '_\\d{4}_') %>% gsub("_", "",.),
                               variable = stringr::str_extract(x, '_[[:lower:]]{2,3}_|_[[:lower:]]{6}_') %>% gsub("_", "",.),
                               value = read.table(x))
) %>% bind_rows() %>% 
  pivot_wider(names_from = variable, values_from = V1) %>%
  mutate(month = as.integer(month),
         year = as.integer(year),
         tas = (tas/10) - 273.15,        #CHELSA's unit = K/10, here converted to C
         tasmin = (tasmin/10) - 273.15,
         tasmax = (tasmax/10) - 273.15,
         pr = pr/100
) %>%
  left_join(., 
            read.csv("data/Draco_locations.csv") %>% rename(locality = `Station.ID.Draco.locality.ID`) ) %>%
  filter(locality %in% c("Hk", "Cr", "Ks", "Ru")) %>%
  dplyr::relocate(any_of(c("Name", "locality", "Longitude", "Latitude", "Altitude"))) %>%
  rowwise() %>%
  mutate(pet_c = as.numeric(SPEI::hargreaves(Tmin = tasmin, Tmax = tasmax, lat = Latitude, Pre = pr, na.rm = T))) %>%
  ## Include filter to the 30-year period
  group_by(locality, month) %>%
  mutate(
    cmi_scaled = scale(cmi),
    pet_scaled = scale(pet_c),
    pr_scaled = scale(pr),
    tas_scaled = scale(tas),
    vpd_scaled = scale(vpd)
  ) %>% ungroup() 


print("Correlation between CHELSA's pet and calculated pet")
cor.test(df$pet, df$pet_c, na.rm = T)


write.csv(df, "data/CHELSA_data.csv", row.names = F)

rm(list = ls())






