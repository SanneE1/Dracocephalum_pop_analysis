## predict future climate using ARIMA. For sensible forecasting I think I need to use the absolute values, rather than the anomalies
## However, these are then easily converted to anomalies using the mean and sd of the OBSERVED data (i.e. chelsa data)


clim_l <- read.csv('data/CHELSA_data.csv') %>%
  filter(locality %in% c("Hk", "Cr", "Ks", "Ru")) %>%
  dplyr::select(locality, month, year, tas_scaled, pr_scaled, pet_scaled) %>%
  arrange(year, month) %>% pivot_longer(c("tas_scaled", "pr_scaled", "pet_scaled"), names_to = "variable", values_to = "value") %>%
  split(., list(.$locality, .$variable) ) %>% 
  lapply(., function(x) 
    ts(x %>% filter(complete.cases(x)) %>% dplyr::select(value), 
       frequency = 12, start = c(x$year[1],x$month[1]))) %>%
  lapply(., auto.arima)


fut_l <- read.csv("data/CHELSA_future_ts_formatted.csv") %>%
  filter(scenario != "historical") %>%
  dplyr::select(locality, model, scenario, month, year, tas_scaled, pr_scaled, pet_scaled) %>%
  arrange(year, month) %>% pivot_longer(c("tas_scaled", "pr_scaled", "pet_scaled"), names_to = "variable", values_to = "value") %>%
  split(., list(.$locality, .$model, .$scenario, .$variable) ) %>% 
  lapply(., function(x) 
    ts(x %>% filter(complete.cases(x)) %>% dplyr::select(value), 
       frequency = 12, start = c(x$year[1],x$month[1]))) %>%
  lapply(., auto.arima)


arima_l <- list(clim_hist_model = clim_l,
                clim_future_model = fut_l)


saveRDS(arima_l,
        file = "results/ARIMA_clim_mods.rds")

rm(list = ls())
