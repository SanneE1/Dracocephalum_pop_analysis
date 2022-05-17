#---------------------------------------------------------------
# formatting of demographic data into long format
#---------------------------------------------------------------
format_raw_excel_data_to_long_format <- function(raw_data_file_location) {
  ## Read in data
  data_sheets <- as.list(c("dataCR", "dataHK", "dataKS", "data RU"))
  
  demo_data_list <- lapply(data_sheets, function(x) 
    readxl::read_excel(raw_data_file_location, sheet = x)
  )
  
  # These cells have * as a value. Should be NA
  demo_data_list[[1]]$S12[which(demo_data_list[[1]]$plant_ID == "CR_765")] <- NA
  demo_data_list[[1]]$S12 <- as.integer(demo_data_list[[1]]$S12)
  demo_data_list[[4]]$l.stem11[which(demo_data_list[[4]]$plant_ID == "RU_90")] <- NA
  demo_data_list[[4]]$l.stem11 <- as.numeric(demo_data_list[[4]]$l.stem11)
  
  # transform excel data into longer format and adding the t1 values
  ### (after checking for unexpected character class in specified columns)
  CR <- excel_to_longer_format(demo_data_list[[1]])
  HK <- excel_to_longer_format(demo_data_list[[2]])
  KS <- excel_to_longer_format(demo_data_list[[3]])
  RU <- excel_to_longer_format(demo_data_list[[4]])
  
  # Merging long format data from all populations
  data <- rbind(
    CR %>% mutate(population = "CR"),
    HK %>% mutate(population = "HK"),
    KS %>% mutate(population = "KS"),
    RU %>% mutate(population = "RU")
  )
  
  rm(list = c("CR", "HK", "KS", "RU", "data_sheets", "demo_data_list"))
  
  # There are cells with "NA" as character, rather than system's NA
  data <- data %>% 
    mutate(stage_t0 = replace(stage_t0, stage_t0 == "NA", NA)) 

  # There are cells with NA in stage_t0 but with recorded flowering or vegetative stems
  data$stage_t0[which(is.na(data$stage_t0) & data$n_fl_stems_t0 > 0)] <- "flow"
  data$stage_t0[which(is.na(data$stage_t0) & data$n_fl_stems_t0 == 0 & data$n_veg_stems_t0 > 0)] <- "veg"
  data$stage_t0[which(is.na(data$stage_t0) & data$n_fl_stems_t0 == 0 & data$n_veg_stems_t0 == 0)] <- "dead"
  
  write.csv(data, "data/Dracocephalum_long_format.csv")
  
  data
  
}






