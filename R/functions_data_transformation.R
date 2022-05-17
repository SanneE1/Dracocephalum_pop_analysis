# Functions to transform the data


transform_raw_data_and_save <- function(raw_data_long_format) {
  # Calculate the average value for the variables with two replicates
  calculate_av_values <- function(data) {
    data %>% 
      mutate(fruit_stem_length_t0 = (fruit_stem_length_1_t0 + fruit_stem_length_2_t0)/2,
             inflor_length_t0 = (inflor_length_1_t0 + inflor_length_2_t0)/2,
             calices_n_t0 = (calices_n_1_t0 + calices_n_2_t0)/2,
             av_seeds_n_t0 = (seeds_n_1_t0 + seeds_n_2_t0)/2) %>%
      select(-c(fruit_stem_length_1_t0, fruit_stem_length_2_t0, 
                inflor_length_1_t0, inflor_length_2_t0,
                calices_n_1_t0, calices_n_2_t0)
             )
  }
  
  ## function to add size t1 to rows
  add_t1 <- function(data) {
    d1 <- data %>% select(-contains("shading"), -contains("Draco"), -"cuscuta_t0") %>%
      mutate(year_t0 = year_t0 - 1)
    names(d1)[3:11] <- gsub(x = names(d1)[3:11], pattern = "_t0", replacement = "_t1")  
    
    left_join(data, d1)
  }
  
  data <- raw_data_long_format %>%
    calculate_av_values(.) %>%
    add_t1(.) %>%
    select(plant_ID, population, year_t0, contains("soil"), "rock", "slope", contains("t0"), contains("t1")) %>%
    mutate(survival_t1 = ifelse(!is.na(stage_t0) & stage_t1 != "dead", 1, 
                                ifelse(stage_t0 == "dead" & stage_t1 == "dead", NA, 0)),
           flower_p_t0 = ifelse(n_fl_stems_t0 > 0, 1, 0),
           flower_p_t1 = ifelse(n_fl_stems_t1 > 0, 1, 0),
           seed_p_t0 = ifelse(av_seeds_n_t0 > 0, 1, 0),
           est_seed_n_t0 = av_seeds_n_t0 * n_fl_stems_t0,
           ln_stems_t0 = log(n_fl_stems_t0 + n_veg_stems_t0),
           ln_stems_t1 = log(n_fl_stems_t1 + n_veg_stems_t1)) %>%
    # select(-c(cuscuta_t0, contains("longest_stem"), contains("Draco"), 
    #           fruit_stem_length_t0, inflor_length_t0, calices_n_t0)) %>%
    select(plant_ID, population, year_t0, contains("soil"), rock, slope,
           contains("stage"), survival_t1, contains("stems_t0"), contains("stems_t1"),
           contains("_p_"), contains("seeds"), est_seed_n_t0, contains("shading"), contains("longest_stem"))
  
  write.csv(data, "data/Dracocephalum_with_vital_rates.csv", row.names = F)
  
  return("data/Dracocephalum_with_vital_rates.csv")
  
}

