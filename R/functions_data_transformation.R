# Functions to transform the data

# Calculate the average value for the variables with two replicates
calculate_av_values <- function(data) {
  data %>% rowwise() %>% 
    mutate(fruit_stem_length_t0 = sum(fruit_stem_length_1_t0, fruit_stem_length_2_t0, na.rm = T)/2,
           inflor_length_t0 = sum(inflor_length_1_t0, inflor_length_2_t0, na.rm = T)/2,
           calices_n_t0 = sum(calices_n_1_t0, calices_n_2_t0, na.rm = T)/2,
           av_seeds_n_t0 = sum(seeds_n_1_t0, seeds_n_2_t0, na.rm = T)/2,
           soil_depth = sum(soil_d1, soil_d2, soil_d3, na.rm = T)/3) %>%
    dplyr::select(-c(fruit_stem_length_1_t0, fruit_stem_length_2_t0, 
              inflor_length_1_t0, inflor_length_2_t0,
              calices_n_1_t0, calices_n_2_t0, 
              soil_d1, soil_d2, soil_d3)
    )
}

## function to add size t1 to rows
add_t1 <- function(data) {
  d1 <- data %>% dplyr::select(-contains("Draco"), -"cuscuta_t0") %>%
    mutate(year_t0 = year_t0 - 1)
  names(d1)[3:18] <- gsub(x = names(d1)[3:18], pattern = "_t0", replacement = "_t1")  
  
  left_join(data, d1)
}

## function to add shading in t-1
add_tm12 <- function(data) {
  d1 <- data %>% dplyr::select(plant_ID, population, year_t0, contains("shading")) %>%
    mutate(year_t0 = year_t0 + 2)
  names(d1)[4:5] <- gsub(x = names(d1)[4:5], pattern = "_t0", replacement = "_tm2")
  names(d1)[6:7] <- gsub(x = names(d1)[6:7], pattern = "_t1", replacement = "_tm1")  
  
  
  left_join(data, d1)

  }



