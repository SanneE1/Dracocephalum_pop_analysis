## Formatting functions


## Function to left join all columns created by longer functions after checking for unexpected character classes
excel_to_longer_format <- function(excel_data) {
  
  # load some internal functions
  ## check for character class in columns that should be numeric/integers
  
  check_unexpected_character_class <- function(data) {
    as.data.frame(lapply(data, class)) %>%
      pivot_longer(everything(), names_to = "column", values_to = "class_type") %>%
      filter(class_type == "character") %>%
      filter(!grepl("plant_ID|plant|mark|old_no|stg|st\\d{2}|note|poz|hl|ska|skl", column))
  }
  
  
  ## Pivot longer functions
  ## functions to transform excell data to longer format
  
  n_fl_stems_longer <- function(excel_data) {
    df <- excel_data %>%
      dplyr::select(plant_ID, contains("fl"), -contains("infl"), -contains("note")) %>%
      pivot_longer(!plant_ID, names_to = "year_t0", values_to = "n_fl_stems_t0") %>%
      mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))
    message( paste("column class:"), typeof(df$n_fl_stems))
    return(df)
  }
  
  n_veg_stems_longer <- function(excel_data) {
    df <- excel_data %>% 
      dplyr::select(plant_ID, contains("veg"), -contains("note")) %>%
      pivot_longer(!plant_ID, names_to = "year_t0", values_to = "n_veg_stems_t0") %>%
      mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))
    message( paste("column class:"), typeof(df$n_veg_stems))
    return(df)
  }
  
  stage_longer <- function(excel_data) {
    df <- excel_data %>%
      dplyr::select(starts_with("st"), -contains("_"), plant_ID, -contains("note")) %>% # As I dedplyr::select everything with _ plant_ID column as to be dplyr::selected after. Change the order, and it won't work!
      pivot_longer(!plant_ID, names_to = "year_t0", values_to = "stage_t0") %>%
      mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))  
    message( paste("column class:"), typeof(df$stage))
    return(df)
  }
  
  herb_shading_longer <- function(excel_data) {
    df <- excel_data %>% 
      dplyr::select(plant_ID, contains("H", ignore.case = F), -contains("note")) %>%
      pivot_longer(!plant_ID, names_to = "year_t0", values_to = "herb_shading_t0") %>%
      mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))
    message( paste("column class:"), typeof(df$herb_shading))
    return(df)
  }
  
  shrub_shading_longer <- function(excel_data) {
    df <- excel_data %>% 
      dplyr::select(plant_ID, contains("S", ignore.case = F), -contains("note")) %>%
      pivot_longer(!plant_ID, names_to = "year_t0", values_to = "shrub_shading_t0") %>%
      mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))
    message( paste("column class:"), typeof(df$shrub_shading))
    return(df)
  }
  
  cuscuta_longer <- function(excel_data) {
    df <- excel_data %>% 
      dplyr::select(plant_ID, contains("cus"), -contains("note")) %>%
      pivot_longer(!plant_ID, names_to = "year_t0", values_to = "cuscuta_t0") %>%
      mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))
    message( paste("column class:"), typeof(df$cuscuta))
    return(df)
  }
  
  stem_length_longer <- function(excel_data) {  ## contains an as.integer statement as there's a * somewhere in RU
    df <- excel_data %>%
      dplyr::select(plant_ID, contains("l.stem"), -contains("note")) %>%
      mutate(across(contains("l.stem"), as.integer)) %>%
      pivot_longer(!plant_ID, names_to = "year_t0", values_to = "longest_stem_t0") %>%
      mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))
    message( paste("column class:"), typeof(df$longest_stem))
    return(df)
  }
  
  draco_longer <- function(excel_data) {
    df <- left_join(
      excel_data %>% dplyr::select(plant_ID, contains("Draco50"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "Draco50_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$")))),
      excel_data %>% dplyr::select(plant_ID, contains("Draco100"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "Draco100_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}$"))))
    )
    message( paste("column class Draco50:"), typeof(df$Draco50))
    message( paste("column class Draco50:"), typeof(df$Draco100))
    return(df)
    
  }
  
  fruit_stem_length_longer <- function(excel_data) {
    df <- left_join(
      excel_data %>% dplyr::select(plant_ID, matches("stem\\d{2}_1"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "fruit_stem_length_1_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}")))),
      excel_data %>% dplyr::select(plant_ID, matches("stem\\d{2}_2"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "fruit_stem_length_2_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}"))))
    )
    message( paste("column class fruit_stem_1:"), typeof(df$fruit_stem_length_1))
    message( paste("column class fruit_stem_2:"), typeof(df$fruit_stem_length_2))
    return(df)
  }
  
  infl_length_longer <- function(excel_data) {
    df <- left_join(
      excel_data %>% dplyr::select(plant_ID, matches("infl\\d{2}_1"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "inflor_length_1_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}")))),
      excel_data %>% dplyr::select(plant_ID, matches("infl\\d{2}_2"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "inflor_length_2_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}"))))
    )
    message( paste("column class inflor_length_1:"), typeof(df$inflor_length_1))
    message( paste("column class inflor_length_2:"), typeof(df$inflor_length_2))
    return(df)
  }
  
  calices_n_longer <- function(excel_data) {
    df <- left_join(
      excel_data %>% dplyr::select(plant_ID, matches("cal\\d{2}_1"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "calices_n_1_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}")))),
      excel_data %>% dplyr::select(plant_ID, matches("cal\\d{2}_2"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "calices_n_2_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}"))))
    )
    message( paste("column class calices_n_1:"), typeof(df$calices_n_1))
    message( paste("column class calices_n_2:"), typeof(df$calices_n_2))
    return(df)
  }
  
  seeds_n_longer <- function(excel_data) {
    df <- left_join(
      excel_data %>% dplyr::select(plant_ID, matches("seeds\\d{2}_1"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "seeds_n_1_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}")))),
      excel_data %>% dplyr::select(plant_ID, matches("seeds\\d{2}_2"), -contains("note")) %>%
        pivot_longer(!plant_ID, names_to = "year_t0", values_to = "seeds_n_2_t0") %>%
        mutate(year_t0 = as.integer(paste0("20", stringr::str_extract(year_t0, "\\d{2}"))))
    )
    message( paste("column class seeds_n_1:"), typeof(df$seeds_n_1))
    message( paste("column class seeds_n_2:"), typeof(df$seeds_n_2))
    return(df)
  }
  
  # check for unexpected character columns in data
  n <- check_unexpected_character_class(excel_data)
  if(nrow(n) > 0) stop(paste("unexpected character columns for", apply( n[,"column"] , 2, paste , collapse = ", ") ))
  
  list_func <- list(n_fl_stems_longer,
                    n_veg_stems_longer,
                    stage_longer,
                    herb_shading_longer,
                    shrub_shading_longer,
                    cuscuta_longer,
                    stem_length_longer,
                    draco_longer,
                    fruit_stem_length_longer,
                    infl_length_longer,
                    calices_n_longer,
                    seeds_n_longer)
  list <- lapply(list_func, function(f) f(excel_data))
  df <- plyr::join_all(list, by = c("plant_ID", "year_t0"), type = "full")
  
  df <- left_join(df,
                  excel_data %>%
                    dplyr::select(plant_ID, soil_d1, soil_d2, soil_d3, rock, slope))
  
  df <- filter(df, rowSums(is.na(df)) != ncol(df)-2)
  df <- filter(df, !is.na(plant_ID) & !is.na(year_t0))
}



