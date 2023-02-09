#---------------------------------------------------------------
# formatting of demographic data into long format
#---------------------------------------------------------------
source("R/functions_data_formatting.R")
source("R/functions_data_transformation.R")

# Location of raw file
raw_data_file_location <- "data/Dracocephalum_03_21.xlsx"

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
raw_data_long_format <- data %>% 
  mutate(stage_t0 = replace(stage_t0, stage_t0 == "NA", NA)) 

# There are cells with NA in stage_t0 but with recorded flowering or vegetative stems
data$stage_t0[which(is.na(data$stage_t0) & data$n_fl_stems_t0 > 0)] <- "flow"
data$stage_t0[which(is.na(data$stage_t0) & data$n_fl_stems_t0 == 0 & data$n_veg_stems_t0 > 0)] <- "veg"
data$stage_t0[which(is.na(data$stage_t0) & data$n_fl_stems_t0 == 0 & data$n_veg_stems_t0 == 0)] <- "dead"

# Transcript error (herb shading = 13 instead of herb = 1, shrub = 3)
data$herb_shading_t0[which(data$plant_ID == "HAK_105" & data$year_t0 == 2017)] <- 1
data$shrub_shading_t0[which(data$plant_ID == "HAK_105" & data$year_t0 == 2017)] <- 3


write.csv(data, "data/Dracocephalum_long_format.csv", row.names = F)

raw_data_long_format <- data

source("analysis/exploratory_plots.R")

# raw_data_long_format <- read.csv("data/Dracocephalum_long_format.csv")

data <- raw_data_long_format %>%
  calculate_av_values(.) %>%
  add_t1(.) %>%
  add_tm1(.) %>%
  add_tm2(.) %>%
  dplyr::select(plant_ID, population, year_t0, contains("soil"), "rock", "slope", contains("tm2"), contains("tm1"), contains("t0"), contains("t1")) %>%
  rowwise() %>%
  mutate(survival_t1 = ifelse(!is.na(stage_t0) & stage_t1 != "dead", 1, 
                              ifelse(stage_t0 == "dead" & stage_t1 == "dead", NA, 0)),
         flower_p_t0 = ifelse(n_fl_stems_t0 > 0, 1, 0),
         flower_p_t1 = ifelse(n_fl_stems_t1 > 0, 1, 0),
         seed_p_t0 = ifelse(av_seeds_n_t0 == 0 | is.na(av_seeds_n_t0), 0, 1),
         est_seed_n_t0 = av_seeds_n_t0 * n_fl_stems_t0,
         ln_stems_t0 = log(n_fl_stems_t0 + n_veg_stems_t0),
         ln_stems_t1 = log(n_fl_stems_t1 + n_veg_stems_t1)) %>%
  dplyr::select(plant_ID, population, year_t0, soil_depth, rock, slope,
         contains("stage"), survival_t1, contains("stems_t0"), contains("stems_t1"),
         contains("_p_"), contains("seeds"), est_seed_n_t0, contains("shading"), contains("longest_stem"))

write.csv(data, "data/Dracocephalum_with_vital_rates.csv", row.names = F)


rm(list = ls())



