## Exploratory plots of data

raw_data_long_format <- read.csv("data/Dracocephalum_long_format.csv")

## Check for outliers using histograms
n_fl_stems_t0 <- ggplot(raw_data_long_format, aes(x=n_fl_stems_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 10000)) + xlim(c(0,45)) +
  facet_wrap(~population) + ggtitle("distribution of observed number of flowering stems")

n_veg_stems_t0 <- ggplot(raw_data_long_format, aes(x=n_veg_stems_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 10000)) + xlim(c(0,45)) +
  facet_wrap(~population)  + ggtitle("distribution of observed number of vegetative stems")

stage_t0 <- ggplot(raw_data_long_format, aes(x=stage_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 10000)) + 
  facet_wrap(~population) + ggtitle("distribution of observed plant stages (with NA = individual not found")

longest_stem_t0 <- ggplot(raw_data_long_format, aes(x=longest_stem_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 1000)) +  
  facet_wrap(~population) + ggtitle("distribution of observed length of the overall longest stem (cm)")

fruit_stem_length <- ggplot(raw_data_long_format %>% 
                              pivot_longer(contains("fruit_stem_length"), names_to = "replicate", values_to = "length"), 
                            aes(x=length, group = replicate, fill = replicate)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 1000)) +  
  facet_wrap(~population)  + ggtitle("distribution of observed Dracocephalum density in a 100cm radius around individual")

inflor_length <- ggplot(raw_data_long_format %>% 
                          pivot_longer(contains("inflor_length"), names_to = "replicate", values_to = "length"), 
                        aes(x=length, group = replicate, fill = replicate)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 500)) + xlim(c(0,125)) + 
  facet_wrap(~population) + ggtitle("distribution of observed length in mm of the inflorescence \n(from the lowest to the highest fruit) at two randomly chosen fruiting stem")

calices_n <- ggplot(raw_data_long_format %>% 
                      pivot_longer(contains("calices_n_"), names_to = "replicate", values_to = "number"), 
                    aes(x=number, group = replicate, fill = replicate)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 100)) + xlim(c(0,30)) +
  facet_wrap(~population) + ggtitle("distribution of observed number of calices at two \n randomly chosen fruiting stems")

seeds_n <- ggplot(raw_data_long_format %>% 
                    pivot_longer(contains("seeds_n_"), names_to = "replicate", values_to = "number"), 
                  aes(x=number, group = replicate, fill = replicate)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 1000)) + xlim(c(0,30)) +  
  facet_wrap(~population)  + ggtitle("distribution of observed number of good seeds (black and hard) \nat two randomly chosen fruiting stems")

herb_shading_t0 <- ggplot(raw_data_long_format, aes(x=herb_shading_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 10000)) + xlim(c(0,10)) +
  facet_wrap(~population)  + ggtitle("distribution of observed herb shading intensity intensity of herb shading \n(0 - no other plants within 20 cm, 10 - very dense herb vegetation - higher than Dracocephalum - completely overgrowing)")

shrub_shading_t0 <- ggplot(raw_data_long_format, aes(x=shrub_shading_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 1000)) + xlim(c(0,10)) + 
  facet_wrap(~population)  + ggtitle("distribution of observed shrub and tree shading intensity \n(0 - no shrubs and trees within 20 cm, 10 - very dense shrub vegetation - higher than Dracocephalum - completely overgrowing)") 

cuscuta_t0 <- ggplot(raw_data_long_format, aes(x=cuscuta_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 100)) + xlim(c(0,10)) + 
  facet_wrap(~population) + ggtitle("distribution of observed intensity of Cuscuta epithymum infestation \n(0 - no Cuscuta, 10 - plant completely infested with Cuscuta)")

Draco50_t0 <- ggplot(raw_data_long_format, aes(x=Draco50_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 1000)) + xlim(c(0,10)) + 
  facet_wrap(~population) + ggtitle("distribution of observed Dracocephalum density in a 50cm radius around individual")

Draco100_t0 <- ggplot(raw_data_long_format, aes(x=Draco100_t0)) +
  geom_bar(position = "dodge", width = 1) + scale_y_log10(limits=c(1, 1000)) + xlim(c(0,10)) + 
  facet_wrap(~population) + ggtitle("distribution of observed Dracocephalum density in a 50cm radius around individual")




## Save as pdf file
pdf(file = "results/exploratory plots/raw_data_observed_distributions.pdf")

print(n_fl_stems_t0)
print(n_veg_stems_t0)
print(stage_t0)
print(longest_stem_t0)
print(fruit_stem_length)
print(inflor_length)
print(calices_n)
print(seeds_n)
print(herb_shading_t0)
print(shrub_shading_t0)
print(cuscuta_t0)
print(Draco50_t0)
print(Draco100_t0)

dev.off()

rm(list = ls())

