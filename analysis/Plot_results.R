library(tidyverse)
library(glmnet)
library(patchwork)
library(gridExtra)

### Script to plot results
rm(list = ls())
set.seed(2)

source("R/functions_plotting.R")

## --------------------------------------------------------------------
## Vital rates
## --------------------------------------------------------------------

surv_mod <- readRDS("results/rds/seasons_surv.rds")
grow_mod <- readRDS("results/rds/seasons_growth.rds")
flowp_mod <- readRDS("results/rds/seasons_flowp.rds")
seedp_mod <- readRDS("results/rds/seasons_seedp.rds")
seedn_mod <- readRDS("results/rds/seasons_seedn.rds")

surv_plot <- plot_coef_clim_shading(surv_mod, save = T, name = "Survival_clim_shading_coef", plot_title = 'Survival') 
grow_plot <- plot_coef_clim_shading(grow_mod, save = T, name = "Growth_clim_shading_coef", plot_title = 'Growth')
flowp_plot <- plot_coef_clim_shading(flowp_mod, save = T, name = "Flower_P_clim_shading_coef", plot_title = 'Flower probability') 
seedp_plot <- plot_coef_clim_shading(seedp_mod, save = T, name = "Seed_P_clim_shading_coef", plot_title = 'Seed probability') 
seedn_plot <- plot_coef_clim_shading(seedn_mod, save = T, name = "Seed_N_clim_shading_coef", plot_title = 'Number of seeds') 


# Code to copy coefficients as table to excel
# 
# df <- data.frame(rownames(coef(surv_mod)),
#                  matrix(coef(surv_mod)),
#                  matrix(coef(grow_mod)),
#                  matrix(coef(flowp_mod)),
#                  matrix(coef(seedp_mod)),
#                  matrix(coef(seedn_mod))) 
# 
# write.table(df,"clipboard",sep="\t")

## --------------------------------------------------------------------
## Extinciton probability
## --------------------------------------------------------------------
source("R/functions_extinction_plots.R")

df <- lapply(list.files("results/ibm/", full.names = T), readRDS) %>% 
  bind_rows() %>%
  mutate(locality = as.factor(locality),
         yr_of_ext = yr_of_ext + 2020,
         yr_of_ext = replace_na(yr_of_ext, 2101),
         scenario = as.character(scenario))

# df$scenario[which(df$model == "No change")] <- "Historical"
df$scenario <- factor(df$scenario, levels = c("rcp45", "rcp85"))

names(df$pop_size) <- paste(df$locality, df$model, df$scenario, c(1:nrow(df)), sep = "_")


ext_pl <- lapply(df$pop_size, function(x) x$plants) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>% tibble::rowid_to_column() %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_plant") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2020) 

ext_sdl <- lapply(df$pop_size, function(x) x$sdl) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sdl") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2020)

ext_sb1 <- lapply(df$pop_size, function(x) x$sb1) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sb1") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2020)

ext_sb2 <- lapply(df$pop_size, function(x) x$sb2) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sb2") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2020)

ext <- ext_pl %>% cbind(.,
                        pop_size_sdl = ext_sdl$pop_size_sdl, 
                        pop_size_sb1 = ext_sb1$pop_size_sb1, 
                        pop_size_sb2 = ext_sb2$pop_size_sb2) %>%
  mutate(below_ext_size = ifelse(pop_size_plant < 10, T, F)) 

rm(ext_pl, ext_sdl, ext_sb1, ext_sb2)


## Actual plotting
# pop_sizes_herb <- pop_dyn_plot(ext %>% filter(Sshading == 0))
# pop_sizes_shrub <- pop_dyn_plot(ext %>% filter(Hshading == 0))
# ext_years_herb <- ext_yr_plot(df %>% filter(Sshading == 0))
# ext_years_shrub <- ext_yr_plot(df %>% filter(Hshading == 0))
# wrap_plots(pop_sizes_herb, guides = "collect")
# wrap_plots(ext_years, guides = "collect", ncol = 4, byrow = F)

pop_size_45_herb <- ggplot(ext %>%
                             filter(scenario == "rcp45" & Sshading == 4)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = Hshading,
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model
                  colour = Hshading)) +
  coord_cartesian(xlim = c(2021, 2100)) + 
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0))+
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous")
pop_size_45_herb

pop_size_45_shrub <- ggplot(ext %>%
                              filter(scenario == "rcp45" & Hshading == 4)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = Sshading,
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model,
                  colour = Sshading)) +
  coord_cartesian(xlim = c(2021, 2100)) + 
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0))+
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub")
pop_size_45_shrub

pop_size_85_herb <- ggplot(ext %>%
                             filter(scenario == "rcp85" & Sshading == 4)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = Hshading,
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model,
                  colour = Hshading)) +
  coord_cartesian(xlim = c(2021, 2100)) + 
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0))+
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous")
pop_size_85_herb

pop_size_85_shrub <- ggplot(ext %>%
                              filter(scenario == "rcp85" & Hshading == 4)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = Sshading,
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model,
                  colour = Sshading)) +
  coord_cartesian(xlim = c(2021, 2100)) + 
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0))+
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub")
pop_size_85_shrub

pop_sizes_45_supp <- pop_size_45_herb + pop_size_45_shrub + plot_layout(guides = "collect") & theme(text = element_text(size = 14))
pop_sizes_45 <- pop_size_45_herb + pop_size_45_shrub + plot_layout(guides = "collect") & coord_cartesian(xlim = c(2021,2050)) & theme(text = element_text(size = 14))

pop_sizes_85_supp <- pop_size_85_herb + pop_size_85_shrub + plot_layout(guides = "collect") & theme(text = element_text(size = 14))
pop_sizes_85 <- pop_size_85_herb + pop_size_85_shrub + plot_layout(guides = "collect") & coord_cartesian(xlim = c(2021,2050)) & theme(text = element_text(size = 14))
pop_sizes_85

ggsave(pop_sizes_45, filename = file.path("results", "IBM_population_sizes_rcp4.5.png"), width = 8, height = 8)
ggsave(pop_sizes_85, filename = file.path("results", "IBM_population_sizes_rcp8.5.png"), width = 8, height = 8)

ggsave(pop_sizes_45_supp, filename = file.path("results", "IBM_population_sizes_rcp4.5_supplement.png"), width = 10, height = 8)
ggsave(pop_sizes_85_supp, filename = file.path("results", "IBM_population_sizes_rcp8.5_supplement.png"), width = 10, height = 8)

## Extinction year distribution --------------------------------------------------------------------

ext_45_herb <- ggplot(df %>%
                        filter(scenario == "rcp45" & Sshading == 4)) +
  geom_boxplot(aes(y = yr_of_ext,
                   # x = model,
                   fill = as.factor(Hshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction") + ggtitle("Herb shading") +
  coord_flip(ylim = c(2022,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  ) 
ext_45_herb

ext_45_shrub <- ggplot(df %>%
                         filter(scenario == "rcp45" & Hshading == 4)) +
  geom_boxplot(aes(y = yr_of_ext,
                   # x = model,
                   fill = as.factor(Sshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Shrub shading") +
  coord_flip(ylim = c(2022,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  ) 
ext_45_shrub

ext_85_herb <- ggplot(df %>%
                        filter(scenario == "rcp85" & Sshading == 4)) +
  geom_boxplot(aes(y = yr_of_ext,
                   # x = model,
                   fill = as.factor(Hshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Herb shading") +
  coord_flip(ylim = c(2022,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  ) 
ext_85_shrub <- ggplot(df %>%
                         filter(scenario == "rcp85" & Hshading == 4)) +
  geom_boxplot(aes(y = yr_of_ext,
                   # x = model,
                   fill = as.factor(Sshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Shrub shading") +
  coord_flip(ylim = c(2022,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  )  

ext_45 <- ext_45_herb / ext_45_shrub + plot_layout(guides = "collect")
ext_85 <- ext_85_herb / ext_85_shrub + plot_layout(guides = "collect")

ggsave(ext_45, filename = file.path("results", "IBM_extinction_year_rcp4.5.png"), width = 4, height = 4)
ggsave(ext_85, filename = file.path("results", "IBM_extinction_year_rcp8.5.png"), width = 4, height = 4)


## --------------------------------------------------------------------
## Visualise shading level distributions
## --------------------------------------------------------------------

df <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
  filter(year_t0 < 2020) %>%
  rowwise() %>%
  mutate(shading_t0 = ifelse(any(!is.na(c(herb_shading_t0, shrub_shading_t0))), 
                             sum(herb_shading_t0, shrub_shading_t0, na.rm = T),
                             NA))%>% 
  group_by(population) %>%
  mutate(mean = mean(shading_t0, na.rm = T))

shading_df <- data.frame(lambda = rep(c(0,2,4,6), each = 1000)) %>%
  rowwise() %>%
  mutate(shading = rpois(1, lambda)) 

shading_dist <- ggplot(df, aes(x = shading_t0)) + 
  geom_histogram(aes(y = ..density..), binwidth = 1) + 
  geom_vline(aes(xintercept = mean, group = population)) +
  facet_grid(rows = vars(population)) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Locality", breaks = NULL, labels = NULL)) + 
  xlab("Shading level") + ylab("") + 
  coord_cartesian(xlim = c(0,17))

shading_mod <- ggplot(shading_df) +
  geom_density(aes(x = shading, y = after_stat(scaled), 
                   colour = as.factor(lambda), 
                   fill = as.factor(lambda)), 
               alpha = 0.5,
               adjust = 2) +
  xlab("Modelled shading level") + ylab("") +
  guides(fill = guide_legend(title = "Shading level"),
         colour = guide_legend(title = "Shading level")) + 
  coord_cartesian(xlim = c(0,17))

plot_shading <- shading_dist + shading_mod + 
  plot_layout(ncol = 1, 
              heights = c(3,1)) +
  plot_annotation(tag_levels = "A")

ggsave("results/shading_distributions.png", plot_shading,
       width = 4, height = 5)


## --------------------------------------------------------------------
## Visualise climate sequences
## --------------------------------------------------------------------
fut_clim <- read.csv("data/CHELSA_future_ts_formatted.csv") %>% 
  filter(complete.cases(.))

fut_spring <- fut_clim %>% 
  filter(month > 2 & month < 6) %>%
  mutate(year_t0 = year - 1, .keep = "unused") %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_spring = mean(pet_scaled, na.rm = T),
            pr_spring = mean(pr_scaled, na.rm = T),
            tas_spring = mean(tas_scaled, na.rm = T)) %>%
  ungroup() 
fut_summer <- fut_clim %>% 
  filter(month > 2 & month < 6) %>%
  rename(year_t0 = year) %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_summer = mean(pet_scaled, na.rm = T),
            pr_summer = mean(pr_scaled, na.rm = T),
            tas_summer = mean(tas_scaled, na.rm = T)) %>%
  ungroup() 
fut_dormant <- fut_clim %>% 
  filter(month < 3 | month > 8) %>%
  mutate(year_t0 = ifelse(month < 3, year - 1, year)) %>%
  group_by(locality, year_t0, model, scenario) %>%
  summarise(pet_dormant = mean(pet_scaled, na.rm = T),
            pr_dormant = mean(pr_scaled, na.rm = T),
            tas_dormant = mean(tas_scaled, na.rm = T)) %>%
  ungroup() 

fut_clim <- left_join(fut_spring, fut_summer)
fut_clim <- left_join(fut_clim, fut_dormant)
rm(fut_summer, fut_spring, fut_dormant)


hist_data <- read.csv("data/CHELSA_data.csv")

hist_spring <- hist_data %>% 
  filter(month > 2 & month < 6) %>%
  mutate(year_t0 = year - 1, .keep = "unused") %>%
  group_by(locality, year_t0) %>%
  summarise(pet_spring = mean(pet_scaled, na.rm = T),
            pr_spring = mean(pr_scaled, na.rm = T),
            tas_spring = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")
hist_summer <- hist_data %>% 
  filter(month > 2 & month < 6) %>%
  rename(year_t0 = year) %>%
  group_by(locality, year_t0) %>%
  summarise(pet_summer = mean(pet_scaled, na.rm = T),
            pr_summer = mean(pr_scaled, na.rm = T),
            tas_summer = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")
hist_dormant <- hist_data %>% 
  filter(month < 3 | month > 8) %>%
  mutate(year_t0 = ifelse(month < 3, year - 1, year)) %>%
  group_by(locality, year_t0) %>%
  summarise(pet_dormant = mean(pet_scaled, na.rm = T),
            pr_dormant = mean(pr_scaled, na.rm = T),
            tas_dormant = mean(tas_scaled, na.rm = T)) %>%
  ungroup() %>%
  mutate(population = toupper(locality), .keep = "unused")

hist_clim <- left_join(hist_spring, hist_summer)
hist_clim <- left_join(hist_clim, hist_dormant)
rm(hist_spring, hist_summer, hist_dormant, hist_data)


a <- left_join(fut_clim, hist_clim) 

