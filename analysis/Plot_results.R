library(tidyverse)
library(glmnet)
library(patchwork)
library(gridExtra)
library(ggh4x)

### Script to plot results
rm(list = ls())
set.seed(2)

# Years of future simulations
start_yr = 2018
end_yr = 2100

source("R/functions_plotting.R")

## --------------------------------------------------------------------
## Vital rates
## --------------------------------------------------------------------

# Code to copy coefficients as table to excel
# 
# surv_mod <- readRDS("results/rds/seasons_surv.rds")
# grow_mod <- readRDS("results/rds/seasons_growth.rds")
# flowp_mod <- readRDS("results/rds/seasons_flowp.rds")
# seedp_mod <- readRDS("results/rds/seasons_seedp.rds")
# seedn_mod <- readRDS("results/rds/seasons_seedn.rds")
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
## Data prep
## --------------------------------------------------------------------
source("R/functions_extinction_plots.R")

df <- lapply(list.files("results/ibm/", full.names = T), readRDS) %>% 
  bind_rows() %>%
  mutate(locality = as.factor(locality),
         yr_of_ext = yr_of_ext + (start_yr - 1),
         yr_of_ext = replace_na(yr_of_ext, end_yr+10),
         scenario = as.character(scenario),
         below_ext_size_2050 = ifelse(yr_of_ext < 2052, T, F))

# df$scenario[which(df$model == "No change")] <- "Historical"
df$scenario <- factor(df$scenario, levels = c("rcp45", "rcp85"))

names(df$pop_size) <- paste(df$locality, df$model, df$scenario, c(1:nrow(df)), sep = "_")


ext_pl <- lapply(df$pop_size, function(x) x$plants) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>% tibble::rowid_to_column() %>%
  pivot_longer(num_range(prefix = "", range = c(1:(end_yr - (start_yr - 1)))), names_to = "time_step", values_to = "pop_size_plant") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + (start_yr - 1)) 

ext_sdl <- lapply(df$pop_size, function(x) x$sdl) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:(end_yr - (start_yr - 1)))), names_to = "time_step", values_to = "pop_size_sdl") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + (start_yr - 1))

ext_sb1 <- lapply(df$pop_size, function(x) x$sb1) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:(end_yr - (start_yr - 1)))), names_to = "time_step", values_to = "pop_size_sb1") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + (start_yr - 1))

ext_sb2 <- lapply(df$pop_size, function(x) x$sb2) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:(end_yr - (start_yr - 1)))), names_to = "time_step", values_to = "pop_size_sb2") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + (start_yr - 1))

ext <- ext_pl %>% cbind(.,
                        pop_size_sdl = ext_sdl$pop_size_sdl, 
                        pop_size_sb1 = ext_sb1$pop_size_sb1, 
                        pop_size_sb2 = ext_sb2$pop_size_sb2) %>%
  mutate(below_ext_size = ifelse(pop_size_plant < 10, T, F)) 

rm(ext_pl, ext_sdl, ext_sb1, ext_sb2)


ext_summary <- ext %>%
  # filter(scenario == "rcp45" & Sshading == 2) %>%
  group_by(Sshading, Hshading, model, scenario, locality, year) %>%
  summarise(
    mean_pop = mean(pop_size_plant, na.rm = TRUE),
    q1 = quantile(pop_size_plant, 0.1, na.rm = TRUE),
    q3 = quantile(pop_size_plant, 0.9, na.rm = TRUE),
    n = sum(!is.na(pop_size_plant)),
    se = sd(pop_size_plant, na.rm = TRUE) / sqrt(n),
    ci_lower = mean_pop - qt(0.975, df = n - 1) * se,
    ci_upper = mean_pop + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )



## --------------------------------------------------------------------
## Population dynamics
## --------------------------------------------------------------------

pop_size_45_herb <- ggplot(ext_summary %>% filter(Sshading == 2 &
                                                    scenario == "rcp45"),
                           aes(x = year, y = mean_pop,
                               colour = Hshading, fill = Hshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  # coord_cartesian(xlim = c(start_yr, 2050)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), cols = vars(model), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



pop_size_45_shrub <- ggplot(ext_summary %>% filter(Hshading == 2 &
                                                     scenario == "rcp45"),
                            aes(x = year, y = mean_pop,
                                colour = Sshading, fill = Sshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  # coord_cartesian(xlim = c(start_yr, 2050)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), cols = vars(model), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



pop_size_85_herb <- ggplot(ext_summary %>% filter(Sshading == 2 &
                                                    scenario == "rcp85"),
                           aes(x = year, y = mean_pop,
                               colour = Hshading, fill = Hshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  # coord_cartesian(xlim = c(start_yr, 2050)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), cols = vars(model), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pop_size_85_shrub <- ggplot(ext_summary %>% filter(Hshading == 2 &
                                                     scenario == "rcp85"),
                            aes(x = year, y = mean_pop,
                                colour = Sshading, fill = Sshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  # coord_cartesian(xlim = c(start_yr, 2050)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), cols = vars(model), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



pop_sizes_45 <- pop_size_45_herb + pop_size_45_shrub + plot_layout(guides = "collect") & theme(text = element_text(size = 14))

pop_sizes_85 <- pop_size_85_herb + pop_size_85_shrub + plot_layout(guides = "collect") & theme(text = element_text(size = 14))

ggsave(pop_sizes_45, filename = file.path("results", "IBM_population_sizes_rcp4.5_per_model.png"), width = 12, height = 6)
ggsave(pop_sizes_85, filename = file.path("results", "IBM_population_sizes_rcp8.5_per_model.png"), width = 12, height = 6)

## Main MS result figure:  --------------------------------------------------------------------

ext_mainMS <- ext %>%
  group_by(Sshading, Hshading, scenario, locality, year) %>%
  summarise(
    mean_pop = mean(pop_size_plant, na.rm = TRUE),
    q1 = quantile(pop_size_plant, 0.1, na.rm = TRUE),
    q3 = quantile(pop_size_plant, 0.9, na.rm = TRUE),
    n = sum(!is.na(pop_size_plant)),
    se = sd(pop_size_plant, na.rm = TRUE) / sqrt(n),
    ci_lower = mean_pop - qt(0.975, df = n - 1) * se,
    ci_upper = mean_pop + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )


fig_45_herb <- ggplot(ext_mainMS %>% filter(Sshading == 2 &
                                                   scenario == "rcp45"),
                           aes(x = year, y = mean_pop,
                               colour = Hshading, fill = Hshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, colour = NA) +
  geom_line(linewidth = 0.8) +
  coord_cartesian(xlim = c(start_yr, end_yr)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous Shading")


fig_45_shrub <- ggplot(ext_mainMS %>% filter(Hshading == 2 &
                                                    scenario == "rcp45"),
                            aes(x = year, y = mean_pop,
                                colour = Sshading, fill = Sshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, colour = NA) +
  geom_line(linewidth = 0.8) +
  coord_cartesian(xlim = c(start_yr, end_yr)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub Shading")

fig_85_herb <- ggplot(ext_mainMS %>% filter(Sshading == 2 &
                                              scenario == "rcp85"),
                      aes(x = year, y = mean_pop,
                          colour = Hshading, fill = Hshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, colour = NA) +
  geom_line(linewidth = 0.8) +
  coord_cartesian(xlim = c(start_yr, end_yr)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous Shading")


fig_85_shrub <- ggplot(ext_mainMS %>% filter(Hshading == 2 &
                                               scenario == "rcp85"),
                       aes(x = year, y = mean_pop,
                           colour = Sshading, fill = Sshading)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, colour = NA) +
  geom_line(linewidth = 0.8) +
  coord_cartesian(xlim = c(start_yr, end_yr)) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub Shading")


fig_45 <- fig_45_herb + fig_45_shrub + plot_layout(guides = "collect") & 
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 300)),
      scale_y_continuous(limits = c(0, 3000)),
      scale_y_continuous(limits = c(0, 70)),
      scale_y_continuous(limits = c(0, 120))
    )) & theme(text = element_text(size = 14))
fig_85_supp <- fig_85_herb + fig_85_shrub + plot_layout(guides = "collect") & 
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 300)),
      scale_y_continuous(limits = c(0, 3000)),
      scale_y_continuous(limits = c(0, 70)),
      scale_y_continuous(limits = c(0, 120))
    )) & theme(text = element_text(size = 14))

fig_45

ggsave(fig_45, filename = file.path("results", "IBM_population_sizes_rcp4.5_2100.png"), width = 8, height = 8)
ggsave(fig_85_supp, filename = file.path("results", "IBM_population_sizes_rcp8.5_2100.png"), width = 8, height = 8)



## --------------------------------------------------------------------
## Seedbank dynamics
## --------------------------------------------------------------------

ext_summary_sb <- ext %>%
  tidyr::pivot_longer(
    cols = c(pop_size_sb1, pop_size_sb2),
    names_to = "seedbank",
    values_to = "pop_size"
  ) %>%
  group_by(Sshading, Hshading, scenario, locality, year, seedbank) %>%
  summarise(
    mean_pop = mean(pop_size, na.rm = TRUE),
    q1 = quantile(pop_size, 0.1, na.rm = TRUE),
    q3 = quantile(pop_size, 0.9, na.rm = TRUE),
    n = sum(!is.na(pop_size)),
    se = sd(pop_size, na.rm = TRUE) / sqrt(n),
    ci_lower = mean_pop - qt(0.975, df = n - 1) * se,
    ci_upper = mean_pop + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

pop_size_45_herb_sb1 <- ggplot(ext_summary_sb %>% filter(Sshading == 2 &
                                                          scenario == "rcp45" &
                                                          seedbank == "pop_size_sb1"),
                              aes(x = year, y = mean_pop,
                                  colour = Hshading, fill = Hshading,
                                  group = interaction(Hshading, seedbank))) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  scale_linetype_discrete(name = "Seed bank") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous seed bank") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pop_size_45_shrub_sb1 <- ggplot(ext_summary_sb %>% filter(Hshading == 2 &
                                                           scenario == "rcp45"&
                                                           seedbank == "pop_size_sb1"),
                               aes(x = year, y = mean_pop,
                                   colour = Sshading, fill = Sshading,
                                   group = interaction(Sshading, seedbank))) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  scale_linetype_discrete(name = "Seed bank") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub seed bank") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pop_size_45_herb_sb2 <- ggplot(ext_summary_sb %>% filter(Sshading == 2 &
                                                          scenario == "rcp45" &
                                                          seedbank == "pop_size_sb2"),
                              aes(x = year, y = mean_pop,
                                  colour = Hshading, fill = Hshading,
                                  group = interaction(Hshading, seedbank))) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  scale_linetype_discrete(name = "Seed bank") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Herbaceous seed bank") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pop_size_45_shrub_sb2 <- ggplot(ext_summary_sb %>% filter(Hshading == 2 &
                                                           scenario == "rcp45"&
                                                           seedbank == "pop_size_sb2"),
                               aes(x = year, y = mean_pop,
                                   colour = Sshading, fill = Sshading,
                                   group = interaction(Sshading, seedbank))) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(aes(yintercept = 10)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2025, 2100, by = 25)) +
  scale_y_continuous(limits = c(0, NA)) +
  viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  scale_linetype_discrete(name = "Seed bank") +
  facet_grid(rows = vars(locality), scales = "free_y") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("Shrub seed bank") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sb1_rcp45 <- pop_size_45_herb_sb1 + pop_size_45_shrub_sb1 + plot_layout(guides = "collect") &
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 50000)),
      scale_y_continuous(limits = c(0, 1500000)),
      scale_y_continuous(limits = c(0, 6000)),
      scale_y_continuous(limits = c(0, 20000))
    )) & theme(text = element_text(size = 14))

sb2_rcp45 <- pop_size_45_herb_sb2 + pop_size_45_shrub_sb2 + plot_layout(guides = "collect") &
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 3000)),
      scale_y_continuous(limits = c(0, 100000)),
      scale_y_continuous(limits = c(0, 800)),
      scale_y_continuous(limits = c(0, 1000))
    )) & theme(text = element_text(size = 14))

ggsave(sb1_rcp45, filename = "results/seedbank_dynamics_sb1.png", width = 6, height = 4)
ggsave(sb2_rcp45, filename = "results/seedbank_dynamics_sb2.png", width = 6, height = 4)



## --------------------------------------------------------------------
## Model validation
## --------------------------------------------------------------------

val_df <- readxl::read_xlsx("data/Population_sizes_validation.xlsx") %>%
  filter(stage == "total") %>%
  dplyr::select(-c(`Haknovec 2`, stage)) %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "obs_size") %>%
  mutate(locality = population |> 
           recode_values(
             "Císařská rokle" ~ "Cr",
             "Haknovec 1" ~ "Hk",
             "Kodská stěna" ~ "Ks",
             "Radotínské údolí" ~ "Ru"
           ))

val_45 <- ggplot(ext_mainMS %>% filter(Sshading == 2 &
                                         Hshading == 2 &
                                         scenario == "rcp45" &
                                         year < 2026),
                 aes(x = year)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.3, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.5, colour = NA) +
  geom_line(aes(y = mean_pop),
            linewidth = 0.8) +
  geom_point(data = val_df %>% filter(year > 2017), 
             aes(x = year, y = obs_size, shape = "Observed data"), colour = "red") + 
  coord_cartesian(xlim = c(start_yr, 2025)) +
  geom_hline(aes(yintercept = 10)) +
  # scale_x_continuous(expand = c(0,0)) +
  # viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  # viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  scale_shape_manual(name = NULL, values = c("Observed data" = 18)) +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("RCP 4.5") +
  theme(text = element_text(size = 14))

val_45_model <- (ggplot(ext_summary %>% filter(Sshading == 2 &
                                                 Hshading == 2 &
                                                 scenario == "rcp45" &
                                                 year < 2026),
                        aes(x = year)) +
                   geom_ribbon(aes(ymin = q1, ymax = q3, group = model), alpha = 0.3, colour = NA) +
                   geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = model), alpha = 0.5, colour = NA) +
                   geom_line(aes(y = mean_pop,
                                 linetype = model),
                             linewidth = 0.8) +
                   geom_point(data = val_df %>% filter(year > 2017), 
                              aes(x = year, y = obs_size, shape = "Observed data"), colour = "red") + 
                   coord_cartesian(xlim = c(start_yr, 2025)) +
                   geom_hline(aes(yintercept = 10)) +
                   # scale_x_continuous(expand = c(0,0)) +
                   # viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
                   # viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
                   scale_shape_manual(name = NULL, values = c("Observed data" = 18)) +
                   facet_grid(rows = vars(locality), scales = "free") +
                   xlab("Simulation year") + ylab("Population size") + ggtitle("RCP 4.5")) +
  theme(text = element_text(size = 14))

val_85 <- ggplot(ext_mainMS %>% filter(Sshading == 2 &
                                         Hshading == 2 &
                                         scenario == "rcp85" &
                                         year < 2026),
                 aes(x = year)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.3, colour = NA) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.5, colour = NA) +
  geom_line(aes(y = mean_pop),
            linewidth = 0.8) +
  geom_point(data = val_df %>% filter(year > 2017), 
             aes(x = year, y = obs_size, shape = "Observed data"), colour = "red") + 
  coord_cartesian(xlim = c(start_yr, 2025)) +
  geom_hline(aes(yintercept = 10)) +
  # scale_x_continuous(expand = c(0,0)) +
  # viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  # viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
  scale_shape_manual(name = NULL, values = c("Observed data" = 18)) +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size") + ggtitle("RCP 8.5") +
  theme(text = element_text(size = 14))

val_85_model <- (ggplot(ext_summary %>% filter(Sshading == 2 &
                                                 Hshading == 2 &
                                                 scenario == "rcp85" &
                                                 year < 2026),
                        aes(x = year)) +
                   geom_ribbon(aes(ymin = q1, ymax = q3, group = model), alpha = 0.3, colour = NA) +
                   geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = model), alpha = 0.5, colour = NA) +
                   geom_line(aes(y = mean_pop,
                                 linetype = model),
                             linewidth = 0.8) +
                   geom_point(data = val_df %>% filter(year > 2017), 
                              aes(x = year, y = obs_size, shape = "Observed data"), colour = "red") + 
                   coord_cartesian(xlim = c(start_yr, 2025)) +
                   geom_hline(aes(yintercept = 10)) +
                   # scale_x_continuous(expand = c(0,0)) +
                   # viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
                   # viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
                   scale_shape_manual(name = NULL, values = c("Observed data" = 18)) +
                   facet_grid(rows = vars(locality), scales = "free") +
                   xlab("Simulation year") + ylab("Population size") + ggtitle("RCP 8.5")) +
  theme(text = element_text(size = 14))


ggsave(val_45 + val_85 + plot_layout(guides = "collect"), 
       filename = file.path("results", "validation_all.png"), width = 8, height = 8)

ggsave(val_45_model + val_85_model + plot_layout(guides = "collect"), 
       filename = file.path("results", "validation_per_model.png"), width = 8, height = 8)




overall_performance <- ext_mainMS %>%
  # filter(Sshading == 2,
  #        Hshading == 2) %>%
  inner_join(
    val_df %>% select(year, locality, obs_size),
    by = c("year", "locality")
  ) %>%
  mutate(
    error = mean_pop - obs_size,
    abs_error = abs(error),
    sq_error = error^2,
    log_error = log1p(mean_pop) - log1p(obs_size)
  ) %>%
  group_by(locality) %>%
  summarise(
    RMSE = sqrt(mean((mean_pop - obs_size)^2)),
    MAE = mean(abs(mean_pop - obs_size)),
    log_RMSE = sqrt(mean(log_error^2)),
    rel_bias = mean((mean_pop - obs_size)/obs_size),
    coverage_80 = mean(obs_size >= q1 & obs_size <= q3),
    coverage_95_CI = mean(obs_size >= ci_lower & obs_size <= ci_upper),
    .groups = "drop"
  )

model_performance <-  ext_summary %>%
  # filter(Sshading == 2,
  #        Hshading == 2) %>%
  inner_join(
    val_df %>% select(year, locality, obs_size),
    by = c("year", "locality")
  ) %>%
  mutate(
    error = mean_pop - obs_size,
    abs_error = abs(error),
    sq_error = error^2,
    log_error = log1p(mean_pop) - log1p(obs_size)
  ) %>%
  group_by(locality, model, scenario) %>%
  summarise(
    RMSE = sqrt(mean((mean_pop - obs_size)^2)),
    MAE = mean(abs(mean_pop - obs_size)),
    log_RMSE = sqrt(mean(log_error^2)),
    rel_bias = mean((mean_pop - obs_size)/obs_size),
    coverage_80_PI = mean(obs_size >= q1 & obs_size <= q3),
    coverage_95_CI = mean(obs_size >= ci_lower & obs_size <= ci_upper),
    .groups = "drop"
  )

print(overall_performance)
print(model_performance)

summary(model_performance$coverage_80_PI)
summary(model_performance$log_RMSE)

## --------------------------------------------------------------------
## Visualise shading level distributions
## --------------------------------------------------------------------

ext_45_herb <- ggplot(df %>%
                        filter(scenario == "rcp45" & Sshading == 2)) +
  geom_boxplot(aes(y = yr_of_ext,
                   x = model,
                   fill = as.factor(Hshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction") + ggtitle("Herb shading") +
  coord_flip(ylim = c(2020,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  ) 
ext_45_herb

ext_45_shrub <- ggplot(df %>%
                         filter(scenario == "rcp45" & Hshading == 2)) +
  geom_boxplot(aes(y = yr_of_ext,
                   x = model,
                   fill = as.factor(Sshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Shrub shading") +
  coord_flip(ylim = c(2022,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  ) 
ext_45_shrub

ext_85_herb <- ggplot(df %>%
                        filter(scenario == "rcp85" & Sshading == 2)) +
  geom_boxplot(aes(y = yr_of_ext,
                   x = model,
                   fill = as.factor(Hshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Herb shading") +
  coord_flip(ylim = c(2022,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  ) 
ext_85_shrub <- ggplot(df %>%
                         filter(scenario == "rcp85" & Hshading == 2)) +
  geom_boxplot(aes(y = yr_of_ext,
                   x = model,
                   fill = as.factor(Sshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Shrub shading") +
  coord_flip(ylim = c(2022,2100)) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  )  

ext_45 <- ext_45_herb / ext_45_shrub + plot_layout(guides = "collect")
ext_85 <- ext_85_herb / ext_85_shrub + plot_layout(guides = "collect")

ggsave(ext_45, filename = file.path("results", "IBM_extinction_year_rcp4.5.png"), width = 6, height = 6)
ggsave(ext_85, filename = file.path("results", "IBM_extinction_year_rcp8.5.png"), width = 6, height = 6)


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

