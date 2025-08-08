library(tidyverse)
library(glmnet)
library(patchwork)

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

plot_coef_clim_shading(surv_mod, save = T, name = "Survival_clim_shading_coef")
plot_coef_clim_shading(grow_mod, save = T, name = "Growth_clim_shading_coef")
plot_coef_clim_shading(flowp_mod, save = T, name = "Flower_P_clim_shading_coef")
plot_coef_clim_shading(seedp_mod, save = T, name = "Seed_P_clim_shading_coef")
plot_coef_clim_shading(seedn_mod, save = T, name = "Seed_N_clim_shading_coef")

## --------------------------------------------------------------------
## Stochastic simulations
## --------------------------------------------------------------------

# files <- list.files("results/stoch_ipms/", full.names = T)
# 
# df <- lapply(files, read.csv) %>% bind_rows
# 
# # mod <- lm(lambda ~ locality + model + scenario + shading + slope, data = df)
# # 
# # drop1(mod, test = "Chisq") %>% 
# #   rownames_to_column("dropped_covariate") %>%
# #   rename(p_value = `Pr(>Chi)`) %>%
# #   mutate(
# #     p_adj_BH = p.adjust(p_value,method="BH",n=96),
# #     p_adj_holm = p.adjust(p_value,method="holm",n=96),
# #     include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.5, "yes", "no")
# #   ) %>% as_tibble()
# 
# 
# df[which(df$time == "hist"),c("model", "scenario")] <- "hist"
# df <- df %>% mutate(model = factor(model, levels = c("hist","ACCESS1",
#                                                      "CESM1", "CMCC", "MIROC5" )))
# 
# 
# stoch_lam45_smooth <- ggplot(na.omit(df %>% filter(scenario %in% c("hist", "rcp45")))) +
#   geom_line(aes(x = shading, y = lambda, linetype = model, colour = model), size = 0.8) +
#   geom_point(aes(x = shading, y = lambda, colour = model)) +
#   facet_wrap(vars(locality), ncol = 4) +
#   scale_colour_manual(name = "Climate model",
#                       values = c("hist" = "#000000",
#                                  "ACCESS1" = "#D55E00",
#                                  "CESM1" = "#0072B2",
#                                  "CMCC" = "#56B4E9",
#                                  "MIROC5" = "#E69F00"),
#                       labels = c("hist" = "Historical",
#                                  "ACCESS1" = "ACCESS1",
#                                  "CESM1" = "CESM1",
#                                  "CMCC" = "CMCC",
#                                  "MIROC5" = "MIROC5")) +
#   scale_linetype_manual(values = c("hist" = 1,
#                                    "ACCESS1" = 2,
#                                    "CESM1" = 3,
#                                    "CMCC" = 4,
#                                    "MIROC5" = 5),
#                         guide = "none") +
#   guides(colour = guide_legend(override.aes = list(linetype = c(1:5)),
#                                title.position = "top", title.hjust = 0.5)) +
#   ylab("Stochastic population growth rate (ln)") + xlab("Shading level") +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal") 
# 
# stoch_lam85_smooth <- ggplot(na.omit(df %>% filter(scenario %in% c("hist", "rcp85")))) +
#   geom_line(aes(x = shading, y = lambda, linetype = model, colour = model), size = 0.8) +
#   geom_point(aes(x = shading, y = lambda, colour = model)) +
#   facet_wrap(vars(locality), ncol = 4) +
#   scale_colour_manual(name = "Model",
#                       values = c("hist" = "#000000",
#                                  "ACCESS1" = "#D55E00",
#                                  "CESM1" = "#0072B2",
#                                  "CMCC" = "#56B4E9",
#                                  "MIROC5" = "#E69F00"),
#                       labels = c("hist" = "Historical",
#                                  "ACCESS1" = "ACCESS1",
#                                  "CESM1" = "CESM1",
#                                  "CMCC" = "CMCC",
#                                  "MIROC5" = "MIROC5")) +
#   scale_linetype_manual(values = c("hist" = 1,
#                                    "ACCESS1" = 2,
#                                    "CESM1" = 3,
#                                    "CMCC" = 4,
#                                    "MIROC5" = 5),
#                         guide = "none") +
#   guides(colour = guide_legend(override.aes = list(linetype = c(1:5)),
#                                title.position = "top", title.hjust = 0.5)) +
#   ylab("Stochastic population growth rate (ln)") + xlab("Shading level")+
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal")
# 
# 
# ggsave(stoch_lam45_smooth, filename = "results/stochastic_lamda_smooth_45.tiff",
#        width = 8, height = 6)
# ggsave(stoch_lam85_smooth, filename = "results/stochastic_lamda_smooth_85.png",
#        width = 8, height = 6)
# 
# 
# ## Visualise ARIMA models for stochastic simulations
# 
# arima_clim <- readRDS("results/rds/ARIMA_clim_mods.rds")
# n_it = 50
# 
# 
# 
# clim_sim_hist <- lapply(arima_clim$clim_hist_model, function(x) {
#   a <- simulate(x, nsim = ((n_it * 12))) %>%
#     ts(., start= c(2023,1), frequency = 12)
#   data.frame(date = zoo::as.yearmon(time(a)),
#              value = a)})
# clim_sim_fut <- lapply(arima_clim$clim_future_model, function(x) {
#   a <- simulate(x, nsim = ((n_it * 12))) %>%
#     ts(., start= c(2023,1), frequency = 12)
#   data.frame(date = zoo::as.yearmon(time(a)),
#              value = a)})
# 
# 
# Cr_pet_a <- ggplot() + 
#   geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.ACCESS1-3.rcp45.pet_scaled`, aes(x = date, y = value, colour = "ACCESS1-3")) +
#   labs(title = "PET")
# Cr_tas_a <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.ACCESS1-3.rcp45.tas_scaled`, aes(x = date, y = value, colour = "ACCESS1-3")) +
#   labs(title = "Average temperature")
# Cr_pr_a <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.ACCESS1-3.rcp45.pr_scaled`, aes(x = date, y = value, colour = "ACCESS1-3")) +
#   labs(title = "Precipitation")
# 
# 
# Cr_pet_m <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) + 
#   geom_line(data = clim_sim_fut$Cr.MIROC5.rcp45.pet_scaled , aes(x = date, y = value, colour = "MIROC5")) 
# Cr_tas_m <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.MIROC5.rcp45.tas_scaled`, aes(x = date, y = value, colour = "MIROC5")) 
# Cr_pr_m <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.MIROC5.rcp45.pr_scaled`, aes(x = date, y = value, colour = "MIROC5")) 
# 
# Cr_pet_c <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.CMCC-CM.rcp45.pet_scaled`, aes(x = date, y = value, colour = "CMCC-CM"))
# Cr_tas_c <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.CMCC-CM.rcp45.tas_scaled`, aes(x = date, y = value, colour = "CMCC-CM"))
# Cr_pr_c <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.CMCC-CM.rcp45.pr_scaled`, aes(x = date, y = value, colour = "CMCC-CM")) 
# 
# Cr_pet_e <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.CESM1-BGC.rcp45.pet_scaled`, aes(x = date, y = value, colour = "CESM1-BGC")) 
# Cr_tas_e <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.CESM1-BGC.rcp45.tas_scaled`, aes(x = date, y = value, colour = "CESM1-BGC")) 
# Cr_pr_e <- ggplot() +
#   geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) +
#   geom_line(data = clim_sim_fut$`Cr.CESM1-BGC.rcp45.pr_scaled`, aes(x = date, y = value, colour = "CESM1-BGC")) 
# 
# col_scale <- scale_colour_manual(name = "ARIMA model",
#                                  values = c("Historic" = "#0033FF",
#                                             "ACCESS1-3" = "#990000",
#                                             "MIROC5" = "#FF0000",
#                                             "CMCC-CM" = "#FF6600",
#                                             "CESM1-BGC" = "#FFCC00"
#                                  ))
# 
# sim_clim <- (Cr_pet_a + Cr_pet_m + Cr_pet_c + Cr_pet_e + 
#                Cr_tas_a + Cr_tas_m + Cr_tas_c +  Cr_tas_e +
#                Cr_pr_a + Cr_pr_m + Cr_pr_c + Cr_pr_e) * col_scale + 
#   plot_layout(ncol = 3, byrow = F, guides = "collect") & theme(legend.position = "bottom")
# 
# 
# ggsave("results/exploratory plots/arima_climate_sim1.tiff", sim_clim,
#        width = 10, height = 7)
# 
# 
# 

## --------------------------------------------------------------------
## Extinciton probability
## --------------------------------------------------------------------
source("R/functions_extinction_plots.R")

df <- lapply(list.files("results/ibm/", full.names = T), readRDS) %>% 
  bind_rows() %>%
  mutate(locality = as.factor(locality),
         yr_of_ext = yr_of_ext + 2022,
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
         year = time_step + 2022) 

ext_sdl <- lapply(df$pop_size, function(x) x$sdl) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sdl") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2022)

ext_sb1 <- lapply(df$pop_size, function(x) x$sb1) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sb1") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2022)

ext_sb2 <- lapply(df$pop_size, function(x) x$sb2) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sb2") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + 2022)

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
                             filter(scenario == "rcp45" & Sshading == 0)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = Hshading,
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model,
                  colour = Hshading)) +
  coord_cartesian(xlim = c(2023, 2070)) + 
  geom_hline(aes(yintercept = 10)) +
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size")
pop_size_45_herb

pop_size_45_shrub <- ggplot(ext %>%
                              filter(scenario == "rcp45" & Hshading == 0)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = Sshading,
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model,
                  colour = Sshading)) +
  coord_cartesian(xlim = c(2023, 2070)) + 
  geom_hline(aes(yintercept = 10)) +
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size")
pop_size_45_shrub

pop_size_85_herb <- ggplot(ext %>%
                             filter(scenario == "rcp85" & Sshading == 0)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = interaction(Hshading, Sshading),
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model,
                  colour = interaction(Hshading, Sshading))) +
  coord_cartesian(xlim = c(2023, 2070)) + 
  geom_hline(aes(yintercept = 10)) +
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size")
pop_size_85_herb

pop_size_85_shrub <- ggplot(ext %>%
                              filter(scenario == "rcp85" & Hshading == 0)) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = interaction(Hshading, Sshading),
                group = rowid), alpha = 0.05) +
  geom_smooth(aes(x = year,
                  y = pop_size_plant,
                  # linetype = model,
                  colour = interaction(Hshading, Sshading))) +
  coord_cartesian(xlim = c(2023, 2070)) + 
  geom_hline(aes(yintercept = 10)) +
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") +
  facet_grid(rows = vars(locality), scales = "free") +
  xlab("Simulation year") + ylab("Population size")
pop_size_85_shrub


pop_size_45_herb + pop_size_45_shrub + plot_layout(guides = "collect")

pop_size_85_herb + pop_size_85_shrub + plot_layout(guides = "collect")





ext_45_herb <- ggplot(df %>%
                         filter(scenario == "rcp45" & Sshading == 0)) +
  geom_boxplot(aes(y = yr_of_ext,
                   # x = model,
                   fill = as.factor(Hshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction") + ggtitle("Herb shading") +
  coord_flip(ylim = c(2022,2070)) + 
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
                              filter(scenario == "rcp45" & Hshading == 0)) +
  geom_boxplot(aes(y = yr_of_ext,
                   # x = model,
                   fill = as.factor(Sshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Shrub shading") +
  coord_flip(ylim = c(2022,2070)) + 
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
                        filter(scenario == "rcp85" & Sshading == 0)) +
  geom_boxplot(aes(y = yr_of_ext,
                   x = model,
                   fill = as.factor(Hshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction") + ggtitle("Herb shading") +
  coord_flip(ylim = c(2022,2070)) + 
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
                             filter(scenario == "rcp85" & Hshading == 0)) +
  geom_boxplot(aes(y = yr_of_ext,
                   x = model,
                   fill = as.factor(Sshading)),
               alpha = 0.5) +
  viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1, name = "Shading level") + 
  ylab("Year of extinction")  + ggtitle("Shrub shading") +
  coord_flip(ylim = c(2022,2070)) + 
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

ext_45_herb / ext_45_shrub + plot_layout(guides = "collect")
ext_85_herb / ext_85_shrub + plot_layout(guides = "collect")


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

shading_df <- data.frame(lambda = rep(c(0,2,4,6), each = 500)) %>%
  rowwise() %>%
  mutate(shading = rpois(1, lambda)) 

shading_dist <- ggplot(df, aes(x = shading_t0)) + 
  geom_histogram(aes(y = ..density..), binwidth = 1) + 
  geom_vline(aes(xintercept = mean, group = population)) +
  facet_grid(rows = vars(population)) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Locality", breaks = NULL, labels = NULL)) + 
  xlab("Shading distribution") + ylab("") + 
  coord_cartesian(xlim = c(0,17))

shading_mod <- ggplot(shading_df) +
  geom_density(aes(x = shading, y = after_stat(scaled), 
                   colour = as.factor(lambda), 
                   fill = as.factor(lambda)), 
               alpha = 0.5,
               adjust = 2) +
  xlab("Modelled shading distribution") + ylab("") +
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


# rcp_scenarios <- c("4.5", "8.5")
# names(rcp_scenarios) <- c("rcp45", "rcp85")
# 
# clim_seq_pr <- ggplot(a) +
#   geom_line(aes(x = date, y = pr_scaled, colour = model )) +
#   geom_line(aes(x = date, y = hist_precip, linetype = "historic"), alpha = 0.8) +
#   facet_grid(rows = vars(model), cols = vars(scenario), labeller = labeller(scenario = rcp_scenarios)) +
#   scale_colour_manual(values = c("CESM1-BGC" = "#0072B2",
#                                  "CMCC-CM" = "#56B4E9",
#                                  "MIROC5" = "#E69F00",
#                                  "ACCESS1-3" = "#D55E00"),
#                       labels = c("CESM1-BGC" = "CESM1",
#                                  "CMCC-CM" = "CMCC",
#                                  "MIROC5" = "MIROC5",
#                                  "ACCESS1-3" = "ACCESS1"),
#                       name = "Predicted Climate \n(2070-2100)") +
#   scale_linetype_manual(values = c("historic" = 1),
#                         labels = c("historic" = "Historic climate \n(1970-2000)"),
#                         name = "",
#                         guide = guide_legend(override.aes = list(color = "black"),
#                                              order = 1)) +
#   scale_y_continuous(sec.axis = sec_axis(~ . , name = "Climate Model", breaks = NULL, labels = NULL)) +
#   scale_x_continuous(sec.axis = sec_axis(~ . , name = "RCP scenario", breaks = NULL, labels = NULL)) +
#   coord_cartesian(xlim = c(2070,2100)) + ylab("Precipitation anomaly") +
#   theme(axis.text.x = element_blank(),
#         axis.title.x.bottom = element_blank(),
#         axis.ticks.x = element_blank())
# 
# 
# clim_seq_tas <- ggplot(a) +
#   geom_line(aes(x = date, y = tas_scaled, colour = model )) +
#   geom_line(aes(x = date, y = hist_tas, linetype = "historic"), alpha = 0.8) +
#   facet_grid(rows = vars(model), cols = vars(scenario), labeller = labeller(scenario = rcp_scenarios)) +
#   scale_colour_manual(values = c("CESM1-BGC" = "#0072B2",
#                                  "CMCC-CM" = "#56B4E9",
#                                  "MIROC5" = "#E69F00",
#                                  "ACCESS1-3" = "#D55E00"),
#                       labels = c("CESM1-BGC" = "CESM1",
#                                  "CMCC-CM" = "CMCC",
#                                  "MIROC5" = "MIROC5",
#                                  "ACCESS1-3" = "ACCESS1"),
#                       name = "Predicted Climate \n(2070-2100)") +
#   scale_linetype_manual(values = c("historic" = 1),
#                         labels = c("historic" = "Historic climate \n(1970-2000)"),
#                         name = "",
#                         guide = guide_legend(override.aes = list(color = "black"),
#                                              order = 1)) +
#   scale_y_continuous(sec.axis = sec_axis(~ . , name = "Climate Model", breaks = NULL, labels = NULL)) +
#   scale_x_continuous(sec.axis = sec_axis(~ . , name = "RCP scenario", breaks = NULL, labels = NULL)) +
#   coord_cartesian(xlim = c(2070,2100)) + ylab("Temperature anomaly") +
#   theme(axis.text.x = element_blank(),
#         axis.title.x.bottom = element_blank(),
#         axis.ticks.x = element_blank())
# 
# 
# clim_seq_pet <- ggplot(a) +
#   geom_line(aes(x = date, y = pet_scaled, colour = model )) +
#   geom_line(aes(x = date, y = hist_pet, linetype = "historic"), alpha = 0.8) +
#   facet_grid(rows = vars(model), cols = vars(scenario), labeller = labeller(scenario = rcp_scenarios)) +
#   scale_colour_manual(values = c("CESM1-BGC" = "#0072B2",
#                                  "CMCC-CM" = "#56B4E9",
#                                  "MIROC5" = "#E69F00",
#                                  "ACCESS1-3" = "#D55E00"),
#                       labels = c("CESM1-BGC" = "CESM1",
#                                  "CMCC-CM" = "CMCC",
#                                  "MIROC5" = "MIROC5",
#                                  "ACCESS1-3" = "ACCESS1"),
#                       name = "Predicted Climate \n(2070-2100)") +
#   scale_linetype_manual(values = c("historic" = 1),
#                         labels = c("historic" = "Historic climate \n(1970-2000)"),
#                         name = "",
#                         guide = guide_legend(override.aes = list(color = "black"),
#                                              order = 1)) +
#   scale_y_continuous(sec.axis = sec_axis(~ . , name = "Climate Model", breaks = NULL, labels = NULL)) +
#   scale_x_continuous(sec.axis = sec_axis(~ . , name = "RCP scenario", breaks = NULL, labels = NULL)) +
#   coord_cartesian(xlim = c(2070,2100)) + ylab("PET anomaly") +
#   theme(axis.text.x = element_blank(),
#         axis.title.x.bottom = element_blank(),
#         axis.ticks.x = element_blank())
# 
# ggsave("results/clim_seq_pr.png", clim_seq_pr,
#        width = 7, height = 9)
# ggsave("results/clim_seq_pet.png", clim_seq_pet,
#        width = 7, height = 9)
# ggsave("results/clim_seq_tas.png", clim_seq_tas,
#        width = 7, height = 9)


## --------------------------------------------------------------------
## Compare future climate mean, sd and autocorrelation
## --------------------------------------------------------------------
# 
# b <- bind_rows(clim_fut, clim_hist %>% 
#                  mutate(scenario = "Historical", model = "Historical") %>%
#                  dplyr::select(-hist_tas) %>%
#                  rename(pr_scaled = hist_precip, 
#                         pet_scaled = hist_pet)) %>% 
#   mutate(model = factor(model, levels = c("Historical", "CESM1-BGC", "CMCC-CM", "MIROC5", "ACCESS1-3")))
# 
# plot_pet_dist <- ggplot(b, aes(x = model, y = pet_scaled, fill = scenario)) + 
#   geom_boxplot()
# 
# 
# plot_pr_dist <- ggplot(b, aes(x = model, y = pr_scaled, fill = scenario)) + 
#   geom_boxplot() 
# 
# ggsave("results/clim_distribution_pet.png", plot_pet_dist,
#        width = 5, height = 5)
# ggsave("results/clim_distribution_pr.png", plot_pr_dist,
#        width = 5, height = 5)
# 
# 
# c <- b %>% 
#   mutate(date = tsibble::yearmonth(date)) %>% 
#   tsibble::as_tsibble(., index = "date", key = c(locality, model, scenario)) 
# 
# c %>%
#   filter(locality == "Cr") %>%
#   group_by(model, scenario) %>% 
#   feasts::ACF(pet_scaled, lag_max = 13) %>% 
#   autoplot() + facet_grid(rows = vars(model), cols = vars(scenario))

## --------------------------------------------------------------------
## Calculate vital rate covariation
## --------------------------------------------------------------------
# 
# FLM_clim_predict <- function(model, precip_vec, pet_vec, 
#                              shading) {
#   lag_vec = c(0:24)
#   ## dummy datalist, won't be using the predictions for n_stems : year_t0 but these need to be provided for predict function
#   new_data <- list(
#     ln_stems_t0 = 2,  
#     population = "CR",
#     tot_shading_t0 = 0,
#     soil_depth = 0,
#     year_t0 = 2016,
#     tot_shading_m = matrix(rep(shading, length(lag_vec)), nrow = 1),
#     lags = matrix(lag_vec, nrow = 1),
#     pr_scaledcovar = matrix(precip_vec, nrow = 1),
#     pet_scaledcovar = matrix(pet_vec, nrow = 1)
#   )
#   
#   pt <- mgcv::predict.gam(model, new_data, type = "terms")
#   return(
#     sum(pt[, grep("scaledcovar", attributes(pt)[[2]][[2]], value = T)])
#   )
# }
# 
# vr_flm <- readRDS("results/rds/VR_FLM.rds")
# 
# cor_df <- data.frame(model = character(),
#                      shading = vector(),
#                      surv = vector(),
#                      growth = vector(),
#                      flower_p = vector(),
#                      abort_p = vector(),
#                      n_seeds = vector())
# 
# d <- b %>% filter(scenario != "rcp85")
# 
# for (cm in unique(d$model)) {
#   
#   clim_vec <- d %>% 
#     filter(model == cm & locality == "Cr") %>% 
#     dplyr::select(date, pr_scaled, pet_scaled) %>% distinct() %>% arrange(date)
#   
#   pr_vec <- list()
#   pet_vec <- list()
#   
#   for(w in c(1:(nrow(clim_vec)-5))) {
#     
#     pr_vec <- append( pr_vec, 
#                       list(matrix(rev(clim_vec$pr_scaled[c(w:(w+24))]), nrow = 1, ncol = (25), byrow = T) )) 
#     pet_vec <- append( pet_vec, 
#                        list(matrix(rev(clim_vec$pet_scaled[c(w:(w+24))]), nrow = 1, ncol = (25), byrow = T) ))
#   }
#   
#   for(s in c(0,2,4,6)){
#     
#     temp_df <- data.frame(
#       model = cm,
#       shading = s,
#       surv = sapply(as.list(c(1:length(pr_vec))), 
#                     function(x) FLM_clim_predict(model = vr_flm$surv, 
#                                                  precip_vec = pr_vec[[x]], 
#                                                  pet_vec = pet_vec[[x]], 
#                                                  shading = s)),
#       growth = sapply(as.list(c(1:length(pr_vec))), 
#                       function(x) FLM_clim_predict(model = vr_flm$growth, 
#                                                    precip_vec = pr_vec[[x]], 
#                                                    pet_vec = pet_vec[[x]], 
#                                                    shading = s)),
#       flower_p = sapply(as.list(c(1:length(pr_vec))), 
#                         function(x) FLM_clim_predict(model = vr_flm$flower_p, 
#                                                      precip_vec = pr_vec[[x]], 
#                                                      pet_vec = pet_vec[[x]], 
#                                                      shading = s)),
#       abort_p = sapply(as.list(c(1:length(pr_vec))), 
#                        function(x) FLM_clim_predict(model = vr_flm$abort_p, 
#                                                     precip_vec = pr_vec[[x]], 
#                                                     pet_vec = pet_vec[[x]], 
#                                                     shading = s)),
#       n_seeds = sapply(as.list(c(1:length(pr_vec))), 
#                        function(x) FLM_clim_predict(model = vr_flm$n_seeds, 
#                                                     precip_vec = pr_vec[[x]], 
#                                                     pet_vec = pet_vec[[x]], 
#                                                     shading = s)))
#     
#     
#     cor_df <- bind_rows(
#       cor_df,
#       temp_df
#     )
#   }
#   
# }
# 
# saveRDS(cor_df, "results/rds/correlation_df.rds")
# covar_df <- readRDS("results/rds/correlation_df.rds")
# 
# get_upper_tri <- function(cormat){
#   cormat[lower.tri(cormat)]<- NA
#   return(cormat)
# }
# 
# 
# cor_df <- split(covar_df, interaction(covar_df$model, covar_df$shading)) %>%
#   lapply(., function(x){ 
#     df <- round(cor(x %>% dplyr::select(-c(model, shading))), 3)
#     df <- get_upper_tri(df)
#     df <- reshape2::melt(df, na.rm = T)
#   }) %>% bind_rows(., .id = "name")
# 
# cor_df$model <- stringr::str_extract(cor_df$name, "^[^\\.]*")
# cor_df$shading <- stringr::str_extract(cor_df$name, "\\d$")
# 
# cor_df <- cor_df %>%
#   mutate(model = factor(model, levels = c("Historical", "CESM1-BGC", "CMCC-CM", "MIROC5", "ACCESS1-3")))
# 
# cor_df_hist <- cor_df %>% filter(model == "Historical") %>%
#   rename(value_hist = value) %>% select(-c(name, model))
# 
# cor_df <- left_join(cor_df, cor_df_hist) %>%
#   rowwise() %>%
#   mutate(cor_diff = value - value_hist)
# 
# plot_cor <- ggplot(cor_df, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="Pearson\nCorrelation") +
#   facet_grid(rows = vars(shading), cols = vars(model)) +
#   # theme_minimal()+ # minimal theme
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1))+
#   coord_fixed() + 
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "bottom",
#     legend.direction = "horizontal")+
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
#                                title.position = "top", title.hjust = 0.5))
# 
# plot_cor_diff <- ggplot(cor_df, aes(as.numeric(Var2), as.numeric(Var1), fill = cor_diff))+
#   geom_tile(color = "white")+
#   geom_text(aes(as.numeric(Var2), as.numeric(Var1), label = value), color = "black", size = 4) +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-0.5,0.5), space = "Lab", 
#                        name="Difference in Pearson Correlation \nbetween with historical climate") +
#   facet_grid(rows = vars(shading), cols = vars(model)) +
#   # theme_minimal()+ # minimal theme
#   # coord_fixed() + 
#   scale_y_continuous(breaks = 1:5,
#                      labels = levels(cor_df$Var2),
#                      sec.axis = sec_axis(~ . , name = "Shading Level", breaks = NULL, labels = NULL)) +
#   scale_x_continuous(breaks = 1:5,
#                      labels = levels(cor_df$Var2),
#                      sec.axis = sec_axis(~ . , name = "Climate Model", breaks = NULL, labels = NULL)) +
#   theme(
#     axis.title.x.bottom = element_blank(),
#     axis.title.y.left = element_blank(),
#     axis.title.x.top = element_text(size= 14),
#     axis.title.y.right = element_text(size= 14),
#     axis.text = element_text(size = 14),
#     axis.text.x = element_text(angle = 45, vjust = 1,
#                                size = 14, hjust = 1),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "bottom",
#     legend.direction = "horizontal")+
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
#                                title.position = "top"))
# 
# plot_cor_diff
# 
# ggsave("results/vr_cor_plot_diff.png", plot_cor_diff,
#        width = 14, height = 10)

