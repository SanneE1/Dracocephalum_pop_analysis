library(tidyverse)
library(patchwork)
library(ggh4x)

### Script to plot population trajectories with alternative extinction thresholds
### (5, 10 [original] and 20 individuals)
rm(list = ls())
set.seed(2)

# Years of future simulations
start_yr = 2018
end_yr = 2100

# Extinction thresholds to show as horizontal lines
thresholds <- data.frame(threshold = factor(c(5, 10, 20)),
                         yintercept = c(5, 10, 20))

## --------------------------------------------------------------------
## Data prep
## --------------------------------------------------------------------

df <- lapply(list.files("results/ibm/", full.names = T), readRDS) %>%
  bind_rows() %>%
  mutate(locality = as.factor(locality),
         scenario = as.character(scenario))

df$scenario <- factor(df$scenario, levels = c("rcp45", "rcp85"))

names(df$pop_size) <- paste(df$locality, df$model, df$scenario, c(1:nrow(df)), sep = "_")

ext <- lapply(df$pop_size, function(x) x$plants) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>% tibble::rowid_to_column() %>%
  pivot_longer(num_range(prefix = "", range = c(1:(end_yr - (start_yr - 1)))), names_to = "time_step", values_to = "pop_size_plant") %>%
  mutate(time_step = as.numeric(time_step),
         Hshading = as.factor(Hshading),
         Sshading = as.factor(Sshading),
         year = time_step + (start_yr - 1))

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

rm(ext)

## --------------------------------------------------------------------
## Population dynamics
## --------------------------------------------------------------------

# Plot trajectories for one scenario and one shading type
plot_trajectory <- function(data, scen, shading_var, title) {

  # Hold the other shading type at level 2
  other_var <- setdiff(c("Hshading", "Sshading"), shading_var)

  data %>%
    filter(scenario == scen & .data[[other_var]] == 2) %>%
    ggplot(aes(x = year, y = mean_pop,
               colour = .data[[shading_var]], fill = .data[[shading_var]])) +
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.1, colour = NA) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(data = thresholds,
               aes(yintercept = yintercept, linetype = threshold),
               colour = "black", linewidth = 0.5) +
    coord_cartesian(xlim = c(start_yr, end_yr)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_linetype_manual(name = "Extinction threshold",
                          values = c("5" = "dotted", "10" = "solid", "20" = "dashed")) +
    viridis::scale_color_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
    viridis::scale_fill_viridis(option = "D", discrete = TRUE, direction = -1, name = "Shading level") +
    facet_grid(rows = vars(locality), scales = "free_y") +
    xlab("Simulation year") + ylab("Population size") + ggtitle(title)
}

# Fixed y-axis limits per locality (same as in plot_results.R)
y_scales <- facetted_pos_scales(
  y = list(
    scale_y_continuous(limits = c(0, 300)),
    scale_y_continuous(limits = c(0, 3000)),
    scale_y_continuous(limits = c(0, 70)),
    scale_y_continuous(limits = c(0, 120))
  ))

fig_45 <- plot_trajectory(ext_mainMS, "rcp45", "Hshading", "Herbaceous Shading") +
  plot_trajectory(ext_mainMS, "rcp45", "Sshading", "Shrub Shading") +
  plot_layout(guides = "collect") & y_scales & theme(text = element_text(size = 14))

fig_85 <- plot_trajectory(ext_mainMS, "rcp85", "Hshading", "Herbaceous Shading") +
  plot_trajectory(ext_mainMS, "rcp85", "Sshading", "Shrub Shading") +
  plot_layout(guides = "collect") & y_scales & theme(text = element_text(size = 14))

fig_45

ggsave(fig_45, filename = file.path("results", "IBM_population_sizes_rcp4.5_2100_thresholds.png"), width = 8, height = 8)
ggsave(fig_85, filename = file.path("results", "IBM_population_sizes_rcp8.5_2100_thresholds.png"), width = 8, height = 8)

