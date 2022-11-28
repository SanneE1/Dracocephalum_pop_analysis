### Script to plot results
rm(list = ls())
set.seed(2)

## --------------------------------------------------------------------
## Vital rates
## --------------------------------------------------------------------
source("R/functions_GAM.R")
lag = 24
vr <- readRDS("results/rds/VR_FLM.rds")

surv <- plot_spline_coeff(best_model = vr$surv,
                  lag = lag,
                  pet = T, shade = T,
                  vital_rate = "Survival")

growth <- plot_spline_coeff(best_model = vr$growth,
                            lag = lag,
                            pr = T, shade = T,
                            vital_rate = "Growth")

flowp <- plot_spline_coeff(best_model = vr$flower_p,
                           lag = lag,
                           pet = T, shade = T,
                           vital_rate = "Flower Probability") 
  
abp <- plot_spline_coeff(best_model = vr$abort_p,
                         lag = lag,
                         pr = T, shade = T,
                         vital_rate = "Seed Probability")

nseeds <- plot_spline_coeff(best_model = vr$n_seeds,
                            lag = lag,
                            pet = T, shade = T,
                            vital_rate = "Number of Seeds")

vr_plot <- surv + growth + flowp + abp + nseeds + guide_area() +
  plot_layout(ncol = 2, byrow = T, guides = "collect") + plot_annotation(tag_levels = "A") & 
  theme(legend.position = "bottom", legend.direction = "vertical",
        text = element_text(size = 12))

ggsave(vr_plot, filename = "results/vitalrate_plot.png", 
       width = 6, height = 6)

## --------------------------------------------------------------------
## Deterministic IPMs
## --------------------------------------------------------------------
det_ipm <- readRDS("results/rds/deterministic_ipm.rds")

ind_lambdas <- ggplot(det_ipm$individual_lambdas %>% mutate(year_t0 = as.integer(year_t0))) + 
  geom_point(aes(x = year_t0, y = lambda, colour = locality), size = 2) +
  scale_x_continuous(breaks = seq(2003,2020, 2), name = NULL) +
  ylab("population growth rate") +
  theme(axis.text.x = element_text(size = 12, angle = 45),
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        legend.position = "bottom")

ggsave(ind_lambdas, filename = "results/individual_lambdas.tiff",
       width = 7, height = 3.5, dpi = 300)


## --------------------------------------------------------------------
## Stochastic simulations
## --------------------------------------------------------------------
df <- read.csv("results/overview_lambda_env_levels.csv") %>%
  mutate(slope = round(slope, digits = 5))

# mod <- lm(lambda ~ locality + model + scenario + shading + slope, data = df)
# 
# drop1(mod, test = "Chisq") %>% 
#   rownames_to_column("dropped_covariate") %>%
#   rename(p_value = `Pr(>Chi)`) %>%
#   mutate(
#     p_adj_BH = p.adjust(p_value,method="BH",n=96),
#     p_adj_holm = p.adjust(p_value,method="holm",n=96),
#     include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.5, "yes", "no")
#   ) %>% as_tibble()


df[which(df$time == "hist"),c("model", "scenario")] <- "hist"
df <- df %>% mutate(model = factor(model, levels = c("hist","ACCESS1",
                         "CESM1", "CMCC", "MIROC5" )))


stoch_lam45_smooth <- ggplot(df %>% filter(scenario %in% c("hist", "rcp45"))) +
  geom_smooth(aes(x = shading, y = lambda, linetype = model, colour = model)) +
  facet_wrap(vars(locality), ncol = 4, scales = "free") +
  scale_colour_manual(name = "Model",
                      values = c("hist" = "blue",
                                 "ACCESS1" = "red",
                                 "CESM1" = "red",
                                 "CMCC" = "red",
                                 "MIROC5" = "red"),
                      labels = c("hist" = "Historical",
                                 "ACCESS1" = "ACCESS1",
                                 "CESM1" = "CESM1",
                                 "CMCC" = "CMCC",
                                 "MIROC5" = "MIROC5")) +
  scale_linetype_manual(values = c("hist" = 1,
                                "ACCESS1" = 2,
                                "CESM1" = 3,
                                "CMCC" = 4,
                                "MIROC5" = 5),
                     guide = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1:5)))) +
  ylab("stochastic population growth rate (log)") + xlab("shading level") +
  theme(legend.position = "bottom")

stoch_lam85_smooth <- ggplot(df %>% filter(scenario %in% c("hist", "rcp85"))) +
  geom_smooth(aes(x = shading, y = lambda, linetype = model, colour = model)) +
  facet_wrap(vars(locality), ncol = 4, scales = "free") +
  scale_colour_manual(name = "Model",
                      values = c("hist" = "blue",
                                 "ACCESS1" = "red",
                                 "CESM1" = "red",
                                 "CMCC" = "red",
                                 "MIROC5" = "red"),
                      labels = c("hist" = "Historical",
                                 "ACCESS1" = "ACCESS1",
                                 "CESM1" = "CESM1",
                                 "CMCC" = "CMCC",
                                 "MIROC5" = "MIROC5")) +
  scale_linetype_manual(values = c("hist" = 1,
                                   "ACCESS1" = 2,
                                   "CESM1" = 3,
                                   "CMCC" = 4,
                                   "MIROC5" = 5),
                        guide = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1:5)))) +
  ylab("stochastic population growth rate (log)") + xlab("shading level")+
  theme(legend.position = "bottom")


ggsave(stoch_lam45_smooth, filename = "results/stochastic_lamda_smooth_45.tiff",
       width = 8, height = 6)
ggsave(stoch_lam85_smooth, filename = "results/stochastic_lamda_smooth_85.tiff",
       width = 8, height = 6)


## Visualise ARIMA models for stochastic simulations

arima_clim <- readRDS("results/rds/ARIMA_clim_mods.rds")
n_it = 50



clim_sim_hist <- lapply(arima_clim$clim_hist_model, function(x) {
  a <- simulate(x, nsim = ((n_it * 12))) %>%
    ts(., start= c(2023,1), frequency = 12)
  data.frame(date = zoo::as.yearmon(time(a)),
             value = a)})
clim_sim_fut <- lapply(arima_clim$clim_future_model, function(x) {
  a <- simulate(x, nsim = ((n_it * 12))) %>%
    ts(., start= c(2023,1), frequency = 12)
  data.frame(date = zoo::as.yearmon(time(a)),
             value = a)})


Cr_pet_a <- ggplot() + 
  geom_line(data = clim_sim_fut$`Cr.ACCESS1-3.rcp45.pet_scaled`, aes(x = date, y = value, colour = "ACCESS1-3")) +
  geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) +
  labs(title = "PET")
Cr_tas_a <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.ACCESS1-3.rcp45.tas_scaled`, aes(x = date, y = value, colour = "ACCESS1-3")) +
  geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) +
  labs(title = "Average temperature")
Cr_pr_a <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.ACCESS1-3.rcp45.pr_scaled`, aes(x = date, y = value, colour = "ACCESS1-3")) +
  geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) +
  labs(title = "Precipitation")


Cr_pet_m <- ggplot() + 
  geom_line(data = clim_sim_fut$Cr.MIROC5.rcp45.pet_scaled , aes(x = date, y = value, colour = "MIROC5")) +
  geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) 
Cr_tas_m <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.MIROC5.rcp45.tas_scaled`, aes(x = date, y = value, colour = "MIROC5")) +
  geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) 
Cr_pr_m <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.MIROC5.rcp45.pr_scaled`, aes(x = date, y = value, colour = "MIROC5")) +
  geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) 

Cr_pet_c <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.CMCC-CM.rcp45.pet_scaled`, aes(x = date, y = value, colour = "CMCC-CM")) +
  geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) 
Cr_tas_c <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.CMCC-CM.rcp45.tas_scaled`, aes(x = date, y = value, colour = "CMCC-CM")) +
  geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) 
Cr_pr_c <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.CMCC-CM.rcp45.pr_scaled`, aes(x = date, y = value, colour = "CMCC-CM")) +
  geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) 

Cr_pet_e <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.CESM1-BGC.rcp45.pet_scaled`, aes(x = date, y = value, colour = "CESM1-BGC")) +
  geom_line(data = clim_sim_hist$Cr.pet_scaled, aes(x = date, y = value, colour = "Historic")) 
Cr_tas_e <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.CESM1-BGC.rcp45.tas_scaled`, aes(x = date, y = value, colour = "CESM1-BGC")) +
  geom_line(data = clim_sim_hist$Cr.tas_scaled, aes(x = date, y = value, colour = "Historic")) 
Cr_pr_e <- ggplot() +
  geom_line(data = clim_sim_fut$`Cr.CESM1-BGC.rcp45.pr_scaled`, aes(x = date, y = value, colour = "CESM1-BGC")) +
  geom_line(data = clim_sim_hist$Cr.pr_scaled, aes(x = date, y = value, colour = "Historic")) 

col_scale <- scale_colour_manual(name = "ARIMA model",
                                 values = c("Historic" = "#0033FF",
                                            "ACCESS1-3" = "#990000",
                                            "MIROC5" = "#FF0000",
                                            "CMCC-CM" = "#FF6600",
                                            "CESM1-BGC" = "#FFCC00"
                                 ))

sim_clim <- (Cr_pet_a + Cr_pet_m + Cr_pet_c + Cr_pet_e + 
                Cr_tas_a + Cr_tas_m + Cr_tas_c +  Cr_tas_e +
                Cr_pr_a + Cr_pr_m + Cr_pr_c + Cr_pr_e) * col_scale + 
  plot_layout(ncol = 3, byrow = F, guides = "collect") & theme(legend.position = "bottom")


ggsave("results/exploratory plots/arima_climate_sim1.tiff", sim_clim,
       width = 10, height = 7)




## --------------------------------------------------------------------
## Extinciton probability
## --------------------------------------------------------------------
df <- readRDS(file = "results/rds/extinction_probability.rds") %>%
  mutate(locality = as.factor(locality))
names(df$pop_size) <- paste(df$locality, df$model, df$scenario, c(1:nrow(df)), sep = "_")

levels(df$scenario) <- c(levels(df$scenario), "historical")
df$scenario[which(df$model == "No change")] <- "historical"

ext_pl <- lapply(df$pop_size, function(x) x$plants) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>% tibble::rowid_to_column() %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_plant") %>%
  mutate(time_step = as.numeric(time_step),
         shading = as.factor(shading),
         year = time_step + 2022) 

ext_sdl <- lapply(df$pop_size, function(x) x$sdl) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sdl") %>%
  mutate(time_step = as.numeric(time_step),
         shading = as.factor(shading),
         year = time_step + 2022)

ext_sb1 <- lapply(df$pop_size, function(x) x$sb1) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sb1") %>%
  mutate(time_step = as.numeric(time_step),
         shading = as.factor(shading),
         year = time_step + 2022)

ext_sb2 <- lapply(df$pop_size, function(x) x$sb2) %>% bind_rows() %>% t() %>%
  cbind(df[,c(1:6)], .) %>%
  pivot_longer(num_range(prefix = "", range = c(1:80)), names_to = "time_step", values_to = "pop_size_sb2") %>%
  mutate(time_step = as.numeric(time_step),
         shading = as.factor(shading),
         year = time_step + 2022)

ext <- ext_pl %>% cbind(.,
                        pop_size_sdl = ext_sdl$pop_size_sdl, 
                        pop_size_sb1 = ext_sb1$pop_size_sb1, 
                        pop_size_sb2 = ext_sb2$pop_size_sb2) %>%
  mutate(below_ext_size = ifelse(pop_size_plant < 10, T, F)) 

  

rm(ext_pl, ext_sdl, ext_sb1, ext_sb2)


pop_size <- ggplot(ext %>%
                     filter(model %in% c("No change", "ACCESS1") &
                              scenario != "rcp85")) +
  geom_line(aes(x = year,
                y = pop_size_plant,
                colour = shading,
                group = rowid), alpha = 0.5) +
  geom_hline(aes(yintercept = 10)) +
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1) +  
  facet_grid(rows = vars(locality), cols = vars(model), scales = "free") +
  xlab("Simulation year") + ylab("Population size")


ggsave(pop_size, filename = "results/plot_population_size_ext_simulation.tiff",
       width = 7, height = 8, units = "in", dpi = 400, type = "cairo")


# 
# mod_all <- lm(yr_of_ext ~ locality + model + scenario + shading + slope, data = df)
# 
# drop1(mod_all, test = "Chisq") %>% 
#   rownames_to_column("dropped_covariate") %>%
#   rename(p_value = `Pr(>Chi)`) %>%
#   mutate(
#     p_adj_BH = p.adjust(p_value,method="BH",n=96),
#     p_adj_holm = p.adjust(p_value,method="holm",n=96),
#     include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.5, "yes", "no")
#   ) %>% as_tibble()
# 
# mod <- lm(yr_of_ext ~ locality + model + shading, data = df)
# summary(mod)



# probability of extinction (0-1) vs time

ext_p <- ext %>% 
  group_by(locality, model, scenario, shading, slope, year) %>%
  summarise(runs = n(),
            extinct = sum(below_ext_size),
            p_extinct = extinct / runs)


ggplot(ext_p) + 
  # geom_line(aes(x = year, 
  #               y = p_extinct, 
  #               colour = shading, 
  #               linetype = model,
  #               group = interaction(model, scenario, slope, shading))) +
  geom_smooth(aes(x = year, 
                y = p_extinct, 
                colour = shading, 
                # linetype = model,
                group = interaction(scenario, slope, shading)), 
              method="glm",
              method.args=list(family="binomial"), 
              se = F) +
  viridis::scale_color_viridis(option = "D", discrete = T, direction = -1) +  
  facet_grid(rows = vars(locality), cols = vars(model))





# Visualise climate sequences --------------------------------------------------------------------
clim <- read.csv("data/CHELSA_future_ts_formatted.csv") %>% 
  mutate(date = zoo::as.yearmon(paste(year, month, sep = "-"))) %>%
  filter(year != 1894)



ggplot(clim %>% filter(scenario != "historical")) + 
  geom_line(aes(x = date, y = pr_scaled, colour = model )) +
  facet_grid(rows = vars(model), cols = vars(scenario))

ggplot(clim %>% filter(scenario != "historical")) + 
  geom_line(aes(x = date, y = tas_scaled, colour = model )) +
  facet_grid(rows = vars(model), cols = vars(scenario))

ggplot(clim %>% filter(scenario != "historical")) + 
  geom_line(aes(x = date, y = pet_scaled, colour = model )) +
  facet_grid(rows = vars(model), cols = vars(scenario))

