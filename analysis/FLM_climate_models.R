### FLM climate model for Growth

VR_climate_FLM <- function(data_for_modeling, transformed_climate, lag = 48) {
  
  # lag = number of months prior to census to include
  
  ## Format dataframe:
  ### response & covariates and add climate variable as lagg since census in columns 
  ### i.e. temp.00 = month of census, temp.01 = temperature one month before census etc.
  ### Looks like for the GAM you need to use a variable name ("temp") and then "." and then add the numbers
  
  demo_data <- read.csv(data_for_modeling) %>%
    mutate(population = factor(population)) %>%
    filter(is.finite(ln_stems_t0)) 
  
  clim_data <- read.csv(transformed_climate) %>%
    climate_wider_for_gam(clim_data = ., 
                          variables = c("Tavg_scaled", "Prec_scaled"), 
                          demo_data = demo_data, 
                          response_t1 = T,
                          lag = lag)
  
  ### Merge the two dataframes and keep only entries for which we have full climate data available
  
  data <- left_join(demo_data, clim_data) %>% 
    filter(complete.cases(Tavg_scaled.0))
  
  # -------------------------------------------------------
  # GAM Models
  # -------------------------------------------------------
  
  
  
  ## -------------------------------------------------------
  ## Survival
  ## -------------------------------------------------------
  
  model_survival <- mgcv::gam(survival_t1 ~ ln_stems_t0 + population +
                                herb_shading_t0 + shrub_shading_t0 +
                                slope +
                                s(year_t0, bs="re") +
                                s(lags, k=lag, by= Tavg_scaledcovar, bs="cs") + 
                                s(lags, k=lag, by= Prec_scaledcovar, bs="cs"),
                              data=data,
                              family = binomial(link = "logit"),
                              method="GCV.Cp", 
                              gamma=1.2)
  
  summary(model_survival)
  
  getBeta.data <- data.frame(lags= c(0:lag), 
                             Tavg_scaledcovar=1, Prec_scaledcovar = 0,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, slope = 14,
                             herb_shading_t0 = 5, shrub_shading_t0 = 5)
  terms.data <- mgcv::predict.gam(model_survival, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_surv_temp <- data.frame(lag=getBeta.data$lags,
                                beta=terms.data$fit[,"s(lags):Tavg_scaledcovar"],
                                se=terms.data$se[,"s(lags):Tavg_scaledcovar"])
  
  getBeta.data <- data.frame(lags=c(0:lag), 
                             Tavg_scaledcovar=0, Prec_scaledcovar = 1,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, slope = 14, 
                             herb_shading_t0 = 5, shrub_shading_t0 = 5)
  terms.data <- mgcv::predict.gam(model_survival, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_surv_prec <- data.frame(lag=getBeta.data$lags,
                                beta=terms.data$fit[,"s(lags):Prec_scaledcovar"],
                                se=terms.data$se[,"s(lags):Prec_scaledcovar"])
  
  
  plot_surv <- ggplot() +
    geom_ribbon(data = betas_surv_temp, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_surv_temp, aes(x = lag * -1, y = beta, colour = "Temperature")) +
    geom_ribbon(data = betas_surv_prec, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_surv_prec, aes(x = lag * -1, y = beta, colour = "Precipitation"))  +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = seq(from = -1 * (lag), to = 0, by = 12),
                       limits = c((lag * -1) - 1, 0),
                       expand = c(0,0)) +
    scale_colour_manual(values = c("red", "blue"),
                        limits = c("Temperature", "Precipitation")) +
    xlab("months") + ylab("climate anomaly coefficient") + ggtitle("Survival") +
    theme(legend.position = "bottom")


  ggsave(plot = plot_surv,
    filename = here::here("results", "surv_temp_precip_spline.png"), 
         width = 5, height = 3, units = "in", dpi = 300, type = "cairo")
  
  ## -------------------------------------------------------
  ## Growth
  ## -------------------------------------------------------

  model_growth <- mgcv::gam(ln_stems_t1 ~ ln_stems_t0 + 
                              herb_shading_t0 +
                              s(year_t0, bs="re") +
                              s(lags, k=lag, by= Tavg_scaledcovar, bs="cs") +
                              s(lags, k=lag, by= Prec_scaledcovar, bs="cs"),
                            data=data %>%
                              filter(survival_t1 == 1),
                            method="GCV.Cp",gamma=1.2)

  summary(model_growth)

  getBeta.data <- data.frame(lags= c(0:lag), # sort(unique(data$lags)),
                             Tavg_scaledcovar=1, Prec_scaledcovar = 0,
                             ln_stems_t0=1, 
                             year_t0 = 2015, slope = 14,
                             herb_shading_t0 = 5)
  terms.data <- mgcv::predict.gam(model_growth, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_growth_temp <- data.frame(lag=getBeta.data$lags,
                                  beta=terms.data$fit[,"s(lags):Tavg_scaledcovar"],
                                  se=terms.data$se[,"s(lags):Tavg_scaledcovar"])

  getBeta.data <- data.frame(lags=c(0:lag),
                             Tavg_scaledcovar=0, Prec_scaledcovar = 1,
                             ln_stems_t0=1,
                             year_t0 = 2015, slope = 14, 
                             herb_shading_t0 = 5)
  terms.data <- mgcv::predict.gam(model_growth, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_growth_prec <- data.frame(lag=getBeta.data$lags,
                                  beta=terms.data$fit[,"s(lags):Prec_scaledcovar"],
                                  se=terms.data$se[,"s(lags):Prec_scaledcovar"])


  plot_growth <- ggplot() +
    geom_ribbon(data = betas_growth_temp,
                aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_growth_temp,
              aes(x = lag * -1, y = beta, colour = "Temperature")) +
    geom_ribbon(data = betas_growth_prec,
                aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_growth_prec,
              aes(x = lag * -1, y = beta, colour = "Precipitation"))  +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = seq(from = -1 * (lag), to = 0, by = 12),
                       limits = c((lag * -1) - 1, 0),
                       expand = c(0,0)) +
      scale_colour_manual(values = c("red", "blue"),
                          limits = c("Temperature", "Precipitation")) +
      xlab("months") + ylab("climate anomaly coefficient") +  ggtitle("Growth") +
      theme(legend.position = "bottom")


  ggsave(plot = plot_growth,
         filename = here::here("results", "growth_temp_precip_spline.png"), 
         width = 5, height = 3, units = "in", dpi = 300, type = "cairo")

  ## -------------------------------------------------------
  ## Flower probability
  ## -------------------------------------------------------

  model_flowerP <- mgcv::gam(flower_p_t0 ~ ln_stems_t0 + population +
                                shrub_shading_t0 + slope + soil_d3 +
                                s(year_t0, bs="re") +
                                s(lags, k=lag, by= Tavg_scaledcovar, bs="cs") +
                                s(lags, k=lag, by= Prec_scaledcovar, bs="cs"),
                              data=data,
                              family = binomial(link = "logit"),
                              method="GCV.Cp",
                              gamma=1.2)

  summary(model_flowerP)

  getBeta.data <- data.frame(lags= c(0:lag), # sort(unique(data$lags)),
                             Tavg_scaledcovar=1, Prec_scaledcovar = 0,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, shrub_shading_t0 = 5,
                             slope = 14, soil_d3 = 4.5)
  terms.data <- mgcv::predict.gam(model_flowerP, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_flow_temp <- data.frame(lag=getBeta.data$lags,
                                beta=terms.data$fit[,"s(lags):Tavg_scaledcovar"],
                                se=terms.data$se[,"s(lags):Tavg_scaledcovar"])

  getBeta.data <- data.frame(lags=c(0:lag),
                             Tavg_scaledcovar=0, Prec_scaledcovar = 1,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, shrub_shading_t0 = 5,
                             slope = 14, soil_d3 = 4.5)
  terms.data <- mgcv::predict.gam(model_flowerP, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_flow_prec <- data.frame(lag=getBeta.data$lags,
                                beta=terms.data$fit[,"s(lags):Prec_scaledcovar"],
                                se=terms.data$se[,"s(lags):Prec_scaledcovar"])


  plot_pflow <- ggplot() +
    geom_ribbon(data = betas_flow_temp, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_flow_temp, aes(x = lag * -1, y = beta, colour = "Temperature")) +
    geom_ribbon(data = betas_flow_prec, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_flow_prec, aes(x = lag * -1, y = beta, colour = "Precipitation"))  +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = seq(from = -1 * (lag), to = 0, by = 12),
                       limits = c((lag * -1) - 1, 0),
                       expand = c(0,0)) +
      scale_colour_manual(values = c("red", "blue"),
                          limits = c("Temperature", "Precipitation")) +
      xlab("months") + ylab("climate anomaly coefficient") + ggtitle("Flower Probability") +
      theme(legend.position = "bottom")


  ggsave(plot = plot_pflow,
         filename = here::here("results", "flowP_temp_precip_spline.png"), 
         width = 5, height = 3, units = "in", dpi = 300, type = "cairo")
  
  ## -------------------------------------------------------
  ## Abortion probability
  ## -------------------------------------------------------
  model_abortionP <- mgcv::gam(seed_p_t0 ~ ln_stems_t0 + population +
                          shrub_shading_t0 + soil_d3 +
                          s(year_t0, bs="re") +
                          s(lags, k=lag, by= Tavg_scaledcovar, bs="cs") +
                          s(lags, k=lag, by= Prec_scaledcovar, bs="cs"),
                        data=data %>% filter(flower_p_t0 == 1),
                        family = binomial(link = "logit"),
                        method="GCV.Cp",
                        gamma=1.2)
  
  
  summary(model_abortionP)
  
  getBeta.data <- data.frame(lags= c(0:lag), # sort(unique(data$lags)),
                             Tavg_scaledcovar=1, Prec_scaledcovar = 0,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, shrub_shading_t0 = 5,
                             soil_d3 = 4.5)
  terms.data <- mgcv::predict.gam(model_abortionP, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_abort_temp <- data.frame(lag=getBeta.data$lags,
                                beta=terms.data$fit[,"s(lags):Tavg_scaledcovar"],
                                se=terms.data$se[,"s(lags):Tavg_scaledcovar"])
  
  getBeta.data <- data.frame(lags=c(0:lag),
                             Tavg_scaledcovar=0, Prec_scaledcovar = 1,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, shrub_shading_t0 = 5,
                             slope = 14, soil_d3 = 4.5)
  terms.data <- mgcv::predict.gam(model_abortionP, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_abort_prec <- data.frame(lag=getBeta.data$lags,
                                beta=terms.data$fit[,"s(lags):Prec_scaledcovar"],
                                se=terms.data$se[,"s(lags):Prec_scaledcovar"])
  
  
  plot_pabort <- ggplot() +
    geom_ribbon(data = betas_abort_temp, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_abort_temp, aes(x = lag * -1, y = beta, colour = "Temperature")) +
    geom_ribbon(data = betas_abort_prec, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_abort_prec, aes(x = lag * -1, y = beta, colour = "Precipitation"))  +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = seq(from = -1 * (lag), to = 0, by = 12),
                       limits = c((lag * -1) - 1, 0),
                       expand = c(0,0)) +
    scale_colour_manual(values = c("red", "blue"),
                        limits = c("Temperature", "Precipitation")) +
    xlab("months") + ylab("climate anomaly coefficient") + ggtitle("Abortion Probability") +
    theme(legend.position = "bottom")
  
  
  ggsave(plot = plot_pabort,
         filename = here::here("results", "abortP_temp_precip_spline.png"), 
         width = 5, height = 3, units = "in", dpi = 300, type = "cairo")
  
  
  ## -------------------------------------------------------
  ## Seed numbers
  ## -------------------------------------------------------

  model_seeds <- mgcv::gam(est_seed_n_t0 ~ ln_stems_t0 + population +
                             shrub_shading_t0 +
                              s(year_t0, bs="re") +
                              s(lags, k=lag, by= Tavg_scaledcovar, bs="cs") +
                              s(lags, k=lag, by= Prec_scaledcovar, bs="cs"),
                            data=data %>%
                              filter(flower_p_t0 == 1 & seed_p_t0 == 1),
                           family = Gamma(link = "log"),
                            method="GCV.Cp",gamma=1.2)

  summary(model_seeds)

  getBeta.data <- data.frame(lags= c(0:lag), # sort(unique(data$lags)),
                             Tavg_scaledcovar=1, Prec_scaledcovar = 0,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, shrub_shading_t0 = 5)
  terms.data <- mgcv::predict.gam(model_seeds, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_seeds_temp <- data.frame(lag=getBeta.data$lags,
                                  beta=terms.data$fit[,"s(lags):Tavg_scaledcovar"],
                                  se=terms.data$se[,"s(lags):Tavg_scaledcovar"])

  getBeta.data <- data.frame(lags=c(0:lag),
                             Tavg_scaledcovar=0, Prec_scaledcovar = 1,
                             ln_stems_t0=1, population = factor("CR"),
                             year_t0 = 2015, shrub_shading_t0 = 5)
  terms.data <- mgcv::predict.gam(model_seeds, newdata=getBeta.data, type="terms",
                                  se=TRUE)
  betas_seeds_prec <- data.frame(lag=getBeta.data$lags,
                                  beta=terms.data$fit[,"s(lags):Prec_scaledcovar"],
                                  se=terms.data$se[,"s(lags):Prec_scaledcovar"])
  plot_seed <- ggplot() +
    geom_ribbon(data = betas_seeds_temp, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_seeds_temp, aes(x = lag * -1, y = beta, colour = "Temperature")) +
    geom_ribbon(data = betas_seeds_prec, aes(x = lag * -1, ymin = beta - se, ymax = beta + se), alpha = 0.2) +
    geom_line(data = betas_seeds_prec, aes(x = lag * -1, y = beta, colour = "Precipitation"))  +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = seq(from = -1 * (lag), to = 0, by = 12),
                       limits = c((lag * -1) - 1, 0),
                       expand = c(0,0)) +
    scale_colour_manual(values = c("red", "blue"),
                        limits = c("Temperature", "Precipitation")) +
    xlab("months") + ylab("climate anomaly coefficient") +  ggtitle("Number of seeds") +
    theme(legend.position = "bottom")


  ggsave(plot = plot_seed,
         filename = here::here("results", "seed_temp_precip_spline.png"), 
         width = 5, height = 3, units = "in", dpi = 300, type = "cairo")
  
 
  ## Summary plot of all 4 FLM

  sum_plot <- plot_surv + plot_growth + plot_pflow + plot_pabort + plot_seed + 
    plot_layout(nrow = 3, guides = "collect") & theme(legend.position = "bottom")

  
  ggsave(plot = sum_plot,
         filename = here::here("results", "summary_temp_precip_spline.png"), 
         width = 8.5, height = 10, units = "in", dpi = 300, type = "cairo")
  




  return(list(surv = model_survival,
              growth = model_growth,
              flower_p = model_flowerP,
              abort_p = model_abortionP,
              n_seeds = model_seeds))

}
