
bootstrap_glmnet_coefs <- function(x, y, alpha, lambda, n_boot = 1000,
                                   family = "gaussian", seed = 123) {
  
  set.seed(seed)
  n <- nrow(x)
  
  boot_coefs <- vector("list", n_boot)
  
  for (i in seq_len(n_boot)) {
    
    idx <- sample.int(n, n, replace = TRUE)
    x_boot <- x[idx, , drop = FALSE]
    y_boot <- y[idx]
    
    fit_boot <- tryCatch(
      glmnet::glmnet(x_boot, y_boot, alpha = alpha, lambda = lambda, family = family),
      error = function(e) NULL
    )
    
    if (!is.null(fit_boot)) {
      boot_coefs[[i]] <- as.matrix(coef(fit_boot))[, 1]
    }
  }
  
  # drop failed fits, bind into a matrix (rows = replicates)
  boot_coefs <- boot_coefs[!sapply(boot_coefs, is.null)]
  coef_matrix <- do.call(rbind, boot_coefs)
  
  return(coef_matrix)
}

plot_coef_clim_shading <- function(model, boot_coefs, save = F, name, plot_title, ci_level = 0.95) {
  
  c_surv <- coef(model)
  alpha_lo <- (1 - ci_level) / 2
  alpha_hi <- 1 - alpha_lo
  
  demo_data <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
    mutate(herb_shading_t0 = scale(herb_shading_t0),
           shrub_shading_t0 = scale(shrub_shading_t0))
  
  herb_center = attr(demo_data$herb_shading_t0, "scaled:center")
  herb_scale = attr(demo_data$herb_shading_t0, "scaled:scale")
  shrub_center = attr(demo_data$shrub_shading_t0, "scaled:center")
  shrub_scale = attr(demo_data$shrub_shading_t0, "scaled:scale")
  
  df <- bind_rows(expand.grid(variable = c("pr", "tas", "pet"),
                              period = c("spring", "summer", "dormant"),
                              Herb_shading = c(0, 2, 4, 6),
                              Shrub_shading = 2),
                  expand.grid(variable = c("pr", "tas", "pet"),
                              period = c("spring", "summer", "dormant"),
                              Herb_shading = 2,
                              Shrub_shading = c(0, 2, 4, 6))) %>%
    rowwise() %>%
    mutate(
      fix_name   = paste(variable, period, sep = "_"),
      herb_name  = rownames(c_surv)[grep(paste0("herb_shading_t\\d:", fix_name), rownames(c_surv))],
      shrub_name = rownames(c_surv)[grep(paste0("shrub_shading_t\\d:", fix_name), rownames(c_surv))],
      
      fix_eff   = c_surv[which(rownames(c_surv) == fix_name)],
      herb_int  = c_surv[which(rownames(c_surv) == herb_name)],
      shrub_int = c_surv[which(rownames(c_surv) == shrub_name)],
      
      # convert raw shading levels to the SCALED values the model was fit on
      Herb_shading_scaled  = (Herb_shading  - herb_center)  / herb_scale,
      Shrub_shading_scaled = (Shrub_shading - shrub_center) / shrub_scale,
      
      coef_value = fix_eff + (herb_int * Herb_shading_scaled) + (shrub_int * Shrub_shading_scaled),
      
      boot_vals = list(
        boot_coefs[, fix_name] +
          boot_coefs[, herb_name]  * Herb_shading_scaled +
          boot_coefs[, shrub_name] * Shrub_shading_scaled
      ),
      ci_lower = quantile(boot_vals, probs = alpha_lo, na.rm = TRUE),
      ci_upper = quantile(boot_vals, probs = alpha_hi, na.rm = TRUE)
    ) %>% ungroup() %>%
    select(-boot_vals, -Herb_shading_scaled, -Shrub_shading_scaled) %>%
    mutate(Herb_shading = as.factor(Herb_shading),   # raw levels retained for plotting/legend
           Shrub_shading = as.factor(Shrub_shading),
           period = ordered(period, levels = c("summer", "dormant", "spring")))
  
  clim_labels = c('pr' = 'Precipitation',
                  'tas' = 'Temperature',
                  'pet' = 'Potential \nEvapotransp.')
  
  herb_coef <- ggplot(df %>% filter(Shrub_shading == 2), 
                      aes(x = period, y = coef_value, colour = Herb_shading, group = Herb_shading)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Herb_shading, group = Herb_shading),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    scale_colour_manual(name = 'Herb shading', 
                        values = RColorBrewer::brewer.pal(5, "Blues")[2:5]) + 
    scale_fill_manual(name = 'Herb shading', 
                      values = RColorBrewer::brewer.pal(5, "Blues")[2:5]) +
    theme_light() + theme(legend.position = "bottom") + guides(colour = guide_legend(title.position = "top")) +
    ylab("Coefficient estimate") + xlab("Time Period") +
    facet_grid(rows = vars(variable), labeller = as_labeller(clim_labels))
  
  shrub_coef <- ggplot(df %>% filter(Herb_shading == 2), 
                       aes(x = period, y = coef_value, colour = Shrub_shading, group = Shrub_shading)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Shrub_shading, group = Shrub_shading),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    scale_colour_manual(name = 'Shrub shading', 
                        values = RColorBrewer::brewer.pal(5, "Greens")[2:5]) + 
    scale_fill_manual(name = 'Shrub shading', 
                      values = RColorBrewer::brewer.pal(5, "Greens")[2:5]) +
    theme_light() + theme(legend.position = "bottom") + guides(colour = guide_legend(title.position = "top")) +
    ylab("Coefficient estimate") + xlab("Time Period") +
    facet_grid(rows = vars(variable), labeller = as_labeller(clim_labels))
  
  combined_coef_plot <- herb_coef + shrub_coef + plot_layout(guides = "collect") + 
    plot_annotation(title = plot_title) & theme(legend.position = 'bottom')
  
  if(save) {
    ggsave(combined_coef_plot, filename = paste0("results/", name, ".png"),
           width = 6, height = 4)
  }
  
  return(combined_coef_plot)
}



build_effect_size_table <- function(model, boot_coefs,
                                    shading_levels = c(0, 2, 4, 6),
                                    ci_level = 0.95) {
  
  c_surv <- coef(model)
  alpha_lo <- (1 - ci_level) / 2
  alpha_hi <- 1 - alpha_lo
  
  demo_data <- read.csv("data/Dracocephalum_with_vital_rates.csv") %>%
    mutate(herb_shading_t0 = scale(herb_shading_t0),
           shrub_shading_t0 = scale(shrub_shading_t0))
  
  herb_center = attr(demo_data$herb_shading_t0, "scaled:center")
  herb_scale = attr(demo_data$herb_shading_t0, "scaled:scale")
  shrub_center = attr(demo_data$shrub_shading_t0, "scaled:center")
  shrub_scale = attr(demo_data$shrub_shading_t0, "scaled:scale")
  vars <- c("pr", "tas", "pet")
  periods <- c("spring", "summer", "dormant")
  
  combos <- expand.grid(variable = vars, period = periods,
                        shading_type = c("Herb", "Shrub"),
                        shading_level = shading_levels,
                        stringsAsFactors = FALSE)
  
  results <- combos %>%
    rowwise() %>%
    mutate(
      fix_name  = paste(variable, period, sep = "_"),
      int_name  = rownames(c_surv)[grep(
        paste0(tolower(shading_type), "_shading_t\\d:", fix_name), rownames(c_surv))],
      
      # scale the shading level using the appropriate herb/shrub scaling params
      shading_scaled = if (shading_type == "Herb") {
        (shading_level - herb_center) / herb_scale
      } else {
        (shading_level - shrub_center) / shrub_scale
      },
      
      fix_eff = c_surv[which(rownames(c_surv) == fix_name)],
      int_eff = if (length(int_name) == 1) c_surv[which(rownames(c_surv) == int_name)] else 0,
      
      coef_value = fix_eff + int_eff * shading_scaled,
      
      boot_vals = list(
        if (length(int_name) == 1) {
          boot_coefs[, fix_name] + boot_coefs[, int_name] * shading_scaled
        } else {
          boot_coefs[, fix_name]
        }
      ),
      ci_lower = quantile(boot_vals, probs = alpha_lo, na.rm = TRUE),
      ci_upper = quantile(boot_vals, probs = alpha_hi, na.rm = TRUE),
      prop_zero_boot = if (length(int_name) == 1) {
        mean(boot_coefs[, int_name] == 0, na.rm = TRUE)
      } else {
        NA_real_
      }
    ) %>%
    ungroup() %>%
    select(-boot_vals) %>%
    mutate(
      ci_excludes_zero = (ci_lower > 0 & ci_upper > 0) | (ci_lower < 0 & ci_upper < 0)
    )
  
  ## flag sign reversals relative to the shading_level == 0 baseline,
  ## within each variable/period/shading_type group
  results <- results %>%
    group_by(variable, period, shading_type) %>%
    mutate(
      baseline_coef = coef_value[shading_level == min(shading_levels)],
      sign_reversal = sign(coef_value) != sign(baseline_coef) & coef_value != 0 & baseline_coef != 0,
      ci_excludes_baseline = baseline_coef < ci_lower | baseline_coef > ci_upper
    ) %>%
    ungroup() %>%
    select(variable, period, shading_type, shading_level,
           coef_value, ci_lower, ci_upper, ci_excludes_baseline, ci_excludes_zero,
           prop_zero_boot, baseline_coef, sign_reversal) %>%
    arrange(shading_type, variable, period, shading_level)
  return(results)
}
