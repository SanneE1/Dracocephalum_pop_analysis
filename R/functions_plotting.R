## Plotting functions
## 

plot_coef_clim_shading <- function(model, save = F, name) {
  
  c_surv <- coef(model)
  
  df <- bind_rows(expand.grid(variable = c("pr", "tas", "pet"),
                              period = c("spring", "summer", "dormant"),
                              Herb_shading = c(0, 2, 4, 6),
                              Shrub_shading = 0),
                  expand.grid(variable = c("pr", "tas", "pet"),
                              period = c("spring", "summer", "dormant"),
                              Herb_shading = 0,
                              Shrub_shading = c(0, 2, 4, 6))) %>%
    rowwise() %>%
    mutate(fix_eff = c_surv[which(rownames(c_surv) == paste(variable, period, sep = "_"))],
           herb_int = c_surv[which(rownames(c_surv) == paste0("herb_shading_t0:", paste(variable, period, sep = "_")))],
           shrub_int = c_surv[which(rownames(c_surv) == paste0("shrub_shading_t0:", paste(variable, period, sep = "_")))],
           coef_value = fix_eff + (herb_int * Herb_shading) + (shrub_int * Shrub_shading)
    ) %>% ungroup() %>%
    mutate(Herb_shading = as.factor(Herb_shading),
           Shrub_shading = as.factor(Shrub_shading),
           period = ordered(period, levels = c("summer", "dormant", "spring")))
  
  clim_labels = c('pr' = 'Precipitation',
                  'tas' = 'Temperature',
                  'pet' = 'Potential \nEvapotranspiration')
  
  herb_coef <- ggplot(df %>% filter(Shrub_shading == 0), 
                      aes(x = period, y = coef_value, colour = Herb_shading, group = Herb_shading)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    scale_colour_manual(name = 'Herb shading', 
                        values = RColorBrewer::brewer.pal(5, "Blues")[2:5]) + 
    theme_light() + theme(legend.position = "bottom") + guides(colour = guide_legend(title.position = "top")) +
    ylab("Coefficient estimate") + xlab("Time Period") +
    facet_grid(rows = vars(variable), labeller = as_labeller(clim_labels))
  
  shrub_coef <- ggplot(df %>% filter(Herb_shading == 0), 
                       aes(x = period, y = coef_value, colour = Shrub_shading, group = Shrub_shading)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    scale_colour_manual(name = 'Shrub shading', 
                        values = RColorBrewer::brewer.pal(5, "Greens")[2:5]) + 
    theme_light() + theme(legend.position = "bottom") + guides(colour = guide_legend(title.position = "top")) +
    ylab("Coefficient estimate") + xlab("Time Period") +
    facet_grid(rows = vars(variable), labeller = as_labeller(clim_labels))
  
  combined_coef_plot <- herb_coef + shrub_coef + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  if(save) {
    
    ggsave(combined_coef_plot, filename = paste0("results/", name, ".png"),
           width = 6, height = 4)
    
  }
  
  return(combined_coef_plot)
  
}