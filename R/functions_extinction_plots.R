
pop_dyn_plot <- function(ext) {
  
  df <- expand.grid(models = c("ACCESS1", "CESM1", "CMCC", "MIROC5", "No change"),
                       population = c("Cr", "Hk", "Ks", "Ru"),
                       stringsAsFactors = F) 
  lims <- data.frame(population = c("Cr", "Hk", "Ks", "Ru"),
                     ylim = c(50,250,35,50))
  pl_df <- left_join(df, lims) %>% split(., sort(as.numeric(rownames(.))))
  
  names(pl_df) <- paste(df$models, df$population)
  
  lapply(pl_df, function(i)
    
    ggplot(ext %>%
             filter(model == i$models &
                      locality == i$population &
                      scenario != "rcp85")) +
      geom_line(aes(x = year,
                    y = pop_size_plant,
                    colour = shading,
                    group = rowid), alpha = 0.5) +
      geom_hline(aes(yintercept = 10)) +
      viridis::scale_color_viridis(option = "D", discrete = T, direction = -1) +  
      ylim(c(0, i$ylim)) +
      xlab("Simulation year") + ylab("Population size")
    
  )
}

ext_yr_plot <- function(df) {
  
  d <- expand.grid(models = c("ACCESS1", "CESM1", "CMCC", "MIROC5", "No change"),
                    population = c("Cr", "Hk", "Ks", "Ru"),
                    stringsAsFactors = F) 
  
  pl_df <- split(d, sort(as.numeric(rownames(d))))
  names(pl_df) <- paste(d$models, d$population)
  
  lapply(pl_df, function(i)
    ggplot(df %>%
           filter(model == i$models &
                    locality == i$population &
                    scenario != "rcp85")) +
    geom_boxplot(aes(y = yr_of_ext,
                     fill = as.factor(shading)),
                 alpha = 0.5) +
    viridis::scale_fill_viridis(option = "D", discrete = T, direction = -1) + 
    ylab("Year of extinction")  +
    coord_flip(ylim = c(2022,2101)) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.background = element_rect(fill='transparent'),
          legend.box.background = element_rect(fill='transparent')
    ) 
  )
  
}
