library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

gcl_results = list.files('Research/SenNetPortalProject/GCLResults/', full.names = T)
scallop_results = list.files('Research/SenNetPortalProject/ScallopResults/', full.names = T)


clusters = paste0('SLC16A7+_cell_',2:4)
samples = 1:3


plots = lapply(clusters, function(cluster){
  all_samples = lapply(samples, function(sample){
    filename = paste0(cluster,"_",sample,'.csv')
    gcl_res = read.csv(paste0('Research/SenNetPortalProject/GCLResults/',filename))
    gcl_res = gcl_res$x
    filename = paste0(cluster,"_",sample,'_Scallop.csv')
    scallop_res = read.csv(paste0('Research/SenNetPortalProject/ScallopResults/',filename))
    scallop_res = scallop_res$freqScore
    
    title = paste0('Sample ',sample,', Cluster ', gsub("_"," ",cluster))
    
    data = data.frame(
      value = c(gcl_res, scallop_res),
      group = factor(c(rep("GCL", length(gcl_res)), rep("Scallop", length(scallop_res))))
    )
    
    data$id = title
    return(data)
  })
  
  all_data = bind_rows(all_samples, .id = "sample")
  
  all_data = all_data %>%
    mutate(
      method = ifelse(str_detect(group, "GCL"), "GCL", "Scallop"),
      sample_label = paste0("Sample ", sample),
      sample_group = factor(
        paste(sample_label, method),
        levels = c(
          "Sample 1 GCL", "Sample 1 Scallop",
          "Sample 2 GCL", "Sample 2 Scallop",
          "Sample 3 GCL", "Sample 3 Scallop"
        )
      )
    )
  
  
  title = gsub('_',' ',cluster)
  
  plot = ggplot(all_data, aes(x = sample_group, y = value, fill = method)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    stat_summary(
      fun = mean,
      fun.min = mean,
      fun.max = mean,
      geom = "crossbar",
      width = 0.3,
      color = "black",
      fatten = 0,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = c("GCL" = "#1f77b4", "Scallop" = "#ff7f0e")) +
    labs(
      title = title,
      x = "Sample and Method",
      y = "Value"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  return(plot)
})

# Remove x-axis for first two plots
plots[[1]] <- plots[[1]] + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")  # remove legend

plots[[2]] <- plots[[2]] + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
# keep legend here

plots[[3]] <- plots[[3]] + 
  theme(legend.position = "none")  # remove legend

# Stack vertically
final_plot <- plots[[1]] / plots[[2]] / plots[[3]] + 
  plot_annotation(
    caption = "Shared X-axis label: Sample and Method"
  )

final_plot
ggsave('Research/SenNetPortalProject/Figure.pdf', final_plot)
ggsave('Research/SenNetPortalProject/Figure.png', final_plot)
