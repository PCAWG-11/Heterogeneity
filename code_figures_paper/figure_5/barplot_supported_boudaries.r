distr_supported_boundaries <- function(supported_transition_points, n_supported_tp_by_type_random_list) {
   # load(paste0(DIR_RESULTS, "supported_boundaries_info.RData"))

   supported <- data.frame(matrix(0, nrow=2,ncol=2))
   colnames(supported) <- c("Real", "Random")
   rownames(supported) <- c("Clonal-subclonal", "Subclonal-subclonal")

   supported["Clonal-subclonal", "Real"] <- supported_transition_points$supported_clonal_subclonal / supported_transition_points$total_clonal_subclonal
   supported["Subclonal-subclonal", "Real"] <- supported_transition_points$supported_subclonal_subclonal / supported_transition_points$total_subclonal_subclonal
   supported <- cbind("bound_type"=rownames(supported), supported)
   
  random_clonal_subclonal_data <- sapply(n_supported_tp_by_type_random_list, function(x) x$supported_clonal_subclonal)
  random_clonal_subclonal_data <- random_clonal_subclonal_data / n_supported_tp_by_type_random_list[[1]]$total_clonal_subclonal

  random_subclonal_subclonal_data <- sapply(n_supported_tp_by_type_random_list, function(x) x$supported_subclonal_subclonal)
  random_subclonal_subclonal_data <- random_subclonal_subclonal_data / n_supported_tp_by_type_random_list[[1]]$total_subclonal_subclonal

  supported["Clonal-subclonal", "Random"] <- mean(random_clonal_subclonal_data)
  supported["Subclonal-subclonal", "Random"] <- mean(random_subclonal_subclonal_data)
   
  err_clonal_subclonal <- c(quantile(random_clonal_subclonal_data)[["25%"]], quantile(random_clonal_subclonal_data)[["75%"]])
  err_subclonal_subclonal <- c(quantile(random_subclonal_subclonal_data)[["25%"]], quantile(random_subclonal_subclonal_data)[["75%"]])

  supported <- melt(supported)
  supported <- cbind(supported, data.frame(matrix(NA, ncol=2)))
  colnames(supported) <-  c("bound_type", "real_random", "value", "ymax", "ymin")
   
   supported[supported$bound_type == "Clonal-subclonal" & supported$real_random == "Random",]$ymin <- err_clonal_subclonal[1]
   supported[supported$bound_type == "Clonal-subclonal" & supported$real_random == "Random",]$ymax <- err_clonal_subclonal[2]
   supported[supported$bound_type == "Subclonal-subclonal" & supported$real_random == "Random",]$ymin <- err_subclonal_subclonal[1]
   supported[supported$bound_type == "Subclonal-subclonal" & supported$real_random == "Random",]$ymax <- err_subclonal_subclonal[2]

  p = ggplot(supported, aes(factor(bound_type), value, fill = real_random, width=0.8)) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    theme_bw() +
    theme(axis.text = element_text(colour="black",size=20,face="plain"),
          axis.title = element_text(colour="black",size=20,face="plain"),
          strip.text = element_text(colour="black",size=20,face="plain"),
          legend.position = "top",
          legend.text=plot_title_theme,
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title=element_blank(),
          plot.margin = margin(0.2,1.5,0.2,0.2,"cm")) +
      xlab("") + ylab("Fraction of boundaries") #+
       # geom_errorbar(aes(ymax = supported$ymax, ymin = supported$ymin), width=.09, position=position_dodge(0.8))

  return(p)
   # suppressWarnings(ggsave(filename = paste0(DIR_RESULTS,"distr_supported_boundaries.pdf"), width = 4, height=4))
}
