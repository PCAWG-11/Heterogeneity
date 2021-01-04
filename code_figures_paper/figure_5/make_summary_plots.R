
make_summary_plots = function(summary_table, annotations, summary_table_file, plot_f_bars_per_cancer_type=F, plot_f_venn_diagrams=T, step = 0.1, start = 0.1) {
  ########### plot 1
  load("additional_stats.RData")
  summary_table$num_snvs = summary_table$num_clonal + summary_table$num_subclonal
  
  stats <- read.csv("Stats.csv", header=T, stringsAsFactors=F)
  # annotate the correct cancer type
  stats$type_pcawg = annotations$histology_detailed[match(stats$tumor_id, annotations$samplename)]
  # take only representative samples
  stats = stats[stats$tumor_id %in% annotations$samplename[annotations$multi_rep],]
  stats_big_tumors = stats[stats$n_tp_outside_transitions > 5,]
  
  not_supported = stats_big_tumors$n_not_supporting_cp_clonal + stats_big_tumors$n_not_supporting_cp_subclonal
  not_supported_compressed = table(not_supported)
  not_supported_compressed <- c(not_supported_compressed[1], ">= 1" = sum(not_supported_compressed[2:length(not_supported_compressed)]))
  not_supported_compressed = data.frame(not_supported_compressed)
  not_supported_compressed$category = row.names(not_supported_compressed)
  not_supported_compressed$category = factor(not_supported_compressed$category, levels=c("0", ">= 1"))
  
  p1 = ggplot(not_supported_compressed) + aes(x=category, y=not_supported_compressed) + 
    geom_bar(stat="identity", fill="grey", colour="black", width=0.6) + 
    xlab("Additional subclones") + ylab("Number of tumours") +
    theme_bw() +
    theme(axis.text = element_text(colour="black",size=20,face="plain"),
          axis.title = element_text(colour="black",size=20,face="plain"),
          strip.text = element_text(colour="black",size=20,face="plain"),
          legend.position = "right",
          legend.text=plot_title_theme,
          legend.title=plot_title_theme,
          plot.title = plot_title_theme,
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(2.47,0.2,0.2,0.2,"cm")) # t, r, b, l
          # plot.margin = margin(0.05,0.5,0.5,0.5,"cm"))
  
  ########### plot 2
  
  stats <- read.csv("Stats.csv", header=T)
  # take only representative samples
  stats = stats[stats$tumor_id %in% annotations$samplename[annotations$multi_rep],]
  load("additional_stats.RData")
  
  summary_table = read.table(summary_table_file, header=T, stringsAsFactors=F)
  subclone_table = summary_table[, c("samplename", "num_subclones")]
  colnames(subclone_table) <- c("tumor_id", "n_clusters")
  
  subclone_table <- subclone_table[subclone_table$tumor_id %in% stats$tumor_id, ]
  
  tab <- c()
  for (cmin in seq(start,0.95,step)) {
    cmax = cmin + step
    d <- sapply(stats[stats$max_change_no_s1_s5 < cmax & stats$max_change_no_s1_s5 > cmin & 
                        stats$n_timesteps > 0,]$tumor_id, toString)
    d <- subclone_table[subclone_table$tumor_id %in% d,]$n_clusters 
    
    max_clusters = 2
    new_item <- data.frame(matrix(0, nrow=1, ncol=(max_clusters + 1)))
    colnames(new_item) <- seq(0,max_clusters)
    
    for (n in seq(0,max_clusters)) {
      if (n < 2) {
        new_item[,toString(n)] <- sum(d == n)
      } else {
        new_item[,toString(n)] <- sum(d >= n)
      }
    }
    
    tab <- rbind(tab, cbind(new_item, total=length(d)))
  }
  rownames(tab) <- paste0((seq(start,0.95,step))*100, "-", (seq(start,0.95,step) + step)*100)
  tab <- tab[tab$total > 0,]
  
  # get the counts
  if (length(which(tab$total > 10)) == 0 ) {
    breakpoint = 1
  } else {
    breakpoint = max(which(tab$total > 10))
  }
  if (breakpoint == 1) {
    col_names <- colnames(tab)
    tab <- data.frame(t(data.frame(apply(tab[breakpoint:nrow(tab),],2,sum))))
    colnames(tab) <- col_names
  } else {
    tab <- rbind(tab[1:(breakpoint-1),], apply(tab[breakpoint:nrow(tab),],2,sum))
  }
  tab_normalized <- t(apply(tab[,-ncol(tab)],1,function(x) x / sum(x)))
  rownames(tab_normalized)[nrow(tab_normalized)] <- paste0(">", (seq(start,0.95,step)[breakpoint])*100)
  
  tab_normalized = as.data.frame(tab_normalized)
  tab_normalized$category = factor(row.names(tab_normalized), levels=row.names(tab_normalized))
  tab_normalized_m = reshape2::melt(tab_normalized, id.vars="category")
  tab_normalized_m$variable = as.character(tab_normalized_m$variable)
  tab_normalized_m$variable[tab_normalized_m$variable==2] = ">=2"
  tab_normalized_m$variable = factor(tab_normalized_m$variable, levels=c("0", "1", ">=2"))
  
  # make the plot
  color_palette = c("#de8f07", "#49a4e4", "#138f60")
  p2 = ggplot(tab_normalized_m) + aes(x=category, y=value, fill=variable) + geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values=color_palette) +
    ylim(0,0.6) + 
    xlab("Maximum signature change") + ylab("Fraction of tumours") +
    theme_bw() +
    theme(axis.text = element_text(colour="black",size=20,face="plain"),
          axis.title = element_text(colour="black",size=20,face="plain"),
          strip.text = element_text(colour="black",size=20,face="plain"),
          legend.position = "top",
          legend.text=plot_title_theme,
          legend.title=plot_title_theme,
          plot.title = plot_title_theme,
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(0.2,1.5,0.2,1.5,"cm")) + # t, r, b, l
  # plot.margin = margin(0.05,0.5,0.5,0.5,"cm"))
  guides(fill=guide_legend(title="Subclones detected", 
                           override.aes=list(size=6), 
                           ncol=3, 
                           keywidth=0.3,
                           keyheight=0.3,
                           default.unit="inch",
                           label.hjust=0.5,
                           label.vjust=0.5,
                           title.hjust=0.5))

  ########### plot 3
  load("supported_boundaries_info.RData")
  p3 = distr_supported_boundaries(supported_transition_points, n_supported_tp_by_type_random_list)
  
  ########### plot 4
  # alternative version of panel f showing cancer types
  load("additional_stats.RData")
  summary_table$num_snvs = summary_table$num_clonal + summary_table$num_subclonal
  
  stats <- read.csv("Stats.csv", header=T, stringsAsFactors=F)
  # annotate the correct cancer type
  stats$type_pcawg = annotations$histology_detailed[match(stats$tumor_id, annotations$samplename)]
  # take only representative samples
  stats = stats[stats$tumor_id %in% annotations$samplename[annotations$multi_rep],]
  
  hist_count = table(stats$type_pcawg)
  hist_to_combine = names(hist_count)[hist_count < 100]
  stats$type_pcawg[stats$type_pcawg %in% hist_to_combine] = "Other"
  
  stats_big_tumors = stats[stats$n_tp_outside_transitions > 5,]
  
  res = data.frame()
  for (histology in unique(stats_big_tumors$type_pcawg)) {
    stats_hist = stats_big_tumors[stats_big_tumors$type_pcawg==histology,]
    not_supported_hist = stats_hist$n_not_supporting_cp_clonal + stats_hist$n_not_supporting_cp_subclonal
    res = rbind(res, data.frame(histology=histology, 
                                equals=sum(not_supported_hist==0)/nrow(stats_hist), 
                                greater=sum(not_supported_hist>0)/nrow(stats_hist)))
  }
  
  res$histology = factor(res$histology, levels=res$histology[order(res$greater, decreasing=F)])
  
  res_m = reshape2::melt(res, id.vars="histology")
  res_m = res_m[res_m$variable=="greater",]
  res_m = res_m[res_m$histology!="Other",]
  colours = c("white", "lightblue")
  names(colours) = c("equal", "greater")
  p4 = ggplot(res_m) + aes(x=histology, y=value, fill=variable) + 
    geom_bar(position="stack", stat="identity", width=0.70, colour="black") + 
    scale_fill_manual(values=colours) + 
    xlab("Cancer type") + ylab("Proportion of samples") + coord_flip(ylim=c(0, 0.6)) + # + scale_y_reverse()
    theme_bw() + 
    theme(axis.text.x=element_text(colour="black",size=16,face="plain", angle=90, vjust=0.5, hjust=1), 
          axis.title.x=element_text(colour="black",size=18,face="plain"), 
          axis.title.y=element_blank(), 
          legend.position="none",
          axis.text.y = element_text(colour="black",size=16,face="plain"),
          strip.text.x = element_text(colour="black",size=18,face="plain"),
          strip.text.y = element_text(colour="black",size=18,face="plain"),
          legend.text=plot_title_theme,
          legend.title=plot_title_theme,
          plot.title = plot_title_theme,
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(0.05,0.5,0.5,0.5,"cm"))
  
  
  if (plot_f_bars_per_cancer_type) {
    
    
    
    temp = data.frame(type_pcawg=stats_big_tumors$type_pcawg, additional_transitions=stats_big_tumors$n_not_supporting_cp_clonal + stats_big_tumors$n_not_supporting_cp_subclonal)
    temp = temp %>% group_by(type_pcawg) %>% summarise(additional_transitions=sum(additional_transitions))
    # get average in samples where at least 1 change detected
    hist_count_with_change = table(stats_big_tumors$type_pcawg[(stats_big_tumors$n_not_supporting_cp_clonal + stats_big_tumors$n_not_supporting_cp_subclonal) > 0])
    temp$mean_per_sample = as.numeric(temp$additional_transitions / hist_count_with_change[match(temp$type_pcawg, names(hist_count_with_change))])
    temp$type_pcawg = factor(temp$type_pcawg, levels=levels(res_m$histology))
    
    p5 = ggplot(temp) + aes(x=type_pcawg, y=mean_per_sample) + 
      geom_bar(stat="identity", width=0.70, fill="orange", colour="black") + 
      scale_y_continuous(breaks=seq(0.0, 1.2, 0.2)) +
      coord_flip() +
      ylab("Mean extra transitions") +
      theme_bw() + 
      theme(axis.text.x=element_text(colour="black",size=16,face="plain", angle=90, vjust=0.5, hjust=1), 
            axis.title.x=element_blank(), 
            axis.title.y=element_blank(), 
            legend.position="left",
            axis.text.y = element_blank(),
            strip.text.x = element_text(colour="black",size=16,face="plain"),
            strip.text.y = element_text(colour="black",size=16,face="plain"),
            legend.text=plot_title_theme,
            legend.title=plot_title_theme,
            plot.title = plot_title_theme,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(0.05,0.5,0.5,0.5,"cm"))
    
    p4 = arrangeGrob(p4, p5, ncol=2, widths=c(11.2/18, 6.8/18))
  } else if (plot_f_venn_diagrams) {
    
    stats <- read.csv("Stats.csv", header=T, stringsAsFactors=F)
    # annotate the correct cancer type
    stats$type_pcawg = annotations$histology_detailed[match(stats$tumor_id, annotations$samplename)]
    # take only representative samples
    stats = stats[stats$tumor_id %in% annotations$samplename[annotations$multi_rep],]
    hist_count = table(stats$type_pcawg)
    
    temp = data.frame(type_pcawg=stats$type_pcawg, additional_transitions=stats$n_not_supporting_cp_clonal + stats$n_not_supporting_cp_subclonal)
    temp = temp %>% group_by(type_pcawg) %>% summarise(additional_transitions=sum(additional_transitions))
    # get average in samples where at least 1 change detected
    hist_count_with_change = table(stats$type_pcawg[(stats$n_not_supporting_cp_clonal + stats$n_not_supporting_cp_subclonal) > 0])
    temp$mean_per_sample = as.numeric(temp$additional_transitions / hist_count[match(temp$type_pcawg, names(hist_count))])
    temp = temp[temp$type_pcawg %in% names(hist_count[hist_count > 24]),]
    temp = temp[order(temp$mean_per_sample, decreasing=T),]
    temp = temp[1:10,]
    temp$type_pcawg = factor(temp$type_pcawg, levels=temp$type_pcawg)
    
    # temp = data.frame(type_pcawg=stats_big_tumors$type_pcawg, additional_transitions=stats_big_tumors$n_not_supporting_cp_clonal + stats_big_tumors$n_not_supporting_cp_subclonal)
    # temp = temp %>% group_by(type_pcawg) %>% summarise(additional_transitions=sum(additional_transitions))
    # # get average in samples where at least 1 change detected
    # hist_count_with_change = table(stats_big_tumors$type_pcawg[(stats_big_tumors$n_not_supporting_cp_clonal + stats_big_tumors$n_not_supporting_cp_subclonal) > 0])
    # temp$mean_per_sample = as.numeric(temp$additional_transitions / hist_count_with_change[match(temp$type_pcawg, names(hist_count_with_change))])
    # # temp$type_pcawg = factor(temp$type_pcawg, levels=levels(res_m$histology))
    
    cancer_types_ordering = as.character(temp$type_pcawg[order(temp$mean_per_sample, decreasing=F)])
    cancer_types_ordering = c("Other", cancer_types_ordering[cancer_types_ordering!="Other"])
    temp$type_pcawg = factor(temp$type_pcawg, levels=cancer_types_ordering)
    
    p5 = ggplot(temp) + aes(x=type_pcawg, y=mean_per_sample) + 
      geom_bar(stat="identity", width=0.70, fill="orange", colour="black") + 
      scale_y_continuous(breaks=seq(0.0, 1, 0.25), limits=c(0,1)) +
      coord_flip() +
      ylab("") +
      theme_bw() + 
      theme(axis.text.x=element_text(colour="black",size=16,face="plain", angle=90, vjust=0.5, hjust=1), 
            axis.title.x=element_text(colour="black",size=18,face="plain"), 
            axis.title.y=element_blank(), 
            legend.position="none",
            axis.text.y = element_text(colour="black",size=16,face="plain"),
            strip.text.x = element_text(colour="black",size=18,face="plain"),
            strip.text.y = element_text(colour="black",size=18,face="plain"),
            legend.text=plot_title_theme,
            legend.title=plot_title_theme,
            plot.title = plot_title_theme,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(0.2,0.5,0.5,1.5,"cm")) +  #t r b l
      ggtitle("Mean extra transitions")
    
    p4 = get_venn_diagrams()
    p4 = arrangeGrob(grobs=list(p4, p5), ncol=2)
  }
    
  # return(arrangeGrob(grobs=list(p3, p2, p1), ncol=3, widths=c(3/8, 2/4, 1/8)))
  # dropped the big bar graph
  return(arrangeGrob(grobs=list(p3, p2, p4), ncol=3, widths=c(1/4,1/4,2/4)))
}

# plot venn diagrams with number of transitions and cluster boundaries
get_venn_diagrams = function() {
  
  # alternative version of panel f showing venn diagrams
  load("additional_stats.RData")
  summary_table$num_snvs = summary_table$num_clonal + summary_table$num_subclonal
  
  stats <- read.csv("Stats.csv", header=T, stringsAsFactors=F)
  stats_reduced <- stats[stats$n_clusters > 1 & stats$n_changepoints > 0,]
  
  library(VennDiagram)
  
  unsupported_changepoints = sum(stats_reduced$n_not_supporting_cp)
  supported_boundaries = (supported_transition_points$supported_clonal_subclonal + supported_transition_points$supported_subclonal_subclonal)/50
  total_boundaries <- (supported_transition_points$total_clonal_subclonal + supported_transition_points$total_subclonal_subclonal)/50
  unsupported_boundaries = total_boundaries-supported_boundaries
  
  grid.newpage()
  p1 = draw.pairwise.venn(unsupported_changepoints+supported_boundaries, 
                          unsupported_boundaries+supported_boundaries, 
                          supported_boundaries, 
                          category = c("Signature transitions", "Cluster boundaries"), 
                          # lty = rep("blank", 2), 
                          fill = c("light blue", "red"),
                          col = rep("black", 2), 
                          lwd = c(4,4),
                          alpha = rep(0.5, 2),
                          cat.fontfamily = rep("Helvetica", 2),
                          cat.cex=rep(1.3,2),
                          cat.pos = c(0, 0), 
                          cat.dist = rep(0.025, 2),
                          fontfamily = rep("Helvetica", 3),
                          cex=1.3,
                          main="test2",
                          sub="test2",
                          main.cex=1.5,
                          main.col="black",
                          main.pos=c(0,0),
                          margin = 0.05)
  r = arrangeGrob(gTree(children=p1), ncol=1)
  return(r)
}



