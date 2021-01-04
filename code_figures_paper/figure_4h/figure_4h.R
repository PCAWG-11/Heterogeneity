### Code for Dentro et al. Figure 4G

library(ggplot2)
library(ggbio)
library(ggrastr)
library(readr)


# read MSAI table and data
eopcplotdf <- read.delim(file = "MSAI_EOPC-DE_samples_regions.tsv")
plotlist <- readRDS(file = "MSAI_EOPC-DE_BAF-LogR_randomised.RDS")

for (idx in 1:nrow(eopcplotdf)) {

  plotdf <- plotlist[[eopcplotdf[idx, "icgc_sample_id"]]]
  
  p.ideo <- Ideogram(genome = "hg19", subchr = paste0("chr", eopcplotdf[idx, "chrom"]), zoom.region = c(eopcplotdf[idx, "start"], eopcplotdf[idx, "end"]))
  # p.ideo
  
  p.baf <- ggplot(data = plotdf, mapping = aes(x = start/1e6)) 
  p.baf <- p.baf + geom_hline(yintercept = median(plotdf$BAF, na.rm = T), color = "#0571b0") + geom_hline(yintercept = 1-median(plotdf$BAF, na.rm = T), color = "#ca0020")
  p.baf <- p.baf + geom_point_rast(mapping = aes(y = BAF), color = "#0571b0", size = 1, stroke = 0, raster.width = 10, raster.height = 3.45) + geom_point_rast(mapping = aes(y = 1-BAF), color = "#ca0020", size = 1, stroke = 0, raster.width = 10, raster.height = 3.45)
  p.baf <- p.baf + theme_minimal() + labs(x = "", y = "BAF") + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(fill = NA))
  p.baf <- p.baf + scale_y_continuous(limits = c(0,1))
  # p.baf
  
  p.logr <- ggplot(data = plotdf, mapping = aes(x = start/1e6)) 
  p.logr <- p.logr + geom_point_rast(mapping = aes(y = LogR), color = "grey40", size = 1, stroke = 0, raster.width = 10, raster.height = 3.45)
  p.logr <- p.logr + theme_minimal() + labs(x = "", y = "LogR") + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(fill = NA))
  p.logr <- p.logr + scale_y_continuous(limits = c(-2,2))
  # p.logr
  
  tks <- tracks(main = "", p.ideo, p.baf, p.logr, heights = c(1, 3, 3))
  # tks
  
  ggsave(filename = paste0("figure_4h_", eopcplotdf[idx, "icgc_donor_id"], "_", eopcplotdf[idx, "icgc_sample_id"], "_chr", eopcplotdf[idx, "chrom"], ".pdf"), plot = tks, width = 9.5, height = 7, units = "in")
  
}
