### Code for Dentro et al. Supplementary Figure 3

library(ggplot2)


plotdf <- read.delim(file = "DPClust_EOPC-DE_results.tsv", as.is = T)

p1 <- ggplot(data = plotdf, mapping = aes(x = tumor_wgs_aliquot_id, y = no.muts, fill = call)) + geom_col()
p1 <- p1 + scale_fill_brewer(type = "qual") + labs( fill = "") + theme_minimal() + facet_wrap(~donorid, drop = T, nrow = 1, scales = "free_x")
p1 <- p1 + labs(x = "Tumor sample", y = "Proportion of mutations") + theme(axis.text.x = element_blank())
# p1
ggsave(filename = paste0("suppfigure_3.pdf"), plot = p1, width = 10, height = 3)


summary(by(data = plotdf, INDICES = plotdf$tumor_wgs_aliquot_id, FUN = function(x) sum(x[x$call == "illusion of clonality", "no.muts"]) / 
             (sum(x[x$call %in% c("clonal", "illusion of clonality"), "no.muts"]))))

