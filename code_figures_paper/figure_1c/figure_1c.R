#
# Creates figure showing agreement obtained by the consensus copynumber pipeline at different levels of confidence
#

library(ggplot2)
library(readr)
library(reshape2)

px2cm = function(px) { return(px*2.54/96) }

dat = readr::read_tsv("summary_stats_agreement.txt")
dat$star_3 = dat$agreement_level_a+dat$agreement_level_b+dat$agreement_level_c
dat$star_2 = dat$agreement_level_d+dat$agreement_level_e+dat$agreement_level_f
dat$star_1 = dat$agreement_level_g+dat$agreement_level_h+dat$agreement_level_i

colrs = c("#009E73", "#56B4E9", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


plot_data = data.frame(samplename=1:nrow(dat))
for (i in 11:13) {
  dat_s = data.frame(samplename=dat[[1]], agreement=rowSums(dat[, 11:i, drop=F]))
  dat_s = dat_s[with(dat_s, order(agreement, decreasing=T)),]
  dat_s$samplename = 1:nrow(dat)
  plot_data = cbind(plot_data, dat_s[,2])
}

colnames(plot_data) = c("tumour", "(Near) complete agreement, clonal", "Strict majority agreement, clonal", "Complete or strict majority vote,\n rounded subclonal")
plot_data_m = reshape2::melt(plot_data, id.vars=c("tumour"))
plot_data_m$variable = factor(plot_data_m$variable, levels=sort(unique(plot_data_m$variable)))

colnames(plot_data_m)[colnames(plot_data_m)=="variable"] = "Consensus"
p = ggplot(plot_data_m) + geom_point(mapping=aes(x=tumour, y=value, colour=Consensus))
p = p + xlab("Tumour") + ylab("Frac. genome agree") + scale_colour_manual(values=colrs) +
  scale_x_continuous(breaks=seq(0, nrow(dat), 500), limits=c(0, nrow(dat))) +
  theme_bw() + theme(axis.title.x=element_text(colour="black",size=22,face="plain"),
                     axis.text.x=element_text(colour="black",size=18,face="plain"),
                     axis.text.y = element_text(colour="black",size=18,face="plain"),  
                     axis.title.y = element_text(colour="black",size=22,face="plain"),
                     strip.text.x = element_text(colour="black",size=24,face="plain"),
                     plot.title = element_text(colour="black",size=24,face="plain"),
                     legend.text=element_text(colour="black",size=18,face="plain"),
                     legend.title=element_text(colour="black",size=20,face="plain"),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=4)))
         
png("copynumber_agreement.png", width=900, height=370)
print(p)
dev.off()

ggsave(filename="copynumber_agreement.pdf", plot=p, width=px2cm(1200), height=px2cm(370), units="cm", useDingbats=FALSE)

# numbers for the paper:
# mean(rowSums(dat[,11:12]))
# median(rowSums(dat[,11:12]))
# sd(rowSums(dat[,11:12]*100))





