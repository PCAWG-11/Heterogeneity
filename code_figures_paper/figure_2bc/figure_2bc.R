#
# Creates plot showing the amount of adjustment to cluster positions and sizes after correcting for the Winner's curse
#

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(scales)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# combine the old and new table, as old contains the original cluster positions and sizes
# the numbers in ccf_correction_batch1_2_2778_samples.tsv (original_cluster_locations) are correct
original_cluster_locations = readr::read_tsv("pcawg_weme2018_final.txt")
dat = readr::read_tsv("pcawg_weme2018_wcc_final.txt")
purity = readr::read_tsv("consensus.20170217.purity.ploidy.txt.gz")
original_cluster_locations$ccf = original_cluster_locations$proportion / purity$purity[match(original_cluster_locations$samplename, purity$samplename)]
dat$ccf = dat$proportion / purity$purity[match(dat$samplename, purity$samplename)]
# dat$consensus_cluster_ccf = dat$consensus_cluster_cp / dat$purity

for (samplename in unique(original_cluster_locations$samplename)) {
  n_muts = sum(dat$n_ssms[dat$samplename==samplename])
  original_cluster_locations$n_ssms[original_cluster_locations$samplename==samplename] = (original_cluster_locations$n_ssms[original_cluster_locations$samplename==samplename] / 100) * n_muts
}

dat = data.frame(aliquot=original_cluster_locations$samplename,
                 obs_cluster_pos=original_cluster_locations$ccf,
                 obs_n_muts=original_cluster_locations$n_ssms,
                 new_cluster_pos=dat$ccf,
                 new_n_muts=dat$n_ssms)

dat$cluster_type = "Subclonal"
dat$cluster_type[dat$obs_cluster_pos==1] = "Clonal"
dat$cluster_type = factor(dat$cluster_type, levels=c("Clonal", "Subclonal"))
dat$cluster_pos_diff = abs(dat$obs_cluster_pos-dat$new_cluster_pos)
dat$cluster_n_muts_diff = abs(dat$obs_n_muts-dat$new_n_muts)
dat$frac_muts_missing = dat$cluster_n_muts_diff/dat$new_n_muts

dat_pos = dat[, c("aliquot", "cluster_pos_diff", "cluster_type")]
colnames(dat_pos) = c("aliquot", "value", "cluster_type")
# dat_size = dat[, c("aliquot", "cluster_n_muts_diff", "cluster_type")]
dat_size = dat[, c("aliquot", "frac_muts_missing", "cluster_type")]
colnames(dat_size) = c("aliquot", "value", "cluster_type")

dat_pos = dat_pos[with(dat_pos, order(dat_pos$value)),]
dat_size = dat_size[with(dat_size, order(dat_size$value)),]

dat_pos$index = 1:nrow(dat_pos)
dat_size$index = 1:nrow(dat_size)

# # making sure the zeroes don't become -1 when log transforming
# dat_size$value[dat_size$value==0] = 0.1


make_plot = function(dat, y_label, plot_title, y_breaks, y_break_labels) {
  mycolours = c("grey", "#000080")
  return(ggplot(as.data.frame(dat)) + aes_string(x="index", y="value", colour="cluster_type") + 
           geom_point(size=1.2) +
           scale_colour_manual(values=mycolours, guide=guide_legend(title = "Cluster type")) + 
           scale_x_continuous(breaks=seq(0, 5000, 1000)) +
           scale_y_continuous(breaks=y_breaks, labels=y_break_labels) +
          theme_bw() +
          theme(axis.text.x = element_text(colour="black",size=18,face="plain"),
          axis.title.x = element_text(colour="black",size=20,face="plain"),
          axis.text.y = element_text(colour="black",size=18,face="plain"), 
          axis.title.y = element_text(colour="black",size=20,face="plain"),
          strip.text.x = element_text(colour="black",size=20,face="plain"),
          strip.text.y = element_text(colour="black",size=20,face="plain"),
          legend.position = "bottom",
          legend.text = element_text(colour="black",size=16,face="plain"),
          legend.title = element_text(colour="black",size=18,face="bold"),
          plot.title = element_text(hjust=0.5, colour="black",size=20,face="plain"),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black")) +
          guides(colour = guide_legend(title = "Cluster type", override.aes = list(size=4))) +
    ylab(y_label) + xlab("Mutation cluster index") + ggtitle(plot_title))
}

math_format0 = function(expr = 10^.x, format = force) {
  quoted <- substitute(expr)
  subs <- function(x) {
    do.call("substitute", list(quoted, list(.x = as.name(x))))
  }
  function(x) {
    x <- format(x)
    c(list(0), lapply(x, subs))
  }
}

trans_breaks0 = function (trans, inv, n = 5, ...) {
  trans <- match.fun(trans)
  inv <- match.fun(inv)
  function(x) {
    c(0, inv(pretty(trans(x), n, ...)))
  }
}

trans_format0 = function (trans, format = scientific_format()) 
{
  if (is.character(trans)) 
    trans <- match.fun(trans)
  function(x) {
    x <- trans(x)
    format(x)
  }
}



# x <- seq(0,2*pi,0.1)
# y <- sin(x)
# 
# df <- data.frame(x, y)
# ggplot(df, aes(x=x, y=y))+
#   geom_point(size=4)+
#   labs(x=expression(Production~rate~" "~mu~moles~NO[3]^{-1}-N~Kg^{-1}),
#        y=expression(Concentration~mg~L^{-1})) 

p1 = make_plot(dat_pos[seq(1, nrow(dat_pos), 1),], "Change in cancer cell fraction", "Mutation cluster position change", seq(0,0.25,0.05), c("0.00", "0.05", "0.10", "0.15", "0.20", "0.25"))
p2 = make_plot(dat_size[seq(1, nrow(dat_size), 1),], "Fraction of SNVs missing", "Mutation cluster size change", seq(0,1,0.25), c("0.00", "0.25", "0.50", "0.75", "1.00")) 
p2 = p2 + ylim(0,1)
# p2 = p2 + scale_y_continuous("Additional SNVs (log10)", trans='log10', breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000), labels=c(0, expression(10^{0}), expression(10^{1}), expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5})), limits=c(-1, 100000))
legend = g_legend(p1)
p = arrangeGrob(
  # grid.arrange(
  arrangeGrob(p1 + theme(legend.position='none'), p2 + theme(legend.position='none'), nrow=1),
  legend, heights=c(13/14, 1/14))

# "^" edited out in Illustrator
to_inch_conversion = 0.393701
ggsave(filename="winners_curse_correction_fracSNVs.pdf", plot=p, height=15*to_inch_conversion, width=30*to_inch_conversion, units="in", useDingbats=FALSE)
png(filename="winners_curse_correction_fracSNVs.png", height=15*to_inch_conversion, width=30*to_inch_conversion, units="in", res=300)
grid.draw(p)
dev.off()

p1 = make_plot(dat_pos[seq(1, nrow(dat_pos), 1),], "Change in CCF", "Mutation cluster position change", seq(0,0.25,0.05), c("0.00", "0.05", "0.10", "0.15", "0.20", "0.25"))
p2 = make_plot(dat_size[seq(1, nrow(dat_size), 1),], "Fraction SNVs missing", "Mutation cluster size change", seq(0,1,0.25), c("0.00", "0.25", "0.50", "0.75", "1.00"))  
p2 = p2 + ylim(0,1)

p = arrangeGrob(
  # grid.arrange(
  arrangeGrob(p1 + theme(legend.position='none'), p2 + theme(legend.position='none'), legend,  
              heights=c(6.5/14, 6.5/14, 1/14), ncol=1))

ggsave(filename="winners_curse_correction_fracSNVs_alt.pdf", plot=p, height=20, width=20, units="cm", useDingbats=FALSE)
png(filename="winners_curse_correction_fracSNVs_alt.png", height=20*to_inch_conversion, width=20*to_inch_conversion, units="in", res=300)
grid.draw(p)
dev.off()

