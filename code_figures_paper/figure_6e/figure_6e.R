#
# Creates figure showing an inventory of drivers that are clinically actionable per cancer type
#

px2cm = function(px) { return(px/100) }

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

createPng = function(p, filename, height, width) {
  png(filename, height=height, width=width)
  print(p)
  dev.off()
}

# Minimum probability of the mutation being subclonal
min_prob_subclonal = 0.5
is_drivers = T
is_coding = T
is_noncoding = T

dat = read.table("TableS3_panorama_driver_mutations_pcawg_clinically_actionable_events.tsv.gz", header=T, stringsAsFactors=F, sep="\t")
# Take only the matches that actually work for that cancer type
dat = dat[dat$match_number != 0,]
# take only SNVs and indels
dat = dat[grepl("MUT", dat$alteration),]
dat$gene = unlist(lapply(1:nrow(dat), function(i) unlist(stringr::str_split(dat$alteration[i], " "))[1]))
dat = dat[dat$effect=="Responsive",]

if (is_drivers) {
  drivers = read.table("TableS3_panorama_driver_mutations_pcawg_annotated_v1.1.tsv.gz", header=T, stringsAsFactors=F, sep="\t")
  drivers$samplename = drivers$sample_id
  # Fetch index of clonal cluster
  clonalid = unlist(lapply(drivers$all_cluster_ccf, function(x) {
    if (grepl(",", x)) {
      ccf_split = as.numeric(unlist(stringr::str_split(x, ",")))
      which(ccf_split < 1.01 & ccf_split > 0.99)
    } else {
      1
    }
  }))
  prob_clonal = unlist(lapply(1:nrow(drivers), function(i) {
    x = drivers$all_assign_probs[i]
    if (grepl(",", x)) {
      as.numeric(unlist(stringr::str_split(x, ",")))[clonalid[i]]
    } else {
      1
    }
  }))
  drivers$prob_clonal = prob_clonal
  drivers$prob_subclonal = 1-prob_clonal
  plot_postfix = "drivers"
  
  # filter for coding status
  if (is_coding & is_noncoding) {
    # keep all
  } else if (is_coding & !is_noncoding) {
    # keep coding drivers only
    drivers = drivers[drivers$category=="coding",]
  } else if (!is_coding & is_noncoding) {
    # keep noncoding only
    drivers = drivers[drivers$category=="noncoding",]
  }
  
} else {
  drivers = read.table("all_annotated_snvIndel_proteinChanging.txt", header=T, stringsAsFactors=F, sep="\t")
  plot_postfix = "proteinalt"
}


############################################################################################
# filter
############################################################################################
summary_table = read.table("summary_table_combined_annotations_v3.txt", header=T, stringsAsFactors=F, sep="\t")

anno = read.table("pcawg.wg11.final_sample_list.txt", header=T, stringsAsFactors=F)
# take representative samples that are powered
anno = anno[anno$multi_rep & anno$power,]

# Take only preferred samples
drivers = drivers[drivers$samplename %in% anno$samplename,]
summary_table = summary_table[summary_table$samplename %in% anno$samplename, ]

drivers$histology_abbreviation = summary_table$histology_abbreviation[match(drivers$sample_id, summary_table$samplename)]

############################################################################################
# annotate
############################################################################################
# Annotate to drivers whether they are actionable
drivers$actionable = F
drivers$effect = NA
drivers$class = NA
drivers$drug = NA
for (i in 1:nrow(drivers)) {
  dat_s = dat[dat$sample==drivers$sample_id[i] & dat$gene==drivers$gene[i],]
  
  m = match(drivers$gene[i], dat_s$gene)
  if (!is.na(m)) {
    drivers$actionable[i] = TRUE
    drivers$effect[i] = paste0(dat_s$effect[m], collapse=";")
    drivers$class[i] = paste0(dat_s$class[m], collapse=";")
    drivers$drug[i] = paste0(dat_s$drug[m], collapse=";")
  }
}
# take only what is actionable
drivers = drivers[drivers$actionable,]

# save a table with the actual numbers
res = drivers %>% group_by(histology_abbreviation, gene) %>% summarise(n=n())
res = dcast(res, formula = histology_abbreviation ~ gene)
res$total = rowSums(res[,2:ncol(res)], na.rm=T)
res[nrow(res)+1,] = c("Total", colSums(res[,2:ncol(res)], na.rm=T))
write.table(res, file="actionable_driver_counts.txt", quote=F, sep="\t", row.names=F)

# deal with samples with duplicates
combined = data.frame()
for (samplename in unique(drivers$sample_id[duplicated(drivers$sample_id)])) {
  # dat_s = dat[dat$sample==samplename,]
  drivers_s = drivers[drivers$sample_id==samplename,]
  
  prob_clonal = prod(drivers_s$prob_clonal)
  prob_subclonal = prod(drivers_s$prob_subclonal)
  prob_both = 1 - (prob_clonal + prob_subclonal)
  dat_s = drivers_s[1,, drop=F]
  dat_s$prob_clonal = prob_clonal
  dat_s$prob_subclonal = prob_subclonal
  dat_s$prob_both = prob_both
  dat_s$gene = "combined"
  if (is_drivers) {
    dat_s[, c("timing", "ccf", "mcn", "mult", "major_cn", "minor_cn", "cluster_ccf")] = NA
  } else {
    dat_s[, c("ccf", "mcn", "mult")] = NA
  }
  combined = rbind(combined, dat_s)
}
drivers = drivers[!drivers$samplename %in% drivers$samplename[duplicated(drivers$sample_id)], ]
drivers$prob_both = 0
drivers = rbind(drivers, combined)

############################################################################################
# summarise the data
############################################################################################

get_counts = function(dat, dat_type, col_name) {
  counts_subclonal = ddply(dat, .(dat$histology_abbreviation), function(x) sum(x[, col_name], na.rm=T))
  colnames(counts_subclonal) = c("histology_abbreviation", "count")
  counts_subclonal$type = dat_type
  return(counts_subclonal)
}

counts_subclonal = get_counts(drivers[, c("sample_id", "histology_abbreviation", "prob_subclonal")], "subclonal", "prob_subclonal")
counts_clonal = get_counts(drivers[, c("samplename", "histology_abbreviation", "prob_clonal")], "clonal", "prob_clonal")
counts_both = get_counts(drivers[, c("samplename", "histology_abbreviation", "prob_both")], "both", "prob_both")

counts_none = data.frame()
for (histology in unique(drivers$histology_abbreviation)) {
  all_hist = summary_table$samplename[summary_table$histology_abbreviation==histology]
  counts_nodriver = sum(!all_hist %in% drivers$sample_id)
  counts_none = rbind(counts_none, data.frame(histology_abbreviation=histology,
                                              count=counts_nodriver,
                                              type="none"))
}

counts_dat = rbind(counts_clonal, counts_subclonal, counts_both, counts_none)
counts_dat$type = factor(counts_dat$type, levels=c("clonal", "both", "subclonal", "none"))
counts_dat$histology_abbreviation = factor(counts_dat$histology_abbreviation, levels=sort(unique(drivers$histology_abbreviation)))
counts_dat$frac = NA
counts_dat$frac_with_target = NA

sel = rep(T, nrow(dat))

# transform counts into fractions per cancer type
hist_sam_counts = table(unique(drivers[sel, c("histology_abbreviation", "samplename")])$histology_abbreviation)
for (hist_type in names(hist_sam_counts)) {
  sel = counts_dat$histology_abbreviation==hist_type
  sel[is.na(sel)] = F
  counts_dat$frac[sel] = counts_dat$count[sel] / hist_sam_counts[[hist_type]]
  counts_dat$frac_with_target[sel] = counts_dat$count[sel] / sum(sel & counts_dat$type!="none")
  counts_dat$frac_with_target[sel & counts_dat$type=="none"] = NA
}

# Obtain the number of events per cancer type, per category
counts_final = data.frame()
for (item in unique(counts_dat$histology_abbreviation)) {
  print(item)
  num_clonal = sum(counts_dat$count[counts_dat$histology_abbreviation==item & counts_dat$type=="clonal"], na.rm=T)
  num_subclonal = sum(counts_dat$count[counts_dat$histology_abbreviation==item & counts_dat$type=="subclonal"], na.rm=T)
  num_both = sum(counts_dat$count[counts_dat$histology_abbreviation==item & counts_dat$type=="both"], na.rm=T)
  num_none = sum(counts_dat$count[counts_dat$histology_abbreviation==item & counts_dat$type=="none"], na.rm=T)
  num_all = num_clonal+num_subclonal+num_both+num_none
  counts_final = rbind(counts_final, data.frame(histology_abbreviation=item, type="clonal", count=num_clonal, frac=num_clonal/num_all, frac_with_target=num_clonal/(num_all-num_none)))
  counts_final = rbind(counts_final, data.frame(histology_abbreviation=item, type="subclonal", count=num_subclonal, frac=num_subclonal/num_all, frac_with_target=num_subclonal/(num_all-num_none)))
  counts_final = rbind(counts_final, data.frame(histology_abbreviation=item, type="both", count=num_both, frac=num_both/num_all, frac_with_target=num_both/(num_all-num_none)))
  counts_final = rbind(counts_final, data.frame(histology_abbreviation=item, type="none", count=num_none, frac=num_none/num_all, frac_with_target=NA))
}

all_clonal = sum(counts_final[counts_final$type=="clonal",]$count)
all_subclonal = sum(counts_final[counts_final$type=="subclonal",]$count)
all_both = sum(counts_final[counts_final$type=="both",]$count)
all_none = sum(counts_final[counts_final$type=="none",]$count)
all_all = sum(all_clonal, all_subclonal, all_both, all_none)

# proportions of total samples
all_subclonal / all_all # 2017 beta release: 0.05146723 2018 v1.1. release: 0.02202909
(all_subclonal+all_both) / all_all # 2017 beta release: 0.1172472 2018 v1.1. release: 0.04654357

# proportions of samples with at least one actionable hit
all_subclonal/(all_all-all_none) # 2017 beta release: 0.09120837 2018 v1.1. release: 0.04504381
(all_subclonal+all_both)/(all_all-all_none) # 2017 beta release: 0.2077813 2018 v1.1. release: 0.09516962
all_clonal/(all_all-all_none) # 2017 beta release: 0.7922187 2018 v1.1. release: 0.9048304

counts_final = rbind(counts_final, data.frame(histology_abbreviation="All", type="clonal", count=all_clonal, frac=all_clonal/all_all, frac_with_target=all_clonal/(all_all-all_none)))
counts_final = rbind(counts_final, data.frame(histology_abbreviation="All", type="subclonal", count=all_subclonal, frac=all_subclonal/all_all, frac_with_target=all_subclonal/(all_all-all_none)))
counts_final = rbind(counts_final, data.frame(histology_abbreviation="All", type="both", count=all_both, frac=all_both/all_all, frac_with_target=all_both/(all_all-all_none)))

histology_order = data.frame()
for (histology in unique(counts_final$histology_abbreviation)) {
  if (histology != "All" & !is.na(histology)) {
    temp = subset(counts_final, counts_final$histology_abbreviation==histology)
    frac = (temp$frac[temp$type=="subclonal"] + temp$frac[temp$type=="both"]) #/ (temp$count[temp$type=="clonal"] + temp$count[temp$type=="both"])
    histology_order = rbind(histology_order, data.frame(histology_abbreviation=histology, frac=frac))
  }
}

counts_final = counts_final[!is.na(counts_final$histology_abbreviation), ]
# remove the none category, as its taken into account by the above calculated fractions already
counts_final = counts_final[counts_final$type!="none",]
counts_final$histology_abbreviation = factor(counts_final$histology_abbreviation, levels=c(as.character(histology_order$histology_abbreviation[order(histology_order$frac)]), "", "All"))
counts_final$type = as.character(counts_final$type)
counts_final$type[counts_final$type=="clonal"] = "Clonal"
counts_final$type[counts_final$type=="subclonal"] = "Subclonal"
counts_final$type[counts_final$type=="both"] = "Both"
counts_final$type = factor(counts_final$type, levels=c("Clonal", "Subclonal", "Both"))

# get counts to plot at the top
hist_sam_counts = table(summary_table$histology_abbreviation)
hist_sam_counts["All"] = sum(hist_sam_counts)


counts_final_plot = counts_final[counts_final$histology_abbreviation %in% c(names(hist_sam_counts)[hist_sam_counts>10], "All"),]
counts_final_plot$histology_abbreviation = factor(counts_final_plot$histology_abbreviation, 
                                                  levels=levels(counts_final_plot$histology_abbreviation)[levels(counts_final_plot$histology_abbreviation) %in% c(as.character(unique(counts_final_plot$histology_abbreviation)), "")])

hist_sam_counts[""] = ""
labels = hist_sam_counts[levels(counts_final_plot$histology_abbreviation)]

# mycolours = c("#e5e5ff", "grey", "#000080") # driver figure colours
mycolours = c("#ef8a62ff", "#67a9cfff", "grey54") # reworked combined driver figure colours
p = ggplot(counts_final_plot) + 
  aes(x=histology_abbreviation, y=frac, fill=type) + 
  geom_bar(position="stack", stat="identity", width=0.75) + 
  scale_fill_manual(values=mycolours, guide=guide_legend(title = "Targetable driver type")) + 
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(breaks=seq(0, 1, 0.1), labels=paste(seq(0, 100, 10), "%", sep=""), limits=c(0,1.04)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black",size=14,face="plain", angle=56, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour="black",size=14,face="plain"), 
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=16,face="plain"),
        legend.position = "bottom",
        legend.text = element_text(colour="black",size=13,face="plain"),
        legend.title = element_text(colour="black",size=14,face="bold"),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab("Cancer Type") + ylab("Percentage of tumours") + #ylim(0,1.04) + 
  annotate("text", x=1:length(labels), y=1.04, label=labels, size=3.5)

createPng(p, paste0("count_per_histology_probsum_altcolour_", plot_postfix, ".png"), 450, 900)
ggsave(paste0("count_per_histology_probsum_altcolour_", plot_postfix, ".pdf"), p, height=14, width=29, units="cm", useDingbats=FALSE)

# figures for the paper
frac_tumours_with_actional_events = sum(counts_final_plot$frac[counts_final_plot$histology_abbreviation=="All"])
tumours_with_onlysubclonal_actionable_events = sum(counts_final_plot$frac[counts_final_plot$histology_abbreviation=="All" & counts_final_plot$type=="Subclonal"])
tumours_with_subclonal_actionable_events = tumours_with_onlysubclonal_actionable_events + sum(counts_final_plot$frac[counts_final_plot$histology_abbreviation=="All" & counts_final_plot$type=="Both"])

frac_tumours_with_onlysubclonal_actionable_events = tumours_with_onlysubclonal_actionable_events / tumours_with_actional_events
frac_tumours_with_subclonal_actionable_events = tumours_with_subclonal_actionable_events / tumours_with_actional_events

print(paste0("Fraction tumours with actionable events: ", frac_tumours_with_actional_events))
print(paste0("Fraction tumours with subclonal actionable events: ", frac_tumours_with_subclonal_actionable_events))
print(paste0("Fraction tumours with only subclonal actionable events: ", frac_tumours_with_onlysubclonal_actionable_events))

