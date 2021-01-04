library(plyr)
library(reshape2)
library(ggplot2)

#######SVs
##side panel
mDrivClon<-read.table("figure_6.sv.input.txt", header = T, sep="\t", stringsAsFactors = F)
sorting.tmp<-mDrivClon[mDrivClon$variable=="subclonal" & mDrivClon$Cancer!="Pancancer",]
organ.order<-sorting.tmp[order(sorting.tmp$value, decreasing = T),"Cancer"]

p4 = ggplot(mDrivClon, aes(y = value, x = Cancer, fill = factor(variable, level=c("subclonal","clonal")))) +
  scale_fill_manual(values=c("#67a9cf","#ef8a62")) + 
  scale_x_discrete(limits=rev(organ.order), name="", position="top")+
  geom_bar(stat = "identity", position = "dodge")+
  scale_y_continuous(name="", limits=c(0,1), position="right")+
  theme(axis.text.x = element_blank(),axis.text.y = element_text(angle = 180, size=12),panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position="")

ggsave(p4,file='figure_6a.svs.side.pdf', width = 10, height = 1,useDingbats=F)

###top panel
propPerSRB<-read.table("proportions_of_drivers_per_SRB.txt", header=T, sep="\t", strip.white = F)

propPerSRB$frac_clonal<-(propPerSRB$clonal/propPerSRB$Sum)
propPerSRB$frac_subclonal<-(propPerSRB$subclonal/propPerSRB$Sum)
gene.order<-propPerSRB[order(propPerSRB$frac_subclonal),"plot_name"]

p5.input<-melt(propPerSRB[c(4:7)], id.vars = c("plot_name","Sum"))

p5 = ggplot(p5.input, aes(y = value, x = plot_name, fill = as.factor(variable)))+
  scale_fill_manual(values=c("#ef8a62", "#67a9cf")) + 
  scale_x_discrete(limits=gene.order, name="", position="top")+
  scale_y_continuous(name="", limits=c(0,1.05), position = "right")+
  geom_bar(stat = "identity", position = "dodge")+
  theme(axis.text.x = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour="black"),
        axis.line.x.bottom=element_line(colour="black"),
        legend.position="")

ggsave(p5,file='figure_6a.svs.top.pdf', width = 10, height = 1,useDingbats=F)

###create centre heatmap
dfBothPlot<-read.table("cancer_type.sv.list.txt", header=T, sep="\t", stringsAsFactors = F)

plot.order<-unique(dfBothPlot[c("Label","gene_label")])
plot.order$Label<-factor(plot.order$Label, levels=gene.order)
plot.order<-plot.order[order(plot.order$Label),"gene_label"]

p6<-ggplot(dfBothPlot, aes(x = gene_label, y = Cancer, fill=clonal))+
  geom_tile()+
  scale_fill_gradient(low="peachpuff", high="#ef8a62", na.value = "white", limits = c(0,1),name="Fraction of cases Clonal") +
  scale_x_discrete(limits=plot.order, position="top", name="")+
  scale_y_discrete(limits=organ.order, name="", position = "right")+
  scale_size_continuous(range = c(1,6), limits=c(0.01,0.25), name="Subclonal")+
  geom_point(colour = "#67a9cf",aes(size=subclonal))+
  theme(axis.text.x = element_text(angle = 90, hjust=0, vjust=-2),panel.background = element_blank(),
        axis.title.y = element_blank(),  
        panel.grid = element_line(colour="grey90", size=0.5),
        axis.text = element_text(size=12))
ggsave(p6, file='figure_6a.svs.main.pdf', width = 13, height = 10, useDingbats=F)

#######SNVs
###side panel
df<-read.table("figure_6.snv.input.txt", header=T, sep="\t", stringsAsFactors = F)

cl_counts<-ddply(df,.(organ), function(z) sum(z$cl_driver))
names(cl_counts)<-c("organ","clonal")
subcl_counts<-ddply(df,.(organ), function(z) sum(z$subcl_driver))
names(subcl_counts)<-c("organ","subclonal")

organ.counts<-data.frame(table(df$organ))
driver.counts<-merge(cl_counts,subcl_counts,all=T)
driver.counts<-merge(organ.counts, driver.counts, all=T, by.y="organ",by.x="Var1")
driver.counts$clonal<-(driver.counts$clonal/driver.counts$Freq)
driver.counts$subclonal<-(driver.counts$subclonal/driver.counts$Freq)

m2<-melt(driver.counts, id.vars = c("Var1","Freq"))
m2[m2$value==0,"value"]<-NA

p2 = ggplot(m2, aes(y = value, x = Var1, fill = as.factor(variable))) +
  scale_fill_manual(values=c("#ef8a62", "#67a9cf")) + 
  scale_x_discrete(limits=organ.order, name="", position="top")+
  geom_bar(stat = "identity", position = "dodge")+
  scale_y_continuous(name="", limits=c(0,1), position="left")+
  theme(axis.text.x = element_blank(),axis.text.y = element_text(angle = 0, size = 12),panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position="")

ggsave(p2,file='figure_6a.snvs.side.pdf', width = 10, height = 1,useDingbats=F)


#####top panel
sample.counts<-read.table("snv.driver.txt", header=T, sep="\t", stringsAsFactors = F)
filtered.counts<-sample.counts[sample.counts$subclonal>3,]

gene.order<-filtered.counts[order(filtered.counts$frac_subclonal),"gene"]

m3<-melt(filtered.counts[c(1,3,7:8)], id.vars = c("gene","num_mutations"))

p3 = ggplot(m3, aes(y = value, x = gene, fill = as.factor(variable)))+
  scale_fill_manual(values=c("#ef8a62", "#67a9cf")) + 
  scale_x_discrete(limits=gene.order, name="", position="top")+
  scale_y_continuous(name="", limits=c(0,1.05))+
  geom_bar(stat = "identity", position = "dodge")+
  theme(axis.text.x = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(colour="black"),
        axis.text.y = element_text(size = 12),
        axis.line.x.bottom=element_line(colour="black"),
        legend.position="")

ggsave(p3,file='figure_6a.snvs.top.pdf', width = 10, height = 1,useDingbats=F)

###create centre heatmap
samples_wDriver<-read.table("cancer_type.gene.list.txt", sep="\t", header = T, stringsAsFactors = F)
samples.filtered<-samples_wDriver[samples_wDriver$gene%in%gene.order,]

samples.filtered<-merge(samples.filtered, sample.counts[c("gene","num_mutations")], all.x=T, by="gene")
samples.filtered$gene_label<-paste0(samples.filtered$gene," (",samples.filtered$num_mutations,")")

plot.order<-unique(samples.filtered[c("gene","gene_label")])
plot.order$gene<-factor(plot.order$gene, levels=gene.order)
plot.order<-plot.order[order(plot.order$gene),"gene_label"]

samples.filtered[samples.filtered$clonal_samples==0,"clonal_samples"]<-NA
samples.filtered[samples.filtered$subclonal_samples==0,"subclonal_samples"]<-NA

p1<-ggplot(samples.filtered, aes(x = gene_label, y = organ, fill=clonal_samples))+
  geom_tile()+
  scale_fill_gradient(low="peachpuff", high="#ef8a62", na.value = "white",limits = c(0,1), name="Fraction of cases Clonal") +
  scale_x_discrete(limits=plot.order, position="top", name="")+
  scale_y_discrete(limits=organ.order, name="")+
  scale_size_continuous(range = c(1,6), limits=c(0.01,0.25), name="Subclonal")+
  geom_point(colour = "#67a9cf",aes(size=subclonal_samples))+
  theme(axis.text.x = element_text(angle = 90, hjust=0, vjust=-2),panel.background = element_blank(),
        axis.title.y = element_blank(),  
        panel.grid = element_line(colour="grey90", size=0.5),
        axis.text = element_text(size=12))
ggsave(p1,file='figure_6a.snvs.main.pdf', width = 13, height = 10.27, useDingbats=F)

