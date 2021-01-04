library(ggplot2)
library(scales)

df<-read.table("nrpcc.summary_table.txt", sep="\t", header=T, stringsAsFactors = F)

myColors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")
names(myColors) <- c(seq(0,7,1), ">7")

###normal overview (exclude PD4120)
p1<-ggplot(data=df,aes(x=(num_clonal+num_subclonal), y=nrpcc, group=as.factor(num_subclones), color=as.factor(num_subclones))) +
  geom_point() +
  scale_colour_manual(values=myColors,name="number of subclones") +
  scale_x_continuous("Total SNVs (log10)",trans = 'log10',breaks = trans_breaks('log10', function(x) 10^x),labels = trans_format('log10', math_format(10^.x)))+
  geom_hline(yintercept = 10, lty=2) + 
  scale_y_continuous("Number of reads per chromosome copy", limits=c(0,50)) +
  annotation_logticks(sides = "b")+
  theme(legend.position="right",panel.background = element_blank(),  axis.line = element_line(colour = "black"),panel.grid.major = element_line(size=0.5, colour="grey95"),panel.grid.minor= element_line(size=0.3, colour="grey95"),axis.title=element_text(size=20),axis.text=element_text(size=18), strip.text.x=element_text(size=20),legend.text = element_text(size = 18),legend.title = element_text(size = 18))
ggsave(plot=p1,height=12,width=15, filename=paste("SuppFig_3a.pdf", sep=""), useDingbats=FALSE)
