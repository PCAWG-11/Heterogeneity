library(ggplot2)
library(scales)

df<-read.table("nrpcc.summary_table.txt", sep="\t", header=T, stringsAsFactors = F)

myColors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")
names(myColors) <- c(seq(0,7,1), ">7")

p1<-ggplot(data=df,aes(x=nrpcc, y=min_ccfs, group=as.factor(num_subclones), color=as.factor(num_subclones))) +
  geom_point() +
  scale_colour_manual(values=myColors,name="number of subclones") +
  scale_x_continuous("Number of reads per chromosome copy", limits=c(0,50))+
  geom_vline(xintercept = 10, lty=2) + 
  scale_y_continuous("Minimum CCF of detected clusters", breaks=seq(0,1,0.1)) +
  theme(legend.position="right",panel.background = element_blank(),  axis.line = element_line(colour = "black"),panel.grid.major = element_line(size=0.5, colour="grey95"),panel.grid.minor= element_line(size=0.3, colour="grey95"),axis.title=element_text(size=20),axis.text=element_text(size=18), strip.text.x=element_text(size=20),legend.text = element_text(size = 18),legend.title = element_text(size = 18))
ggsave(plot=p1,height=12,width=15, filename=paste("SuppFig_3b.pdf", sep=""), useDingbats=FALSE)
