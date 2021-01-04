library(ggplot2)
library(RColorBrewer)

wgdStatusSeparation<-function(wgd_anno){
  
  p1<-ggplot(data=wgd_anno, aes(y=ploidy, x=hom, colour=wgd_status))+
    geom_point()+
    scale_colour_manual(values=c("#404040","#E69F00"),labels=c("No WGD","WGD"),name="")+
    scale_y_continuous(name="Ploidy")+
    scale_x_continuous(name="Fraction of genome with LOH")+
    geom_abline(intercept = 2.9, slope = -2, lty=2, colour="#009E73")+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank(),legend.position="bottom",axis.title=element_text(size=22),axis.text=element_text(size=20),panel.grid.major.y = element_line(size=0.5, colour="grey95"),
          panel.grid.major.x = element_line(size=0.5, colour="grey95"), strip.text.x=element_text(size=18),strip.text.y=element_text(size=22),legend.text = element_text(size = 20),legend.title = element_text(size = 20))
  ggsave(plot=p1,height=10,width=12, filename=paste("1c.wgdSeparation.pdf", sep=""), useDingbats=FALSE)

}

input<-read.table("wgd.status.txt", stringsAsFactors = F, sep = "\t", header=T)
wgdStatusSeparation(input)

