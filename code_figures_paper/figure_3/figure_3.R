library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(plyr)
library(ggforce)
library(scales)

source("pcawg.colour.palette.R")

numSubclonesBarplot<-function(power, order){
  all_plot<-power
  
  power<-power[(power$organ!="Other"),]
  
  colorsSub <- c("#E69F00", brewer.pal(4,"Blues")[2:4])
  names(colorsSub) <- c(seq(0,2),"3+")
  
  power$plotSubs<-ifelse(power$num_subclones>=3, "3+", power$num_subclones)
  dat<-prop.table(t(table(power$plotSubs,power$organ)),1)
  all.dat<-prop.table(t(table(power$plotSubs)),1)
  row.names(all.dat)<-"All"
  dat<-rbind(dat, rep.int(0,4))
  dat<-rbind(dat,all.dat)
  
  dat.melt<-melt(dat)
  names(dat.melt)<-c("organ","variable","value")
  dat.melt$value<-as.numeric(dat.melt$value)
  dat.melt$plot<-"num_subclones"
  
  facet_names <- c(`num_subclones` = "# subclones")
  
  sample_counts<-table(power[,"organ"])
  all_counts<-sum(table(all_plot[,"organ"]))
  
  names(all_counts)<-"All"
  final_counts = c(sample_counts, " ", all_counts)
  final_counts<-final_counts[order]
  
  plot1 <- ggplot(data=dat.melt,aes(x=organ,y=value, fill=as.factor(variable))) + geom_bar(width=.9, position="dodge", stat="identity")+ scale_fill_manual(values=colorsSub,name="Number of subclones")+
    scale_y_continuous(name="Fraction", limits=c(0,1.1),breaks=c(0,0.25,0.5,0.75,1))+scale_x_discrete(name="", limits=order, drop=F)+
    theme(panel.background = element_blank(),legend.position="bottom",axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x=element_text(size=12),strip.text.y=element_text(size=16),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.text = element_text(size = 14),legend.title = element_text(size = 14))+annotate("text", x=order, y=rep(1.05, length(order)), label=final_counts)+facet_grid(plot ~ ., labeller = as_labeller(facet_names))
  ggsave(plot=plot1,height=5,width=15, filename=paste("figure_3a.pdf", sep=""), useDingbats=FALSE)
}

makeScurvePanels<-function(power, order){
  allPlot<-power
  power<-power[(power$organ!="Other"),]
  
  power$frac_subclonal = ifelse(is.na(power$frac_clonal_postWCC),0,1-power$frac_clonal_postWCC)
  meltM<-cbind(power[c("organ", "samplename","frac_subclonal")], plot="SNVs")
  names(meltM)<-c("organ","variable","value","plot")
  
  allPlot$frac_subclonal = ifelse(is.na(allPlot$frac_clonal_postWCC),0,1-allPlot$frac_clonal_postWCC)
  meltM.all<-cbind(allPlot[c("organ", "samplename","frac_subclonal")], plot="SNVs")
  names(meltM.all)<-c("organ","variable","value","plot")
  
  counts<-data.frame(table(meltM$organ))
  longestDs<-max(counts$Freq)
  for(i in 1:nrow(counts)){
    row<-counts[i,]
    if(row$Freq<longestDs){
      diff<-(longestDs-row$Freq)
      for(j in 1:floor(diff/2)){
        meltM<-rbind(meltM, c(as.character(row$Var1),stringi::stri_rand_strings(1, 10) ,-1,"SNVs"))
      }
      for(j in 1:ceiling(diff/2)){
        meltM<-rbind(meltM, c(as.character(row$Var1),stringi::stri_rand_strings(1, 10) ,2,"SNVs"))
      }
    }
  }
  for(k in 1:longestDs){
    meltM<-rbind(meltM, c(" ",stringi::stri_rand_strings(1, 10) ,-1,"SNVs"))
    meltM<-rbind(meltM, c("All",stringi::stri_rand_strings(1, 10) ,-1,"SNVs"))
  }
  
  meltM<-meltM[order(match(meltM$organ,order),meltM$value),]
  sample.order<-meltM$variable
  
  dup<-meltM.all
  dup$variable<-paste(dup$variable, "a", sep="_")
  dup$value<--1
  meltM.all<-rbind(meltM.all,dup)
  dup$variable<-paste(dup$variable, "b", sep="_")
  dup$value<-2
  meltM.all<-rbind(meltM.all,dup)
  
  meltM.all<-meltM.all[order(meltM.all$value),]
  sample.order.all<-meltM.all$variable
  
  power$frac_sub_indel = (1-power$frac_clonal_indels)
  meltI<-cbind(power[c("organ", "samplename","frac_sub_indel")], plot="Indels")
  names(meltI)<-c("organ","variable","value","plot")

  buffer<-meltM[meltM$value<0 | meltM$value>1,]
  buffer$plot<-"Indels"
  meltI<-rbind(meltI, buffer)
  
  allPlot$frac_sub_indel = (1-allPlot$frac_clonal_indels)
  meltI.all<-cbind(allPlot[c("organ", "samplename","frac_sub_indel")], plot="Indels")
  names(meltI.all)<-c("organ","variable","value","plot")
  
  plotData<-power
  plotData$frac_aber_subclonal = ifelse((plotData$clonal.arm.consensus+plotData$subclonal.arm.consensus)>9,(plotData$subclonal.arm.consensus/(plotData$clonal.arm.consensus+plotData$subclonal.arm.consensus)),2)
  plotData$frac_aber_clonal = 1-plotData$frac_aber_subclonal
  meltC<-cbind(plotData[c("organ", "samplename","frac_aber_subclonal")], plot="CNAs")
  names(meltC)<-c("organ","variable","value","plot")
  
  plotData.all<-allPlot
  plotData.all$frac_aber_subclonal = ifelse((plotData.all$clonal.arm.consensus+plotData.all$subclonal.arm.consensus)>9,(plotData.all$subclonal.arm.consensus/(plotData.all$clonal.arm.consensus+plotData.all$subclonal.arm.consensus)),2)
  plotData.all$frac_aber_clonal = 1-plotData.all$frac_aber_subclonal
  meltC.all<-cbind(plotData.all[c("organ", "samplename","frac_aber_subclonal")], plot="CNAs")
  names(meltC.all)<-c("organ","variable","value","plot")
  
  plotData<-power
  plotData$frac_sub_sv = ifelse((plotData$num_clonal_svs+plotData$num_subclonal_svs)>9, (1-plotData$frac_clonal_svs),2)
  meltS<-cbind(plotData[c("organ", "samplename","frac_sub_sv")], plot="SVs")
  names(meltS)<-c("organ","variable","value","plot")
  
  plotData.all<-allPlot
  plotData.all$frac_sub_sv = ifelse((plotData.all$num_clonal_svs+plotData.all$num_subclonal_svs)>9, (1-plotData.all$frac_clonal_svs),2)
  meltS.all<-cbind(plotData.all[c("organ", "samplename","frac_sub_sv")], plot="SVs")
  names(meltS.all)<-c("organ","variable","value","plot")
  
  final<-rbind(meltM, meltI,meltS, meltC)
  final$colour<-ifelse((final$organ=="Kidney-RCC.clearcell" | final$organ=="Kidney-RCC.papillary"),"Kidney-RCC",ifelse((final$organ=="Skin-Melanoma.acral" | final$organ=="Skin-Melanoma.cutaneous"),"Skin-Melanoma",final$organ))
  final$pch<-ifelse((final$organ=="Kidney-RCC.papillary" | final$organ=="Skin-Melanoma.acral"), 24,21)
  
  final.all<-rbind(meltM.all, meltI.all,meltS.all, meltC.all)
  
  organ.colours.borders<-organ.colours
  organ.colours.borders[names(organ.colours.borders)=="Lung-AdenoCA"]<-"black"
  organ.colours.borders[names(organ.colours.borders)=="Lung-SCC"]<-"black"
  
  from.x<-seq(floor(longestDs/3),nrow(meltM),by=longestDs)
  to.x<-seq(ceiling(longestDs/3*2),nrow(meltM),by=longestDs)
  line.df<-data.frame(from.x,to.x)
  
  line.df.muts<-line.df
  line.df.muts$plot<-"SNVs"
  line.df.muts$organ<-unique(meltM$organ)
  
  medians<-ddply(meltM, .(organ), function(z) median(as.numeric(z[as.numeric(z$value)!=-1 & as.numeric(z$value)!=2,"value"]), na.rm=T))
  names(medians)<-c("organ","median")
  line.df.muts<-merge(line.df.muts, medians, all.x = T)
  line.df.muts[line.df.muts$organ=="All","median"]<-median((1-allPlot$frac_clonal), na.rm=T)
  
  line.df.ind<-line.df
  line.df.ind$plot<-"Indels"
  line.df.ind$organ<-unique(meltM$organ)
  
  medians<-ddply(meltI, .(organ), function(z) median(as.numeric(z[as.numeric(z$value)!=-1 & as.numeric(z$value)!=2,"value"]), na.rm=T))
  names(medians)<-c("organ","median")
  line.df.ind<-merge(line.df.ind, medians, all.x = T)
  line.df.ind[line.df.ind$organ=="All","median"]<-median((1-allPlot$frac_clonal_indels), na.rm=T)
  
  line.df.sv<-line.df
  line.df.sv$plot<-"SVs"
  line.df.sv$organ<-unique(meltM$organ)
  
  medians<-ddply(meltS, .(organ), function(z) median(as.numeric(z[as.numeric(z$value)!=-1 & as.numeric(z$value)!=2,"value"]), na.rm=T))
  names(medians)<-c("organ","median")
  line.df.sv<-merge(line.df.sv, medians, all.x = T)
  line.df.sv[line.df.sv$organ=="All","median"]<-median(1-(allPlot[(allPlot$num_clonal_svs+allPlot$num_subclonal_svs)>9,]$frac_clonal_svs), na.rm=T)
  
  line.df.cna<-line.df
  line.df.cna$plot<-"CNAs"
  line.df.cna$organ<-unique(meltM$organ)
  
  medians<-ddply(meltC, .(organ), function(z) median(as.numeric(z[as.numeric(z$value)!=-1 & as.numeric(z$value)!=2,"value"]), na.rm=T))
  names(medians)<-c("organ","median")
  line.df.cna<-merge(line.df.cna, medians, all.x = T)
  
  sub.df<-allPlot[(allPlot$clonal.arm.consensus+allPlot$subclonal.arm.consensus)>9,]
  line.df.cna[line.df.cna$organ=="All","median"]<-median((sub.df$subclonal.arm.consensus/(sub.df$subclonal.arm.consensus+sub.df$clonal.arm.consensus)), na.rm=T)
  
  line.df.all<-rbind(line.df.muts,line.df.ind,line.df.sv,line.df.cna)
  
  final$plot = factor(final$plot, levels=c('SNVs','Indels','SVs','CNAs'))
  line.df.all$plot<-factor(line.df.all$plot, levels=c('SNVs','Indels','SVs','CNAs'))
  
  p1<-ggplot(data=final,aes(x=variable, y=as.numeric(value))) +
    geom_point(aes(fill = colour, shape = pch), colour="black", stroke=0.05)+
    scale_shape_identity()+
    scale_x_discrete(name="", limit=sample.order, expand=c(-0.5,0))+
    scale_fill_manual(values=organ.colours,name="")+
    scale_y_continuous(name="Fraction subclonal", limits = c(0,1.1),breaks=c(0,0.25,0.5,0.75,1))+
    geom_segment(aes(x = from.x, y = median, xend = to.x, yend = median), line.df.all, size=1)+
    theme(legend.position="",axis.title=element_text(size=16),axis.text=element_text(size=14),
          strip.text.x=element_blank(),
          strip.text.y=element_text(size=16),panel.grid = element_blank(),panel.grid.major.y = element_line(size=0.5, colour="grey95"),
          axis.text.x = element_blank(),
          legend.text = element_text(angle=45, size = 14),
          panel.background = element_blank(), legend.title = element_text(size = 14),axis.ticks=element_blank())+
    facet_grid(plot ~ .)
  
  p.all<-ggplot(data=final.all,aes(x=variable, y=as.numeric(value))) +
    geom_point(colour = "slategrey",pch=20, stroke=0.05)+
    scale_x_discrete(name="", limit=sample.order.all)+
    scale_y_continuous(name="Fraction subclonal", limits = c(0,1.1),breaks=c(0,0.25,0.5,0.75,1))+
    theme(legend.position="",axis.title=element_blank(),axis.text=element_blank(),
          strip.text.x=element_blank(),strip.text.y=element_blank(),panel.grid = element_blank(),panel.grid.major.y = element_line(size=0.5, colour="grey95"),axis.text.x = element_blank(),legend.text = element_text(angle=45, size = 14),
          panel.background = element_blank(), legend.title = element_text(size = 14),axis.ticks=element_blank())+
    facet_grid(plot ~ .)
  
  ggsave(plot=p1,height=6.5,width=15, filename="figure_3b-e.pdf", useDingbats=FALSE)
  ggsave(plot=p.all,height=6.5,width=1, filename="figure_3b-e.all.pdf", useDingbats=FALSE)
}

mutationBurdenViolin<-function(power, sample.order){
  power$num_muts<-(power$num_clonal+power$num_subclonal)
  
  all_plot<-power
  
  power<-power[(power$organ!="Other"),]
  
  meltN<-cbind(power[c("organ", "samplename","num_muts")], plot="num_muts")
  names(meltN)<-c("organ","variable","value","plot")
  facet_names <- c(`num_muts` = "# mutations")
  
  dup<-cbind(all_plot[c("samplename","num_muts")], plot="num_muts")
  dup[,"organ"]<-"All"
  names(dup)<-c("variable","value","plot","organ")
  
  meltN<-rbind(meltN,dup)
  meltN$colour<-ifelse((meltN$organ=="Kidney-RCC.clearcell" | meltN$organ=="Kidney-RCC.papillary"),"Kidney-RCC",ifelse((meltN$organ=="Skin-Melanoma.acral" | meltN$organ=="Skin-Melanoma.cutaneous"),"Skin-Melanoma",meltN$organ))
  meltN$pch<-ifelse((meltN$organ=="Kidney-RCC.papillary" | meltN$organ=="Skin-Melanoma.acral"), 24,21)
  
  plot1 <- ggplot(data=meltN,aes(x=organ,y=value))+
    geom_violin(scale =  "width",adjust = .5)+
    geom_sina(aes(fill = colour, shape = pch),colour="black", stroke=0.05,size = 2, scale=F, maxwidth=0.7)+
    scale_shape_identity()+
    scale_y_continuous("Total SNVs",trans = 'log10',breaks = trans_breaks('log10', function(x) 10^x),labels = trans_format('log10', math_format(10^.x)))+
    scale_x_discrete(name="", limits=order,drop=F)+
    scale_fill_manual(values=organ.colours,name="")+
    theme(panel.background = element_blank(),legend.position="",axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x = element_blank(),panel.grid = element_blank(),panel.grid.major.y = element_line(size=0.5, colour="grey95"),
          strip.text.x=element_text(size=12),strip.text.y=element_text(size=16),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.text = element_text(size = 14),legend.title = element_text(size = 14))+
    facet_grid(plot ~ ., labeller = as_labeller(facet_names))
  ggsave(plot=plot1,height=1.8,width=15, filename=paste("figure_3f.pdf", sep=""), useDingbats=FALSE)
}

fracAberratedPanel<-function(power, order){
  power$frac_abberated<-(1-power$noCNA.consensus)
  all_plot<-power
  
  power<-power[(power$organ!="Other"),]
  
  meltA<-cbind(power[c("organ", "samplename","frac_abberated")], plot="Genome aberrated")
  names(meltA)<-c("organ","variable","value","plot")
  
  dup<-cbind(all_plot[c("samplename","frac_abberated")], plot="Genome aberrated")
  dup[,"organ"]<-"All"
  names(dup)<-c("variable","value","plot","organ")
  
  meltA<-rbind(meltA,dup)
  
  meltA$colour<-ifelse((meltA$organ=="Kidney-RCC.clearcell" | meltA$organ=="Kidney-RCC.papillary"),"Kidney-RCC",ifelse((meltA$organ=="Skin-Melanoma.acral" | meltA$organ=="Skin-Melanoma.cutaneous"),"Skin-Melanoma",meltA$organ))
  meltA$pch<-ifelse((meltA$organ=="Kidney-RCC.papillary" | meltA$organ=="Skin-Melanoma.acral"), 24,21)
  
  p1 <- ggplot(data=meltA,aes(x=organ,y=value))+
    geom_violin(scale =  "width",adjust = .5)+
    geom_sina(aes(fill = colour, shape = pch),colour="black", stroke=0.05,size = 2, scale=F, maxwidth=0.7)+
    scale_shape_identity()+
    scale_y_continuous("Fraction",breaks = c(0,0.25,0.5,0.75,1))+
    scale_x_discrete(name="", limits=order,drop=F)+
    scale_fill_manual(values=organ.colours,name="")+
    theme(panel.background = element_blank(),legend.position="",axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x = element_blank(),panel.grid = element_blank(),panel.grid.major.y = element_line(size=0.5, colour="grey95"),
          strip.text.x=element_text(size=12),strip.text.y=element_text(size=12),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.text = element_text(size = 14),legend.title = element_text(size = 14))+
    facet_grid(plot ~ .)
    ggsave(plot=p1,height=1.8,width=15, filename=paste("figure_3g.pdf", sep=""), useDingbats=FALSE)
}

wgdHeatmap<-function(power, order){
  
  prop.wgd<-data.frame(prop.table(table(power$organ, power$wgd_status),1))
  prop.wgd<-prop.wgd[prop.wgd$Var2=="wgd",]
  prop.wgd<-rbind(prop.wgd,data.frame(Var1=" ",Var2="wgd",Freq=0))
  all.wgd<-data.frame(prop.table(table(power$wgd_status)))
  prop.wgd<-rbind(prop.wgd,data.frame(Var1="All",Var2="wgd",Freq=all.wgd[all.wgd$Var1=="wgd","Freq"]))
  
  gg <- ggplot(prop.wgd, aes(x=Var1,y=1, fill=Freq))+
    geom_tile(color="white", size=0.1)+
    theme(legend.position="bottom",axis.title=element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),axis.text=element_blank(),strip.text.x=element_blank(),strip.text.y=element_blank(),axis.text.x = element_blank(),legend.text = element_text(angle=45, size = 12),legend.title = element_text(size = 14))+
    scale_x_discrete(limits=order)+
    scale_fill_gradient(low="white",high="red", name="% WGD", limits=c(0, 1))
  
  ggsave(plot=gg,height=1.58,width=15, filename="figure_3h.pdf", useDingbats=FALSE)
  
}

powerHeatmap<-function(power, order){
  
  power<-power[(power$organ!="Other"),]
  
  mean.power<-ddply(power,.(organ),function(z) mean(z$nrpcc))
  prop.power<-rbind(mean.power,data.frame(organ=c("","All"),V1=c(0,mean(power$nrpcc))))
  
  gg <- ggplot(prop.power, aes(x=organ,y=1, fill=V1))+
    geom_tile(color="white", size=0.1)+
    theme(legend.position="bottom",axis.title=element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),axis.text=element_blank(),strip.text.x=element_blank(),strip.text.y=element_blank(),axis.text.x = element_blank(),legend.text = element_text(angle=45, size = 12),legend.title = element_text(size = 14))+
    scale_x_discrete(limits=order)+
    scale_fill_gradient(low="white",high="darkgreen", name="Mean NRPCC", limits=c(10, 30))
  
  ggsave(plot=gg,height=1.58,width=15, filename="figure_3i.pdf", useDingbats=FALSE)
  
}

df<-read.table("figure3.input.txt", header=T, sep="\t", stringsAsFactors = F)

organ.colours<-pcawg.colour.palette(x=c(unique(df$histology_abbreviation)," ","All"),scheme="tumour.subtype")
names(organ.colours)<-c(unique(df$histology_abbreviation)," ", "All")

order<-c(rev(names(sort(tapply(df$frac_clonal_postWCC, df$organ, median, na.rm=TRUE ))))," ","All")
order<-order[order!="Other"]

numSubclonesBarplot(df, order)
makeScurvePanels(df, order)
mutationBurdenViolin(df, order)
fracAberratedPanel(df, order)
wgdHeatmap(df, order)
powerHeatmap(df, order)

