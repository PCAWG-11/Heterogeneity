library(plyr)
library(GGally)

source("pcawg.colour.palette.R")

subclonalCorrPlot<-function(power){
  power$frac_aber_subclonal = ifelse((power$clonal.arm.consensus+power$subclonal.arm.consensus)>4,(power$subclonal.arm.consensus/(power$clonal.arm.consensus+power$subclonal.arm.consensus)),NA)
  power$frac_sub_sv = ifelse((power$num_clonal_svs+power$num_subclonal_sv)>4, (1-power$frac_clonal_svs),NA)

  power$frac_subclonal<-(1-power$frac_clonal_postWCC)
  power$frac_subclonal_indels<-(1-power$frac_clonal_indels)
  power$mut.burden<-(power$num_clonal+power$num_subclonal)
  power$frac_abberated<-(1-power$noCNA.consensus)
  
  power<-power[power$organ!="Other",]
  
  p.values<-ddply(power, .(organ), function(z) cor.test(z$frac_subclonal,z$mut.burden)$p.value)
  p.values$adjusted<-p.adjust(p.values$V1, method="BH")
  p.values$adjusted<-round(p.values$adjusted,3)
  p.values$V1<-round(p.values$V1,3)
  
  p.values<-ddply(power, .(organ), function(z) cor.test(z$num_subclones,z$mut.burden, method = "spearman")$p.value)
  p.values$adjusted<-p.adjust(p.values$V1)
  p.values$adjusted<-round(p.values$adjusted,3)
  p.values$V1<-round(p.values$V1,3)
  
  p.values<-ddply(power[complete.cases(power$frac_aber_subclonal, power$frac_abberated),], .(organ), function(z) ifelse(nrow(z)>10, cor.test(z$frac_aber_subclonal,z$frac_abberated)$p.value, NA))
  p.values$adjusted<-p.adjust(p.values$V1)
  
  p.values<-ddply(power[complete.cases(power$frac_sub_sv, power$frac_abberated),], .(organ), function(z) ifelse(nrow(z)>10, cor.test(z$frac_sub_sv,z$frac_abberated)$p.value, NA))
  p.values$adjusted<-p.adjust(p.values$V1)
  
  power$colour<-ifelse((power$organ=="Kidney-RCC.clearcell" | power$organ=="Kidney-RCC.papillary"),"Kidney-RCC",ifelse((power$organ=="Skin-Melanoma.acral" | power$organ=="Skin-Melanoma.cutaneous"),"Skin-Melanoma",power$organ))
  power$pch<-ifelse((power$organ=="Kidney-RCC.papillary" | power$organ=="Skin-Melanoma.acral"), 24,21)
  power$lineCol<-ifelse(grepl("Lung",power$colour), "Skin-Melanoma", power$colour)
  
  plotData<-power[c("samplename","colour","pch", "frac_subclonal","frac_subclonal_indels","frac_sub_sv","frac_aber_subclonal","lineCol")]
  names(plotData)<-c("samplename","colour","pch","SNVs","Indels","SVs","CNAs","lineCol")
  
  ggally_mysmooth <- function(data, mapping, ...){
    ggplot(data = data, mapping=mapping) +
      geom_density(aes(y=..scaled..),fill="lightgrey", alpha=0.7)
  }
  
  pm<-ggpairs(plotData, columns=4:7, lower= "blank", upper=list(continuous = "points", combo = "dot_no_facet"), ggplot2::aes(fill=colour,colour=lineCol, shape=pch, stroke=0.05),diag = "blank")
  for(i in 1:pm$nrow) {
    for(j in 1:pm$ncol){
      pm[i,j] <- pm[i,j] + 
        scale_fill_manual(values=organ.colours) + scale_colour_manual(values=organ.colours) + scale_shape_identity()+theme(legend.position="",axis.title=element_text(size=16),axis.text=element_text(size=10),
                                                                                                                           strip.text=element_text(size=16),panel.grid = element_blank(),panel.grid.major= element_line(size=0.5, colour="grey95"),legend.text = element_text(size = 14),
                                                                                                                           panel.background = element_blank(), legend.title = element_text(size = 14))+
        scale_x_continuous(limits=c(0,1))+
        scale_y_continuous(limits=c(0,1))
    }
  }
  ggsave(plot=pm,height=10,width=10, filename="SuppFigure_S4.upper.pdf", useDingbats=FALSE)
  
  pm<-ggpairs(plotData, columns=4:7, lower= "blank", upper="blank", ggplot2::aes(shape=pch, stroke=0.05),diag = list(continuous = ggally_mysmooth))
  for(i in 1:pm$nrow) {
    for(j in 1:pm$ncol){
      pm[i,j] <- pm[i,j] + 
        scale_fill_manual(values=organ.colours) + scale_colour_manual(values=organ.colours) + scale_shape_identity()+theme(legend.position="",axis.title=element_text(size=16),axis.text=element_text(size=10),
                                                                                                                           strip.text=element_text(size=16),panel.grid = element_blank(),panel.grid.major= element_line(size=0.5, colour="grey95"),legend.text = element_text(size = 14),
                                                                                                                           panel.background = element_blank(), legend.title = element_text(size = 14))+
        scale_x_continuous(limits=c(0,1))+
        scale_y_continuous(limits=c(0,1))
    }
  }
  ggsave(plot=pm,height=10,width=10, filename="SuppFigure_S4.diag.pdf", useDingbats=FALSE)
  
  plotData<-power[c("organ","colour","pch", "frac_subclonal","frac_subclonal_indels","frac_sub_sv","frac_aber_subclonal")]
  names(plotData)<-c("organ","colour","pch","SNVs","Indels","SVs","CNAs")
  maps<-unique(power[c("organ","colour","pch")])
  
  print(cor.test(plotData$SNVs,plotData$Indels, method = "pearson"))
  print(cor.test(plotData$SNVs,plotData$SVs, method = "pearson"))
  print(cor.test(plotData$SNVs,plotData$CNAs, method = "pearson"))
  print(cor.test(plotData$Indels,plotData$SVs, method = "pearson"))
  print(cor.test(plotData$Indels,plotData$CNAs, method = "pearson"))
  print(cor.test(plotData$SVs,plotData$CNAs, method = "pearson"))
  
  snvs.summ<-ddply(plotData,.(organ), function(z) data.frame(median=median(z$SNVs, na.rm = T),high=quantile(z$SNVs, na.rm = T)[4],low=quantile(z$SNVs, na.rm = T)[2]))
  indels.summ<-ddply(plotData,.(organ), function(z) data.frame(median=median(z$Indels, na.rm = T),high=quantile(z$Indels, na.rm = T)[4],low=quantile(z$Indels, na.rm = T)[2]))
  cnas.summ<-ddply(plotData,.(organ), function(z) data.frame(median=median(z$CNAs, na.rm = T),high=quantile(z$CNAs, na.rm = T)[4],low=quantile(z$CNAs, na.rm = T)[2]))
  svs.summ<-ddply(plotData,.(organ), function(z) data.frame(median=median(z$SVs, na.rm = T),high=quantile(z$SVs, na.rm = T)[4],low=quantile(z$SVs, na.rm = T)[2]))
  
  df.vec<-list(SNVs=snvs.summ,Indels=indels.summ,SVs=svs.summ,CNAs=cnas.summ)
  for(i in 1:(length(df.vec)-1)){
    first.df<-df.vec[[i]]
    mut1<-names(df.vec)[i]
    for(j in (i+1):length(df.vec)){
      second.df<-df.vec[[j]]
      mut2<-names(df.vec)[j]
      comb<-merge(first.df, second.df, by="organ", all=T)
      comb<-merge(comb, maps, by="organ", all=T)
      p1<-ggplot(comb, aes(x=median.x, y=median.y,fill=colour, shape=pch))+geom_point(size=2.5, stroke=0.05)+
        scale_fill_manual(values=organ.colours) + scale_shape_identity()+theme(legend.position="",axis.title=element_blank(),axis.text=element_blank(),
                                                                               strip.text=element_text(size=16),panel.grid = element_blank(),panel.grid.major= element_line(size=0.5, colour="grey95"),legend.text = element_text(size = 14),
                                                                               panel.background = element_blank(), legend.title = element_text(size = 14), axis.ticks = element_blank())+
        scale_x_continuous(limits=c(0,1))+
        scale_y_continuous(limits=c(0,1))+
        geom_errorbar(aes(ymin=low.y, ymax=high.y, colour=colour), width=0)+
        geom_errorbarh(aes(xmin=low.x, xmax=high.x, colour=colour))+
        scale_colour_manual(values=organ.colours)
      
      print(ggsave(plot=p1,height=2.5,width=2.5, filename=paste(mut1,"_",mut2,"SuppFig_4.lower.pdf",sep=""), useDingbats=FALSE))
    }
  }
}

df<-read.table("input.figureS4.txt", header=T, sep="\t", stringsAsFactors = F)

organ.colours<-pcawg.colour.palette(x=c(unique(df$histology_abbreviation)," ","All"),scheme="tumour.subtype")
names(organ.colours)<-c(unique(df$histology_abbreviation)," ", "All")

subclonalCorrPlot(df)


