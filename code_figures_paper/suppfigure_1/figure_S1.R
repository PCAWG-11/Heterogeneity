library(scales)
library(tidyverse)
library(readxl)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.95/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5,0.1, length(na.omit(as.data.frame(cbind(x,y)))[,1]) ,cex=2)
}

##get supplementary data for Aran et al (2015) [10.1038/ncomms9971]
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx", destfile = "aran.purity.xlsx")
published.values <- read_excel("aran.purity.xlsx", skip=2, na="NaN")
published.values <- published.values[c(1:7)]

consensus.purity<-read.table("TCGAsamples.consensus_purity.txt", sep="\t", header=T, stringsAsFactors = F)
merged.purity<-merge(published.values, consensus.purity[c(2,3)], by.x="Sample ID", by.y="submitted_specimen_id")

names(merged.purity)[names(merged.purity)=="IHC"]<-"H&E staining"
names(merged.purity)[names(merged.purity)=="consensus_purity"]<-"PCAWG consensus purity"

pdf("figure_S1.pdf", width=10, height=8)
pairs(merged.purity[,c(3:8)],
      lower.panel=function(...) {smoothScatter(..., nrpoints=Inf, add=T)}, upper.panel=panel.cor,
      pch=".", xlim=c(0,1),ylim=c(0,1),cex.axis=1.5)
dev.off()
