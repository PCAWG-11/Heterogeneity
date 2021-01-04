#######################################################
## Load results from dNdS runs
load("aGB.032019.Rda")
names(gbDNDS) <- gsub(".032019.txt.gz","",names(gbDNDS))
lDndsout <- lD <- gbDNDS
#######################################################
## Load colours from PCAWG colourPalette
source("colourPalette.R")
pcp <- pcawg.colour.palette(scheme="all",return.scheme=T)
cols <- c(pcp$coding.snv$colours[pcp$coding.snv$levels=="nonsynonymous"],
          pcp$coding.snv$colours[pcp$coding.snv$levels=="stopgain"],
          pcp$coding.snv$colours[pcp$coding.snv$levels=="splicing"],
          "grey15",
          "grey0")
##########################################################



#######################################################
addIntervals <- function(xx,y1,y2,xd=.3,col="grey",LWD=3)
{
    segments(xx,y1,xx,y2,lwd=LWD,col=col)
    segments(xx-xd,y1,xx+xd,y1,lwd=LWD,col=col)
    segments(xx-xd,y2,xx+xd,y2,lwd=LWD,col=col)
}

plotDNDS <- function(bp,
                     bpD,
                     bpU,
                     names.lSt,
                     allns,
                     cols=c("grey60","grey45", "grey30","grey15","grey0"),
                     keep,
                     space=NULL,
                     ylim=c(0,5.5))
{
    trans <- .9
    coeff <- .88
    layout(mat=cbind(c(1,1,1),c(3,2,4)),widths=c(65,10),heights=c(1,1,1))
    par(mar=c(3,4.5,1,1))
    kk <- barplot(bp,
                  border=NA,
                  ylim=ylim,
                  beside=T,space=space,
                  ylab="dN/dS (566 cancer genes)",
                  cex.lab=1.2,
                  cex.axis=1.5,
                  col=rgb(0,0,0,0))
    grid()
    kk <- barplot(bp,
                  border=NA,
                  beside=T,space=space,
                  yaxt="n",
                  col=rgb(0,0,0,0),
                  add=T)
    abline(h=1,lty=2,col=rgb(.1,.1,.1,.5),lwd=1.6)
    addIntervals(as.vector(kk),bpD,bpU,col=cols)
    points(as.vector(kk),bp,col=cols,pch=19,cex=1.4)
    xx <- as.vector(kk)
    inter <- diff(xx)[1]/2
    axis(side=1,
         at=xx[seq(3,length(kk),length.out=length(names.lSt))]-inter,
         cex.axis=.6,
         names.lSt,
         tick=F)
    par(mar=c(0,0,0,0))
    plot(0,0,col=rgb(0,0,0,0),
         xlim=c(0,1),
         ylim=c(0,1),
         xlab="",ylab="",
         frame=F,xaxt="n",yaxt="n")
    legend("left",
           cex=.7,pch=15,
           box.col=rgb(0,0,0,0),
           legend=c("missense","nonsense","splice site" ,"truncating", "all")[keep],
           text.col=rgb(.5,.5,.5),
           col=cols)
}
##########################################################


##########################################################
## Keep results from running on all samples clonal vs. subclonal
print(names(lD)[1:2])
lDndsout <- list("Clonal"=lD[[1]],
                 "Sublonal"=lD[[2]])
#######################################################
## Plot results
pdf("Figure6b.dNdS.COSMIC84.032019.pdf",width=4.53,height=2.45)
whichtoplot <- 1:length(lDndsout)
keep <- c(1,2,3,5)
plotDNDS(bp=unlist(lapply(whichtoplot,function(x) lDndsout[[x]]$mle[keep])),
         bpD=unlist(lapply(whichtoplot,function(x) lDndsout[[x]]$cilow[keep])),
         bpU=unlist(lapply(whichtoplot,function(x) lDndsout[[x]]$cihigh[keep])),
         names.lSt=names(lDndsout),
         allns=c("",""),
         cols=cols[keep],
         keep=keep,
         space=unlist(lapply(length(whichtoplot),function(x) c(1,rep(0,length(keep)-1)))),
         ylim=c(0,3))
dev.off()
#########################################################







