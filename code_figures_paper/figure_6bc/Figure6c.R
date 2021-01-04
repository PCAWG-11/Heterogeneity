#######################################################
## Load results from dNdS runs and reformat
load("aGB.032019.Rda")
names(gbDNDS) <- gsub(".032019.txt.gz","",names(gbDNDS))
out <- gbDNDS
nms <- names(out)
remRCC <- grepl("Kidney-RCC",nms) & (!grepl("clearcell",nms) & !grepl("papillary",nms))
keep <- which(grepl("subclonal",nms) & !grepl("all",nms) &!remRCC)
dnds <- lapply(keep,function(x) out[[x]])
names(dnds) <- gsub("subclonal","",nms[keep])
keepC <- which(!grepl("subclonal",nms) & !grepl("all",nms)&!remRCC)
dndsC <- lapply(keepC,function(x) out[[x]])
names(dndsC) <- gsub("clonal","",nms[keepC])
#######################################################
## Load colours from PCAWG colourPalette
source("colourPalette.R")
pcp <- pcawg.colour.palette(scheme="all",return.scheme=T)
cols <- c(pcp$coding.snv$colours[pcp$coding.snv$levels=="nonsynonymous"],
          pcp$coding.snv$colours[pcp$coding.snv$levels=="stopgain"],
          pcp$coding.snv$colours[pcp$coding.snv$levels=="splicing"],
          "grey15",
          "grey0")
names(cols) <- c("wnon","wnon","wspl","wall")
##########################################################


##########################################################
## Load cohort data and annotation for which sample to keep (i.e. that is powered and representative if multi-sample)
summ <- read.table("summary_table_combined_annotations_v4.txt",header=T,sep="\t")
summ$histology_abbreviation <- summ$histology_detailed
samp2Inc <- read.table("pcawg.wg11.final_sample_list.txt",
                       sep="\t",
                       header=T)
samp2Inc <- samp2Inc[rowSums(samp2Inc[,3:4])==2,1]
keep <- summ$samplename%in%as.character(samp2Inc) & summ$is_preferred
sum(keep)
allcounts <- table(summ[keep,"histology_abbreviation"])
##########################################################


##########################################################
## Count samples per cancer types and re-order by cohort size
Ns <- allcounts[names(dnds)]
ord <- order(Ns,decreasing=T)
neword <- names(dnds)[ord]
dnds <- lapply(neword,function(x) dnds[[x]])
names(dnds) <- neword
dndsC <- lapply(neword,function(x) dndsC[[x]])
names(dndsC) <- neword
##########################################################


##########################################################
## Plotting functions
plotIma<-function(ima,xlab="",ylab="",...)
  {
    plot(1:2, type='n',xlab=xlab,ylab=ylab, xaxt="n",yaxt="n", frame.plot=FALSE, ...)
    rasterImage(ima, 1, 1, 2, 2,interpolate=F)
  }

plotIma <- function(im)
{
    plot(0,0,
         col=rgb(0,0,0,0),
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="",
         frame=F,
         xlim=c(1,2),
         ylim=c(1,2))
    xpas <- 1/nrow(im[,,1])
    ypas <- 1/ncol(im[,,1])
    for(i in 1:nrow(im[,,1]))
    {
        for(j in 1:ncol(im[,,1]))
        {
            polygon(1+c((j-1)*ypas,j*ypas,j*ypas,(j-1)*ypas),
                    2-c((i-1)*xpas,(i-1)*xpas,i*xpas,i*xpas),
                    col=rgb(im[i,j,1],im[i,j,2],im[i,j,3]),
                    border=NA)
        }
    }
}

plotNames <- function(dnds, mut=c(1,2,3,5),CEX=.8)
{
    pasY <- 1/length(dnds)/2
    plot(0,0,col=rgb(0,0,0,0),xlab="",ylab="",yaxt="n",xaxt="n",frame=F,xlim=c(0,1),ylim=c(1,2))
    text(rep(.7,length(dnds)),
         seq(1+pasY,2-pasY,length.out=length(dnds)),
         paste0(names(dnds))[length(dnds):1],pos=2,cex=CEX)
    text(rep(.85,length(dnds)),
         seq(1+pasY,2-pasY,length.out=length(dnds)),
         paste0("(",allcounts[names(dnds)][length(dnds):1],")"),cex=CEX)
}

plotHeatmap <- function(dnds, mut=c(1,2,3,5))
{
    nms <- names(dnds)
    dnds <- lapply(dnds,function(x)
    {
        kk <- try(x[mut,],silent=T)
        if(inherits(kk,"try-error")) NA
        kk
    })
    names(dnds) <- nms
    heatmap <- t(sapply(dnds,function(x)
    {
        kk <- try(x[,"cilow"]>1,silent=T)
        if(inherits(kk,"try-error")) return(c(F,F,F,F))
        kk
    }))
    im <- array(1,dim=c(dim(heatmap),3))
    for(i in 1:4) im[,i,1][heatmap[,i]] <- col2rgb(cols)[1,i]/255
    for(i in 1:4) im[,i,2][heatmap[,i]] <- col2rgb(cols)[2,i]/255
    for(i in 1:4) im[,i,3][heatmap[,i]] <- col2rgb(cols)[3,i]/255
    pasY <- 1/length(dnds)/2
    pasX <- 1/length(mut)/2
    im[,,1] <- im[,,1]
    im[,,2] <- im[,,2]
    im[,,3] <- im[,,3]
    par(mar=c(0,0,0,2))
    plotIma(im)
}

plotAll <- function(dnds,dndsC,COLSEG=rgb(.3,.3,.3,.2))
{
    layout(mat=cbind(c(1,1,1,1,1,1,1),c(2,2,2,2,2,2,2),c(3,3,3,3,3,3,3),c(3,3,2,2,2,4,4)+2),
           widths=c(5.6,2,2,2.5),heights=c(1,1,1,1,1,1,1))
    par(mar=c(0,3,0,1))
    plotNames(dnds,CEX=.4)
    par(mar=c(0,0,0,1))
    plotHeatmap(dnds)
    segments(1,1,1,2,col=COLSEG)
    segments(2,1,2,2,col=COLSEG)
    plotHeatmap(dndsC)
    segments(1,1,1,2,col=COLSEG)
    segments(2,1,2,2,col=COLSEG)
    par(mar=c(0,0,0,0))
    plot(0,0,col=rgb(0,0,0,0),
         xlim=c(0,1),
         ylim=c(0,1),
         xlab="",ylab="",
         frame=F,xaxt="n",yaxt="n")
    legend("left",
           cex=.8,
           pch=15,
           box.col=rgb(0,0,0,0),
           legend=c("dNdS>1","   missense","   nonsense","   splice site" , "   all"),
           text.col=rgb(.5,.5,.5),
           col=c(rgb(0,0,0,0),cols))
}
##########################################################



##########################################################
pdf("dnds.heatmap.COSMIC84.032019.pdf",
    width=3.2,
    height=3.8)
plotAll(dndsC,dnds)
dev.off()
##########################################################


##########################################################
## Two extra tests on dN/dS signal vs. cohort size
##########################################################
## dN/dS needs to be powered to detect selection
## running on low amount of samples/low amount of mutations
## decreases the power to detect selection
##########################################################



##########################################################
## some cancer types didnt have enough mutations to run dN/dS and thus no result data was output
## reading the result data for those cancer types will fail
## this function returns 0 if x leads to an error, i.e. in this case it assumes no positive dN/dS if no data
mytry <- function(x)
{
    kk <- try(x,silent=T)
    if(inherits(kk,"try-error")) return(0)
    kk
}
##########################################################


##########################################################
## Mann-Whitney test for difference in cohort size when there is a positive dN/dS signal
vs <- sapply(dnds,function(x) mytry(max(x[,"mle"])))
ns <- allcounts[names(dnds)]
gs <- split(ns,vs>1)  ##p-value = 0.01648
wilcox.test(gs[[1]],gs[[2]])
##########################################################


##########################################################
## Mann-Whitney test for difference in cohort size when there is a positive & significant dN/dS signal
vs <- sapply(dndsC,function(x) mytry(max(x[,"cilow"])))
ns <- allcounts[names(dndsC)]
gs <- split(ns,vs>1)  ##p-value = 0.0002193
wilcox.test(gs[[1]],gs[[2]])
##########################################################
