###############################################################
## Load colour library
library(RColorBrewer)
###############################################################
## Load precomputed pairwise scores/similarities
load("pcawg_similarities_before_MTimer_no_wcc.Rdata")
load("simclone1000_similarities_with_truth_before_MTimer_noWCC.Rdata")
###############################################################
## assign truth and simulation scores/similarities
tru <- pcawg.overall.table[-c(18),-c(18)]
simT <- sim <- simclone1000.overall.table
###############################################################
## Plotting parameters
lwdSEG <- 1.5
COL1 <- rgb(.7,.2,.3)
COL2 <- rgb(.2,.2,.6)
###############################################################



###############################################################
## -------------- functions     ---------------------------- ##
###############################################################
getColnames <- function(methods)
{
    cols <- col2rgb(RColorBrewer:::brewer.pal(12,"Paired"))/255
    N1 <- 4## C
    N2 <- 10 ## I
    N3 <- 6## R
    sapply(methods,function(x) if(x%in%c("CICC","CSR","WM"))
                                   return(rgb(cols[1,N1],cols[2,N1],cols[3,N1]))
                               else if(grepl("rd_",x))
                                   return(rgb(cols[1,N3],cols[2,N3],cols[3,N3]))
                               else
                                   return(rgb(cols[1,N2],cols[2,N2],cols[3,N2])))
}

getColnames <- function(methods)
{
    cols <- cbind(c(124,190,243),c(214,173,254),c(255,168,113))/255
    cols <- (cbind(c(200,200,200),c(150,150,150),c(55,55,55))/255)[,3:1][,c(2,1,3)]
    N1 <- 1
    N2 <- 2
    N3 <- 3
    sapply(methods,function(x) if(x%in%c("CICC","CSR","WM"))
                                   return(rgb(cols[1,N1],cols[2,N1],cols[3,N1]))
                               else if(grepl("rd_",x))
                                   return(rgb(cols[1,N3],cols[2,N3],cols[3,N3]))
                               else if(x=="Truth")
                                   return(rgb(0,0,0))
                               else
                                   return(rgb(cols[1,N2],cols[2,N2],cols[3,N2])))
}

getCol1 <- function(x)
{
    vec <- c(sapply(seq(0,1,length.out=100),function(x) rgb(.7,.1,.3,(1-x)*.9+.1)),
             sapply(seq(0,1,length.out=100),function(x) rgb(.3,.1,.7,.1+x*.9)))
    vec[round(x*199)+1]
}

getCol3 <- function(x)
{
    vec <- c(sapply(seq(0,1,length.out=100),function(x) rgb(.7,.4,.3,(1-x)*.9+.1)),
             sapply(seq(0,1,length.out=100),function(x) rgb(.3,.4,.7,.1+x*.9)))
    vec[round(x*199)+1]
}

###############################################################
## this is what a pre-cancerous code looks like:
getCol2 <- getCol1
getCol1 <- getCol3
getCol3 <- getCol2
getCol2 <- getCol1
###############################################################

getTEXT <- function(x)
{
    if(x%in%c("CICC","Weme","CSR")) return("C")
    if(grepl("rd",x)) return("R")
    else return("I")
}

getCOLTEXT <- function(x)
{
    cols <- RColorBrewer:::brewer.pal(12,"Paired")
    cols <- paste0(cols,"BB")
    N1 <- 4## C
    N2 <- 10 ## I
    N3 <- 6## R
    if(x%in%c("CICC","Weme","CSR")) return(cols[N1])
    if(grepl("rd",x)) return(cols[N3])
    else return(cols[N2])
}

getCOLTEXT <- function(x)
{
    cols <- cbind(c(124,190,243),c(214,173,254),c(255,168,113))/255
    cols <- (cbind(c(200,200,200),c(150,150,150),c(55,55,55))/255)[,3:1][,c(2,1,3)]
    N1 <- 1
    N2 <- 2
    N3 <- 3
    paste0(sapply(x,function(x) if(x%in%c("CICC","CSR","Weme"))
                                   return(rgb(cols[1,N1],cols[2,N1],cols[3,N1]))
                               else if(grepl("rd_",x))
                                   return(rgb(cols[1,N3],cols[2,N3],cols[3,N3]))
                               else if(grepl("Truth",x))
                                   return(rgb(0,0,0))
                               else
                                   return(rgb(cols[1,N2],cols[2,N2],cols[3,N2]))),"BB")
}

addBorder <- function(i,j,xpas,ypas,N=18,
                      col1=COL1,
                      col2=COL2,
                      lwds=lwdSEG)
{
    if(i==1 & j==1)
    {
        segments((i-1)*xpas,j*ypas,i*xpas,j*ypas,col=col1,lwd=lwds)
        segments((i)*xpas,j*ypas,i*xpas,(j-1)*ypas,col=col2,lwd=lwds)
    }
    else if(i==j)
    {
        segments((i-1)*xpas,j*ypas,i*xpas,j*ypas,col=col1,lwd=lwds)
        segments((i-1)*xpas,(j-1)*ypas,(i-1)*xpas,(j)*ypas,col=col1,lwd=lwds)
        segments((i)*xpas,j*ypas,i*xpas,(j-1)*ypas,col=col2,lwd=lwds)
        segments((i-1)*xpas,(j-1)*ypas,i*xpas,(j-1)*ypas,col=col2,lwd=lwds)
    }
    else
    {
        if(i==1)
        {
            segments((i-1)*xpas,(j-1)*ypas,(i-1)*xpas,(j)*ypas,col=col1,lwd=lwds)
        }
        if(i==N & j!=N)
        {
            segments((i)*xpas,j*ypas,i*xpas,(j-1)*ypas,col=col2,lwd=lwds)
        }
        if(j==2)
        {
            segments((i-1)*xpas,(j-1)*ypas,i*xpas,(j-1)*ypas,col=col2,lwd=lwds)
        }
        if(j==N & i!=N)
        {
            segments((i-1)*xpas,j*ypas,i*xpas,j*ypas,col=col1,lwd=lwds)
        }
    }
}

myMap <- function(x,y,cexTEXT=.9)
{
    plot(0,0,xlim=c(0,1),ylim=c(0,1),
         xlab="",ylab="",
         xaxt="n",yaxt="n",
         col=rgb(0,0,0,0),frame=F)
    xpas <- 1/ncol(x)
    ypas <- 1/nrow(x)
    if(F){
        x <- x-min(x,na.rm=T)
        x <- x/max(x,na.rm=T)}
    ord <- 1:nrow(x)
    x <- x[ord,ord]
    bord <- 0
    for(i in 1:nrow(x))
    {
        for(j in 1:ncol(x))
        {
            if(i==1) polygon(c((i-1)*xpas,(i-1)*xpas,(i)*xpas-bord,(i)*xpas-bord),
                            c((j-1)*ypas,(j)*ypas-bord,(j)*ypas-bord,(j-1)*ypas),
                            col=getCol2(x[i,j]),border=NA)
            else if(j>i)
                polygon(c((i-1)*xpas,(i-1)*xpas,(i)*xpas-bord,(i)*xpas-bord),
                        c((j-1)*ypas,(j)*ypas-bord,(j)*ypas-bord,(j-1)*ypas),
                        col=getCol1(x[i,j]),border=NA)
            else polygon(c((i-1)*xpas,(i-1)*xpas,(i)*xpas-bord,(i)*xpas-bord),
                         c((j-1)*ypas,(j)*ypas-bord,(j)*ypas-bord,(j-1)*ypas),
                         col=getCol3(x[i,j]),border=NA)
            if(i==j)
            {
                bb <- xpas/20
                polygon(c((i-1)*xpas+bb,(i-1)*xpas+bb,(i)*xpas-bb,(i)*xpas-bb),
                        c((j-1)*ypas+bb,(j)*ypas-bb,(j)*ypas-bb,(j-1)*ypas+bb),
                        col=getCOLTEXT(colnames(x)[i]),border=NA)
                if(F) text((i-1)*xpas+xpas/2,
                (j-1)*ypas+ypas/2,
                getTEXT(colnames(x)[i]),
                cex=cexTEXT,
                col=getCOLTEXT(colnames(x)[i]))
            }
        }
    }
}

plotHeatmap <- function(tru,sim)
{
    all <- rbind(rep(NA,18),cbind(rep(NA,17),tru))
    rownames(all) <- colnames(all) <- colnames(sim)
    for(i in 1:nrow(sim))
        for(j in i:ncol(sim))
        {
            all[i,j] <- sim[i,j]
        }
    layout(mat=cbind(c(1,1,1,1,1,1,1),c(5,2,2,6,3,3,4)),
           widths=c(4,1.3),
           heights=c(1,1,1,1,1,1,1))
    par(mar=c(1,0,1.3,0))
    myMap(all,NA)
    NN <- 18
    par(mar=c(0,0,0,0))
    plot(0,0,
         col=rgb(0,0,0,0),
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="",
         frame=F,
         xlim=c(0,1),
         ylim=c(0,1))
    legend("left",
           legend=c("Truth", "Consensus","Individual","Random"),
           col=c(getCOLTEXT("Truth"),
                 getCOLTEXT("CICC"),
                 getCOLTEXT("Phylogic"),
                 getCOLTEXT("rd_")),
           pch=15,##c("C","I","R"),
           cex=.78,box.col=rgb(0,0,0,0))
    par(mar=c(0,0,0,0))
    plot(0,0,
         col=rgb(0,0,0,0),
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="",
         frame=F,
         xlim=c(0,1),
         ylim=c(0,1))
    xx1 <- 0
    xx2 <- .2
    ss <- seq(0,1,length.out=50)
    scale <- .6
    for(i in 2:length(ss))
    {
        polygon(c(xx1,xx1,xx2,xx2),
                c(ss[i-1],ss[i],ss[i],ss[i-1])*scale,
                col=getCol1(ss[i-1]),border=NA)
    }
    for(i in 2:length(ss))
    {
        polygon(c(xx1,xx1,xx2,xx2)+.3,
                c(ss[i-1],ss[i],ss[i],ss[i-1])*scale,
                col=getCol3(ss[i-1]),border=NA)
    }
    text(.15,0,"no similarity",cex=.68,pos=4)
    text(.15,1*scale,"high similarity",cex=.68,pos=4)
}

minmaxNorm <- function(tru)
{
    tru <- tru-min(tru,na.rm=T)
    tru <- tru/max(tru,na.rm=T)
    tru
}
minmaxNorm.truth <- function(simT,N=18)
{
    simT[N,] <- minmaxNorm(simT[N,])
    simT[,N] <- minmaxNorm(simT[,N])
    simT
}
###############################################################




###############################################################
set.seed(1234) ##not used
###############################################################
## order methods
ORD <- c(18,15,16,17,
         11:1,
         14:12)
###############################################################
## normalise scores separately
tru <- minmaxNorm(tru)
simT <- minmaxNorm(simT)
simT <- minmaxNorm.truth(simT)
###############################################################
## plot heatmap and colour key
pdf("heatmap.similarities.pdf",
    width=3.12,
    height=2.6)
plotHeatmap(tru[ORD[-c(1)],ORD[-c(1)]],simT[ORD,ORD])
dev.off()
###############################################################
## The rest is just having fun in inkscape...
###############################################################
