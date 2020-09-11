# Functions
library(tripack)


voronoiFilter <- function(dat, s) {
  subset <- dat
  dropped <- vector()
  for (i in 1:(nrow(dat)-s)) {
    v <- voronoi.mosaic(x=subset[,'x'], y=subset[,'y'], duplicate='error')
    info <- cells(v)
    areas <- unlist(lapply(info, function(x) x$area))
    smallest <- which(areas == min(areas, na.rm=TRUE))
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}

voronoiFilter3D <- function(dat, s) {
  subset <- dat
  dropped <- vector()
  for (i in 1:(nrow(dat)-s)) {
    v1 <- voronoi.mosaic(x=subset[,'x'], y=subset[,'y'], duplicate='error')
    info1 <- cells(v1)
    areas1 <- unlist(lapply(info1, function(x) x$area))

    v2 <- voronoi.mosaic(x=subset[,'x'], y=subset[,'z'], duplicate='error')
    info2 <- cells(v2)
    areas2 <- unlist(lapply(info2, function(x) x$area))
    
    v3 <- voronoi.mosaic(x=subset[,'y'], y=subset[,'z'], duplicate='error')
    info3 <- cells(v3)
    areas3 <- unlist(lapply(info3, function(x) x$area))
    
    areas <- apply(cbind(areas1, areas2, areas3), 1, min, na.rm=TRUE)
    
    smallest <- which(areas == min(areas, na.rm=TRUE))
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], dat[,'z'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], subset[smallest,'z'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}

# 

<<<<<<< HEAD






diaz<-function(ll, pc12){
# DIAZ PLOT
################ KERNEL DENSITY ESTIMATION ##############################
library("ks")
library("vegan")

H <- Hpi(x=pc12)      # optimal bandwidth estimation
est<- kde(x=pc12, H=H, compute.cont=TRUE)    # kernel density estimation

# set contour probabilities for drawing contour levels
cl<-contourLevels(est, prob=c(0.5, 0.05, 0.001), approx=TRUE)

fit<-envfit(pc12, dat[,cats]) # use envfit for drawing arrows
fit2<-fit$vectors$arrows*-1 # drawing line segments opposite arrows
fit2[1,]<-fit2[1,]*3 # leafarea
fit2[2,]<-fit2[2,]*3 # leaf N
fit2[3,]<-fit2[3,]*3 # LMA
fit2[4,]<-fit2[4,]*3 # height
fit2[5,]<-fit2[5,]*3 # seed mass
 fit2[6,]<-fit2[6,]*3 # stem density
 fit2[7,]<-fit2[7,]*3 # stem density

par(mar=c(4,4,2,2))
plot(est, cont=seq(1,100,by=1), display="filled.contour2", add=FALSE, ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5),las=1) 
plot(est,abs.cont=cl[1], labels=c(0.5),labcex=0.75, add=TRUE, lwd=0.75, col="grey30")
plot(est,abs.cont=cl[2], labels=c(0.95),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
plot(est,abs.cont=cl[3], labels=c(0.99),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
points( pc12[,], pch=16, cex=0.25, col="black") 
plot(fit, cex=0.90, col=1, labels=list(vectors = c("GR","CW","MCS","SD","CH","SAV","IS")))
segments(0,0, fit2[,1], fit2[,2], col=1, lty=2, lwd=1)
mtext("PC1", cex=0.75, side=1, line=0.5, adj=1)
mtext("PC2", cex=0.75, side=2, line=0.5, at=5.3) #, las=2)
=======
voronoiFilterAb <- function(dat, s) {
  subset <- dat
  dropped <- vector()
  for (i in 1:(nrow(dat)-s)) {
    v <- voronoi.mosaic(x=subset[,'x'], y=subset[,'y'], duplicate='error')
    info <- cells(v)
    areas <- unlist(lapply(info, function(x) x$area))
    
    # Incorporating relative abundance (b/w 0 and 1) simply by multiplying with area
    areas <- areas * subset$abund
    
    # Find smallest value, no longer just the voronoi area
    smallest <- which(areas == min(areas, na.rm=TRUE))
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
>>>>>>> 8c242c5076b173c2b5683a5483a5896beee0aef0
}
