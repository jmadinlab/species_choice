
# -- mike wd -- # You really got to stop doing this! 
#setwd("/Users/mikemcwilliam/Documents/PostDoc/species_choice")

source("data_prep.R")

# Anyalsis
library(tripack)
library(rgl)

source("R/functions.R")

# points selection

n <- 100 # total species
s <- 20  # species selected
dat <- data.frame(x=rnorm(n), y=rnorm(n), z=rnorm(n))
plot(y ~ x, dat, col="grey")

dat2D <- voronoiFilter(dat, s)
dat3D <- voronoiFilter3D(dat, s)

points(y ~ x, dat2D, col="blue", pch=20, cex=0.6)
points(y ~ x, dat3D, col="red", pch=3)

pp <- scatterplot3d(dat, angle=55, phi=20)
pp$points3d(dat2D, col="blue", pch=20, cex=0.6)
pp$points3d(dat3D, col="red", pch=3)

plot3d( x = dat, type="s", radius=0.1, col=rgb(0,0,0,0.1))
spheres3d(x = dat2D, col = "red", radius = 0.2, alpha=0.3)
spheres3d(x = dat3D, col = "blue", radius = 0.2, alpha=0.3)

# random
#points(y ~ x, dat[sample(1:n, s, replace=FALSE),], col="red", pch=20, cex=0.6)

# most evenly spread










dat<-read.csv("data/data.csv")
dat$Abundance.GBR<-factor(dat$Abundance.GBR, levels=c("rare","uncommon", "common"))


# pnas traits
cats<-c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca<-prcomp(dat[,cats], center=T, scale=T)
points<-data.frame(species=dat$species, pc1=pca$x[,1],pc2=pca$x[,2],pc3=pca$x[,3],pc4=pca$x[,4])
#biplot(pca)





# DIAZ PLOT
library("vegan")
library("ks")

# use princomp for PCA for being consistent with scaling of scores in Sandra's analysis 
prin<-princomp((dat[,cats]), cor = TRUE, scores = TRUE)
pc12<-prin$scores[,1:2]
ll<-prin$loadings

#pc12[,1]<-pc12[,1]*-1
#ll[,1]<-ll[,1]*-1

# # check for consistent results
# pc12[,1]<-pc12[,1]*-1
# #plot(pc12[,1]*-1, pca12[,1])
# plot(pc12[,1], dat[,18])
# cor(pc12[,1], dat[,18])

################ KERNEL DENSITY ESTIMATION ##############################
H <- Hpi(x=pc12)      # optimal bandwidth estimation
est<- kde(x=pc12, H=H, compute.cont=TRUE)     # kernel density estimation

# set contour probabilities for drawing contour levels
cl<-contourLevels(est, prob=c(0.5, 0.05, 0.001), approx=TRUE)

fit<-envfit(pc12, dat[,cats]) # use envfit for drawing arrows, can be also done using trait loadings
fit2<-fit$vectors$arrows*-1 # drawing line segments in arrow opposites direction for pretty layout

fit2<-fit$vectors$arrows*-1
fit2[1,]<-fit2[1,]*3 # leafarea
fit2[2,]<-fit2[2,]*3 # leaf N
fit2[3,]<-fit2[3,]*3 # LMA
fit2[4,]<-fit2[4,]*3 # height
fit2[5,]<-fit2[5,]*3 # seed mass
 fit2[6,]<-fit2[6,]*3 # stem density
 fit2[7,]<-fit2[7,]*3 # stem density


#fit$vectors

#pdf("_PCA_BLOB_FINAL.pdf", width=8, height=6)
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
#dev.off()







# 2D


points<-data.frame(x=jitter(prin$scores[,1], amount=0.25), y=jitter(prin$scores[,2], amount=0.25))
n <- nrow(points) # total species
s <- 20  # species selected

spread<-voronoiFilter(points, s)

#plot(y ~ x, points, col="grey")
# random
#points(y ~ x, points[sample(1:n, s, replace=FALSE),], col="red", pch=20, cex=0.6)
# most evenly spread
#points(y ~ x, spread, col="blue", pch=20, cex=0.6)

spread$species<-dat$species[as.numeric(rownames(spread))]
spread$abundance<-dat$Abundance.GBR[as.numeric(rownames(spread))]

ggplot()+
#geom_segment(data=spread, aes(x=x,y=y,xend=0,yend=0),size=0.15)+
geom_point(data=points, aes(x,y), shape=21, col="grey", size=2)+
geom_point(data=spread, aes(x,y, col=abundance))+
geom_text(data=spread, aes(x,y, label=species), size=2, position=position_nudge(x=0.7, y=0.1))+
scale_colour_manual(values=c("red","orange","green"))+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),legend.title=element_blank(), legend.position=c(0.85,0.9), legend.key.size=unit(2, "mm"))




















