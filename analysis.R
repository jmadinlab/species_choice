# Species choice

# gbr <- read.csv("data/data_20201005.csv", as.is=TRUE)
# dim(gbr)
# library("ggplot2")
# library("cowplot")
# library("FD")
library("vegan")
library("tripack")
library("ks")
# library("car")
library("rgl")
library("spatstat")
library("sp")

source("R/functions.R")
source("R/data_prep.R")

# points
##############################
dat<-read.csv("data/data.csv")
dat$range <- dat$Range.size / max(dat$Range.size)
dat$abund[dat$Abundance.GBR=="common"] <- 1
dat$abund[dat$Abundance.GBR=="uncommon"] <- 0.25
dat$abund[dat$Abundance.GBR=="rare"] <- 0.1
dat$bleach <- 1 - (dat$BRI / 100)

dat2 <- dat[c("species", "range", "abund", "bleach", "restore", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")]
dat2 <- na.omit(dat2)
dim(dat2)


cats<-c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat2[cats], cor=TRUE, scores=TRUE)
# biplot(pca, col=c("grey", "red"), cex=0.6)

space <- data.frame(x=jitter(pca$scores[,1], amount=0.1), y=jitter(pca$scores[,2], amount=0.1), z=jitter(pca$scores[,3], amount=0.1))

dat2 <- cbind(dat2, space)
dim(dat2)

head(dat2)
dat2$species_n <- paste0("Species ", 1:nrow(dat2))

write.csv(dat["species"], "output/species.csv", row.names = FALSE)

# Density in trait space
##############################
H12 <- Hpi(x=space[,1:2])      # optimal bandwidth estimation
H13 <- Hpi(x=space[,c(1, 3)])      # optimal bandwidth estimation
H23 <- Hpi(x=space[,2:3])      # optimal bandwidth estimation
est12 <- kde(x=space[,1:2], H=H12, eval.points=space[,1:2], compute.cont=TRUE)    # kernel density estimation
est13 <- kde(x=space[,c(1, 3)], H=H13, eval.points=space[,c(1, 3)], compute.cont=TRUE)
est23 <- kde(x=space[,2:3], H=H23, eval.points=space[,2:3], compute.cont=TRUE)

dat2$density <- apply(cbind(est12$estimate, est13$estimate, est23$estimate), 1, mean)
dat2$density <- dat2$density / max(dat2$density, na.rm=TRUE)

# trait vectors 
##############################
fit <- envfit(dat2[c("x", "y")], dat2[,cats]) # use envfit for drawing arrows
fit2 <- fit$vectors$arrows*-3 # drawing line segments opposite arrows

png("figs/fig_1.png", width = 4, height = 7, units = 'in', res = 300)

par(mfrow=c(2, 1), mar=c(1.5, 1, 1, 1), oma=c(3, 3, 0, 0))
plot(y ~ x, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
plot(fit, cex=0.90, col=1, labels=list(vectors = c("GR","CW","MCS","SD","CH","SAV","R")))
segments(0,0, fit2[,1], fit2[,2], col=1, lty=2, lwd=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("A", line=-1, adj=0)
## TRAIT METHODS
s <- 20
dat2Dv <- dd2D(dat2, s, vars="", voronoi=TRUE, trait=TRUE)
plot(y ~ x, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(y ~ x, dat2Dv, col="red", pch=1, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
dat2Dn <- dd2D(dat2, s, vars="", voronoi=FALSE, trait=TRUE)
points(y ~ x, dat2Dn, col="black", pch=20, cex=1)
legend("topright", pch=c(20, 1), legend=c("Nearest neighbor", "Voronoi area"), bty="n", col=c("black", "red"), cex=0.8)
arrows(0, 0, -3, 2.5, length=0.1)
text(-3, 2.7, "Cementors", cex=0.8)
arrows(0, 0, 4, -0.5, length=0.1)
text(3.95, -0.2, "Fillers", cex=0.8)
arrows(0, 0, -1.5, -2, length=0.1)
text(-1.5, -2.2, "Builders", cex=0.8)
mtext("B", line=-1, adj=0)
mtext("PC1", 1, line=1, outer=TRUE)
mtext("PC2", 2, line=1, outer=TRUE)

dev.off()

# DESIRABILITY
##############################

png("figs/fig_2.png", width = 3, height = 8, units = 'in', res = 300)

par(mfrow=c(3, 1), mar=c(1.5, 1, 1, 1), oma=c(3, 3, 0, 0))

dat2Drab <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=TRUE)
dat2Drabf <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=FALSE)

plot(y ~ x, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(y ~ x, dat2Drabf, col="black", pch=20, cex=1)
points(y ~ x, dat2Drab, col="blue", pch=1, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("A", line=-1, adj=0)

pts <- chull(dat2Drabf$x, dat2Drabf$y)
pts <- c(pts, pts[1])
polygon(dat2Drabf$x[pts], dat2Drabf$y[pts], col=NA, border="black", lty=3)

x1 <- dat2Drab$x
y1 <- dat2Drab$y
hpts <- chull(x = x1, y = y1)
hpts <- c(hpts, hpts[1])
xy.coords <- cbind(x1, y1)
chull.coords <- xy.coords[hpts,]
chull.poly <- Polygon(chull.coords, hole=F)
chull.poly@area

pts <- chull(dat2Drab$x, dat2Drab$y)
pts <- c(pts, pts[1])
polygon(dat2Drab$x[pts], dat2Drab$y[pts], col=NA, border="blue", lty=3)


legend("topright", pch=c(20, 1), legend=c("Winners", "Winners & function"), bty="n", col=c("black", "blue"), cex=0.8)

dat2Drabr <- dd2D(dat2, s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE)

plot(y ~ x, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(y ~ x, dat2Drabr, col="green", pch=1, cex=1)
points(y ~ x, dat2Drab, col="black", pch=20, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("B", line=-1, adj=0)

legend("topright", pch=c(20, 1), legend=c("Winners & function", "Add restorability"), bty="n", col=c("black", "green"), cex=0.8)

dat2Drabro <- dd2D(dat2, s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE, oppo=TRUE)

plot(y ~ x, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(y ~ x, dat2Drabro, col="red", pch=1, cex=1)
points(y ~ x, dat2Drabr, col="black", pch=20, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("C", line=-1, adj=0)

legend("topright", pch=c(20, 1), legend=c("Winners, function & restorability", "Lossers, function & restorability"), bty="n", col=c("black", "red"), cex=0.8)

mtext("PC1", 1, line=1, outer=TRUE)
mtext("PC2", 2, line=1, outer=TRUE)

dev.off()

####

###

temp <- data.frame(
  Winners=dat2Drabf$species_n,
  Winners.Function=dat2Drab$species_n,
  Winners.Function.Restore=dat2Drabr$species_n
)

temp2 <- data.frame(
  Winners=dat2Drabf$value,
  Winners.Function=dat2Drab$value,
  Winners.Function.Restore=dat2Drabr$value
)

spp <- as.character(unique(unlist(temp)))

mat <- data.frame(matrix(NA, length(spp), ncol(temp)))
row.names(mat) <- spp
names(mat) <- names(temp)
mat2 <- mat

for (j in names(temp)) {
  mat[j] <- spp %in% temp[[j]]
  mat2[spp %in% temp[[j]],j] <- temp2[[j]]
}

mat2[,1] <- (mat2[,1] - min(mat2[,1], na.rm=TRUE)) / max(mat2[,1], na.rm=TRUE)
mat2[,2] <- (mat2[,2] - min(mat2[,2], na.rm=TRUE)) / max(mat2[,2], na.rm=TRUE)
mat2[,3] <- (mat2[,3] - min(mat2[,3], na.rm=TRUE)) / max(mat2[,3], na.rm=TRUE)


png("figs/fig_3.png", width = 5, height = 12, units = 'in', res = 300)

par(mar=c(1, 8, 13, 8))
image(t(as.matrix(mat2)), axes=FALSE, col=hcl.colors(20))
axis(2, at=seq(0, 1, (1/(nrow(mat)-1))), labels=spp, las=2, cex.axis=1, font=3)
axis(3, at=seq(0, 1, (1/(ncol(mat)-1))), labels=c("Winners", "Winners, function", "Winners, function, restorability"), cex.axis=1, las=2)
axis(4, at=seq(0, 1, (1/(nrow(mat)-1)))[mat[["Winners.Function.Restore"]]], labels=spp[mat[["Winners.Function.Restore"]]], las=2, cex.axis=1, font=3)

dev.off()

#### OPPO

dat2Drao <- dd2D(dat2, s, vars=c("abund", "range"), voronoi=FALSE, trait=TRUE, oppo=TRUE)
dat2Drabfo <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=FALSE, oppo=TRUE)
dat2Drabo <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=TRUE, oppo=TRUE)
dat2Drabro <- dd2D(dat2, s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE, oppo=TRUE)


temp <- data.frame(
  Loosers=dat2Drabfo$species_n,
  Loosers.Function=dat2Drabo$species_n,
  Loosers.Function.Restore=dat2Drabro$species_n
)

spp <- as.character(unique(unlist(temp)))

mat <- data.frame(matrix(NA, length(spp), ncol(temp)))
row.names(mat) <- spp
names(mat) <- names(temp)

for (j in names(temp)) {
  mat[j] <- spp %in% temp[[j]]
}


png("figs/fig_4.png", width = 5, height = 12, units = 'in', res = 300)

par(mar=c(1, 8, 13, 8))
image(t(as.matrix(mat)), axes=FALSE, col=c("white", "black"))
axis(2, at=seq(0, 1, (1/(nrow(mat)-1))), labels=spp, las=2, cex.axis=1, font=3)
axis(3, at=seq(0, 1, (1/(ncol(mat)-1))), labels=c("Losers", "Losers, function", "Losers, function, restorability"), cex.axis=1, las=2)
axis(4, at=seq(0, 1, (1/(nrow(mat)-1)))[mat[["Loosers.Function.Restore"]]], labels=spp[mat[["Loosers.Function.Restore"]]], las=2, cex.axis=1, font=3)

dev.off()

### NUMBERS

spp <- sort(as.character(dat2$species_n))
store <- data.frame(spp=spp)

# dd2D(dat2, s=1, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE)
  
for (s in c(25, 20, 15, 10:1)) {
  num <- as.character(dd2D(dat2, s=s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE)$species_n)
  
  store <- cbind(store, data.frame(ss=spp %in% num))
  
}

row.names(store) <- spp
store <- store[,-1]
store <- store[apply(store, 1, sum) > 0,]
store <- store[order(apply(store, 1, sum)),]

png("figs/fig_5.png", width = 5.5, height = 9, units = 'in', res = 300)

par(mar=c(1, 8, 4, 1))
image(t(as.matrix(store)), axes=FALSE, col=c("white", "black"))
axis(2, at=seq(0, 1, (1/(nrow(store)-1))), labels=row.names(store), las=2, cex.axis=1, font=3)
axis(3, at=seq(0, 1, (1/(length(c(25, 20, 15, 10:1))-1))), labels=c(25, 20, 15, 10:1), cex.axis=1, las=2)
mtext("Number of species", side=3, line=2)

dev.off()


###

seg<-cbind(dat3D[,c("x","y","z")], x0=rep(0,nrow(dat3D)), y0=rep(0,nrow(dat3D)),z0=rep(0,nrow(dat3D)))
head(seg)
  #writeWebGL()
plot3d(dat$x, dat$y, dat$z, type="s", xlab="PC1", ylab="PC2", zlab="PC3", size=0.5, box=FALSE, axes=F, col.panel = "black",)
spheres3d(dat3D$x, dat3D$y, dat3D$z,  col="red", radius=0.01)
axes3d(c("x", "y", "z"), col="white")
segments3d(x=as.vector(t(seg[,c("x","x0")])),y=as.vector(t(seg[,c("y","y0")])), z=as.vector(t(seg[,c("z","z0")])), col='red', alpha=0.5)
           bg3d("slategrey")

# 3D -weighted
##############################
dat3Dw <- voronoiFilter3D.w(dat, s)

png("figs/figure3.png", width = 5, height = 5.5, units = 'in', res = 300)

plot(y ~ x, dat, col="grey", cex=1, xlab="PC1", ylab="PC2")
points(y ~ x, dat3D, col="blue", pch=20, cex=0.6)
 points(y ~ x, dat3Dw, col="red", pch=5, cex=1.2)
legend("topleft", pch=c(0, 20, 5), col=c("grey", "blue", "red"), bty="n", legend=c("species", "even spread alone", "even spread & quality"), cex=0.6)

dev.off()


seg<-cbind(dat3Dw[,c("x","y","z")], xend=rep(0,nrow(dat3D)), yend=rep(0,nrow(dat3Dw)), zend=rep(0,nrow(dat3Dw)))
head(seg)
#writeWebGL()
plot3d(dat$x, dat$y, dat$z, type="s", xlab="PC1", ylab="PC2", zlab="PC3", size=0.5, box=FALSE, axes=F, col.panel = "black",)
spheres3d(dat3Dw$x, dat3Dw$y, dat3Dw$z,  col="red", radius=0.01)
axes3d(c("x", "y", "z"), col="white")
segments3d(x=as.vector(t(seg[,c("x","xend")])),y=as.vector(t(seg[,c("y","yend")])), z=as.vector(t(seg[,c("z","zend")])), col='red', alpha=0.5)
bg3d("slategrey")

# 3D -weighted/density
##############################
dat3Dwd <- voronoiFilter3D.w.d(dat, s)

png("figs/figure4.png", width = 5, height = 5.5, units = 'in', res = 300)

plot(y ~ x, dat, col="grey", cex=1, xlab="PC1", ylab="PC2")
points(y ~ x, dat3D, col="blue", pch=20, cex=0.6)
points(y ~ x, dat3Dw, col="red", pch=5, cex=1.2)
points(y ~ x, dat3Dwd, col="green", pch=3, cex=1.5)
legend("topleft", pch=c(0, 20, 5), col=c("grey", "blue", "red", "green"), bty="n", legend=c("species", "even spread alone", "even spread & quality", "even spread, quality & density"), cex=0.6)

dev.off()




# do species have a range of tolerances? 

ggplot()+
geom_histogram(data=dat, aes(x=BRI), bins=15, fill="white", col="grey")+
geom_point(data=dat3Dw, aes(x=BRI, y=1) , fill="red", shape=25)+
#scale_x_log10()+ 
ggtitle("bleaching response index of selected species (red)")







# could the species build a reef? - clipperton
# cemnters are flat, builders have large volume, fillers have branches. 
# all are relatively large.
 
colnames(dat)
bfc<-data.frame(rbind(
apply(dat[dat$cat_colonydiameter>3 & dat$cat_colonyheight>3 & dat$cat_SA_vol<4, c("x","y")], 2, mean), # builders
apply(dat[dat$cat_colonydiameter>3 & dat$cat_spacesize>3, c("x","y")], 2, mean), # fillers
apply(dat[dat$cat_colonydiameter>3 & dat$cat_colonyheight==1, c("x","y")], 2, mean))) # cementers
bfc$type<-c("builders", "fillers", "cementers")

clip<-dat[dat$clipperton==1,]

png("figs/figure5.png", width = 5, height = 5.5, units = 'in', res = 300)

ggplot()+
stat_density2d(data=dat, aes(x,y), size=0.2, bins=8)+
geom_point(data=dat, aes(x,y), shape=21, col="grey", size=0.2)+
geom_point(data=rbind(cbind(dat3Dw[,c("x","y")], type="Selected species"),cbind(dat[dat$clipperton==1,c("x","y")], type="Clipperton species")), aes(x,y, col=type, shape=type))+
geom_point(data=bfc, aes(x,y), size=1, shape=4, stroke=2)+
geom_text(data=bfc, aes(x+0.05,y+0.02, label=type), size=3)+
theme_bw()+theme(panel.grid.minor=element_blank(),legend.position=c(0.85, 0.9), panel.grid.major=element_blank(), legend.title=element_blank())+
labs(x="PC1",y="PC2")+
scale_colour_manual(values=c("black","red"), name="l")+
scale_shape_manual(values=c(16,21), name="l")

dev.off()




























# mike old code


# Find nearest neighbours in each node...
library("FNN")
library("reshape2")

n=7
#gower<-gowdis(dat[,cats])
knn<-get.knn(dat[,c("x","y")],k=n)$nn.index
near<-NULL
for(x in 1:n){
col<-dat$species[knn[,x]]
near<-cbind(near, x=col)}
colnames(near)<-paste("n", seq(1:n), sep="")
near<-as.data.frame(near)
near$species <- dat$species
near<-near[near$species %in% dat3Dw$species,]
near$group<-c(1:nrow(near))
near<-melt(near, id.var=c("group", "species"), value.name="near")
near[,c("x1","y1")]<-dat[match(near$species, dat$species),c("x","y")]
near[,c("x2","y2")]<-dat[match(near$near, dat$species), c("x","y")]
head(near)

# remove duplicates
near<-near[!near$near %in% near$species,]
near<-near[!duplicated(near$near),]

ggplot()+
geom_point(data=dat, aes(x,y), shape=21, col="grey", size=2)+
geom_point(data=dat3Dw, aes(x,y))+
geom_segment(data=near, aes(x=x1, y=y1, xend=x2, yend=y2), size=0.1)+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),legend.title=element_blank(), legend.position=c(0.85,0.9), legend.key.size=unit(2, "mm"))+
labs(x="principal component 1", y="principal component 2")



