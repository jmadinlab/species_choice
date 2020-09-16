# Species choice

library("ggplot2")
library("cowplot")
library("FD")
library("vegan")
library("tripack")
library("ks")
library("car")
library("rgl")

source("R/functions.R")
source("R/data_prep.R")

# points
##############################
dat<-read.csv("data/data.csv")
dat$Abundance.GBR <- factor(dat$Abundance.GBR, levels=c("rare","uncommon", "common"))

cats<-c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pcoa <- pcoa(gowdis(dat[,cats]))

space <- data.frame(x=jitter(pcoa$vectors[,1], amount=0.01), y=-jitter(pcoa$vectors[,2], amount=0.01), z=jitter(pcoa$vectors[,3], amount=0.01))
dat <- cbind(dat, space)

write.csv(dat["species"], "output/species.csv", row.names = FALSE)

# Density in trait space
##############################
H12 <- Hpi(x=space[,1:2])      # optimal bandwidth estimation
H13 <- Hpi(x=space[,c(1, 3)])      # optimal bandwidth estimation
H23 <- Hpi(x=space[,2:3])      # optimal bandwidth estimation
est12 <- kde(x=space[,1:2], H=H12, eval.points=space[,1:2], compute.cont=TRUE)    # kernel density estimation
est13 <- kde(x=space[,c(1, 3)], H=H13, eval.points=space[,c(1, 3)], compute.cont=TRUE)
est23 <- kde(x=space[,2:3], H=H23, eval.points=space[,2:3], compute.cont=TRUE)

dat$density <- apply(cbind(est12$estimate, est13$estimate, est23$estimate), 1, mean)
dat$density <- dat$density / max(dat$density, na.rm=TRUE)

hist(dat$density)
plot(space[,1:2], cex=dat$density*10)

# trait vectors 
##############################
prin<-princomp((dat[,cats]), cor = TRUE, scores = TRUE)
pc12<-prin$scores[,1:2]
#ll<-prin$loadings
fit<-envfit(pc12, dat[,cats]) # use envfit for drawing arrows
fit2<-fit$vectors$arrows*-3 # drawing line segments opposite arrows

plot( pc12[,], pch=16, cex=0.5, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5),las=1) 
plot(fit, cex=0.90, col=1, labels=list(vectors = c("GR","CW","MCS","SD","CH","SAV","R")))
segments(0,0, fit2[,1], fit2[,2], col=1, lty=2, lwd=1)
mtext("PC1", cex=0.75, side=1, line=0.5, adj=1)
mtext("PC2", cex=0.75, side=2, line=0.5, at=3.5) 

# DESIRABILITY
##############################

# Range size
dat$Range.size <- dat$Range.size / max(dat$Range.size)

# Abundance
dat$abund <- NA
dat$abund[dat$Abundance.GBR=="common"] <- 1
dat$abund[dat$Abundance.GBR=="uncommon"] <- 0.25
dat$abund[dat$Abundance.GBR=="rare"] <- 0.1

# Bleaching Index
# miz <- read.csv("data/mizerek/338_2018_1702_MOESM1_ESM.csv", as.is=TRUE)
# miz <- miz[c("Revised.species.name", "BI")]
# names(miz) <- c("species", "BI")
# dat <- merge(dat, miz, all.x=TRUE)
# hist(miz$BI)
dat$BI <- runif(nrow(dat))

# Genus age
dat$Genus.fossil.age[is.na(dat$Genus.fossil.age)] <- mean(dat$Genus.fossil.age, na.rm=TRUE)
dat$genus_age <- dat$Genus.fossil.age / max(dat$Genus.fossil.age, na.rm=TRUE)
sum(is.na(dat$Genus.fossil.age))
#hist(dat$genus_age)

# 3D
##############################
s <- 20
dat3D <- voronoiFilter3D(dat, s)

png("figs/figure1.png", width = 5, height = 5.5, units = 'in', res = 300)
plot(y ~ x, dat, col="grey", cex=1, xlab="PC1", ylab="PC2")
points(y ~ x, dat3D, col="blue", pch=20, cex=0.6)
legend("topleft", pch=c(0, 20), col=c("grey", "blue"), bty="n", legend=c("all species", "widely spread species"), cex=0.6)
dev.off()


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



