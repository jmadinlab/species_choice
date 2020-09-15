# Species choice

library("ggplot2")
library("cowplot")
library("FD")
library("vegan")
library("tripack")
source("data_prep.R")
library("ks")
library("car")
library("rgl")

source("R/functions.R")

# points
##############################
dat<-read.csv("data/data.csv")
dat$Abundance.GBR <- factor(dat$Abundance.GBR, levels=c("rare","uncommon", "common"))

cats<-c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pcoa <- pcoa(gowdis(dat[,cats]))

space <- data.frame(x=jitter(pcoa$vectors[,1], amount=0.01), y=jitter(pcoa$vectors[,2], amount=0.01), z=jitter(pcoa$vectors[,3], amount=0.01))
dat <- cbind(dat, space)

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


seg<-cbind(dat3D[,c("x","y","z")], data.frame(xend=rep(0,nrow(dat3D)),yend=rep(0,nrow(dat3D)),zend=rep(0,nrow(dat3D))))
head(seg)
  #writeWebGL()
plot3d(dat$x, dat$y, dat$z, type="s", xlab="PC1", ylab="PC2", zlab="PC3", size=0.5, box=FALSE, axes=F, col.panel = "black",)
spheres3d(dat3D$x, dat3D$y, dat3D$z,  col="red", radius=0.01)
axes3d(c("x", "y", "z"), col="white")
segments3d(x=as.vector(t(seg[,c("x","xend")])),y=as.vector(t(seg[,c("y","yend")])), z=as.vector(t(seg[,c("z","zend")])), col='red', alpha=0.5)
           bg3d("slategrey")

# 3D -weighted
##############################
dat3Dd <- voronoiFilter3DDi(dat, s)

png("figs/figure3.png", width = 5, height = 5.5, units = 'in', res = 300)

plot(y ~ x, dat, col="grey", cex=1, xlab="PC1", ylab="PC2")
points(y ~ x, dat3D, col="blue", pch=20, cex=0.6)
 points(y ~ x, dat3Dd, col="red", pch=5, cex=1.2)
legend("topleft", pch=c(0, 20, 5), col=c("grey", "blue", "red"), bty="n", legend=c("species", "even spread alone", "even spread & quality"), cex=0.6)

dev.off()


seg<-cbind(dat3D[,c("x","y","z")], data.frame(xend=rep(0,nrow(dat3D_Di)),yend=rep(0,nrow(dat3D_Di)),zend=rep(0,nrow(dat3D_Di))))
head(seg)
#writeWebGL()
plot3d(dat$x, dat$y, dat$z, type="s", xlab="PC1", ylab="PC2", zlab="PC3", size=0.5, box=FALSE, axes=F, col.panel = "black",)
spheres3d(dat3D_Di$x, dat3D_Di$y, dat3D_Di$z,  col="red", radius=0.01)
axes3d(c("x", "y", "z"), col="white")
segments3d(x=as.vector(t(seg[,c("x","xend")])),y=as.vector(t(seg[,c("y","yend")])), z=as.vector(t(seg[,c("z","zend")])), col='red', alpha=0.5)
bg3d("slategrey")

# 3D -weighted
##############################
dat3Ddd <- voronoiFilter3DDi(dat, s)

png("figs/figure4.png", width = 5, height = 5.5, units = 'in', res = 300)

plot(y ~ x, dat, col="grey", cex=1, xlab="PC1", ylab="PC2")
points(y ~ x, dat3D, col="blue", pch=20, cex=0.6)
points(y ~ x, dat3Dd, col="red", pch=5, cex=1.2)
points(y ~ x, dat3Ddd, col="green", pch=3, cex=1.5)
legend("topleft", pch=c(0, 20, 5), col=c("grey", "blue", "red", "green"), bty="n", legend=c("species", "even spread alone", "even spread & quality", "even spread, quality & density"), cex=0.6)

dev.off()




# do species have a range of tolerances? 

ggplot()+
geom_histogram(data=dat, aes(x=BRI), bins=15, fill="white", col="grey")+
geom_point(data=dat3D_Di, aes(x=BRI, y=1) , fill="red", shape=25)+
#scale_x_log10()+
ggtitle("bleaching response index of selected species (red)")






















# mike old code



# select traits
dat<-read.csv("data/data.csv")
dat$Abundance.GBR<-factor(dat$Abundance.GBR, levels=c("rare","uncommon", "common"))
colnames(dat)

cats<-c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")

library(rgl)
pp <- scatterplot3d(dat, angle=55, phi=20)
pp$points3d(dat2D, col="blue", pch=20, cex=0.6)
pp$points3d(dat3D, col="red", pch=3)

plot3d( x = dat, type="s", radius=0.1, col=rgb(0,0,0,0.1))
spheres3d(x = dat2D, col = "red", radius = 0.2, alpha=0.3)
spheres3d(x = dat3D, col = "blue", radius = 0.2, alpha=0.3)

# generate components

# prcomp
#pca<-prcomp(dat[,cats], center=T, scale=T)
#biplot(pca)

# princomp
prin<-princomp((dat[,cats]), cor = TRUE, scores = TRUE)
pc12<-prin$scores[,1:2]
ll<-prin$loadings
#space<-data.frame(x=jitter(prin$scores[,1], amount=0.25), y=jitter(prin$scores[,2], amount=0.25), z=jitter(prin$scores[,3], amount=0.25))

# pcoa
pcoa<-pcoa(gowdis(dat[,cats]))
space<-data.frame(x=jitter(pcoa$vectors[,1], amount=0.01), y=jitter(pcoa$vectors[,2], amount=0.01), z=jitter(pcoa$vectors[,3], amount=0.01))














# trait vectors 



# princomp
prin<-princomp((dat[,cats]), cor = TRUE, scores = TRUE)
pc12<-prin$scores[,1:2]
ll<-prin$loadings
#space<-data.frame(x=jitter(prin$scores[,1], amount=0.25), y=jitter(prin$scores[,2], amount=0.25), z=jitter(prin$scores[,3], amount=0.25))

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
plot( pc12[,], pch=16, cex=0.25, col="grey", add=FALSE, ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5),las=1) 
#plot(est, cont=seq(1,100,by=1), display="filled.contour2", add=FALSE, ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5),las=1) 
#plot(est,abs.cont=cl[1], labels=c(0.5),labcex=0.75, add=TRUE, lwd=0.75, col="grey30")
#plot(est,abs.cont=cl[2], labels=c(0.95),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
#plot(est,abs.cont=cl[3], labels=c(0.99),labcex=0.75, add=TRUE, lwd=0.5, col="grey60")
#points( pc12[,], pch=16, cex=0.25, col="black") 
plot(fit, cex=0.90, col=1, labels=list(vectors = c("GR","CW","MCS","SD","CH","SAV","R")))
segments(0,0, fit2[,1], fit2[,2], col=1, lty=2, lwd=1)
mtext("PC1", cex=0.75, side=1, line=0.5, adj=1)
mtext("PC2", cex=0.75, side=2, line=0.5, at=3.5) #, las=2)







# 2D spread....
points<-cbind(space, abun=dat$Abundance.GBR, species=dat$species)

use<-points[,c("x","y")]
#use<-points[points$abun=="common",c("x","y")] # only common species

n <- nrow(use) # total species
s <- 40  # species selected
spread<-voronoiFilter(use, s)
spread$species<-points$species[as.numeric(rownames(spread))]

ggplot()+
geom_point(data=points, aes(x,y), shape=21, col="grey", size=2)+
geom_point(data=spread, aes(x,y))+
geom_text(data=spread, aes(x,y, label=species), size=2, position=position_nudge(x=diff(range(spread$x))*0.07, y=diff(range(spread$y))*0.01))+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),legend.title=element_blank(), legend.position=c(0.85,0.9), legend.key.size=unit(2, "mm"))+
annotation_custom(ggplotGrob(
ggplot()+
geom_segment(data=spread, aes(x=x,y=y,xend=0,yend=0),size=0.25)+
geom_point(data=points, aes(x,y), shape=21, col="grey", size=0.15)+
geom_point(data=spread, aes(x,y), size=0.5)+theme_void()), 
xmin = 3, xmax = 5, ymin = 2, ymax = 3)+
labs(x="principal component 1", y="principal component 2")














# Find nearest neighbours in each node...
library("FNN")
library("reshape2")

n=7
#gower<-gowdis(dat[,cats])
knn<-get.knn(points[,c("x","y")],k=n)$nn.index
near<-NULL
for(x in 1:n){
col<-dat$species[knn[,x]]
near<-cbind(near, x=col)}
colnames(near)<-paste("n", seq(1:n), sep="")
near<-as.data.frame(near)
near$species <- dat$species
near<-near[near$species %in% spread$species,]
near$group<-c(1:nrow(near))
near<-melt(near, id.var=c("group", "species"), value.name="near")
near[,c("x1","y1")]<-points[match(near$species, points$species),c("x","y")]
near[,c("x2","y2")]<-points[match(near$near, points$species), c("x","y")]
head(near)

# remove duplicates
near<-near[!near$near %in% near$species,]
near<-near[!duplicated(near$near),]

ggplot()+
geom_point(data=points, aes(x,y), shape=21, col="grey", size=2)+
geom_point(data=spread, aes(x,y))+
geom_segment(data=near, aes(x=x1, y=y1, xend=x2, yend=y2), size=0.1)+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),legend.title=element_blank(), legend.position=c(0.85,0.9), legend.key.size=unit(2, "mm"))+
annotation_custom(ggplotGrob(
ggplot()+
geom_segment(data=spread, aes(x=x,y=y,xend=0,yend=0),size=0.25)+
geom_point(data=points, aes(x,y), shape=21, col="grey", size=0.15)+
geom_point(data=spread, aes(x,y), size=0.5)+theme_void()), 
xmin = 3, xmax = 5, ymin = 2, ymax = 3)+
labs(x="principal component 1", y="principal component 2")











# Desirability index (out of 100)
###########################
# 50 points for common (categorical - ctb)
# 20 points for large size (ranked)
# 15 points for large range (ranked)
# 15 points for old (top 75% quartile)
# restorability ?
colnames(dat)

# abundance 
dat$D.abun<-ifelse(dat$Abundance.GBR=="common", 50, ifelse(dat$Abundance.GBR=="uncommon", 10, 0))

# colony size
dat$D.size<-dat$cat_colonydiameter/max(dat$cat_colonydiameter)*20

# range size
#qR<-quantile(dat$Range.size,0.25)
#dat$D.rang<-ifelse(dat$Range.size>qR, 15, 0)
#hist(dat$Range.size)
#abline(v=qR)
dat$D.range<-dat$Range.size/max(dat$Range.size)*15

# genus age 
#hist(log10(dat$Species.age.phylogeny))
#hist(log10(dat$Genus.fossil.age))
#qA<-quantile(dat$Genus.fossil.age,0.25, na.rm=T)
#dat$D.age<-ifelse(dat$Genus.fossil.age>qA, 15, 0)
dat$D.age<-dat$Genus.fossil.age/max(dat$Genus.fossil.age, na.rm=T)*10
dat$D.age[is.na(dat$D.age)]<-15 # benefit of doubt
#Species.age.phylogeny
dat$D.age2<-dat$Species.age.phylogeny/max(dat$Species.age.phylogeny, na.rm=T)*5

dat$D.index<-dat$D.abun+dat$D.size+dat$D.rang+dat$D.age+dat$D.age2


# Most desirable in each group

des<-unique(melt(near[,c("group","species","near")], id.var="group", value.name="species"))
des<-des[order(as.numeric(des$group)),]
des$D.index<-dat$D.index[match(des$species, dat$species)]

d.agg<-aggregate(D.index~group, des, max)
d.max<-merge(d.agg, des)
d.max[order(as.numeric(d.max$group)),]
nrow(d.max)

points$D.index<-dat$D.index[match(points$species, dat$species)]
select<-points[points$species %in% d.max$species,]

ggplot()+
geom_point(data=points, aes(x,y, size=D.index), shape=21, col="grey")+
geom_point(data=select, aes(x,y))+
geom_text(data=select, aes(x,y, label=species), size=2, position=position_nudge(x=diff(range(spread$x))*0.07, y=diff(range(spread$y))*0.01))+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),legend.title=element_blank(), legend.position=c(0.85,0.9), legend.key.size=unit(2, "mm"))+
labs(x="principal component 1", y="principal component 2")+
scale_radius(range=c(0,3))

























# self thinning... 
# find closest pairs, remove lowest D.index
head(select)
rownames(select)<-select$species
knn2<-as.matrix(dist(select[,c("x","y")]))
knn2[upper.tri(knn2, diag = TRUE)]<-NA
knn2<-melt(knn2)
knn2<-knn2[!is.na(knn2$value),]

pairs<-knn2[order(knn2$value),][c(1:15),]
pairs$D1<-dat$D.index[match(pairs$Var1, dat$species)]
pairs$D2<-dat$D.index[match(pairs$Var2, dat$species)]
pairs

min(c(pairs$D1, pairs$D2))
































points$cols<-as.factor(ifelse(points$species %in% spread$species, "red","white"))
points$labs<-ifelse(points$cols=="red",points$species, "")
sub<-points[points$cols=="red",]
seg<-cbind(sub, data.frame(xend=rep(0,nrow(sub)),yend=rep(0,nrow(sub)),zend=rep(0,nrow(sub))))
head(points)

  #writeWebGL()
plot3d(points$x, points$y, points$z, group=points$group, type="s",   
xlab="PC1", ylab="PC2", zlab="PC3", col=points$cols, size=0.5, box=FALSE, axes=F, col.panel = "black",)
axes3d(c("x", "y", "z"), col="white", ntick=3)
segments3d(x=as.vector(t(seg[,c("x","xend")])),y=as.vector(t(seg[,c("y","yend")])), z=as.vector(t(seg[,c("z","zend")])), col='red', alpha=0.1)
           bg3d("slategrey") 
  text3d(points$x, points$y, points$z, texts=points$labs, cex=0.5,col="red", adj=c(1.1, 1.1))





   


library("plotly")
           
     
     fig <- plot_ly(points, x = ~x, y = ~y, z = ~z, 
     color = ~cols, colors = c('red', 'grey'), 
     marker = list(symbol = 'circle'),
     size = 0.25,
     showgrid = FALSE,
     hovertemplate = ~paste(species))
    # text = ~paste(species))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(
       xaxis = list(title = 'PC1', showspikes=FALSE, showgrid = FALSE),
       yaxis = list(title = 'PC2', showspikes=FALSE, showgrid = FALSE),
       zaxis = list(title = 'PC3', showspikes=FALSE, showgrid = FALSE), 
                     paper_bgcolor = 'grey',
                     plot_bgcolor = "grey"))
     
fig
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
           
           
           