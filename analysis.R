
# -- mike wd -- #
#setwd("/Users/mikemcwilliam/Documents/PostDoc/species_choice")

library("ggplot2")
library("cowplot")
source("data_prep.R")


# Anyalsis
library(tripack)
source("R/functions.R")
library("FD")

# points selection
#n <- 100 # total species
#s <- 20  # species selected
#dat <- data.frame(x=rnorm(n), y=rnorm(n))
#plot(y ~ x, dat, col="grey")
# random
#points(y ~ x, dat[sample(1:n, s, replace=FALSE),], col="red", pch=20, cex=0.6)
# most evenly spread
#points(y ~ x, voronoiFilter(dat, s), col="blue", pch=20, cex=0.6)


# select traits
dat<-read.csv("data/data.csv")
dat$Abundance.GBR<-factor(dat$Abundance.GBR, levels=c("rare","uncommon", "common"))
colnames(dat)

cats<-c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")



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



































# 3D interactive plots

library("car")
library("rgl")

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
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
           
           
           