# Species choice

# LIBRARIES
library("vegan")
library("tripack")
library("ks")
library("spatstat")
library("sp")
library("psych")

# FUNCTIONS
source("R/functions.R")

# DATA
source("R/data_prep.R")
dat<-read.csv("data/data.csv")

# TRAIT SPACE
cats<-c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
# biplot(pca, col=c("grey", "red"), cex=0.6)

space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))

dat2 <- cbind(dat, space)
dim(dat2)
write.csv(dat2["species"], "output/species.csv", row.names = FALSE)

# CORRELATIONS AMONG VARS
png("figs/fig_0.png", width = 5, height = 5, units = 'in', res = 300)
pairs.panels(dat2[c("abund", "range", "bleach", "restore", "PC1", "PC2")], hist.col="lightgrey")
dev.off()

# Density in trait space
##############################
# H12 <- Hpi(x=space[,1:2])      # optimal bandwidth estimation
# H13 <- Hpi(x=space[,c(1, 3)])      # optimal bandwidth estimation
# H23 <- Hpi(x=space[,2:3])      # optimal bandwidth estimation
# est12 <- kde(x=space[,1:2], H=H12, eval.points=space[,1:2], compute.cont=TRUE)    # kernel density estimation
# est13 <- kde(x=space[,c(1, 3)], H=H13, eval.points=space[,c(1, 3)], compute.cont=TRUE)
# est23 <- kde(x=space[,2:3], H=H23, eval.points=space[,2:3], compute.cont=TRUE)
# 
# dat2$density <- apply(cbind(est12$estimate, est13$estimate, est23$estimate), 1, mean)
# dat2$density <- dat2$density / max(dat2$density, na.rm=TRUE)

# trait vectors 
##############################
fit <- envfit(dat2[c("PC1", "PC2")], dat2[,cats]) # use envfit for drawing arrows
fit2 <- fit$vectors$arrows*-3 # drawing line segments opposite arrows

s <- 20
dat2Dv <- dd2D(dat2, s, vars="", voronoi=TRUE, trait=TRUE)
dat2Dn <- dd2D(dat2, s, vars="", voronoi=FALSE, trait=TRUE)

png("figs/fig_1.png", width = 4, height = 7, units = 'in', res = 300)

par(mfrow=c(2, 1), mar=c(1.5, 1, 1, 1), oma=c(3, 3, 0, 0))
plot(PC2 ~ PC1, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
plot(fit, cex=0.90, col=1, labels=list(vectors = c("GR","CW","MCS","SD","CH","SAV","R")))
# segments(0,0, fit2[,1], fit2[,2], col=1, lty=2, lwd=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("A", line=-1, adj=0)

# TRAIT METHODS
plot(PC2 ~ PC1, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(PC2 ~ PC1, dat2Dv, col="black", pch=4, cex=1)
points(PC2 ~ PC1, dat2Dn, col="black", pch=1, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
legend("topright", pch=c(1, 4), legend=c("Nearest neighbor", "Voronoi area"), bty="n", col=c("black", "black"), cex=0.8)
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
dat2Drab <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=TRUE)
dat2Drabf <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=FALSE)
dat2Drabr <- dd2D(dat2, s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE)
dat2Drabro <- dd2D(dat2, s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE, oppo=TRUE)

png("figs/fig_2.png", width = 3, height = 8, units = 'in', res = 300)

par(mfrow=c(3, 1), mar=c(1.5, 1, 1, 1), oma=c(3, 3, 0, 0))

plot(PC2 ~ PC1, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(PC2 ~ PC1, dat2Drabf, col="black", pch=4, cex=1)
points(PC2 ~ PC1, dat2Drab, col="black", pch=1, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("A", line=-1, adj=0)

pts <- chull(dat2Drabf$PC1, dat2Drabf$PC2)
pts <- c(pts, pts[1])
polygon(dat2Drabf$PC1[pts], dat2Drabf$PC2[pts], col=NA, border="black", lty=2)

x1 <- dat2Drabf$PC1
y1 <- dat2Drabf$PC2
hpts <- chull(x = x1, y = y1)
hpts <- c(hpts, hpts[1])
xy.coords <- cbind(x1, y1)
chull.coords <- xy.coords[hpts,]
chull.poly <- Polygon(chull.coords, hole=F)
chull.poly@area

pts <- chull(dat2Drab$PC1, dat2Drab$PC2)
pts <- c(pts, pts[1])
polygon(dat2Drab$PC1[pts], dat2Drab$PC2[pts], col=NA, border="grey", lty=1)


legend("topright", pch=c(4, 1), legend=c("Resistance", "Resistance & function"), bty="n", col=c("black", "black"), cex=0.8)


plot(PC2 ~ PC1, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(PC2 ~ PC1, dat2Drabr, col="black", pch=3, cex=1)
points(PC2 ~ PC1, dat2Drab, col="black", pch=1, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("B", line=-1, adj=0)

legend("topright", pch=c(1, 3), legend=c("Resistance & function", "Resistance, function & restoration"), bty="n", col=c("black", "black"), cex=0.8)


plot(PC2 ~ PC1, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
points(PC2 ~ PC1, dat2Drabro, col="red", pch=1, cex=1)
points(PC2 ~ PC1, dat2Drabr, col="black", pch=3, cex=1)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("C", line=-1, adj=0)

legend("topright", pch=c(20, 1), legend=c("Resistance, function & restoration", "Vulnerability, function & restoration"), bty="n", col=c("black", "red"), cex=0.8)

mtext("PC1", 1, line=1, outer=TRUE)
mtext("PC2", 2, line=1, outer=TRUE)

dev.off()

####

###

Winners <- dat2Drabf$value
names(Winners) <- dat2Drabf$species_n
Winners.Function <- dat2Drab$value
names(Winners.Function) <- dat2Drab$species_n
Winners.Function.Restore <- dat2Drabr$value
names(Winners.Function.Restore) <- dat2Drabr$species_n

spp <- as.character(unique(c(names(Winners.Function.Restore), names(Winners.Function), names(Winners))))

mat <- data.frame(matrix(NA, length(spp), 3))
row.names(mat) <- spp
names(mat) <- c("Winners", "Winners.Function", "Winners.Function.Restore")
mat2 <- mat

mat["Winners"] <- spp %in% names(Winners)
mat2[spp %in% names(Winners),"Winners"] <- Winners

mat["Winners.Function"] <- spp %in% names(Winners.Function)
mat2[spp %in% names(Winners.Function),"Winners.Function"] <- Winners.Function

mat["Winners.Function.Restore"] <- spp %in% names(Winners.Function.Restore)
mat2[spp %in% names(Winners.Function.Restore),"Winners.Function.Restore"] <- Winners.Function.Restore

# mat2[is.na(mat2)] <- 0

mat2 <- mat2[order(mat2$Winners.Function.Restore, decreasing=FALSE),]
mat2 <- mat2[order(mat2$Winners.Function, decreasing=FALSE),]
mat2 <- mat2[order(mat2$Winners, decreasing=FALSE),]

png("figs/fig_3.png", width = 5, height = 12, units = 'in', res = 300)

par(mar=c(1, 8, 13, 8))
image(t(as.matrix(mat2)), axes=FALSE, col=heat.colors(20, rev = TRUE))

text(rep(0, nrow(mat2)), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,1]), 2)), cex=0.7)

text(rep(0.5, nrow(mat2)), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,2]), 2)), cex=0.7)
text(rep(1, nrow(mat2)), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,3]), 2)), cex=0.7)

axis(2, at=seq(0, 1, (1/(nrow(mat2)-1))), labels=spp, las=2, cex.axis=1, font=3)
axis(3, at=seq(0, 1, (1/(ncol(mat2)-1))), labels=c("Resistance", "Resistance & function", "Resistance, function & restoration"), cex.axis=1, las=2)
axis(4, at=seq(0, 1, (1/(nrow(mat2)-1)))[!is.na(mat2[["Winners.Function.Restore"]])], labels=spp[!is.na(mat2[["Winners.Function.Restore"]])], las=2, cex.axis=1, font=3)

dev.off()

#### OPPO

dat2Drao <- dd2D(dat2, s, vars=c("abund", "range"), voronoi=FALSE, trait=TRUE, oppo=TRUE)
dat2Drabfo <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=FALSE, oppo=TRUE)
dat2Drabo <- dd2D(dat2, s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=TRUE, oppo=TRUE)
dat2Drabro <- dd2D(dat2, s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE, oppo=TRUE)


Losers <- dat2Drabfo$value
names(Losers) <- dat2Drabfo$species_n

Losers.Function <- dat2Drabo$value
names(Losers.Function) <- dat2Drabo$species_n

Losers.Function.Restore <- dat2Drabro$value
names(Losers.Function.Restore) <- dat2Drabro$species_n

spp <- as.character(unique(c(names(Losers.Function.Restore), names(Losers.Function), names(Losers))))

mat <- data.frame(matrix(NA, length(spp), 3))
row.names(mat) <- spp
names(mat) <- c("Losers", "Losers.Function", "Losers.Function.Restore")
mat2 <- mat

mat["Losers"] <- spp %in% names(Losers)
mat2[spp %in% names(Losers),"Losers"] <- Losers

mat["Losers.Function"] <- spp %in% names(Losers.Function)
mat2[spp %in% names(Losers.Function),"Losers.Function"] <- Losers.Function

mat["Losers.Function.Restore"] <- spp %in% names(Losers.Function.Restore)
mat2[spp %in% names(Losers.Function.Restore),"Losers.Function.Restore"] <- Losers.Function.Restore

# mat2[is.na(mat2)] <- 0

mat2 <- mat2[order(mat2$Losers.Function.Restore, decreasing=FALSE),]
mat2 <- mat2[order(mat2$Losers.Function, decreasing=FALSE),]
mat2 <- mat2[order(mat2$Losers, decreasing=FALSE),]


png("figs/fig_4.png", width = 5, height = 12, units = 'in', res = 300)

par(mar=c(1, 8, 13, 8))
image(t(as.matrix(mat2)), axes=FALSE, col=heat.colors(20, rev = TRUE))

text(rep(0, nrow(mat2)), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,1]), 2)), cex=0.7)
text(rep(0.5, nrow(mat2)), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,2]), 2)), cex=0.7)
text(rep(1, nrow(mat2)), seq(0, 1, (1/(nrow(mat2)-1))), (round(as.matrix(mat2[,3]), 2)), cex=0.7)

axis(2, at=seq(0, 1, (1/(nrow(mat2)-1))), labels=spp, las=2, cex.axis=1, font=3)
axis(3, at=seq(0, 1, (1/(ncol(mat2)-1))), labels=c("Vulnerability", "Vulnerability & function", "Vulnerability, function & restoration"), cex.axis=1, las=2)
axis(4, at=seq(0, 1, (1/(nrow(mat2)-1)))[!is.na(mat2[["Losers.Function.Restore"]])], labels=spp[!is.na(mat2[["Losers.Function.Restore"]])], las=2, cex.axis=1, font=3)

dev.off()

### NUMBERS
library(raster)
xmn <- -4
xmx <- 4
ymn <- -3
ymx <- 3

spp <- sort(as.character(dat2$species_n))
store <- data.frame(spp=spp)
ns <- seq(361, 3, -1)
# ns <- c(100, 50, 25, 15, 10, 7, 5, 4, 3)
store_hull <- data.frame() 

# png("figs/fig_5.png", width = 7, height = 8, units = 'in', res = 300)
# par(mfrow=c(4, 3), mar=c(0,0,0,0), oma=c(0,0,0,0))

for (s in ns) {
  ch_f <- dd2D(dat2, s=s, vars=c(), voronoi=FALSE, trait=TRUE)
  ch_e <- dd2D(dat2, s=s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=FALSE)
  ch_fe <- dd2D(dat2, s=s, vars=c("abund", "range", "bleach"), voronoi=FALSE, trait=TRUE)
  ch_fer <- dd2D(dat2, s=s, vars=c("abund", "range", "bleach", "restore"), voronoi=FALSE, trait=TRUE)
  
  x1 <- ch_f$PC1
  y1 <- ch_f$PC2
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln5_f <- length(unique(cells))
  cellm5_f <- mean(table(cells))
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln10_f <- length(unique(cells))
  cellm10_f <- mean(table(cells))
  pts_f <- chull(x = x1, y = y1)
  pts_f <- c(pts_f, pts_f[1])
  xy.coords <- cbind(x1, y1)
  chull.coords <- xy.coords[pts_f,]
  chull.poly <- Polygon(chull.coords, hole=F)
  area_f <- chull.poly@area
  spread_f <- mean(nndist(x1, y1))

  x1 <- ch_e$PC1
  y1 <- ch_e$PC2
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln5_e <- length(unique(cells))
  cellm5_e <- mean(table(cells))
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln10_e <- length(unique(cells))
  cellm10_e <- mean(table(cells))
  pts_e <- chull(x = x1, y = y1)
  pts_e <- c(pts_e, pts_e[1])
  xy.coords <- cbind(x1, y1)
  chull.coords <- xy.coords[pts_e,]
  chull.poly <- Polygon(chull.coords, hole=F)
  area_e <- chull.poly@area
  spread_e <- mean(nndist(x1, y1))
  
  x1 <- ch_fe$PC1
  y1 <- ch_fe$PC2
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln5_fe <- length(unique(cells))
  cellm5_fe <- mean(table(cells))
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln10_fe <- length(unique(cells))
  cellm10_fe <- mean(table(cells))
  pts_fe <- chull(x = x1, y = y1)
  pts_fe <- c(pts_fe, pts_fe[1])
  xy.coords <- cbind(x1, y1)
  chull.coords <- xy.coords[pts_fe,]
  chull.poly <- Polygon(chull.coords, hole=F)
  area_fe <- chull.poly@area
  spread_fe <- mean(nndist(x1, y1))

  x1 <- ch_fer$PC1
  y1 <- ch_fer$PC2
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln5_fer <- length(unique(cells))
  cellm5_fer <- mean(table(cells))
  r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)
  cells <- cellFromXY(r, cbind(x1, y1))
  celln10_fer <- length(unique(cells))
  cellm10_fer <- mean(table(cells))
  pts_fer <- chull(x = x1, y = y1)
  pts_fer <- c(pts_fer, pts_fer[1])
  xy.coords <- cbind(x1, y1)
  chull.coords <- xy.coords[pts_fer,]
  chull.poly <- Polygon(chull.coords, hole=F)
  area_fer <- chull.poly@area
  spread_fer <- mean(nndist(x1, y1))
  
  area_cli <- NA
  spread_cli <- NA
  celln5_cli <- NA
  cellm5_cli <- NA
  celln10_cli <- NA
  cellm10_cli <- NA
  if (s == 7) {
    ch_cli <- dat2[dat2$clipperton==1,]
    x1 <- ch_cli$PC1
    y1 <- ch_cli$PC2
    r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
    cells <- cellFromXY(r, cbind(x1, y1))
    celln5_cli <- length(unique(cells))
    cellm5_cli <- mean(table(cells))
    r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)
    cells <- cellFromXY(r, cbind(x1, y1))
    celln10_cli <- length(unique(cells))
    cellm10_cli <- mean(table(cells))
    pts_cli <- chull(x = x1, y = y1)
    pts_cli <- c(pts_cli, pts_cli[1])
    xy.coords <- cbind(x1, y1)
    chull.coords <- xy.coords[pts_cli,]
    chull.poly <- Polygon(chull.coords, hole=F)
    area_cli <- chull.poly@area
    spread_cli <- mean(nndist(x1, y1))
  }

  # if (s %in% c(100, 75, 50, 25, 20, 15, 12, 10, 7, 6, 4, 2)) {
  #   plot(PC2 ~ PC1, dat2, pch=20, cex=1, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-2.178684, 2.541263), xlim=c(-3.116689, 3.894455), axes=FALSE) 
  #   # points(PC2 ~ PC1, ch1, col="black", pch=4, cex=1)
  #   polygon(ch_f$PC1[pts_f], ch_f$PC2[pts_f], col=NA, border=rgb(0,0,1,0.5), lwd=3)
  #   polygon(ch_e$PC1[pts_e], ch_e$PC2[pts_e], col=NA, border=rgb(1,0,0,0.5), lwd=3)
  #   polygon(ch_fe$PC1[pts_fe], ch_fe$PC2[pts_fe], col=NA, border=rgb(0,1,0,0.5), lwd=3)
  #   polygon(ch_fer$PC1[pts_fer], ch_fer$PC2[pts_fer], col=NA, border=rgb(1,1,0,0.5), lwd=3)
  #   if (s == 7) {
  #     polygon(ch_cli$PC1[pts_cli], ch_cli$PC2[pts_cli], col=NA, border=rgb(0,0,0,0.5), lwd=3)
  #     arrows(0, 0, -3, 2.5, length=0.1)
  #     text(-3, 2.7, "Cementors", cex=0.8)
  #     arrows(0, 0, 4, -0.5, length=0.1)
  #     text(3.95, -0.2, "Fillers", cex=0.8)
  #     arrows(0, 0, -1.5, -2, length=0.1)
  #     text(-1.5, -2.2, "Builders", cex=0.8)
  #   }
  #   mtext(paste0("n=", s), line=-2, adj=0, cex=2)
  # }
  
  rarea_store <- c()
  rspread_store <- c()
  rcelln5_store <- c()
  rcelln10_store <- c()
  rcellm5_store <- c()
  rcellm10_store <- c()
  for (i in 1:1000) {
    samp <- dat2[c("PC1", "PC2")][sample(1:nrow(dat2), size=s, replace=F),]

    x1 <- samp$PC1
    y1 <- samp$PC2
    r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
    cells5 <- cellFromXY(r, cbind(x1, y1))
    r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)
    cells10 <- cellFromXY(r, cbind(x1, y1))
    pts <- chull(x = x1, y = y1)
    pts <- c(pts, pts[1])
    xy.coords <- cbind(x1, y1)
    chull.coords <- xy.coords[pts,]
    chull.poly <- Polygon(chull.coords, hole=F)
    rarea_store <- c(rarea_store, chull.poly@area)
    rspread_store <- c(rspread_store, mean(nndist(x1, y1)))
    rcelln5_store <- c(rcelln5_store, length(unique(cells5)))
    rcelln10_store <- c(rcelln10_store, length(unique(cells10)))
    rcellm5_store <- c(rcellm5_store, mean(table(cells5)))
    rcellm10_store <- c(rcellm10_store, mean(table(cells10)))
  }
  rarea_store <- sort(rarea_store)
  rspread_store <- sort(rspread_store)
  rcelln5_store <- sort(rcelln5_store)
  rcelln10_store <- sort(rcelln10_store)
  rcellm5_store <- sort(rcellm5_store)
  rcellm10_store <- sort(rcellm10_store)
  
  store_hull <- rbind(store_hull, data.frame(n=s, 
                                             area_f=area_f, area_e=area_e, area_fe=area_fe, area_fer=area_fer, area_cli=area_cli, area_rand=rarea_store[500], area_rand_lo=rarea_store[25], area_rand_hi=rarea_store[975], 
                                             spread_f=spread_f, spread_e=spread_e, spread_fe=spread_fe, spread_fer=spread_fer, spread_cli=spread_cli, spread_rand=rspread_store[500], spread_rand_lo=rspread_store[50], spread_rand_hi=rspread_store[950], 
                                             celln5_f=celln5_f, celln5_e=celln5_e, celln5_fe=celln5_fe, celln5_fer=celln5_fer, celln5_cli=celln5_cli, celln5_rand=rcelln5_store[500], celln5_rand_lo=rcelln5_store[50], celln5_rand_hi=rcelln5_store[950],
                                             cellm5_f=cellm5_f, cellm5_e=cellm5_e, cellm5_fe=cellm5_fe, cellm5_fer=cellm5_fer, cellm5_cli=cellm5_cli, cellm5_rand=rcellm5_store[500], cellm5_rand_lo=rcellm5_store[50], cellm5_rand_hi=rcellm5_store[950],
                                             celln10_f=celln10_f, celln10_e=celln10_e, celln10_fe=celln10_fe, celln10_fer=celln10_fer, celln10_cli=celln10_cli, celln10_rand=rcelln10_store[500], celln10_rand_lo=rcelln10_store[50], celln10_rand_hi=rcelln10_store[950],
                                             cellm10_f=cellm10_f, cellm10_e=cellm10_e, cellm10_fe=cellm10_fe, cellm10_fer=cellm10_fer, cellm10_cli=cellm10_cli, cellm10_rand=rcellm10_store[500], cellm10_rand_lo=rcellm10_store[50], cellm10_rand_hi=rcellm10_store[950]
  ))

}
# dev.off()

write.csv(store_hull, "output/occupancy.csv", row.names=FALSE)

## PLOTS

n5_max <- max(store_hull$celln5_f)
m5_max <- max(store_hull$cellm5_f)
n10_max <- max(store_hull$celln10_f)
m10_max <- max(store_hull$cellm10_f)

png("figs/fig_6a.png", width = 4, height = 6, units = 'in', res = 300)

par(mfrow=c(2, 1), oma=c(5, 2, 0, 0), mar=c(2, 2, 2, 1))

plot(sqrt(store_hull$n), (store_hull$celln5_rand/n5_max), type="l", xlab="", ylab="", axes=FALSE, lty=2, xlim=c(0, 10), ylim=c(0, 1))
polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c((store_hull$celln5_rand_lo/n5_max), rev((store_hull$celln5_rand_hi/n5_max))), border=NA, col=rgb(0,0,0,0.2))

mtext("A. 5x5 occupancy grid", 3, 0, adj=0)

axis(2, las=2)
axis(1, at=1:10, labels=NA, las=2)

lines(sqrt(store_hull$n), (store_hull$celln5_e/n5_max), col="red")
lines(sqrt(store_hull$n), (store_hull$celln5_f/n5_max), col="black")
lines(sqrt(store_hull$n), (store_hull$celln5_fe/n5_max), col="blue")
lines(sqrt(store_hull$n), (store_hull$celln5_fer/n5_max), col="green")

points(sqrt(7), (store_hull$celln5_cli[!is.na(store_hull$celln5_cli)]/n5_max), col="black", pch=8)



## 10

plot(sqrt(store_hull$n), (store_hull$celln10_rand/n10_max), type="l", xlab="", ylab="", axes=FALSE, lty=2, xlim=c(0, 10), ylim=c(0, 1))
polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c((store_hull$celln10_rand_lo/n10_max), rev((store_hull$celln10_rand_hi/n10_max))), border=NA, col=rgb(0,0,0,0.2))

mtext("B. 10x10 occupancy grid", 3, 0, adj=0)

axis(2, las=2)
axis(1, at=1:10, labels=c(1:10)^2, las=2)

lines(sqrt(store_hull$n), (store_hull$celln10_e/n10_max), col="red")
lines(sqrt(store_hull$n), (store_hull$celln10_f/n10_max), col="black")
lines(sqrt(store_hull$n), (store_hull$celln10_fe/n10_max), col="blue")
lines(sqrt(store_hull$n), (store_hull$celln10_fer/n10_max), col="green")

points(sqrt(7), (store_hull$celln10_cli[!is.na(store_hull$celln10_cli)]/n10_max), col="black", pch=8)

legend("bottomright", legend=c("Random", "Least vulnerable only", "Function & least vulnerable", "Function, least vulnerable & restore", "Function only", "Clipperton Atoll"), lty=c(2, 1, 1, 1, 1, NA), pch=c(NA, NA, NA, NA, NA, 8), col=c("black", "red", "blue", "green", "black", "black"), bty="n", cex=0.5, xpd=TRUE)

mtext("Number of species", 1, 1, outer=TRUE)
mtext("Proportion occupancy", 2, 1, outer=TRUE)


# REDUNDANCY
# 
# plot(sqrt(store_hull$n), (store_hull$cellm5_rand), type="l", xlab="", ylab="Mean redundancy", axes=FALSE, lty=2, xlim=c(0, 10), ylim=c(1, 6))
# polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c((store_hull$cellm5_rand_lo), rev((store_hull$cellm5_rand_hi))), border=NA, col=rgb(0,0,0,0.2))
# 
# mtext("C. 5x5", 3, 0, adj=0)
# 
# axis(2, las=2)
# axis(1, at=1:10, labels=c(1:10)^2, las=2)
# 
# lines(sqrt(store_hull$n), (store_hull$cellm5_e), col="red")
# lines(sqrt(store_hull$n), (store_hull$cellm5_f), col="black")
# lines(sqrt(store_hull$n), (store_hull$cellm5_fe), col="blue")
# lines(sqrt(store_hull$n), (store_hull$cellm5_fer), col="green")
# 
# points(sqrt(7), (store_hull$cellm5_cli[!is.na(store_hull$cellm5_cli)]), col="black", pch=8)
# 
# # 10
# 
# plot(sqrt(store_hull$n), (store_hull$cellm10_rand), type="l", xlab="", ylab="", axes=FALSE, lty=2, xlim=c(0, 10), ylim=c(1, 6))
# polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c((store_hull$cellm10_rand_lo), rev((store_hull$cellm10_rand_hi))), border=NA, col=rgb(0,0,0,0.2))
# 
# mtext("D. 10x10", 3, 0, adj=0)
# 
# axis(2, las=2)
# axis(1, at=1:10, labels=c(1:10)^2, las=2)
# 
# lines(sqrt(store_hull$n), (store_hull$cellm10_e), col="red")
# lines(sqrt(store_hull$n), (store_hull$cellm10_f), col="black")
# lines(sqrt(store_hull$n), (store_hull$cellm10_fe), col="blue")
# lines(sqrt(store_hull$n), (store_hull$cellm10_fer), col="green")
# 
# points(sqrt(7), (store_hull$cellm10_cli[!is.na(store_hull$cellm10_cli)]), col="black", pch=8)


dev.off()
### BARS

thresh <- 0.8

n5 <- c(
store_hull$n[which.min(abs((store_hull$celln5_rand / n5_max) - thresh))],
store_hull$n[which.min(abs((store_hull$celln5_e / n5_max) - thresh))],
store_hull$n[which.min(abs((store_hull$celln5_fe / n5_max) - thresh))],
store_hull$n[which.min(abs((store_hull$celln5_fer / n5_max) - thresh))],
store_hull$n[which.min(abs((store_hull$celln5_f / n5_max) - thresh))]
)

n10 <- c(
  store_hull$n[which.min(abs((store_hull$celln10_rand / n10_max) - thresh))],
  store_hull$n[which.min(abs((store_hull$celln10_e / n10_max) - thresh))],
  store_hull$n[which.min(abs((store_hull$celln10_fe / n10_max) - thresh))],
  store_hull$n[which.min(abs((store_hull$celln10_fer / n10_max) - thresh))],
  store_hull$n[which.min(abs((store_hull$celln10_f / n10_max) - thresh))]
)

png("figs/fig_6_n5.png", width = 5, height = 4, units = 'in', res = 300)
  bp <- barplot(sqrt(n5), col=c("grey","red", "blue", "green", "black"), lty=c(1,1,1,1,2), border=NA, las=2, ylim=c(0, 15), axes=FALSE)
  arrows(bp[1], sqrt(store_hull$n[which.min(abs((store_hull$celln5_rand_lo / n5_max) - thresh))]), bp[1], sqrt(store_hull$n[which.min(abs((store_hull$celln5_rand_hi / n5_max) - thresh))]), angle=90, code=3)
  axis(2, at=c(0, 2, 4, 6, 8, 10, 12, 14), labels=c(0, 2, 4, 6, 8, 10, 12, 14)^2, las=2, cex.axis=1.7)
  mtext("Species at 0.8 occupancy", 3, 0, cex=2)
dev.off()

png("figs/fig_6_n10.png", width = 5, height = 4, units = 'in', res = 300)
  barplot(sqrt(n10), col=c("grey","red", "blue", "green", "black"), lty=c(1,1,1,1,2), border=NA, las=2, ylim=c(0, 15), axes=FALSE)
  arrows(bp[1], sqrt(store_hull$n[which.min(abs((store_hull$celln10_rand_lo / n10_max) - thresh))]), bp[1], sqrt(store_hull$n[which.min(abs((store_hull$celln10_rand_hi / n10_max) - thresh))]), angle=90, code=3)
  axis(2, at=c(0, 2, 4, 6, 8, 10, 12, 14), labels=c(0, 2, 4, 6, 8, 10, 12, 14)^2, las=2, cex.axis=1.7)
  mtext("Species at 0.8 occupancy", 3, 0, cex=2)
dev.off()

### REDUN

thresh <- 2

m5 <- c(
  store_hull$n[which.min(abs((store_hull$cellm5_rand) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm5_e) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm5_fe) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm5_fer) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm5_f) - thresh))]
)

m10 <- c(
  store_hull$n[which.min(abs((store_hull$cellm10_rand) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm10_e) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm10_fe) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm10_fer) - thresh))],
  store_hull$n[which.min(abs((store_hull$cellm10_f) - thresh))]
)

png("figs/fig_6_m5.png", width = 5, height = 4, units = 'in', res = 300)
  barplot(m5, col=c("grey","red", "blue", "green", "black"), lty=c(1,1,1,1,2), border=NA, las=2, cex.axis=1.5)
  mtext("Species at mean redundancy of 2", 3, 1, cex=1.5)
dev.off()

png("figs/fig_6_m10.png", width = 5, height = 4, units = 'in', res = 300)
  barplot(m10, col=c("grey","red", "blue", "green", "black"), lty=c(1,1,1,1,2), border=NA, las=2, cex.axis=1.5)
dev.off()

### SUBS

png("figs/fig_6_sub5.png", width = 3.5, height = 4, units = 'in', res = 300)

plot(PC2 ~ PC1, dat2, pch=20, cex=0.8, ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
x1 <- dat2$PC1
y1 <- dat2$PC2
r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 5, ncol = 5)
xx <- seq(xmn, xmx, length=6)[1:5]
yy <- seq(ymx, ymn, length=6)[1:5]
xxd <- diff(xx)[1]/2
yyd <- diff(yy)[1]/2
xx <- xx + xxd
yy <- yy + yyd
xxv <- rep(xx, length(yy))
yyv <- rep(yy, each=length(xx))
cells <- unique(cellFromXY(r, cbind(x1, y1)))
rect(xxv[cells] - xxd, yyv[cells] - yyd, xxv[cells] + xxd, yyv[cells] + yyd, border="grey")

dev.off()

png("figs/fig_6_sub10.png", width = 3.5, height = 4, units = 'in', res = 300)

plot(PC2 ~ PC1, dat2, pch=20, cex=0.8, ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE) 
x1 <- dat2$PC1
y1 <- dat2$PC2
r <- raster(extent(xmn, xmx, ymn, ymx), nrow = 10, ncol = 10)
xx <- seq(xmn, xmx, length=11)[1:10]
yy <- seq(ymx, ymn, length=11)[1:10]
xxd <- diff(xx)[1]/2
yyd <- diff(yy)[1]/2
xx <- xx + xxd
yy <- yy + yyd
xxv <- rep(xx, length(yy))
yyv <- rep(yy, each=length(xx))
cells <- unique(cellFromXY(r, cbind(x1, y1)))
rect(xxv[cells] - xxd, yyv[cells] - yyd, xxv[cells] + xxd, yyv[cells] + yyd, border="grey")

dev.off()








png("figs/fig_6b.png", width = 5, height = 5, units = 'in', res = 300)

plot(sqrt(store_hull$n), (store_hull$cellm_rand), type="l", xlab="Number of species", ylab="Functional redundancy", axes=FALSE, lty=2, xlim=c(1, 10))
polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c((store_hull$cellm_rand_lo), rev((store_hull$cellm_rand_hi))), border=NA, col=rgb(0,0,0,0.2))

axis(2, las=2)
axis(1, at=1:20, labels=c(1:20)^2)

lines(sqrt(store_hull$n), (store_hull$cellm_e), col="red")
lines(sqrt(store_hull$n), (store_hull$cellm_f), col="black")
lines(sqrt(store_hull$n), (store_hull$cellm_fe), col="blue")
lines(sqrt(store_hull$n), (store_hull$cellm_fer), col="green")

points(sqrt(7), (store_hull$cellm_cli[!is.na(store_hull$area_cli)]), col="black", pch=8)

legend("topleft", legend=c("Random", "Function only", "Least vulnerable only", "Function & least vulnerable", "Function, least vulnerable & restore", "Clipperton Atoll"), lty=c(2, 1, 1, 1, 1, NA), pch=c(NA, NA, NA, NA, NA, 8), col=c("black", "black", "red", "blue", "green", "black"), bty="n", cex=0.8)

dev.off()






dd <- diff(rev((store_hull$spread_fe * store_hull$n)))
nn <- 4:100
plot(nn, dd)
abline(h=0, lty=2)
spl <- predict(loess(dd ~ nn, span=0.5))
lines(spl, col="red")

dd <- diff(rev(store_hull$area_f / store_hull$area_rand))
nn <- 4:100
plot(nn, dd)
abline(h=0, lty=2)
spl <- predict(loess(dd ~ nn, span=0.5))
lines(spl, col="red")



# plot(sqrt(store_hull$n), (store_hull$newmet2)/max(store_hull$newmet2), type="l", ylim=c(0, 1), col="red")
# lines(sqrt(store_hull$n), (store_hull$newmet1)/max(store_hull$newmet2), col="blue")
# lines(sqrt(store_hull$n), (store_hull$newmet3)/max(store_hull$newmet2), col="green")
# lines(sqrt(store_hull$n), (store_hull$newmet5)/max(store_hull$newmet2), col="orange")
# lines(sqrt(store_hull$n), (store_hull$newm)/max(store_hull$newmet2), col="black")
# polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c((store_hull$newm_lo)/max(store_hull$newmet2), rev((store_hull$newm_hi)/max(store_hull$newmet2))), border=NA, col=rgb(0,0,0,0.2))


plot(sqrt(store_hull$n), (store_hull$spread_f / store_hull$spread_rand), type="l", col="black", axes=FALSE, xlab="Number of species", ylab="Relative to best functional spread")
axis(2, las=2)
axis(1, at=1:20, labels=c(1:20)^2, las=2)

plot(sqrt(store_hull$n), (store_hull$spread_area_f^2 / store_hull$spread_area_rand), type="l", col="black", axes=FALSE, xlab="Number of species", ylab="Relative to best functional spread")
axis(2, las=2)
axis(1, at=1:20, labels=c(1:20)^2, las=2)

plot(sqrt(store_hull$n), (store_hull$area_f / store_hull$area_rand), type="l", col="black", axes=FALSE, xlab="Number of species", ylab="Relative to best functional spread")
axis(2, las=2)
axis(1, at=1:20, labels=c(1:20)^2, las=2)




abline(h=0, lty=2)
polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c(log10(store_hull$newm_lo / store_hull$newmet2), rev(log10(store_hull$newm_hi / store_hull$newmet2))), border=NA, col=rgb(0,0,0,0.1))

lines(sqrt(store_hull$n), log10(store_hull$spread_area_fe / store_hull$spread_area_f), col="blue")
lines(sqrt(store_hull$n), log10(store_hull$spread_area_e / store_hull$spread_area_f), col="green")
lines(sqrt(store_hull$n), log10(store_hull$spread_area_fer / store_hull$spread_area_f), col="red")

legend("bottomright", legend=c("Least vulnerable only", "Least vulnerable & function", "Least vulnerable, function & restore"), lty=c(1, 1, 1), col=c("green", "blue", "red"), bty="n", cex=0.8)



##########



plot(sqrt(store_hull$n), (store_hull$newmet2), type="l", xlab="Number of species", ylab="Proportion trait space", col="red", axes=FALSE)
axis(2, las=2)
axis(1, at=1:10, labels=c(1:10)^2)

lines(sqrt(store_hull$n), (store_hull$rand)/max(store_hull$hullarea2), col="black")
polygon(c(sqrt(store_hull$n), rev(sqrt(store_hull$n))), c((store_hull$rand_lo)/max(store_hull$hullarea2), rev((store_hull$rand_hi)/max(store_hull$hullarea2))), border=NA, col=rgb(0,0,0,0.2))

lines(sqrt(store_hull$n), (store_hull$hullarea1)/max(store_hull$hullarea2), col="blue")


lines(sqrt(store_hull$n), (store_hull$hullarea3)/max(store_hull$hullarea2), col="green")

legend("bottomright", legend=c("Function random", "Function only", "Least vulnerable only", "Function & least vulnerable"), lty=c(1, 1, 1, 1), col=c("black", "red", "green","blue"), bty="n", cex=0.8)






store <- store[,-1]
store <- as.data.frame(sapply(store, as.numeric))
row.names(store) <- spp
store <- store[!apply(store, 1, function(x) all(is.na(x))),]
store <- store[order(apply(store, 1, function(x) sum(!is.na(x)))),]

png("figs/fig_5.png", width = 5.5, height = 9, units = 'in', res = 300)

par(mar=c(1, 8, 4, 1))
image(t(as.matrix(store)), axes=FALSE, col=heat.colors(20))
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



