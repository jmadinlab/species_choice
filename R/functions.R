# Functions

# Original, area only
dd2D <- function(dat, s, vars, voronoi=FALSE, trait=TRUE, oppo=FALSE, se=FALSE) {

  if (se) {
    dat$PC1r <- dat$PC1 + rnorm(length(dat$PC1), 0, sd(dat$PC1) / sqrt(length(dat$PC1)))
    dat$PC2r <- dat$PC2 + rnorm(length(dat$PC2), 0, sd(dat$PC2) / sqrt(length(dat$PC2)))
  } else {
    dat$PC1r <- dat$PC1
    dat$PC2r <- dat$PC2
  }
  
  subset <- dat
  dropped <- vector()
  for (i in 1:(nrow(dat)-s)) {
    
    if (trait) {
      if (voronoi) {
        v <- voronoi.mosaic(x=subset[,'PC1r'], y=subset[,'PC2r'], duplicate='error')
        info <- cells(v)
        areas <- unlist(lapply(info, function(x) x$area))
        areas[is.na(areas)] <- mean(areas, na.rm=TRUE)
      } else {
        # areas <- apply(nndist(subset$xr, subset$yr, k=1:3), 1, mean)
        areas <- nndist(subset$PC1r, subset$PC2r, k=1)
        areas <- areas / max(areas)
      }
    } else {
      areas <- rep(1, length(subset$PC2r))
    }
    
    if ("range" %in% vars) {
      if (oppo) {
        areas <- areas * (1 - subset$range)
      } else {
        areas <- areas * subset$range
      }
    }
    if ("lat_poleward" %in% vars) {
      if (oppo) {
        areas <- areas * (1 - subset$lat_poleward)
      } else {
        areas <- areas * subset$lat_poleward
      }
    }
    if ("abund" %in% vars) {
      if (oppo) {
        areas <- areas * (1 - subset$abund)
      } else {
        areas <- areas * subset$abund
      }
    }
    if ("age" %in% vars) {
      if (oppo) {
        areas <- areas * (1 - subset$age)
      } else {
        areas <- areas * subset$age
      }
    }
    if ("bleach" %in% vars) {
      if (oppo) {
        areas <- areas * (1 - subset$bleach)
      } else {
        areas <- areas * (subset$bleach)
      }
    }
    if ("restore" %in% vars) {
      areas <- areas * (subset$restore)
    }
    if ("density" %in% vars) {
      areas <- areas * subset$density
    }
    
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'PC1r'], dat[,'PC2r'], sep='_') == paste(subset[smallest,'PC1r'], subset[smallest,'PC2r'], sep='_')))
    subset <- subset[-smallest,]
    areas <- areas[-smallest]
  }
  temp2 <- dat[-dropped,]
  temp2$value <- areas
  return(temp2)
}

# 2D and building in other factors.
voronoiFilter2D.w.d <- function(dat, s) {
  subset <- dat
  dropped <- vector()
  for (i in 1:(nrow(dat)-s)) {
    v <- voronoi.mosaic(x=subset[,'x'], y=subset[,'y'], duplicate='error')
    info <- cells(v)
    areas <- unlist(lapply(info, function(x) x$area))
    
    # Incorporating relative abundance (b/w 0 and 1) simply by multiplying with area
    areas <- areas * subset$Range.size
    areas <- areas * subset$abund
    #areas <- areas * subset$BI
    areas <- areas * subset$genus_age
    areas <- areas * subset$density
    
    # Find smallest value, no longer just the voronoi area
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}



# 3D 
# 
# voronoiFilter3D <- function(dat, s) {
#   subset <- dat
#   dropped <- vector()
#   for (i in 1:(nrow(dat)-s)) {
#     v1 <- voronoi.mosaic(x=subset[,'x'], y=subset[,'y'], duplicate='error')
#     info1 <- cells(v1)
#     areas1 <- unlist(lapply(info1, function(x) x$area))
# 
#     v2 <- voronoi.mosaic(x=subset[,'x'], y=subset[,'z'], duplicate='error')
#     info2 <- cells(v2)
#     areas2 <- unlist(lapply(info2, function(x) x$area))
#     
#     v3 <- voronoi.mosaic(x=subset[,'y'], y=subset[,'z'], duplicate='error')
#     info3 <- cells(v3)
#     areas3 <- unlist(lapply(info3, function(x) x$area))
#     
#     areas <- apply(cbind(areas1, areas2, areas3), 1, min, na.rm=TRUE)
#     
#     smallest <- which(areas == min(areas, na.rm=TRUE))[1]
#     dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], dat[,'z'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], subset[smallest,'z'], sep='_')))
#     subset <- subset[-smallest,]
#   }
#   return(dat[-dropped,])
# }
# 


# 3D weighted
voronoiFilter3D <- function(dat, s, vars="") {
  subset <- dat2
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
    
    if ("range" %in% vars) {
      areas <- areas * subset$range
    }
    if ("abund" %in% vars) {
      areas <- areas * subset$abund
    }
    if ("age" %in% vars) {
      areas <- areas * subset$age
    }
    if ("bi" %in% vars) {
      areas <- areas * subset$bi
    }
    if ("density" %in% vars) {
      areas <- areas * subset$density
    }
    
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], dat[,'z'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], subset[smallest,'z'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}

# 
# 
# # 3D weighted with density
# 
# voronoiFilter3D.w.d <- function(dat, s) {
#   subset <- dat
#   dropped <- vector()
#   for (i in 1:(nrow(dat)-s)) {
#     v1 <- voronoi.mosaic(x=subset[,'x'], y=subset[,'y'], duplicate='error')
#     info1 <- cells(v1)
#     areas1 <- unlist(lapply(info1, function(x) x$area))
# 
#     v2 <- voronoi.mosaic(x=subset[,'x'], y=subset[,'z'], duplicate='error')
#     info2 <- cells(v2)
#     areas2 <- unlist(lapply(info2, function(x) x$area))
#     
#     v3 <- voronoi.mosaic(x=subset[,'y'], y=subset[,'z'], duplicate='error')
#     info3 <- cells(v3)
#     areas3 <- unlist(lapply(info3, function(x) x$area))
#     
#     areas <- apply(cbind(areas1, areas2, areas3), 1, min, na.rm=TRUE)
#     
#     # Incorporating relative abundance (b/w 0 and 1) simply by multiplying with area
#     areas <- areas * subset$range
#     areas <- areas * subset$abund
#     #areas <- areas * subset$BI
#     areas <- areas * subset$age
#     areas <- areas * 1/subset$trait_dens
#     
#     smallest <- which(areas == min(areas, na.rm=TRUE))[1]
#     dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], dat[,'z'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], subset[smallest,'z'], sep='_')))
#     subset <- subset[-smallest,]
#   }
#   return(dat[-dropped,])
# }
# 


