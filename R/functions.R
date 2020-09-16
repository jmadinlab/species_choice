# Functions

# Original, area only
voronoiFilter2D <- function(dat, s) {
  subset <- dat
  dropped <- vector()
  for (i in 1:(nrow(dat)-s)) {
    v <- voronoi.mosaic(x=subset[,'x'], y=subset[,'y'], duplicate='error')
    info <- cells(v)
    areas <- unlist(lapply(info, function(x) x$area))
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
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
    
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], dat[,'z'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], subset[smallest,'z'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}



# 3D weighted
voronoiFilter3D.w <- function(dat, s) {
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
    
    # Incorporating relative abundance (b/w 0 and 1) simply by multiplying with area
    areas <- areas * subset$Range.size
    areas <- areas * subset$abund
    #areas <- areas * subset$BI
    areas <- areas * subset$genus_age
    
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], dat[,'z'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], subset[smallest,'z'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}



# 3D weighted with density

voronoiFilter3D.w.d <- function(dat, s) {
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
    
    # Incorporating relative abundance (b/w 0 and 1) simply by multiplying with area
    areas <- areas * subset$Range.size
    areas <- areas * subset$abund
    #areas <- areas * subset$BI
    areas <- areas * subset$genus_age
    areas <- areas * subset$density
    
    smallest <- which(areas == min(areas, na.rm=TRUE))[1]
    dropped <- c(dropped, which(paste(dat[,'x'], dat[,'y'], dat[,'z'], sep='_') == paste(subset[smallest,'x'], subset[smallest,'y'], subset[smallest,'z'], sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}

















































