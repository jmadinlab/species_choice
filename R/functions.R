# Functions

voronoiFilter <- function(dat, s) {
  subset <- dat
  dropped <- vector()
  for (i in 1:(n-s)) {
    v <- voronoi.mosaic(x=subset[,'x'],y=subset[,'y'],duplicate='error')
    info <- cells(v)
    areas <- unlist(lapply(info,function(x) x$area))
    smallest <- which(areas == min(areas,na.rm=TRUE))
    dropped <- c(dropped,which(paste(dat[,'x'],dat[,'y'],sep='_') == paste(subset[smallest,'x'],subset[smallest,'y'],sep='_')))
    subset <- subset[-smallest,]
  }
  return(dat[-dropped,])
}





