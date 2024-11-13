
extract.rast.multip <- function (ss, rast, mult = 1, cov.name = "val.1", func = median) {
  
  for (i in 1:length(ss)) {
    tmpS <- ss[[i]][, c("X", "Y")] * mult
    id <- factor(1:nrow(tmpS))
    tmpR <- cbind(tmpS, id)
    aa <- rasterFromXYZ(tmpR, crs = projection(rast))
    bb <- rasterToPolygons(aa)
    
    for (j in 1:length(rasters)){
      
      r1 <- extract(rast[[j]], bb, fun = func)
      ss[[i]][, cov.name[[j]]] <- r1 
              }
  }
  return(ss)
}
