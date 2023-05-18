
library(raster)
library(rgdal)
library(rgeos)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Density")
d1 <- raster("Densisty_all_1.tif")
d2 <- raster("Densisty_all_2.tif")
d3 <- raster("Densisty_all_3.tif")
d4 <- raster("Densisty_all_4.tif")
d5 <- raster("Densisty_all_5.tif")
coordinates(d2) == coordinates(d5)

plot(d1)
length(values(d1))
nrow(coordinates(d1))

# Scale
X_mean <- mean(coordinates(d1)[,1])
X_sd <- sd(coordinates(d1)[,1])
X_sc <- (coordinates(d1)[,1] - X_mean) / X_sd

Y_mean <- mean(coordinates(d1)[,2])
Y_sd <- sd(coordinates(d1)[,2])
Y_sc <- (coordinates(d1)[,2] - Y_mean) / Y_sd

df <- data.frame(X_sc = X_sc, Y_sc = Y_sc, 
                 y1 = values(d1), y2 = values(d2), y3 = values(d3),
                 y4 = values(d4), y5 = values(d5))
df2 <- df[complete.cases(df),c(1:2)]
coordinates(df2) <- df2[,c(1:2)]
plot(df2)


save(df, file = "Density_scaledCoords.RData")

load("Density_scaledCoords.RData")

#####

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol")
core <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol/corepol.shp", 
        layer = "corepol")
plot(core)

pol_sc <- list()
centroid_sc <- list()

for (i in 1:4){
  
  X2_sc <- (core@polygons[[1]]@Polygons[[i]]@coords[,1] - X_mean)/X_sd
  Y2_sc <- (core@polygons[[1]]@Polygons[[i]]@coords[,2] - Y_mean)/Y_sd
  
  coord_sc <- cbind(X2_sc, Y2_sc)
  
  Pl <- Polygon(coord_sc)
  ID <- i
  Pls <- Polygons(list(Pl), ID=ID)
  SPls <- SpatialPolygons(list(Pls))
  df <- data.frame(value=1, row.names=ID)
  SPDF <- SpatialPolygonsDataFrame(SPls, df)
  pol_sc[[i]] <- SPDF
  centroid_sc[[i]] <- gCentroid(SPDF,byid=TRUE)
}

pols_sc <- do.call(bind, pol_sc) 
centroids_sc <- do.call(bind, centroid_sc)

# Initial location estimated by Deon model
init <- data.frame(X_sc = 0.37, Y_sc = 0.18)
coordinates(init) <- init[1,1:2]

plot(df2, pch = 21, cex = 0.2, axes = TRUE)
plot(pols_sc, add = TRUE)
points(centroids_sc, col = "red", pch = 19)
points(init, col = "blue", pch = 19)



plot(pol_sc[[1]])
points(centroid_sc[[1]])

coordinates(centroid_sc[[1]])
