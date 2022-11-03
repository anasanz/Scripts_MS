

# MAPS FOR BEAR WORKSHOP

library(sf)
library(raster)


map <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2.shp") %>%
  st_transform(map, crs = 32631) # WGS84 31N
plot(map$geometry)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
forest <- raster("forest_hrbear.tif")

# Study area
setwd("D:/MargSalas/Oso/SCR/Workshop Pyrenees")
pdf(file = "study_area.pdf")
plot(forest)
plot(map$geometry, add = TRUE)
dev.off()

# Zoon study area

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_final_1719")
load("tdf2019.RData")
X <- tdf2019[,c(2,3)] # Extent of the study area

x <- rep(seq(min(X[,1]), max(X[,1]), by = 15000), 5)
y <- rep(seq(min(X[,2]), max(X[,2]), by = 15000), each = 16)
X2 <- cbind(x,y) #Generate example traps

setwd("D:/MargSalas/Oso/SCR/Workshop Pyrenees")
pdf(file = "1.study_area_zoom.pdf")
plot(forest, xlim = c(min(X[,1]), max(X[,1])), ylim = c(min(X[,2]), max(X[,2])))
plot(map$geometry, add = TRUE)
points(X2, pch = 4)
dev.off()

# Zoom study area detections

setwd("D:/MargSalas/Oso/SCR/Workshop Pyrenees")
pdf(file = "2.study_area_zoom_det1.pdf")
plot(forest, xlim = c(min(X[,1]), max(X[,1])), ylim = c(min(X[,2]), max(X[,2])))
plot(map$geometry, add = TRUE)
points(X2, pch = 4)
points(X2[c(11,12,27),] , pch = 7, col = "red")
dev.off()

setwd("D:/MargSalas/Oso/SCR/Workshop Pyrenees")
pdf(file = "3.study_area_zoom_det1.pdf")
plot(forest, xlim = c(min(X[,1]), max(X[,1])), ylim = c(min(X[,2]), max(X[,2])))
plot(map$geometry, add = TRUE)
points(X2, pch = 4)
points(X2[c(26,42),] , pch = 7, col = "red")
dev.off()

# Zoom study area detections + buffer

xmin <- min(X[,1]) - 25000
ymin <- min(X[,2]) - 25000
xmax <- max(X[,1]) + 25000
ymax <- max(X[,2]) + 25000

setwd("D:/MargSalas/Oso/SCR/Workshop Pyrenees")
pdf(file = "4.study_area_zoom_det1.pdf")
plot(forest, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
plot(map$geometry, add = TRUE)
points(X2, pch = 4)
dev.off()


