
rm(list = ls())

library(raster)
library(rgeos)
library(rgdal)
library(viridis)
library(dichromat)
library(ggplot2)
library(ggmap)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

# Load political map study area
sa <- readOGR("D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2_WGS84_31N_all")
eur <- readOGR("D:/MargSalas/Oso/Datos/GIS/Countries", "esp_fr_2")
eur <- spTransform(eur, crs(sa))

sa_t <- spTransform(sa, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

# Load State space
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf@proj4string <- CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs")
Xbuf_t <- spTransform(Xbuf, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
e <- extent(Xbuf)

# Load traps

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

load("tdf2017_effort.RData")
load("tdf2018_effort.RData")
load("tdf2019_effort.RData")
load("tdf2020_effort.RData")
load("tdf2021_effort.RData")

tdf_all <- rbind(tdf2017[,1:3], tdf2018[,1:3], tdf2019[,1:3],
                 tdf2020[,1:3], tdf2021[,1:3]) # Join to define state space
rownames(tdf_all) <- 1:nrow(tdf_all)
X <- tdf_all[,c(2,3)]
colnames(X) <- c('x', 'y')
coordinates(X) <- X[,c(1:2)]
X@proj4string <- CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs")
X_t <- spTransform(X, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
X_t@data$x <- coordinates(X_t)[,1]
X_t@data$y <- coordinates(X_t)[,2]


# Load satellite map study area
myMap <- get_map(location = bbox(sa_t), source = "stamen", maptype = "terrain", crop = FALSE) # Basemap

# PLOT STUDY AREA

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots")
#pdf("ss_sat.pdf",7,7)
p1 <- ggmap(myMap) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") + 
  geom_polygon(data = Xbuf_t, aes(x = long, y = lat), colour = "black", fill = "transparent", size = 0.8)
p1
#dev.off()

# PLOT STUDY AREA + TRAPS

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots")
#pdf("ss_sat_traps.pdf",7,7)

p2 <- ggmap(myMap) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") + 
  geom_polygon(data = Xbuf_t, aes(x = long, y = lat), colour = "black", fill = "transparent", size = 0.8) +
  geom_point(data = X_t@data, aes(x = x , y = y), colour = adjustcolor("darkred", alpha = 0.5), size = 0.5)
p2

#dev.off()


#----   GET A RASTER FOR EACH  HABITAT ---- 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
distcore <- raster("logDistcore_hrbear.tif")
habitat.r.distcore <- crop(distcore, e) # Crop it to extent of state-space
distcoreMask <- rasterize(Xbuf, habitat.r.distcore, mask = TRUE)

forest <- raster("forest_hrbear.tif")
habitat.r.forest <- crop(forest, e) # Crop it to extent of state-space
forestMask <- rasterize(Xbuf, habitat.r.forest, mask = TRUE)

dem <- raster("dem_hrbear.tif")
habitat.r.dem <- crop(dem, e) # Crop it to extent of state-space
demMask <- rasterize(Xbuf, habitat.r.dem, mask = TRUE)


setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Capas_David")
rug <- raster("RUGOSITAT_REPR_PYR_3.tif")
rug <- projectRaster(rug, crs = crs(distcore), method="bilinear")
habitat.r.rug <- crop(rug, e) # Crop it to extent of state-space
rugMask <- rasterize(Xbuf, habitat.r.rug, mask = TRUE)
rugMask <- resample(rugMask, distcoreMask, method = 'bilinear')




setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots")
#pdf("ss_traps.pdf",7,7)

plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2))
points(X, col = adjustcolor("darkred", alpha.f = 0.2), pch = 19)
plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2), add = TRUE)

#dev.off()


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots")
#pdf("ss_distCore.pdf",7,7)

plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2))
plot(distcoreMask, legend = FALSE, col = viridis(50), add = TRUE)
plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2), add = TRUE)

dev.off()


#plot(distcoreMask, legend = FALSE, col = dichromat(viridis(50, direction = -1), type = "tritan"), add = TRUE)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots")
pdf("ss_roguhness.pdf",7,7)

plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2))
plot(distcoreMask, legend = FALSE, col = viridis(50), add = TRUE)
plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2), add = TRUE)

dev.off()


#---- PLOT SEX PROPORTION ---- 

Xbuf2 <- readOGR("Buffer_8500_traps.shp") # Load sampling area (where we estimate abundance)


plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2))

