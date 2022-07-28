

## -------------------------------------------------
##             VARIABLES for PRE-ANALYSIS
## ------------------------------------------------- 

rm(list = ls())

library(raster)
library(terra)
library(rgdal)


## -------------------------------------------------
##                 Forest variable
##              CORINE land use cover 2018
## ------------------------------------------------- 


setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Corine_Landcover_2018_100m")
cor <- raster("Corine_WGS84_clip.tif")

lcor <- layerize(cor, classes=c(23,24,25), bylayer=TRUE, suffix='numbers')

# The layer 23 is broad-leaved forest
# The layer 24 is coniferous forest
# The layer 25 is mixed forest

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Corine_Landcover_2018_100m")

l1 <- raster(lcor, layer = 1)
#writeRaster(l1, filename = 'l1', format = 'GTiff')

l2 <- raster(lcor, layer = 2)
#writeRaster(l2, filename = 'l2', format = 'GTiff')

l3 <- raster(lcor, layer = 3)
#writeRaster(l3, filename = 'l3', format = 'GTiff')

l<- l1+l2+l3
# writeRaster(l, filename = 'forest', format = 'GTiff')


# Join with forest category from Andorra

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Corine_Landcover_2018_100m")
cor <- raster("forest_WGS8431N.tif") # Procesed forest layer
values(cor)[is.na(values(cor))] <- 0

and <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Corine_Landcover_2018_100m/Andorra/Habitats_2012_6", layer = "habitats2013_nivell4_WGS8431N")
and <- and[which(and@data$Txt_niv4 == "Boscos"), ]
and$Txt_niv4 <- as.factor(and$Txt_niv4)

r1 <- raster(resolution=res(cor), extent(cor))
crs(r1) <- "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs"
ras.and <- rasterize(and, r1, field = "Txt_niv4")
values(ras.and)[is.na(values(ras.and))] <- 0

all <- mosaic(cor,ras.and, fun = max)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Corine_Landcover_2018_100m")
writeRaster(all, filename = 'forestALL_WGS8431N', format = 'GTiff')

## -------------------------------------------------
##                 Roads variable
## ------------------------------------------------- 

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads")
roads <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = "allroads_clip_WGS84_31N")

## ---- Subset ----
# Categories:
# roads1-> Asphalted major roads + highway links: code starting by 511, 513
# roads2 -> Asphalted minor roads: 512+5141+5143
# roads3 -> Tracks: 514, except 5141, 5143, 5147 (track_grade 5, hardly visible, hay que esta categoría y con la 5146)
# roads4 -> Combines roads2 and roads3
# roads5 -> Paths: 515 + 5147 (track_grade5)

roads1 <- roads[grep(pattern = '511|513', roads$code), ]
roads2 <- roads[grep(pattern = '512|5141|5143', roads$code), ]

roads3 <- roads[grep(pattern = '514', roads$code), ]
roads3 <- roads3[-which(roads3$code %in% c("5147","5143","5141")), ]

roads4 <- rbind(roads2,roads3)
roads5 <- roads[grep(pattern = '515|5147', roads$code), ]

#writeOGR(roads1, dsn = "D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = 'roads1', driver = "ESRI Shapefile")
#writeOGR(roads2, dsn = "D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = 'roads2', driver = "ESRI Shapefile")
#writeOGR(roads3, dsn = "D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = 'roads3', driver = "ESRI Shapefile")
#writeOGR(roads4, dsn = "D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = 'roads4', driver = "ESRI Shapefile")
#writeOGR(roads5, dsn = "D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = 'roads5', driver = "ESRI Shapefile")

roads <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = "allroads_clip_WGS84_31N")

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads")
roads1 <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/OpenStreetMap_roads", layer = "roads1")

# Create roads raster from dem
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/EU_DEM")
dem <- raster("dem_WGS84_31N_Clip.tif")
rs <- dem 
rs.aggregate <- aggregate(rs, fact=40)
res(rs.aggregate)
rs[] <- 1:ncell(rs)

plot(rs)

# Intersect lines with raster "polygons" and add length to new lines segments
rsp <- rasterToPolygons(rs)
rp <- intersect(roads, rsp)
rp$length <- gLength(rp, byid=TRUE) / 1000
x <- tapply(rp$length, rp$layer, sum)
r <- raster(rs)
r[as.integer(names(x))] <- x

head(roads1@data)


## -------------------------------------------------
##                 Altitude
## ------------------------------------------------- 

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/EU_DEM")
dem <- raster("dem_WGS84_31N_Clip.tif")

## ---- Calculate slope and roughness ----

# Since the resolution of the state space is quite low (5x5km), and I calculate
# a moving window at the bear hr size, average altitude or slope might show
# high variability, and roughness might be a better indicator of density at this scale

slope <- terrain(dem,'slope', unit = 'degrees', neighbors = 8, filename = 'slope')
writeRaster(slope, filename = 'clip_slope', format = 'GTiff')

roughness <- terrain(dem, 'roughness', neighbors = 8, filename = 'roughness')
writeRaster(roughness, filename = 'clip_roughness', format ='GTiff')


## -------------------------------------------------
##  Moving window of habitat variables (home range size bear)
## ------------------------------------------------- 

## ---- Calculate radius circular home range (A=pi*r2) ----

# According to maelis the area was 200km2 in a preliminary analysis
r <- sqrt(200/pi) # 8km = 8000 m

## ---- Load layers ----

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Corine_Landcover_2018_100m")
cor <- raster("forestALL_WGS8431N.tif")

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/EU_DEM")
dem <- raster("dem_WGS84_31N_Clip.tif")
slope <- raster("clip_slope.tif")
rough <- raster("clip_roughness.tif")

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe")
distcore <- raster("dist_corekernel.tif")
logDistcore <- raster("dist_corekernel_log.tif")

dens1km <- raster("dens1km.tif")
dens200m <- raster("dens200m.tif")
dens200m_preST <- raster("dens200m_preStudy1016.tif")



## ---- Function ----

hrscale <- function(data = data) {
  
  factor <- 5000/res(data)[1] # Determine the factor to use in aggregate

  dat5x5 <- terra::aggregate(data, fact = factor)
  
  xy <- xyFromCell(dat5x5, 1:ncell(dat5x5)) # Takes the center coordinates of the pixels
  
  ## Extract mean values from buffer (radious = radious hr size)
  
  values <-extract(dat5x5, xy, method = 'simple',
                   buffer = 8000,
                   small = TRUE, fun = mean, na.rm = TRUE,
                   df = TRUE, factors = TRUE, sp = TRUE)
  
  dat2 <- dat5x5 # New raster to store
  values(dat2) <- values[,2] # Just substitutes the values in the same order
  
  return(dat2)
  
}

## ---- Forest ----

cor2 <- hrscale(data = cor)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
writeRaster(cor2, filename = "forest_hrbear.tif")

## ---- Altitude ----

dem2 <- hrscale(data = dem)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
writeRaster(dem2, filename = "dem_hrbear.tif")

## ---- Slope ----

slope2 <- hrscale(data = slope)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
writeRaster(slope2, filename = "slope_hrbear.tif")

## ---- Roughness ----

rough2 <- hrscale(data = rough)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
writeRaster(rough2, filename = "rough_hrbear.tif")

## ---- Distance to core ----

core2 <- hrscale(data = distcore)
logcore2 <- hrscale(data = logDistcore)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
writeRaster(core2, filename = "distcore_hrbear.tif")
writeRaster(logcore2, filename = "logDistcore_hrbear.tif")

## ---- Density of observations ----

obsDens1km <- hrscale(data = dens1km)
obsDens200m <- hrscale(data = dens200m)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
writeRaster(obsDens1km, filename = "obsDens1km_hrbear.tif")
writeRaster(obsDens200m, filename = "obsDens200m_hrbear.tif")

# Al final al rescalar esta variable al tamaño del hr size me parece que da igual 1km que 200m

##`Density of observations during the period of systematic sampling before our study period

obsDens200m_preST <- hrscale(data = dens200m_preST)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
writeRaster(obsDens200m_preST, filename = "obsDens200m_preST_hrbear.tif")

plot(dens1km)
plot(dens200m)
plot(dens200m_preST)


