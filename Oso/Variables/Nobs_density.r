
## -------------------------------------------------
##                  Nobs <- DENSITY
## ------------------------------------------------- 
# Covariate that represents the density of observations of all kinds in space

rm(list = ls())

library(tidyverse)
library(sf)
library(sp)
library(lubridate)
library(raster)

## ---- 1. Clean observations ----

map <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2.shp")

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)  %>% 
  filter(Confirmed_Individual != "Indetermined") %>%
  mutate(date = as.Date(Date, format = "%d/%m/%Y"),
         an = year(date),
         month = month(date)) %>%
  st_as_sf(coords = c("x_long","y_lat"), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Sample the locations from Radiotracking to reduce their weight

os_track <- os %>%
  filter(Method == "Radiotracking")

# Radiotracking database reduced (for each individual, 1 location per month and year)
os_track_reduced <- os_track %>% 
  group_by(Confirmed_Individual,month,an) %>%
  sample_n(1)

os <- os %>% # Database without radiotracking
  filter(Method != "Radiotracking")

#Join
os_final <- rbind(os,os_track_reduced)

# Remove the observations that are repeated in two different methods (e.g., hair and camera)
# - By removing the observations that are duplicated for an individual within the same day and coord
os_dup <- os_final[,c(2,22,8,9)] 
os_final2 <- os_final[-which(duplicated(os_dup)),]

# - Next step could be removing the observations within a buffer of 10 m within a day?

## ---- 2. Store in raster ----

# Function to rasterize point density
pointcount <- function(r, pts){
  # make a raster of zeroes like the input
  r2 <- r
  r2[] <- 0
  # get the cell index for each point and make a table:
  counts <- table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] <- counts
  return(r2)
}

#### A. Data whole study period (1996-2021)

# Transform coordinate system os to fit
os_final2 <- os_final2 %>% st_transform(crs = "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs")
os_final2 <- as(os_final2, "Spatial")

# Create empty raster with extent study area rasterize polygon core area
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/EU_DEM")
dem <- raster("dem_WGS84_31N_Clip.tif")
rs <- dem 
rs.aggregate1km <- aggregate(rs, fact = 40) # 1000m resolution
values(rs.aggregate1km) <- NA # Raster to calculate 

rs.aggregate200m <- aggregate(rs, fact = 8) # 200m resolution, similar to forest
values(rs.aggregate200m) <- NA # Raster to calculate 

# Calculate number of points per pixel at 2 resolutions
dens1km <- pointcount(rs.aggregate1km, os_final2)
dens200m <- pointcount(rs.aggregate200m, os_final2)

# Check and visualize

plot(dens1km, ext = bbox(os_final2))

#e <- drawExtent()

plot(dens1km, ext = e)
points(os_final2)

plot(dens200m, ext = e)
points(os_final2)

# Me parece que es mejor el de 200m de resolución y luego escalar

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe")
writeRaster(dens1km, filename = 'dens1km', format = 'GTiff')
writeRaster(dens200m, filename = 'dens200m', format = 'GTiff')

#### B. Data systematic sampling before our study period (2010-2016)

os_final2 <- os_final[-which(duplicated(os_dup)),]

# Transform coordinate system os to fit
os_final2 <- os_final2 %>% 
  st_transform(crs = "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs") %>% 
  filter (Year %in% c(2010:2016))
os_final2 <- as(os_final2, "Spatial")

# Create empty raster with extent study area rasterize polygon core area
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/EU_DEM")
dem <- raster("dem_WGS84_31N_Clip.tif")
rs <- dem 
rs.aggregate200m <- aggregate(rs, fact = 8) # 200m resolution, similar to forest
raster::values(rs.aggregate200m) <- NA # Raster to calculate 

# Calculate number of points per pixel 
dens200m <- pointcount(rs.aggregate200m, os_final2)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe")
writeRaster(dens200m, filename = 'dens200m_preStudy1016', format = 'GTiff')

