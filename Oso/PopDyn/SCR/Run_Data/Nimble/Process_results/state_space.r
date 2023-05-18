## -------------------------------------------------
##             Save state space without buffer
##        AIM: Estimate abundance only in trap array
## ------------------------------------------------- 
rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)
library(raster)
library(rgeos)
library(oSCR)
library(terra)
library(sp)
library(dplyr)
library(parallel)
library(sf)
library(lubridate)



setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
#source('dbinomLocal_normalBear.R')
source('dbinomLocal_normalBear_rbinom2.R')

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

#---- 1. LOAD THE DETECTION DATA ---- 

load("edf1721.RData")

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

load("tdf2017_effort.RData")
load("tdf2018_effort.RData")
load("tdf2019_effort.RData")
load("tdf2020_effort.RData")
load("tdf2021_effort.RData")

tdf_all <- rbind(tdf2017[,1:3], tdf2018[,1:3], tdf2019[,1:3],
                 tdf2020[,1:3], tdf2021[,1:3]) # Join to define state space
rownames(tdf_all) <- 1:nrow(tdf_all)


# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Nere", "Goiat")), ]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
#----   2.1 GET TRAPS---- 
X <- tdf_all[,c(2,3)]
colnames(X) <- c('x', 'y')
J <- dim(X)[1]

#----   2.2 DEFINE STATE SPACE EXTENT ---- 
# State space coordinates
# Buffer: 25000 (used by Maelis, also ~3*sigma in pre-analysis where sig = 6640)
xmin <- min(X[,1]) - 25000
ymin <- min(X[,2]) - 25000
xmax <- max(X[,1]) + 25000
ymax <- max(X[,2]) + 25000
e <- as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons") # Extent of state space

#----   2.3 GET A RASTER FOR THE HABITAT ---- 
# USE A FOREST RASTER TO GET A BASIS FOR THE HABITAT RASTER
# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
distcore <- raster("logDistcore_hrbear.tif")
# Crop it to extent of state-space
habitat.r <- crop(distcore, e) 
plot(habitat.r)

#----   2.4 DEFINE THE BUFFER AREA, DECIDE WHICH AREA IS CONSIDERED NOT BUFFER ---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

pid <- sapply(slot(Xbuf, "polygons"), function(x) slot(x, "ID"))  # All this is only to convert to spdf and save it takes ages otherwise to make each time
p.df <- data.frame( ID=1:length(Xbuf), row.names = pid) 
Xbuf <- SpatialPolygonsDataFrame(Xbuf, p.df)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
writeOGR(Xbuf, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "Buffer_statespace", driver = "ESRI Shapefile")   

# Convex hull
hpts <- chull(coordinates(Xpoints))
hpts <- c(hpts, hpts[1])

# Smaller buffer
Xbuf2 <- gBuffer(Xpoints, width = 8500)

pid <- sapply(slot(Xbuf2, "polygons"), function(x) slot(x, "ID"))  # All this is only to convert to spdf and save it takes ages otherwise to make each time
p.df <- data.frame( ID=1:length(Xbuf2), row.names = pid) 
Xbuf2 <- SpatialPolygonsDataFrame(Xbuf2, p.df)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
writeOGR(Xbuf2, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "Buffer_8500_traps", driver = "ESRI Shapefile")      


# Buffer using all observations of all types
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)  %>% 
  filter(Confirmed_Individual != "Indetermined") %>%
  mutate(date = as.Date(Date, format = "%d/%m/%Y"),
         an = year(date),
         month = month(date))

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
      coordinates(os_final) <- os_final[,c(11,12)]
      proj4string(os_final) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
      os_final <- spTransform(os_final, CRS = CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs"))
      
      setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
      writeOGR(os_final, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "allobs_forbuffer", driver = "ESRI Shapefile")  
      
      
      #osbuf <- gBuffer(os_final, width = 8500)
      
      pid <- sapply(slot(osbuf, "polygons"), function(x) slot(x, "ID"))   # All this is only to convert to spdf and save it takes ages otherwise to make each time
      p.df <- data.frame( ID=1:length(osbuf), row.names = pid) 
      p.df <- data.frame( ID = "buf")
      osbuf <- SpatialPolygonsDataFrame(osbuf, p.df)

      setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
      writeOGR(osbuf, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "Buffer_8500_allobs", driver = "ESRI Shapefile")      

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
osbuf <- readOGR("Buffer_8500_allobs.shp")

# Load all to check

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
osbuf <- readOGR("Buffer_8500_allobs.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

# Plots 
plot(e)
plot(Xbuf, col = "lightblue", add = TRUE)
plot(full_buf, add = TRUE)
points(Xpoints)
lines(X[hpts, ], col = "red")
plot(Xbuf2, col = adjustcolor("pink", alpha = 0.5), add = TRUE)
points(os_final, pch = 19, col = "violet")
plot(osbuf, col = adjustcolor("green", alpha = 0.5), add = TRUE)


# It seems that the whole range is more reasonable (osbuf), but need to adjust it
# 1- Join it with the buffer of the traps (is similar but it will be give a more reasonable area)
#     -> Anyway the model knows that in the buffer of the traps there are no individuals, so it shouldn't change
# 2- Clip it with state space (outside the state space, the model assumes there can't be individuals)
# 3 - Remove super outliers


proj4string(Xbuf2) <- proj4string(osbuf)
full_buf <- rbind(osbuf,Xbuf2)
full_buf <- aggregate(full_buf, dissolve = T)
full_buf_cropped <- crop(full_buf, Xbuf)

pid <- sapply(slot(full_buf_cropped, "polygons"), function(x) slot(x, "ID"))   # All this is only to convert to spdf and save it takes ages otherwise to make each time
p.df <- data.frame( ID=1:length(full_buf_cropped), row.names = pid) 
p.df <- data.frame( ID = "buf")
full_buf_cropped <- SpatialPolygonsDataFrame(full_buf_cropped, p.df)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
writeOGR(full_buf_cropped, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "Buffer_8500_allobs_cropped", driver = "ESRI Shapefile")  


# Load the layer of the buffers without outliers (removed manually arcgis)
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
osbuf_noOut <- readOGR("allobs_forbuffer_noOutliers_dis.shp")
osbuf_noOut <- aggregate(osbuf_noOut, dissolve = T)

full_buf_noOut <- rbind(osbuf_noOut,Xbuf2) # Join with buffer traps
full_buf_noOut <- aggregate(full_buf_noOut, dissolve = T)
full_buf_noOut_cropped <- crop(full_buf_noOut, Xbuf)

pid <- sapply(slot(full_buf_noOut_cropped, "polygons"), function(x) slot(x, "ID"))   # All this is only to convert to spdf and save it takes ages otherwise to make each time
p.df <- data.frame( ID=1:length(full_buf_noOut_cropped), row.names = pid) 
p.df <- data.frame( ID = "buf")
full_buf_noOut_cropped <- SpatialPolygonsDataFrame(full_buf_noOut_cropped, p.df)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
writeOGR(full_buf_noOut_cropped, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "Buffer_8000_allobs_cropped_noOut", driver = "ESRI Shapefile")  


# Plots 
plot(e)
plot(Xbuf, col = "lightblue", add = TRUE)
plot(full_buf_noOut_cropped, col = adjustcolor("green", alpha = 0.5), add = TRUE)
plot(Xbuf2, col = "pink", add = TRUE)

plot(full_buf_cropped, add = TRUE)
points(Xpoints)
plot(Xbuf2, col = adjustcolor("pink", alpha = 0.5), add = TRUE)
plot(osbuf, col = adjustcolor("green", alpha = 0.5), add = TRUE)

###############################################################P
# Aproximation of what Maelis could have done to estimate abundance (buffer size)

# Join opportunistic observations to trap array

# 1. Get depredations from 2010 to 2020

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)  %>% 
  filter(Confirmed_Individual != "Indetermined") %>%
  filter(Year == 2010 | Year == 2011 | Year == 2012 | Year == 2013 | Year == 2014 | Year == 2015 | Year == 2016 | Year == 2017 | Year == 2018 | Year == 2019 | Year == 2020)
  
os <- os %>%
  filter(Obs_type == "Predation" | Remarks == "1 oveja muerta" | Remarks == "Livestock predation. 2 sheep"
         | Remarks == " Livestock predation (1 sheep). Info completada con Taula General" 
         | Remarks == "De la depredació del dia 23-24/06/2012 a Campaus"
         | Remarks == "Livestock predation (Method is \"prospection suite à dégât\")" 
         | Remarks == "De l'atac de 14-07"
         | Remarks == "Livestock predation (Method is \"prospection suite à dégât\"). The other cub (Auberta) dies in 2014"
         | Remarks == "Caramelles with cubs. Livestock predation (Method is \"prospection suite à dégât\")"
         | Remarks == "ID not found: S23-SLO11. Livestock predation (Method is \"prospection suite à dégât\")"
         | Remarks == "S23-SLO11 + Esmolet. Livestock predation (Method is \"prospection suite à dégât\")"
         | Remarks == "Livestock predation (Method is \"prospection suite à dégât\"). Plusieurs crottes avec différentes dates estimées : 25 au 28/07, 26 au 29/07, 18 au 19/07, 23 au 24/07."
         | Remarks == "depredació a 1 ovella de Toniquet" 
         | Remarks == "ID France: S29-SLO6. Livestock predation (Method is \"prospection suite à dégât\")"
         | Remarks == "Livestock predation (Method is \"prospection suite à dégât\"). Sheep"
         | Remarks == "Indeterminat. Ataque ovejas en Baqueira. Qualitat 100%."
         | Remarks == "Livestock predation. Identified in database Seguiment França Ossos 2016-2020.xls" 
         | Remarks == "Livestock predation (Method is \"prospection suite à dégât\"). Identified in database Seguiment França Ossos 2016-2020.xls" 
         | Remarks == "2 ovejas asfixiadas en el corral de un ataque anterior"
         | Remarks == "Livestock predation. 5 hives. Attributed to Goiat" 
         | Remarks == "ID France: S29-SLO6. Livestock predation" 
         | Remarks == "ID France: S29-SLO6.Livestock predation (Method is \"prospection suite à dégât\")"
         | Remarks == "Mont Atac colmenas" 
         | Remarks == "1 Oveja" 
         | Remarks == "2 Ovejas y 1 Cordero" 
         | Remarks == "Sant Joan  1 oveja 1 cordero 1 cabra" 
         | Remarks == "ID France: New18_11. Livestock predation (Method is \"prospection suite à dégât\")" 
         | Remarks == "Fecha Noche oso. 1 sheep injured" 
         | Remarks == "Fecha Noche oso. 14 beehives" 
         | Remarks == "ID France: New 18-15. Livestock predation (Method is \"prospection suite à dégât\")" 
         | Remarks == "ID France: New 18-14. Livestock predation (Method is \"prospection suite à dégât\")" 
         | Remarks == "Fecha Noche oso. 1 goat" 
         | Remarks == "Fecha Noche oso. 1 sheep, 1 ram" 
         | Remarks == "Fecha Noche oso. 1 sheep" 
         | Remarks == "Fecha Noche oso. 2 beehives"
         | Remarks == "Fecha Noche oso. 2 cattle"
         | Remarks == "ID France: New 19-04. Livestock predation (Method is \"prospection suite à dégât\")" 
         | Remarks == "ID France: New 18-18. Livestock predation (Method is \"prospection suite à dégât\")" 
         | Remarks == "Fecha Noche oso. 14 beehiuves" 
         | Remarks == "Fecha Noche oso. 11 beehives" 
         | Remarks == "Fecha recogida dato. 1 beehive"  
         | Remarks == "Fecha Noche oso. 5 beehives"  
         | Remarks == "Fecha Noche oso. 4 beehives"  
         | Remarks == "Fecha Noche oso. 13 beehives"  
         | Remarks == "Fecha Noche oso. 12 beehives"
         | Remarks == "Fecha recogida dato. 3 beehives"
         | Remarks ==  "Fecha recogida dato. 1 sheep" 
         | Remarks == "Livestock predation (Method is \"prospection suite à dégât\"). Cubs died at some point in 2019"  
         | Remarks == "ID France: New 18-06. Livestock predation (Method is \"prospection suite à dégât\")"  
         | Remarks == "Fecha Noche oso. 2 sheep" 
         | Remarks == "New 18-06. Livestock predation (Method is \"prospection suite à dégât\")"  
         | Remarks == "Fecha Noche oso. 1 sheep, 1 lamb" 
         | Remarks == "2 CABRAS MUERTAS"  
         | Remarks == "CORDERA"  
         | Remarks ==  "1 MARDANO+1OVEJA+1CORDERA+1OVEJA NO ENCONTRADA+5 DESAPARECIDAS" 
         | Remarks == "1Ovins" 
         | Remarks == "3Ovins"
         | Remarks ==  "2Ovins" 
         | Remarks ==  "4Ovins"  
         | Remarks ==  "1 CABRA"  
         | Remarks ==  "8Ovins"   
         | Remarks ==  "8 COLMENAS"  
         | Remarks ==  "2 COLMENAS" 
         | Remarks ==   "15Ovins" 
         )
           
coordinates(os) <- os[,c(11,12)]
proj4string(os) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
os <- spTransform(os, CRS = CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs"))

coordinates(Xpoints) <- Xpoints[,c(1,2)]
proj4string(Xpoints) <- proj4string(os)

os@data <- data.frame(ID = 1:length(os))
Xpoints@data <- data.frame(ID = 1:length(Xpoints))

traps_opport <- rbind(os, Xpoints)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
writeOGR(traps_opport, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "allobs_trapsOpport", driver = "ESRI Shapefile")  

osMAE <- readOGR("Min_Pol_Maelis_clip.shp")

Xbuf3 <- gBuffer(traps_opport, width = 10)

plot(e)
plot(Xbuf, col = "lightblue", add = TRUE)
plot(full_buf_noOut_cropped, col = adjustcolor("green", alpha = 0.5), add = TRUE)
plot(Xbuf2, col = "pink", add = TRUE)
plot(osMAE, col = adjustcolor("red", alpha = 0.5), add = TRUE)


unique(os$Remarks)

os <- os %>%
  filter(Obs_type %in% c("Radiotracking")


# SAVE THE HABITAT COORDINATES UNSCALED

distcoreMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

G <- coordinates(distcoreMask)[!is.na(distcoreMask[]),]
colnames(G) <- c("x","y")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
save(G, file = "habcoord.RData") 


