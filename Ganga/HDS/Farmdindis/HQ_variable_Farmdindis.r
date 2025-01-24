## -------------------------------------------------
##            HABITAT QUALITY VARIABLE
## ------------------------------------------------- 


rm(list=ls())

library(dplyr)
library(stringr)
library(rgdal)
library(sf)
library(raster)
library(rgeos)
library(mapview)
library(spatialEco)

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
dat <- read.csv(file = "Data_HDS_Farmdindis.csv")

## ---- Choose analyzed transects and assign habitat quality ----

tr <- readOGR("D:/MargSalas/Ganga/Data/GIS", "Trans_Ganga_FarmDindis") 

proj4string(tr)
# We take all transects where species is present or absent (all categories)
# Abs-Pres classification (David Giralt):
#   0: Transectos donde nunca se ha visto ganga (2010-2022)
#   1: Transectos donde se han visto al menos 1 año en el periodo 2010-2022
#   2: Transectos donde se ha visto al menos 1 año entre 2010-2016

## ---- Extract zone of habitat quality where each transect is located: buffer around transects ----

trbuf1 <- gBuffer(tr, byid = TRUE, width = 500, capStyle = "FLAT")
trbuf2 <- gBuffer(tr, byid = TRUE, width = 500, capStyle = "SQUARE")

# I want a combination of this two but it is very challenging...

mcp <- function(x, percentile=95){
  
  centroid <- sf::st_centroid(sf::st_union(x))
  dist <- as.numeric(sf::st_distance(x, centroid))
  within_percentile_range <- dist <= quantile(dist, percentile/100)
  x_filter <- st_union(x[within_percentile_range,])
  st_convex_hull(x_filter)
  
} # Function to calculate mcp of st objects

pol_list <- list()
for (i in 1:length(trbuf1)){
  pol <- trbuf1[i, ]
  poltr <- tr[i,]
  pol <- st_as_sf(pol)
  poltr <- st_as_sf(poltr)
  
  pointspol <- extract.vertices(pol)
  pointstr <- extract.vertices(poltr)
  
  all <- rbind(pointspol, pointstr)
  
  h <- mcp(all)
  
  pol_list[[i]] <- as(h,"Spatial")
  pol_list[[i]]$ID <-i
}
all_pol <- do.call(rbind, pol_list)


# I finally get the maximum polygon of trbuf1 and all_pol

pol_list2 <- list()

for (i in 1:length(trbuf1)){
  pol <- trbuf1[i, ]
  poltr <- all_pol[i,]
  pol <- st_as_sf(pol)
  poltr <- st_as_sf(poltr)
  
  single_sf <- dplyr::bind_rows(pol, poltr)
  dissolve_sf <- st_union(single_sf)
  
  pol_list2[[i]] <- as(dissolve_sf,"Spatial")
  pol_list2[[i]]$ID <-i
}
all_pol2 <- do.call(rbind, pol_list2)


#mapview(trbuf1) +
#  mapview(trbuf2) +
#  mapview(all_pol) +
#  mapview(tr) + 
#  mapview(all_pol) + 
#  mapview(all_pol2)

# Not perfect but I keep all_pol2. I save it to improve it manually
#writeOGR(all_pol2, "Buffer_Farmdindis", dsn = "D:/MargSalas/Ganga/Data/FarmdindisDS/GIS", driver = "ESRI Shapefile")

all_pol2 <- readOGR("D:/MargSalas/Ganga/Data/GIS", "BufferImproved_Farmdindis") 

## -------------------------------------------------
##                  SDM 2011
##   Extract habitat quality of transects 2010-2015
## ------------------------------------------------- 

pa1 <- readOGR("D:/MargSalas/Ganga/Data/GIS", "zonesGanga2011")
pa1$ZONA <- as.numeric(pa1$ZONA)
# Zonas optimas = 1 (high); Buenas = 2 (medium); Adecuadas: 3 (low)

mapview(all_pol2) + mapview(pa1, zcol = "ZONA")
mapview(raster_pa1)

## ---- 1 Rasterize and extract % habitat quality in buffers ----
# Rasterize
ext <- extent(pa1)
r <- raster(pa1, res = 25)
raster_pa1 <- rasterize(pa1, r, field = "ZONA")

habBuf <- extract(raster_pa1, all_pol2, na.rm = TRUE)

df_habBuf <- data.frame(matrix(NA, nrow = length(habBuf), ncol = 5))
colnames(df_habBuf) <- c("TransectID", "hq0", "hq1", "hq2", "hq3")
df_habBuf$TransectID <- tr@data$Codi

for (i in 1:length(habBuf)){
  df_habBuf[i,2] <- length(which(is.na(habBuf[[i]])))/length(habBuf[[i]])
  df_habBuf[i,3] <- length(which(habBuf[[i]] == 1))/length(habBuf[[i]])
  df_habBuf[i,4] <- length(which(habBuf[[i]] == 2))/length(habBuf[[i]])
  df_habBuf[i,5] <- length(which(habBuf[[i]] == 3))/length(habBuf[[i]])
}


rowSums(df_habBuf[,c(2:5)])
df_habBuf$habBuf_assign_BUF <- NA

for (i in 1:nrow(df_habBuf)){
  df_habBuf$habBuf_assign_BUF[i] <- colnames(df_habBuf)[which(df_habBuf[i,] == apply(df_habBuf[,c(2:5)], 1, max)[i])]
}

mapview(raster_pa1) + 
  mapview(tr) + 
  mapview(all_pol2)

df_habBuf2011 <- df_habBuf
## ---- 2. Estimate TRANSECT HABITAT QUALITY variable ----

# Do it from the buffer, it is more realistic about the area sampled by the transect
# Create "Weighted habitat quality"/transect
# Like this you take into account the proportion of each habitat in the covariate for lambda
# And estimate a relationship between lambda and habitat quality

# Re-arrange data frame so that it goes from lowest (hq0) to highest (hq3) quality
df_habBuf2 <- df_habBuf
colnames(df_habBuf2) <- c("transectID", "hq0", "hq3", "hq2", "hq1")
df_habBuf2 <- df_habBuf2[,c(1,2,5,4,3)]

# Estimate variable
df_habBuf2$WeightedQuality2011 <- NA

for (i in 1:nrow(df_habBuf2)){
  df_habBuf2$WeightedQuality2011[i] <- df_habBuf2$hq0[i]*0 + df_habBuf2$hq1[i]*1 + df_habBuf2$hq2[i]*2 + df_habBuf2$hq3[i]*3
}

summary(df_habBuf2$WeightedQuality2011)

# It goes from 0 habitat quality to 3 of maximum habitat quality (if there would be a 100% of habitat 3: 3*1 = 3)
# Explore how what is the habitat quality associated to the observations

# Join with observation data
dat2 <- left_join(dat, df_habBuf2, by = "transectID")

## ---- 3. Explore how observations are related to habitat quality ----

dat2 <- dat2[-which(dat2$Count < 1), ]

# Relation detections - hab.quality
par(mfrow = c(1,1))
hist(dat2$WeightedQuality2011[which(dat2$Year %in% c(2010:2015))], main = "Detections - HQ (2010-2015)")
hist(dat2$WeightedQuality2011, main = "Detections - HQ (all years)")

lowhq10_15 <- dat2[which(dat2$Year %in% c(2010:2015) & dat2$WeightedQuality2011 < 0.5), ]

# First years (more fitted to 2011 SDM)

year <- 2010:2015

par(mfrow = c(2,3))
for (t in 1:length(year)){
  dy <- dat2[which(dat2$Year %in% year[t]), ]
  hist(dy$WeightedQuality, main = paste("Detections - HQ (", year[t], ")", sep = ""), breaks = 10)
}

# Save
setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
write.csv(df_habBuf2, file = "HQvariable2011_Farmdindis.csv")

# Save removing transects that are hq = 0
#df_habBuf3 <- df_habBuf2[which(df_habBuf2$WeightedQuality != 0), ]
#
#setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
#write.csv(df_habBuf3, file = "HQvariable2.csv")

## ---- 4. Area of each habitat quality ----

library(rgeos)

area1 <- as.vector(gArea(pa1,byid = TRUE)) #Area
pa1@data <- cbind(pa1@data,area1)
df1 <- as.data.frame(pa1@data)

area_zona1 <- aggregate(df1$area1, by = list(df1$ZONA), FUN = sum)
area_zona_HA1 <- cbind(area_zona1, area_zona1$x/10000) # En HA

## CAREFUL!! IN MY ANALYSIS 1 IS LOWEST AND 3 IS HIGHEST, CHANGE CAT NAMES So THAT IT IS THIS WAY
area_zona_HA1$Group.1[1] <- 3
area_zona_HA1$Group.1[3] <- 1

area_zona_HA1 <- arrange(area_zona_HA1, Group.1)


## -------------------------------------------------
##                  SDM 2021
##   Extract habitat quality of transects 2016-2021
## -------------------------------------------------

pa2 <- readOGR("D:/MargSalas/Ganga/Data/GIS", "zonesGanga2021") 

## ---- 1. Rasterize and extract % habitat quality in buffers ----

# Rasterize
ext <- extent(pa2)
r <- raster(pa2, res = 25)
raster_pa2 <- rasterize(pa2, r, field = "X1")

habBuf <- extract(raster_pa2, all_pol2, na.rm = TRUE)

df_habBuf <- data.frame(matrix(NA, nrow = length(habBuf), ncol = 5))
colnames(df_habBuf) <- c("TransectID", "hq0", "hq1", "hq2", "hq3")
df_habBuf$TransectID <- tr@data$Codi

for (i in 1:length(habBuf)){
  df_habBuf[i,2] <- length(which(is.na(habBuf[[i]])))/length(habBuf[[i]])
  df_habBuf[i,3] <- length(which(habBuf[[i]] == 1))/length(habBuf[[i]])
  df_habBuf[i,4] <- length(which(habBuf[[i]] == 2))/length(habBuf[[i]])
  df_habBuf[i,5] <- length(which(habBuf[[i]] == 3))/length(habBuf[[i]])
}
df_habBuf[22,c(2:5)] <- c(1,0,0,0)
df_habBuf[23,c(2:5)] <- c(1,0,0,0)

rowSums(df_habBuf[,c(2:5)])
df_habBuf$habBuf_assign_BUF <- NA

for (i in 1:nrow(df_habBuf)){
  df_habBuf$habBuf_assign_BUF[i] <- colnames(df_habBuf)[which(df_habBuf[i,] == apply(df_habBuf[,c(2:5)], 1, max)[i])]
}

mapview(raster_pa2) + 
  mapview(tr) + 
  mapview(all_pol2)


## ---- 2. Estimate TRANSECT HABITAT QUALITY variable ----
# Do it from the buffer, it is more realistic about the area sampled by the transect
# Create "Weighted habitat quality"/transect
# Like this you take into account the proportion of each habitat in the covariate for lambda
# And estimate a relationship between lambda and habitat quality

# Re-arrange data frame so that it goes from lowest (hq0) to highest (hq3) quality
df_habBuf2 <- df_habBuf
colnames(df_habBuf2) <- c("transectID", "hq0", "hq3", "hq2", "hq1")
df_habBuf2 <- df_habBuf2[,c(1,2,5,4,3)]

# Estimate variable
df_habBuf2$WeightedQuality2021 <- NA

for (i in 1:nrow(df_habBuf2)){
  df_habBuf2$WeightedQuality2021[i] <- df_habBuf2$hq0[i]*0 + df_habBuf2$hq1[i]*1 + df_habBuf2$hq2[i]*2 + df_habBuf2$hq3[i]*3
}

summary(df_habBuf2$WeightedQuality2021)

# It goes from 0 habitat quality to 3 of maximum habitat quality (if there would be a 100% of habitat 3: 3*1 = 3)
# Explore how what is the habitat quality associated to the observations

# Join with observation data
dat2 <- left_join(dat, df_habBuf2, by = "transectID")
#dat2 <- left_join(dat2, df_habBuf2, by = "transectID")

dat2$WeightedQuality2011 == dat2$WeightedQuality2021
## ---- 3. Explore how observations are related to habitat quality ----

dat2 <- dat2[-which(dat2$Count < 1), ]

# Relation detections - hab.quality

hist(dat2$WeightedQuality2021, main = "Detections - HQ (all years)")

hist(dat2$WeightedQuality2021[which(dat2$Year %in% c(2015:2022))], main = "Detections - HQ (2015-2022)")
lowhq15_20 <- dat2[which(dat2$Year %in% c(2015:2022) & dat2$WeightedQuality2021 < 0.5), ]


# Last 6 years (more fitted to updated SDM)

year <- 2016:2022

par(mfrow = c(2,4))
for (t in 1:length(year)){
  dy <- dat2[which(dat2$Year %in% year[t]), ]
  hist(dy$WeightedQuality, main = paste("Detections - HQ (", year[t], ")", sep = ""), breaks = 10)
}

# It looks almost quadratic

# Save
setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
write.csv(df_habBuf2, file = "HQvariable2021_Farmdindis.csv")

## Save removing transects that are hq = 0
#df_habBuf3 <- df_habBuf2[which(df_habBuf2$WeightedQuality != 0), ]

#setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
#write.csv(df_habBuf3, file = "HQvariable2.csv")

## ---- 4. Area of each habitat quality ----

pa2 <- readOGR("D:/MargSalas/Ganga/Data/GIS", "zonesGanga2021") 

library(rgeos)

area2 <- as.vector(gArea(pa2,byid = TRUE)) #Area
pa2@data <- cbind(pa2@data,area2)
df2 <- as.data.frame(pa2@data)

area_zona2 <- aggregate(df2$area2, by = list(df2$X1), FUN = sum)
area_zona_HA2 <- cbind(area_zona2, area_zona2$x/10000) # En HA

## CAREFUL!! IN MY ANALYSIS 1 IS LOWEST AND 3 IS HIGHEST, CHANGE CAT NAMES So THAT IT IS THIS WAY
area_zona_HA2$Group.1[1] <- 3
area_zona_HA2$Group.1[3] <- 1

area_zona_HA2 <- arrange(area_zona_HA2, Group.1)
write.csv(area_zona_HA2, file = "HQ_area.csv")

## ---- 5. Add column "HQ" to categorical ----

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")

hq_dat2011 <- read.csv(file = "HQvariable2011_Farmdindis.csv") 
hq_dat2021 <- read.csv(file = "HQvariable2021_Farmdindis.csv") 

hq_dat2011$hq2011 <- apply(hq_dat2011[,3:6], 1, which.max)-1 
hq_dat2021$hq2012 <- apply(hq_dat2021[,3:6], 1, which.max)-1

write.csv(hq_dat2011, file = "HQvariable2011_Farmdindis.csv")
write.csv(hq_dat2021, file = "HQvariable2021_Farmdindis.csv")


      