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

tr <- readOGR("D:/MargSalas/Ganga/Data/FarmdindisDS/GIS", "Trans_Ganga_FarmDindis") 

# 1- We take all transects where species is present or absent (all categories)
# Abs-Pres classification (David Giralt):
#   0: Transectos donde nunca se ha visto ganga (2010-2022)
#   1: Transectos donde se han visto al menos 1 año en el periodo 2010-2022
#   2: Transectos donde se ha visto al menos 1 año entre 2010-2016


# 1 Extract zone of habitat quality where each transect is located ----

## ---- 1.1. Buffer around transects ----

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
#writeOGR(all_pol2, "Buffer", dsn = "D:/MargSalas/Ganga/Data/FarmdindisDS/GIS", driver = "ESRI Shapefile")

all_pol2 <- readOGR("D:/MargSalas/Ganga/Data/FarmdindisDS/GIS", "BufferImproved") 


## ---- 1.2. Rasterize and extract % habitat quality in buffers ----

pa <- readOGR("D:/MargSalas/Ganga/Data/FarmdindisDS/GIS", "zonesGanga") 
# Rasterize
ext <- extent(pa)
r <- raster(pa, res = 25)
raster_pa <- rasterize(pa, r, field = "X1")

habBuf <- extract(raster_pa, all_pol2, na.rm = TRUE)

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

mapview(raster_pa) + 
  mapview(tr) + 
  mapview(all_pol2)

## ---- 2.3. Rasterize and extract % habitat quality in line transect ----

trline <- gBuffer(tr, byid = TRUE, width = 1, capStyle = "FLAT") # paths are 2 meters wide

habLine <- extract(raster_pa, trline, na.rm = TRUE)

df_habLine <- data.frame(matrix(NA, nrow = length(habLine), ncol = 5))
colnames(df_habLine) <- c("TransectID", "hq0", "hq1", "hq2", "hq3")
df_habLine$TransectID <- tr@data$Codi

for (i in 1:length(habLine)){
  df_habLine[i,2] <- length(which(is.na(habLine[[i]])))/length(habLine[[i]])
  df_habLine[i,3] <- length(which(habLine[[i]] == 1))/length(habLine[[i]])
  df_habLine[i,4] <- length(which(habLine[[i]] == 2))/length(habLine[[i]])
  df_habLine[i,5] <- length(which(habLine[[i]] == 3))/length(habLine[[i]])
}
df_habLine[21,c(2:5)] <- c(1,0,0,0)
df_habLine[22,c(2:5)] <- c(1,0,0,0)
df_habLine[23,c(2:5)] <- c(1,0,0,0)

rowSums(df_habLine[,c(2:5)])
df_habLine$habLine_assign_LINE <- NA

for (i in 1:nrow(df_habLine)){
  df_habLine$habLine_assign_LINE[i] <- colnames(df_habLine)[which(df_habLine[i,] == apply(df_habLine[,c(2:5)], 1, max)[i])]
}

mapview(raster_pa) + 
  mapview(tr) + 
  mapview(trline) +
  mapview(all_pol2)

## ---- 1.3. Join both datasets ----

df_hab <- cbind(df_habBuf, df_habLine[,c(2:6)]) 

## ---- 2. Estimate TRANSECT HABITAT QUALITY variable ----
# Do it from the buffer, it is more realistic about the area sampled by the transect
# Create "Weighted habitat quality"/transect
# Like this you take into account the proportion of each habitat in the covariate for lambda
# And estimate a relationship between lambda and habitat quality

# Re-arrange data frame so that it goes from lowest (hq0) to highest (hq3) quality
df_habBuf2 <- df_habBuf
colnames(df_habBuf2) <- c("TransectID", "hq0", "hq3", "hq2", "hq1")
df_habBuf2 <- df_habBuf2[,c(1,2,5,4,3)]

# Estimate variable
df_habBuf2$WeightedQuality <- NA

for (i in 1:nrow(df_habBuf2)){
  df_habBuf2$WeightedQuality[i] <- df_habBuf2$hq0[i]*0 + df_habBuf$hq1[i]*1 + df_habBuf$hq2[i]*2 + df_habBuf$hq3[i]*3
}

summary(df_habBuf2$WeightedQuality)

# Save
setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
write.csv(df_habBuf2, file = "HQvariable.csv")
