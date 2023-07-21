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

setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
dat <- read.csv(file = "Data_HDS_Specific.csv")

## ---- Choose analyzed transects and assign habitat quality ----

tr <- readOGR("D:/MargSalas/Ganga/Data/GIS", "Trans_Ganga_500_Specific")
tr2 <- readOGR("D:/MargSalas/Ganga/Data/GIS", "Trans_Ganga_FarmDindis") 

tr <- spTransform(tr, CRSobj = proj4string(tr2))

## ---- Extract zone of habitat quality where each transect is located: buffer around transects ----

trbuf1 <- gBuffer(tr, byid = TRUE, width = 400, capStyle = "FLAT") # We take 400 because is the truncation distance of specific transects
trbuf2 <- gBuffer(tr, byid = TRUE, width = 400, capStyle = "SQUARE")

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


mapview(trbuf1) +
  mapview(trbuf2) +
  mapview(all_pol) +
  mapview(tr) + 
  mapview(all_pol) + 
  mapview(all_pol2)

# This one doesn't need improvement, so I keep it as it is and I don't improve it manually (as I did with Farmdindis)
#writeOGR(all_pol2, "Buffer_Specific", dsn = "D:/MargSalas/Ganga/Data/GIS", driver = "ESRI Shapefile")

all_pol2 <- readOGR("D:/MargSalas/Ganga/Data/GIS", "Buffer_Specific") 

## -------------------------------------------------
##                  SDM 2021
##   Extract habitat quality of transects (take obviously 2021)
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
df_habBuf$TransectID <- tr@data$CODI

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

## ---- 3. Explore how observations are related to habitat quality ----

dat2 <- dat2[-which(dat2$Count < 1), ]

# Relation detections - hab.quality

hist(dat2$WeightedQuality2021, main = "Detections - HQ (all years)")

# Save
setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
write.csv(df_habBuf2, file = "HQvariable2021_Specific.csv")

## Save removing transects that are hq = 0
#df_habBuf3 <- df_habBuf2[which(df_habBuf2$WeightedQuality != 0), ]

#setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
#write.csv(df_habBuf3, file = "HQvariable2.csv")

## ---- 4. Area of each habitat quality (same in farmdindis and specific) ----

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

