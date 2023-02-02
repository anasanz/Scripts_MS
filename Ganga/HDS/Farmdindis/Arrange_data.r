
## -------------------------------------------------
##      Arrange data DS Farmdindis 2010-2022
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

dat <- read.csv("Transectes_Farmdindis_Ganga_TOTS anys.csv", sep = ";")

colnames(dat)[which(colnames(dat) == "ï..Codi_seca")] <- "Region.Label"
colnames(dat)[which(colnames(dat) == "anys")] <- "Year"
colnames(dat)[which(colnames(dat) == "Hora_inici")] <- "Start_time"
colnames(dat)[which(colnames(dat) == "Hora_final")] <- "End_time"
colnames(dat)[which(colnames(dat) == "Nuvolositat")] <- "Clouds"
colnames(dat)[which(colnames(dat) == "Temperatura")] <- "Temp"
colnames(dat)[which(colnames(dat) == "codiEspecie")] <- "Species"
colnames(dat)[which(colnames(dat) == "Num")] <- "Count" 
dat$Effort <- 500

## ---- Create variable transectID, than matches with the code of the GIS layers (i.e., two digits: 09) ----

#1. Add a 0 before the transect number
for (i in 1:nrow(dat)){ 
  dat$Num_transecte[i] <- paste(0,dat$Num_transecte[i], sep = "")
}

#2. Keep only the last 2 digits (or 3 in the case of the transects that contain 100)

for (i in 1:nrow(dat)) { 
  cent <- substr(dat$Num_transecte[i], 4,4)
  cent <- as.numeric(cent) # NA if it doesnt have 4 digits
  if(is.na(cent)) { # if is NA (has 3 digits)
    dat$Num_transecte[i] <- str_sub(dat$Num_transecte[i], start = -2) # Keep the last 2
  } else { dat$Num_transecte[i] <- str_sub(dat$Num_transecte[i], start = -3)} # Otherwise, keep the last 3
}
# Create variable by pasting it
for (i in 1:nrow(dat)){ 
  dat$transectID[i] <- paste(dat$Region.Label[i],dat$Num_transecte[i], sep = "")
}

## ---- Transect-Year variable ----

for (i in 1:nrow(dat)){ 
  dat$T_Y[i] <- paste(dat$transectID[i],dat$Year[i], sep = "_")
}

## ---- Add info of distance bands ----

# Banda 1: 0-25
# Banda 2: 25-50
# Banda 3: 50-100
# Banda 4: 100-200
# Banda 5: 200-500

dat$distance <- NA # Medium point of each bin

for (i in 1:nrow(dat)){
  if (!is.na(dat$Banda[i])){
    if (dat$Banda[i] == 1) {dat$distance[i] = 12.5}
    else if (dat$Banda[i] == 2) {dat$distance[i] = 37.5}
    else if (dat$Banda[i] == 3) {dat$distance[i] = 75}
    else if (dat$Banda[i] == 4) {dat$distance[i] = 150}
    else if (dat$Banda[i] == 5) {dat$distance[i] = 350}
  }}

## ---- Explore detection curves ----

# All years

hist(dat$distance[which(!is.na(dat$distance))], breaks = c(0,25,50,100,200,500),
     main = "PTALC all years", col = "grey", freq = FALSE) 

years <- unique(dat$Year)
for (t in 1:length(years)){
  hist(dat$distance[which(dat$Year == years[t] & !is.na(dat$distance))], breaks = c(0,25,50,100,200,500),
       main = paste("PTALC", years[t]), col = "grey", freq = FALSE) 
}

## ---- Choose analyzed transects and assign habitat quality ----

tr <- readOGR("D:/MargSalas/Ganga/Data/FarmdindisDS", "Trans_Ganga_FarmDindis") 

# 1- We take all transects where species is present or absent (all categories)
# Abs-Pres classification (David Giralt):
#   0: Transectos donde nunca se ha visto ganga (2010-2022)
#   1: Transectos donde se han visto al menos 1 año en el periodo 2010-2022
#   2: Transectos donde se ha visto al menos 1 año entre 2010-2016


# 2. Extract zone of habitat quality where each transect is located ----

## ---- 2.1. Buffer around transects ----

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


mapview(trbuf1) +
  mapview(trbuf2) +
  mapview(all_pol) +
  mapview(tr) + 
  mapview(all_pol) + 
  mapview(all_pol2)

# Not perfect but I keep all_pol2, possibly change it manually

## ---- 2.2. Rasterize and extract % habitat quality in buffers ----

pa <- readOGR("D:/MargSalas/Ganga/Data/FarmdindisDS", "zonesGanga") 
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

## ---- 2.3. Join both datasets ----

df_hab <- cbind(df_habBuf, df_habLine[,c(2:6)]) 

datcounts <- dat[which(dat$Count != 0 & dat$Count != -1), ] # Dataset with only the counts


# Test how many transects we would keep if we would delete those in non-habitat
df_hab_subset <- df_hab[which(df_hab$hq0 != 1), ]

dat1 <- datcounts[which(datcounts$transectID %in% df_hab_subset$TransectID), ] # Keeping obs from those transects
nrow(datcounts) - nrow(dat1) # We would lose 9 obs

# Test how many transects we would keep if we would delete those with majority of non-habitat
onlyHabBuf <- df_hab[which(df_hab$habBuf_assign_BUF != "hq0"), ]
dat2 <- datcounts[which(datcounts$transectID %in% onlyHabBuf$TransectID), ]
nrow(datcounts) - nrow(dat2) # Only 32 obs

hist(dat2$distance[which(!is.na(dat2$distance))], breaks = c(0,25,50,100,200,500),
     main = "PTALC all years", col = "grey", freq = FALSE) 

for (t in 1:length(years)){
  
  hist(dat$distance[which(dat$Year == years[t] & !is.na(dat$distance))], breaks = c(0,25,50,100,200,500),
       main = paste("PTALC", years[t]), col = "grey", freq = FALSE) 
  
  hist(dat2$distance[which(dat2$Year == years[t] & !is.na(dat2$distance))], breaks = c(0,25,50,100,200,500),
       main = paste("PTALC", years[t]), col = "grey", freq = FALSE) 
}

#### For the detection curve it only looks worse for 2016 and 2021