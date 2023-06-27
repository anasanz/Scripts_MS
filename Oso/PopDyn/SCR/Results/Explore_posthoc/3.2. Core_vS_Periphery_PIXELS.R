## -------------------------------------------------
##           SPATIAL NUCLEUS vS PERIPHERY
##        SECOND TRIAL: PIXELS DISTANCE TO CORE 
## ------------------------------------------------- 

# POST-HOC ANALYSIS

# Identify the spatial nucleus (core) of the population and the periphery, to estimate the age structure in both areas
# The age structure will be estimated as the number of individuals of each age class in nucleus and periphery

# The core and periphery need to occupy the sampling area (not the state space), so I need to check how to define it
# - 1º trial: Buffer of female average dispersal distance around the core polygons (origin) for the distance to core covariate
# - 2º trial: Classification with distance to core covariate (e.g., take pixels with highest and lowest Distance to Core values)

rm(list = ls())

# library(rgdal)
# library(dplyr)
library(rgeos)
# library(MCMCvis)
library(rgdal)
library(raster)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
logDistcore <- raster("logDistcore_hrbear.tif")


# Load core polygon
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol")
corepol <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol/corepol.shp", "corepol")
corepol@data$ID <- 1
corepol@data <- corepol@data[2]

# Load state space
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp") # Load sampling area (where we estimate abundance)
proj4string(Xbuf2) <- proj4string(corepol)
Xbuf2@data$ID <- 1

plot(logDistcore)
plot(Xbuf, add = TRUE)
plot(Xbuf2, border = "white", add = TRUE)
plot(corepol, border = "red", add = TRUE)

## ---- 1. Create polygons from DistCore Variable ----

xyDistCore <- data.frame(xyFromCell(logDistcore, 1:ncell(logDistcore))) # Extract central xy coordinates from logDistcore raster
coordinates(xyDistCore) <- xyDistCore[ ,c("x", "y")]
proj4string(xyDistCore) <- proj4string(logDistcore)

pbuf <- over(xyDistCore,Xbuf2) # Central cell points that overlap with sampling buffer

plot(logDistcore)
points(xyDistCore[which(pbuf == 1),])

logDistCore.buf <- logDistcore[which(pbuf == 1)] # Values of DistCore within Buffer
q1 <- summary(logDistCore.buf)[2] # Quantiles to know which values subset to create core and peripery polygon
q2 <- quantile(logDistCore.buf, 0.4)
q3 <- summary(logDistCore.buf)[3]
q4 <- quantile(logDistCore.buf, 0.6)
q5 <- summary(logDistCore.buf)[5]

## ---- 1.1. Polygon definition from MEDIAN ----

plot(logDistcore)
points(xyDistCore[which(values(logDistcore) < q3 & pbuf == 1),])
points(xyDistCore[which(values(logDistcore) > q3 & pbuf == 1),], col = "red")

pCorePeri <- logDistcore # New raster with values according to new polygons
pCorePeri[] <- NA
pCorePeri[which(values(logDistcore) < q3 & pbuf == 1)] <- 1
pCorePeri[which(values(logDistcore) > q3 & pbuf == 1)] <- 2

pCoreMedian <- rasterToPolygons(pCorePeri, fun=function(x){x == 1}) # Convert raster to polygons
pPerMedian <- rasterToPolygons(pCorePeri, fun=function(x){x == 2})

pCoreMedian <- SpatialPolygonsDataFrame(aggregate(pCoreMedian,dissolve=T), data = data.frame(ID = 1)) # Aggregate in 1 polygon and convert to spdf
pPerMedian <- SpatialPolygonsDataFrame(aggregate(pPerMedian,dissolve=T), data = data.frame(ID = 1))


## ---- 1.2. Polygon definition from 4th and 6th DECILES ----

plot(logDistcore)
points(xyDistCore[which(values(logDistcore) < q2 & pbuf == 1),])
points(xyDistCore[which(values(logDistcore) > q4 & pbuf == 1),], col = "red")

pCorePeri <- logDistcore # New raster with values according to new polygons
pCorePeri[] <- NA
pCorePeri[which(values(logDistcore) < q2 & pbuf == 1)] <- 1
pCorePeri[which(values(logDistcore) > q4 & pbuf == 1)] <- 2

pCoreDeciles <- rasterToPolygons(pCorePeri, fun=function(x){x == 1}) # Convert raster to polygons
pPerDeciles <- rasterToPolygons(pCorePeri, fun=function(x){x == 2})

pCoreDeciles <- SpatialPolygonsDataFrame(aggregate(pCoreDeciles,dissolve=T), data = data.frame(ID = 1)) # Aggregate in 1 polygon and convert to spdf
pPerDeciles <- SpatialPolygonsDataFrame(aggregate(pPerDeciles,dissolve=T), data = data.frame(ID = 1))


par(mfrow = c(1,2))

plot(Xbuf)
plot(logDistcore, add = TRUE) # Check
plot(pCoreMedian, add = TRUE)
plot(pPerMedian, col = "red", add = TRUE)
mtext("Median", side = 3, line = -3, cex = 2)

plot(Xbuf)
plot(logDistcore, add = TRUE) # Check
plot(pCoreDeciles, add = TRUE)
plot(pPerDeciles, col = "red", add = TRUE)
mtext("4th and 6th Deciles", side = 3, line = -3, cex = 2)



## ---- 2. Estimation of age structure in nuc vS periphery, different definitions from - to + restrictive ----

library(nimbleSCR) # Load here to avoid problems with raster

# Load data

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("myResults_3-3.1_param.RData")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams") # Load results
load("myResults_3-3.1_sxy.RData")

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y') # Unscale the sxy coordinates
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled

## Proportion of each age class relative to the zone (core or periphery) ALL YEARS TOGETHER 
# Function for estimating age structure in two polygons: nucleous and periphery

age_st_allYears <- function(nucleus, periphery){
  
  # 1. Total abundance in each zone
  
  zones <- c(nucleus, periphery)
  
  # Matrix to store abundance in the buffer each iteration and ALL years
  NIn_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(myResultsSXYZ$sims.list$z[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        NIn_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }
    }
  }
  
  NIn_trapBuf_AllYears <- apply(NIn_trapBuf, 3, rowSums)
  
  # 2. Abundance of each age category in each zone
  
  ## CUBS
  ZZcubs <- myResultsSXYZ$sims.list$z
  ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  
  cub_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZcubs[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        cub_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  cub_trapBuf_AllYears <- apply(cub_trapBuf, 3, rowSums)
  
  ## SUBADULTS 
  
  ZZsub <- myResultsSXYZ$sims.list$z
  ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  
  subad_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZsub[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        subad_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  subad_trapBuf_AllYears <- apply(subad_trapBuf, 3, rowSums)
  
  
  ## ADULTS
  
  ZZad <- myResultsSXYZ$sims.list$z
  ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  
  ad_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        ad_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ad_trapBuf_AllYears <- apply(ad_trapBuf, 3, rowSums)
  
  
  ## ADULT FEMALES
  
  ZZadFEM <- myResultsSXYZ$sims.list$z
  ZZadFEM[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  ZZadFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  adFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZadFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        adFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  adFEM_trapBuf_AllYears <- apply(adFEM_trapBuf, 3, rowSums)
  
  
  ## ADULT MALES
  
  ZZadMAL <- myResultsSXYZ$sims.list$z
  ZZadMAL[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  ZZadMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  adMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:2){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZadMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        adMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  adMAL_trapBuf_AllYears <- apply(adMAL_trapBuf, 3, rowSums)
  
  
  # 3. Estimate abundance proportion of each age class 
  
  prop_AllYears <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], 5, length(zones))) # Dim1 = iterations, Dim2 = c(cub, subad, ad, adFEM, adMAL), Dim3 = Core/Periphery
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      prop_AllYears[ite, 1, s] <- cub_trapBuf_AllYears[ite,s]/NIn_trapBuf_AllYears[ite,s]
      prop_AllYears[ite, 2, s] <- subad_trapBuf_AllYears[ite,s]/NIn_trapBuf_AllYears[ite,s]
      prop_AllYears[ite, 3, s] <- ad_trapBuf_AllYears[ite,s]/NIn_trapBuf_AllYears[ite,s]
      prop_AllYears[ite, 4, s] <- adFEM_trapBuf_AllYears[ite,s]/NIn_trapBuf_AllYears[ite,s]
      prop_AllYears[ite, 5, s] <- adMAL_trapBuf_AllYears[ite,s]/NIn_trapBuf_AllYears[ite,s]
    } } 
  
  results <- list(prop_AllYears = prop_AllYears, N = NIn_trapBuf_AllYears, cub = cub_trapBuf_AllYears, subad = subad_trapBuf_AllYears, ad = ad_trapBuf_AllYears, adfem = adFEM_trapBuf_AllYears, admal = adMAL_trapBuf_AllYears)
  return(results)
  
}

# Run function

age_peri1_allY <- age_st_allYears(nucleus = pCoreMedian, periphery = pPerMedian)
age_peri2_allY <- age_st_allYears(nucleus = pCoreDeciles, periphery = pPerDeciles)

# Estimate Mean and SD

mean_prop1_all <- round(apply(age_peri1_allY$prop, 3, colMeans),3) #### pMedian
colnames(mean_prop1_all) <- c("Core", "Periphery")
rownames(mean_prop1_all) <- c("cub", "subad", "ad", "adFEM", "adMALE")

mean_prop2_all <- round(apply(age_peri2_allY$prop, 3, colMeans),3) #### pDeciles
colnames(mean_prop2_all) <- c("Core", "Periphery")
rownames(mean_prop2_all) <- c("cub", "subad", "ad", "adFEM", "adMALE")

