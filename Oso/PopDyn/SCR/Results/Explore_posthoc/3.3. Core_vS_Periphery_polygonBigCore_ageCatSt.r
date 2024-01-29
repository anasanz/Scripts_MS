## ------------------------------------------------------------------------------------
##                      SPATIAL NUCLEUS vS PERIPHERY
##  Definition BIG CORE: Polygon female dispersal distance around original core
## ------------------------------------------------------------------------------------

# POST-HOC ANALYSIS

# Identify the spatial nucleus (core) of the population and the periphery, to estimate the age structure in both areas
# The age structure will be estimated as the number of individuals of each age class in nucleus and periphery


rm(list = ls())

## ---- Load stuff ----

library(rgdal)
library(dplyr)
library(rgeos)
library(MCMCvis)
library(rgdal)
library(nimbleSCR)
library(raster)
library(sf)

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# Load Dispersal distances
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
d <- read.csv("disp_distance.csv", header = TRUE, row.names = NULL, sep = ",")
mean_dist <- d %>% # Calculate dispersal distance by sex
  group_by(Sex) %>%
  summarise(
    count = n(),
    mean = mean(Distance, na.rm = TRUE),
    sd = sd(Distance, na.rm = TRUE),
    se = stderr(Distance, na.rm = TRUE)
  )

distfem <- mean_dist$mean[1]

# Load core polygon
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol")
corepol <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol/corepol.shp", "corepol")
corepol@data$ID <- 1
corepol@data <- corepol@data[2]

# Load state space
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
#Xbuf2 <- readOGR("Buffer_8500_traps.shp") # Load sampling area (where we estimate abundance)
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array
proj4string(Xbuf2) <- proj4string(corepol)

# Load distCore (ONLY FOR VISUALIZATION PURPOSES HERE)
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
logDistcore <- raster("logDistcore_hrbear.tif")

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

## ---- 1. Define core and periphery (big core polygon) ----

# Buffer of disp. distances arounf core area
nuc <- gBuffer(corepol, width = distfem, byid = TRUE)

# Intersection Disp distance Nucleous and Core (BIG CORE) vS Trap buffer 
# Nucleous: Polygon used to define core area in DistCore variable, buffered by female dispersal distance and cropped so that it falls inside buffer (nuc2)
# Periphery: Sampling buffer (extracting the nuc): Xbuf2_peri2
nuc2 <- crop(nuc, Xbuf2)
Xbuf2_peri2 <- gDifference(Xbuf2, nuc2)
p.df <- data.frame( ID=1 )
Xbuf2_peri2 <- SpatialPolygonsDataFrame(Xbuf2_peri2, p.df) 

plot(Xbuf)
plot(Xbuf2, add = TRUE)
plot(Xbuf2_peri2, col = adjustcolor("green", alpha = 0.3), add = TRUE)
plot(nuc2, col = adjustcolor("red", alpha = 0.3), add = TRUE)


setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
#writeOGR(nuc2, "D:/MargSalas/Oso/Datos/GIS/Countries/nuc.shp", layer = "nuc", driver = "ESRI Shapefile")
#writeOGR(Xbuf2_peri2, "D:/MargSalas/Oso/Datos/GIS/Countries/per.shp", layer = "per", driver = "ESRI Shapefile")

# Load again as st
nuc <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/nuc.shp")
per <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/per.shp")

st_area(nuc)/st_area(per) # Core is almost half of periphery
area <- c(st_area(nuc), st_area(per))

area <- units::set_units(area, km^2)

## ---- 2. Estimation of age structure in nuc vS periphery per sex and year ----
## ---- 2.1. Abundance ----

age_st_sex <- function(nucleus, periphery){
  
  # 1. Total abundance in each zone
  
  zones <- c(nucleus, periphery)
  
  # Matrix to store abundance in the buffer each iteration and year
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
  
  
  # 2. Abundance of each age category and sex in each zone
  
  ## FEMALE ALL
  
  ZZFEM <- myResultsSXYZ$sims.list$z
  ZZFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  allFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        allFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## FEMALE CUBS
  ZZcubsFEM <- myResultsSXYZ$sims.list$z
  ZZcubsFEM[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZcubsFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  
  cubFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZcubsFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        cubFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}

  
  ## FEMALE SUBADULTS 
  
  ZZsubFEM <- myResultsSXYZ$sims.list$z
  ZZsubFEM[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  ZZsubFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  subadFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZsubFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        subadFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## FEMALE ADULTS
  
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
  
  ## MALE ALL
  ZZMAL <- myResultsSXYZ$sims.list$z
  ZZMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  allMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        allMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE CUBS
  
  ZZcubsMAL <- myResultsSXYZ$sims.list$z
  ZZcubsMAL[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZcubsMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  cubMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZcubsMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        cubMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE SUBADULTS
  
  ZZsubMAL <- myResultsSXYZ$sims.list$z
  ZZsubMAL[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  ZZsubMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  subadMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZsubMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        subadMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}

  
  ## MALE ADULTS
  
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
  
  
  # 3. Estimate abundance proportion of each age class 
  
  prop_years <- list()
  prop_years_zone <- list()
  
  prop <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], 8, dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # Dim1 = iterations, Dim2 = c(cub, subad, ad, adFEM, adMAL), Dim3 = Years, Dim4 = Core/Periphery
  
  for(s in 1:length(zones)){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
      for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
        prop[ite, 1, t, s] <- cubFEM_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        prop[ite, 2, t, s] <- cubMAL_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        prop[ite, 3, t, s] <- subadFEM_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        prop[ite, 4, t, s] <- subadMAL_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        prop[ite, 5, t, s] <- adFEM_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        prop[ite, 6, t, s] <- adMAL_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        prop[ite, 7, t, s] <- allFEM_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        prop[ite, 8, t, s] <- allMAL_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
        
      } } }
  
  results <- list(prop = prop, N = NIn_trapBuf, cubFEM = cubFEM_trapBuf, cubMAL = cubMAL_trapBuf, subadFEM = subadFEM_trapBuf, subadMAL = subadMAL_trapBuf, adFEM = adFEM_trapBuf, adMAL = adMAL_trapBuf, allFEM = allFEM_trapBuf, allMAL = allMAL_trapBuf) 
                  #cub = cub_trapBuf, subad = subad_trapBuf, ad = ad_trapBuf, adfem = adFEM_trapBuf, admal = adMAL_trapBuf)
  return(results)
  
}
ageSt_bear <- age_st_sex(nucleus = nuc2, periphery = Xbuf2_peri2)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
save(ageSt_bear, file = "Nbuffer_newSize.RData")

## ---- 2.2. Density ----

agesex_DENS <- function(nucleus, periphery, area){
  
  # 1. Total abundance in each zone
  
  zones <- c(nucleus, periphery)
  
  # Matrix to store abundance in the buffer each iteration and year
  NIn_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(myResultsSXYZ$sims.list$z[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        NIn_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
        
      }
    }
  }
  
  
  # 2. Abundance of each age category and sex in each zone
  
  ## FEMALE ALL
  
  ZZFEM <- myResultsSXYZ$sims.list$z
  ZZFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  allFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        allFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## FEMALE CUBS
  ZZcubsFEM <- myResultsSXYZ$sims.list$z
  ZZcubsFEM[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZcubsFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  
  cubFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZcubsFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        cubFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## FEMALE SUBADULTS 
  
  ZZsubFEM <- myResultsSXYZ$sims.list$z
  ZZsubFEM[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  ZZsubFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  subadFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZsubFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        subadFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## FEMALE ADULTS
  
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
        
        adFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE ALL
  ZZMAL <- myResultsSXYZ$sims.list$z
  ZZMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  allMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        allMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE CUBS
  
  ZZcubsMAL <- myResultsSXYZ$sims.list$z
  ZZcubsMAL[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZcubsMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  cubMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZcubsMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        cubMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE SUBADULTS
  
  ZZsubMAL <- myResultsSXYZ$sims.list$z
  ZZsubMAL[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  ZZsubMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  subadMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZsubMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        subadMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## MALE ADULTS
  
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
        
        adMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T)/area[s] # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  results <- list(totDens = NIn_trapBuf, cubFEM = cubFEM_trapBuf, cubMAL = cubMAL_trapBuf, subadFEM = subadFEM_trapBuf, subadMAL = subadMAL_trapBuf, adFEM = adFEM_trapBuf, adMAL = adMAL_trapBuf, allFEM = allFEM_trapBuf, allMAL = allMAL_trapBuf) 
  #cub = cub_trapBuf, subad = subad_trapBuf, ad = ad_trapBuf, adfem = adFEM_trapBuf, admal = adMAL_trapBuf)
  return(results)
  
}
agesex_DENS_bear <- agesex_DENS(nucleus = nuc2, periphery = Xbuf2_peri2, area = area)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
save(agesex_DENS_bear, file = "Densbuffer_newSize.RData")

## ---- 3. Plots ----

rm(list = ls())

library(nimbleSCR)
library(nimble)
library(rgdal)

source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")
source("D:/MargSalas/Scripts_MS/Functions/PlotViolinsHoriz.r")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("Nbuffer_newSize.RData") # Obtained from script: 3.3. Core_vS_Periphery_polygonBigCore_ageCatSt.r

## ---- 3.1. Absolute number/Ageclass/year/sex ----

#colZn2 <- c("#8073ac", "#e08214")
colZn4 <- c("#9970ab", "#a6dba0")
#colZn5 <- c("#d6604d", "#4393c3")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("Core_Periphery_AgeSex.pdf", 5, 5)

par(mfrow = c(4,2),
    mar = c(1,1,1.3,1),
    oma = c(2,3,2,2))

ageSt_list <- ageSt_bear[-c(1,2)]
nAgeSex <- length(ageSt_list)

zones <- dim(ageSt_list[[1]])[3]
offset <- seq(0,0.25,length.out = zones)
at <- c(1,2,3,4,5)

leg <- c(FALSE, TRUE, rep(FALSE, 6))
tit1 <- c("Females", "Males", rep(FALSE, 6))
tit2 <- c("Cubs", FALSE, "Subadults", FALSE, "Adults", FALSE, "All", FALSE)


scaleY_plots <- c(rep(max(ageSt_list[[1]], ageSt_list[[2]]),2), rep(max(ageSt_list[[3]], ageSt_list[[4]], na.rm = TRUE),2),
                 rep(max(ageSt_list[[5]], ageSt_list[[6]]),2), rep(max(ageSt_list[[7]], ageSt_list[[8]]),2))

for (n in 1:nAgeSex){
  
  if(all(complete.cases(ageSt_list[[n]])) == FALSE){ # delete rows with NA if there are
    delete <- unique(c(which(!complete.cases(ageSt_list[[n]][,,1])), which(!complete.cases(ageSt_list[[n]][,,2]))))
    ageSt_list[[n]] <- ageSt_list[[n]][-delete,,]
  }
  

  plot(1, ylim = c(0,scaleY_plots[n]+1), 
       xlim = c(0.5, at[5] + max(offset) + 0.5)  , 
       type ="n", 
       yaxt="n", 
       xaxt="n", 
       xlab = " ", ylab = "", main = " ",
       cex.axis = 0.8, axes = FALSE)
  
  axis(1, c(1:ncol(ageSt_list[[n]])), labels = c(2017:2021), 
       at = at + max(offset)/2, las = 1, cex.axis = 0.75, pos = 0, lwd.ticks = 0.2, lwd = 0.2 )
  axis(2, cex.axis = 0.75, pos = 0.75, lwd.ticks = 0.2, lwd = 0.2)
  
  
  m <- data.frame(matrix(NA, nrow = 5, ncol = 2))     
  for(s in 1:zones){
    for (i in 1:5){
      plot.violins3(list(ageSt_list[[n]][ ,i,s]),
                    x = i,
                    at = at[i]+offset[s],
                    violin.width = 0.05,
                    plot.ci = 0.85,
                    col = colZn4[s],
                    add = T,
                    alpha = 0.8,
                    scale.width = FALSE,
                    border.col = colZn4[s],
                    horizontal = FALSE)
      m[i,s] <- median(ageSt_list[[n]][ ,i,s]) 
    }
  }
  
  for(s in 1:zones){
    lo <- loess(m[,s]~ c(at+offset[s]))
    lines(predict(lo)~ c(at+offset[s]), col=colZn4[s], lwd=1.5)
  }
  if(leg[n] == TRUE){
    legend("topright", inset = c(+0.1, -0.05), legend = c("Core", "Periphery"), fill = colZn4, border = NA, bty = "n", horiz = FALSE, cex = 0.9)
  }
  if(tit1[n] != FALSE){
    mtext(tit[n], side = 3, cex = 0.8, line = 1)
  }
  if(tit2[n] != FALSE){
    mtext(tit2[n], side = 2, cex = 0.75, line = 1.5)
  }
}

dev.off()

## ---- 3.2. Absolute number/Ageclass/ALL years/sex (Function Cyril) ----

ageSt_list <- ageSt_bear[-c(1,2)]
nAgeSex <- length(ageSt_list)
zones <- dim(ageSt_list[[1]])[3]

for(n in 1:length(ageSt_list)){
  if(all(complete.cases(ageSt_list[[n]])) == FALSE){ # delete rows with NA if there are
    delete <- unique(c(which(!complete.cases(ageSt_list[[n]][,,1])), which(!complete.cases(ageSt_list[[n]][,,2]))))
    ageSt_list[[n]] <- ageSt_list[[n]][-delete,,]
  }
}

ageSt_bear_F <- ageSt_list[c(1,3,5,7)]
ageSt_bear_M <- ageSt_list[c(2,4,6,8)]

lim.fem <- max(unlist(ageSt_bear_F), na.rm = TRUE)
lim.mal <- max(unlist(ageSt_bear_M), na.rm = TRUE)+5 # To make it more equal both sides
colZn4 <- c("#9970ab", "#a6dba0")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("Core_Periphery_AgeSex_Allyears.pdf", 7, 5)

par(mfrow = c(1,1), 
    mar = c(3, 3, 3, 6),
    xpd=TRUE)

plot(-1000, xlim=c(-lim.mal,lim.fem), ylim=c(0,(length(ageSt_bear_F)+1)),
     type ="n", 
     yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = " ",
     cex.axis = 0.8, axes = FALSE)

axis(1, seq(-25,+30, by = 5), labels = abs(seq(-25,+30, by = 5)), 
     at = seq(-25,+30, by = 5), las = 1, cex.axis = 0.75, pos = 0, lwd.ticks = 0.2, lwd = 0.2)


#abline(v=0)
segments(x0=0,y0=0,x1=0,y1=5)


for(i in 1:length(ageSt_bear_F)){
  
  
  plot.violinsHoriz(dat.list = list(apply(ageSt_bear_F[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T,horizontal = T,col="#9970ab",border.col = "#9970ab", alpha = 0.5)
  
  plot.violinsHoriz(dat.list =list(-apply(ageSt_bear_M[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T,horizontal = T,col="#9970ab",border.col = "#9970ab", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(apply(ageSt_bear_F[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T,horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(-apply(ageSt_bear_M[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T,horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
  
} 

text(-15,5, "Males")
text(15,5, "Females")
mtext(c("Cubs", "Subadults", "Adults", "All"), side = 2, at = c(1,2,3,4), las = 1, line = -4)
legend("topright", inset=c(-0.1,0), legend = c("Core", "Periphery"), fill = colZn4, border = NA)

dev.off()

## ---- 3.3. Proportion/Ageclass/sex ----

prop <- list(ageSt_bear$prop[,1,,], ageSt_bear$prop[,2,,], ageSt_bear$prop[,3,,], ageSt_bear$prop[,4,,],
             ageSt_bear$prop[,5,,], ageSt_bear$prop[,6,,], ageSt_bear$prop[,7,,], ageSt_bear$prop[,8,,])

for(n in 1:length(prop)){
  if(all(complete.cases(prop[[n]])) == FALSE){ # delete rows with NA if there are
    delete <- unique(c(which(!complete.cases(prop[[n]][,,1])), which(!complete.cases(prop[[n]][,,2]))))
    prop[[n]] <- prop[[n]][-delete,,]
  }
}

prop_F <- prop[c(1,3,5,7)]
prop_M <- prop[c(2,4,6,8)]
colZn4 <- c("#9970ab", "#a6dba0")


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("Core_Periphery_AgeSex_Allyears_PROP.pdf", 7, 7)

par(mfrow = c(1,1), 
    mar = c(3, 3, 3, 6),
    xpd=TRUE)

plot(-0.8, xlim=c(-0.8,1), ylim=c(0,(length(prop_F)+1)),
     type ="n", 
     yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = " ",
     cex.axis = 0.8, axes = FALSE)

axis(1, seq(-0.6,0.8, by = 0.2), labels = abs(seq(-0.6,0.8, by = 0.2)), 
     at = seq(-0.6,0.8, by = 0.2), las = 1, cex.axis = 0.75, pos = 0, lwd.ticks = 0.2, lwd = 0.2)


for(i in 1:length(prop_F)){
  
  polygon(x = c(0, median(apply(prop_F[[i]], c(1,3), mean)[,1]), median(apply(prop_F[[i]], c(1,3), mean)[,1]), 0),
          y = c(i,i,i+0.4,i+0.4), col = adjustcolor("#9970ab", alpha.f = 0.3), border = "#9970ab")

  plot.violinsHoriz(dat.list = list(apply(prop_F[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T, violin.width = 0.2, horizontal = T,col="#9970ab", border.col = "#9970ab", alpha = 0.5)
  

  polygon(x = c(0, -median(apply(prop_M[[i]], c(1,3), mean)[,1]), -median(apply(prop_M[[i]], c(1,3), mean)[,1]), 0),
          y = c(i,i,i+0.4,i+0.4), col = adjustcolor("#9970ab", alpha.f = 0.3), border = "#9970ab")
  
  plot.violinsHoriz(dat.list =list(-apply(prop_M[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T, violin.width = 0.2, horizontal = T,col="#9970ab",border.col = "#9970ab", alpha = 0.5)
  
  
  polygon(x = c(0, median(apply(prop_F[[i]], c(1,3), mean)[,2]), median(apply(prop_F[[i]], c(1,3), mean)[,2]), 0),
          y = c(i,i,i-0.4,i-0.4), col = adjustcolor("#a6dba0", alpha.f = 0.3), border = "#a6dba0")
  
  plot.violinsHoriz(dat.list =list(apply(prop_F[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T, violin.width = 0.2, horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
  
  
  polygon(x = c(0,-median(apply(prop_M[[i]], c(1,3), mean)[,2]), -median(apply(prop_M[[i]], c(1,3), mean)[,2]), 0),
          y = c(i,i,i-0.4,i-0.4), col = adjustcolor("#a6dba0", alpha.f = 0.3), border = "#a6dba0")
  
  plot.violinsHoriz(dat.list =list(-apply(prop_M[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T, violin.width = 0.1, horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
}

segments(x0=0,y0=0,x1=0,y1=5)

text(-0.4,5, "Males")
text(0.4,5, "Females")
mtext(c("Cubs", "Subadults", "Adults", "All"), side = 2, at = c(1,2,3,4), las = 1, line = -2)
legend("topright", inset=c(-0.2,0), legend = c("Core", "Periphery"), fill = colZn4, border = NA)

dev.off()

## ---- 3.4. Proportion/Ageclass/sex no bars ----

prop <- list(ageSt_bear$prop[,1,,], ageSt_bear$prop[,2,,], ageSt_bear$prop[,3,,], ageSt_bear$prop[,4,,],
             ageSt_bear$prop[,5,,], ageSt_bear$prop[,6,,], ageSt_bear$prop[,7,,], ageSt_bear$prop[,8,,])

for(n in 1:length(prop)){
  if(all(complete.cases(prop[[n]])) == FALSE){ # delete rows with NA if there are
    delete <- unique(c(which(!complete.cases(prop[[n]][,,1])), which(!complete.cases(prop[[n]][,,2]))))
    prop[[n]] <- prop[[n]][-delete,,]
  }
}
prop_F <- prop[c(7,5,3,1)]
prop_M <- prop[c(8,6,4,2)]
colZn4 <- c("#9970ab", "#a6dba0")


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("Core_Periphery_AgeSex_Allyears_PROP_nobars.pdf", 6, 3)

par(mfrow = c(1,1), 
    mar = c(3, 3, 3, 6),
    xpd=TRUE)

plot(-0.8, xlim=c(-0.8,1), ylim=c(0,(length(prop_F)+1)),
     type ="n", 
     yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = " ",
     cex.axis = 0.8, axes = FALSE)

axis(1, seq(-0.6,0.8, by = 0.2), labels = abs(seq(-0.6,0.8, by = 0.2)), 
     at = seq(-0.6,0.8, by = 0.2), las = 1, cex.axis = 1.2, pos = 0, lwd.ticks = 0.2, lwd = 0.2)


for(i in 1:length(prop_F)){
  
  
  plot.violinsHoriz(dat.list = list(apply(prop_F[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#9970ab", border.col = "#9970ab", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(-apply(prop_M[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#9970ab",border.col = "#9970ab", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(apply(prop_F[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(-apply(prop_M[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
}

segments(x0=0,y0=0,x1=0,y1=5)

text(-0.4,5, "Males")
text(0.4,5, "Females")
mtext(c("All", "Adults", "Subadults", "Cubs"), side = 2, at = c(1,2,3,4), las = 1, line = -2)
legend("topright", inset=c(-0.2,0), legend = c("Core", "Periphery"), fill = colZn4, border = NA)

dev.off()

## ---- 4. Study area and location of core/periphery ----

# Load europe
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = crs(sa))
and <- sa[which(sa$NAME_0 == "Andorra"),]
eur2 <- st_union(and,eur)

# Plot
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots")
pdf("sa_coreper2.pdf",7,7)

plot(st_geometry(eur2), col = "khaki", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000))
plot(Xbuf2, col = "#9970ab", border = "#9970ab", add = TRUE)
plot(Xbuf2_peri2, col = "#a6dba0", border = "#a6dba0", add = TRUE)
plot(st_geometry(eur2), border = adjustcolor("white", alpha.f = 0.3),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)

dev.off()
