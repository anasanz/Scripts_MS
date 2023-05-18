
## -------------------------------------------------
##           SPATIAL NUCLEUS vS PERIPHERY
##              FIRST TRIAL: POLYGONS 
## ------------------------------------------------- 

# POST-HOC ANALYSIS

# Identify the spatial nucleus (core) of the population and the periphery, to estimate the age structure in both areas
# The age structure will be estimated as the number of individuals of each age class in nucleus and periphery

# The core and periphery need to occupy the sampling area (not the state space), so I need to check how to define it
# - 1º trial: Buffer of female average dispersal distance around the core polygons (origin) for the distance to core covariate
# - 2º trial: Classification with distance to core covariate (e.g., take pixels with highest and lowest Distance to Core values)

rm(list = ls())

library(rgdal)
library(dplyr)
library(rgeos)
library(MCMCvis)
library(rgdal)
library(nimbleSCR)
library(raster)

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# Load Dispersal distances
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
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
Xbuf2 <- readOGR("Buffer_8500_traps.shp") # Load sampling area (where we estimate abundance)
proj4string(Xbuf2) <- proj4string(corepol)

# Load distCore (ONLY FOR VISUALIZATION PURPOSES HERE)
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
logDistcore <- raster("logDistcore_hrbear.tif")

## ---- 1. Buffer of disp. distances arounf core area----
nuc <- gBuffer(corepol, width = distfem, byid = TRUE)

plot(Xbuf)
plot(Xbuf2, col = "grey", border = NA, add = TRUE)
plot(corepol, col = "red", add = TRUE)
plot(nuc, border = "red", add = TRUE)

# To compare core area vS periphery, I need to create the periphery polygon, extracting the core area.
# I will test the age structure defining core vS periphery in 2 different ways

## ---- 1.1. Core (SMALL CORE) vS Trap buffer ----

# Nucleous: Polygon used to define core area in DistCore variable (corepol)
# Periphery: Sampling buffer (extracting the corepol): Xbuf2_peri1

Xbuf2_peri1 <- gDifference(Xbuf2, corepol)

## ---- 1.2. Intersection Disp distance Nucleous and Core (BIG CORE) vS Trap buffer ----

# Nucleous: Polygon used to define core area in DistCore variable, buffered by female dispersal distance and cropped so that it falls inside buffer (nuc2)
# Periphery: Sampling buffer (extracting the nuc): Xbuf2_peri2

nuc2 <- crop(nuc, Xbuf2)
Xbuf2_peri2 <- gDifference(Xbuf2, nuc2)

####### Plot density values in core vS periphery (code below): density_scale_checkBUFFERS.pdf ######
# Here I see:
# - Using small core is probably too restrictive
# - The occidental nucleus does not have reproductions, therefore will push the proportion of cubs in the core down (as it has mostly adults).
#   If we don't be consider it a reproductive core, the proportions might be more realistic.
#   or maybe we can give the results with both (Occidental + Central and Only central)

## ---- 1.3. Intersection Disp distance Nucleous and Core, only the central nucleus (BIG CENTRAL CORE) vS Trap buffer ----

# Nucleous: Polygon used to define core area in DistCore variable, buffered by female dispersal distance and cropped so that it falls inside buffer, select only central/reproductive core (nuc3)
# Periphery: Sampling buffer (extracting the nuc): Xbuf2_peri2

coords <- nuc2@polygons[[1]]@Polygons[[1]]@coords
poly1 <- sp::Polygon(cbind(coords[,1],coords[,2]))
firstPoly <- sp::Polygons(list(poly1), ID = 1)
nuc3 <- sp::SpatialPolygons(list(firstPoly))

Xbuf2_peri3 <- gDifference(Xbuf2, nuc3)


par(mfrow = c(1,3))

plot(Xbuf)
plot(logDistcore, add = TRUE) # Check
plot(corepol, add = TRUE)
plot(Xbuf2_peri1, col = "red", add = TRUE)
mtext("Small core", side = 3, line = -3, cex = 2)

plot(Xbuf)
plot(logDistcore, add = TRUE) # Check
plot(nuc2, add = TRUE)
plot(Xbuf2_peri2, col = "red", add = TRUE)
mtext("Big core", side = 3, line = -3, cex = 2)

plot(Xbuf)
plot(logDistcore, add = TRUE) # Check
plot(nuc3, add = TRUE)
plot(Xbuf2_peri3, col = "red", add = TRUE)
mtext("Big central core", side = 3, line = -3, cex = 2)



## ---- 2. Estimation of age structure in nuc vS periphery, different definitions from - to + restrictive ----

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

## ---- 2.1. Proportion of each age class relative to the zone (core or periphery) PER YEAR ----

# I think it makes more sense than in the whole polygon, otherwise it is hard to compare

# Function for estimating age structure in two polygons: nucleous and periphery

age_st <- function(nucleus, periphery){
  
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
  
  colMeans(NIn_trapBuf[,,1]) 
  colMeans(NIn_trapBuf[,,2])
  
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
  
  colMeans(cub_trapBuf[,,1])
  colMeans(cub_trapBuf[,,2])
  
  
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
  
  colMeans(subad_trapBuf[,,1])
  colMeans(subad_trapBuf[,,2])
  
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
  
  colMeans(adFEM_trapBuf[,,1])
  colMeans(adFEM_trapBuf[,,2])
  
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
  
  colMeans(adMAL_trapBuf[,,1])
  colMeans(adMAL_trapBuf[,,2])
  
  # 3. Estimate abundance proportion of each age class 
  
  prop_years <- list()
  prop_years_zone <- list()
  
  prop <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], 5, dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # Dim1 = iterations, Dim2 = c(cub, subad, ad, adFEM, adMAL), Dim3 = Years, Dim4 = Core/Periphery
  
  for(s in 1:length(zones)){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      prop[ite, 1, t, s] <- cub_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
      prop[ite, 2, t, s] <- subad_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
      prop[ite, 3, t, s] <- ad_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
      prop[ite, 4, t, s] <- adFEM_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
      prop[ite, 5, t, s] <- adMAL_trapBuf[ite,t,s]/NIn_trapBuf[ite,t,s]
    } } }
  
  results <- list(prop = prop, N = NIn_trapBuf, cub = cub_trapBuf, subad = subad_trapBuf, ad = ad_trapBuf, adfem = adFEM_trapBuf, admal = adMAL_trapBuf)
  return(results)
  
}

# Run function

age_peri1 <- age_st(nucleus = corepol, periphery = Xbuf2_peri1)
age_peri2 <- age_st(nucleus = nuc2, periphery = Xbuf2_peri2)
age_peri3 <- age_st(nucleus = nuc3, periphery = Xbuf2_peri3)


# Estimate Mean and SD
#### SMALL CORE AND PERIPHERY
mean_prop1 <- array(NA, c(dim(myResultsSXYZ$sims.list$z)[3], 5, 2)) # Dim1 = years, Dim2 = Ages, Dim3 = Zone 

for(s in 1:2){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    mean_prop1[t,,s] <- colMeans(age_peri1$prop[,,t,s])
  }
}

#### BIG CORE AND PERIPHERY 
mean_prop2 <- array(NA, c(dim(myResultsSXYZ$sims.list$z)[3], 5, 2)) # Dim1 = years, Dim2 = Ages, Dim3 = Zone 

for(s in 1:2){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    mean_prop2[t,,s] <- colMeans(age_peri2$prop[,,t,s])
  }
}

 
#### BIG CENTRAL CORE AND PERIPHERY 
mean_prop3 <- array(NA, c(dim(myResultsSXYZ$sims.list$z)[3], 5, 2)) # Dim1 = years, Dim2 = Ages, Dim3 = Zone 

for(s in 1:2){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    mean_prop3[t,,s] <- colMeans(age_peri3$prop[,,t,s])
  }
}


## ---- 2.2. Proportion of each age class relative to the zone (core or periphery) ALL YEARS TOGETHER ----

# I think it makes more sense than in the whole polygon, otherwise it is hard to compare

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

age_peri1_allY <- age_st_allYears(nucleus = corepol, periphery = Xbuf2_peri1)
age_peri2_allY <- age_st_allYears(nucleus = nuc2, periphery = Xbuf2_peri2)
age_peri3_allY <- age_st_allYears(nucleus = nuc3, periphery = Xbuf2_peri3)

# Estimate Mean and SD

mean_prop1_all <- round(apply(age_peri1_allY$prop, 3, colMeans),3) #### SMALL CORE AND PERIPHERY
colnames(mean_prop1_all) <- c("Core", "Periphery")
rownames(mean_prop1_all) <- c("cub", "subad", "ad", "adFEM", "adMALE")

mean_prop2_all <- round(apply(age_peri2_allY$prop, 3, colMeans),3)#### BIG CORE AND PERIPHERY 
colnames(mean_prop2_all) <- c("Core", "Periphery")
rownames(mean_prop2_all) <- c("cub", "subad", "ad", "adFEM", "adMALE")

mean_prop3_all <- round(apply(age_peri3_allY$prop, 3, colMeans),3) #### BIG CENTRAL CORE AND PERIPHERY 
colnames(mean_prop3_all) <- c("Core", "Periphery")
rownames(mean_prop3_all) <- c("cub", "subad", "ad", "adFEM", "adMALE")

## ---- Results and conclusion ----

#### SMALL CORE AND PERIPHERY
mean_prop1
mean_prop1_all

#### BIG CORE AND PERIPHERY 
mean_prop2
mean_prop2_all

#### BIG CENTRAL CORE AND PERIPHERY 
mean_prop3
mean_prop3_all

#### CONCLUSIONS: 
# - The number of adults is very high in the core, so the proportion of cubs is not really higher in the core.
# - The more realistic proportion is mean_prop2, because it takes a larger area at the core and therefore integrates all cubs
# - mean_prop3 seemed good but it is not (its similar), I think because the density of both adults and cubs is lower in the occidental nucleus.

## ---- Check absolute number in Big core (it seems the most realistic) ----

dim(age_peri2$cub)

cubs <- round(apply(age_peri2$cub, 3, colMeans),3)
colnames(cubs) <- c("Core", "Periphery")
rownames(cubs) <- c("2017", "2018", "2019", "2020", "2021")

subad <- round(apply(age_peri2$subad, 3, colMeans),3) 
ad <- round(apply(age_peri2$ad, 3, colMeans),3) 
adFEM <- round(apply(age_peri2$adfem, 3, colMeans),3) 
adMAL <- round(apply(age_peri2$admal, 3, colMeans),3)


## -------------------------------------------------
##    APPENDIX OF SCRIPT: Plot density values each age category with all buffers on top
## ------------------------------------------------- 

#rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(rgdal)
library(raster)
library(viridis)

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")

# Load functions
#sourceCpp("GetSpaceUse_PD.cpp")
sourceCpp("GetDensity_PD.cpp")
source("getDensityInput.R")

# Load buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

# Conver to rasters

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
r <- raster("logDistcore_hrbear.tif")
raster::values(r) <- NA
r <- crop(r,Xbuf) 

Xbuf_raster <- rasterize(Xbuf, r)
raster::values(Xbuf_raster)[which(is.na(raster::values(Xbuf_raster)))] <- 0 # Raster of 0 and 1

Xbuf2_raster <- rasterize(Xbuf2, r)
raster::values(Xbuf2_raster)[which(raster::values(Xbuf2_raster) == 1)] <- 2
raster::values(Xbuf2_raster)[which(is.na(raster::values(Xbuf2_raster)))] <- 0 # Raster of 0 and 2

ras <- overlay(Xbuf_raster, Xbuf2_raster, fun = max)
f <- rasterize(Xbuf, ras, mask = TRUE)
#values(f)[which(values(f) == 1)] <- 0
#values(f)[which(values(f) == 2)] <- 1

# Load posterior distribution
library(nimbleSCR) # Load nimbleSCR here, otherwise it gets in conflict with raster package, weird

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

# Load original habitat coordinates
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

#give dim names to your posteriors sxy
dim(myResultsSXYZ$sims.list$sxy) # YOURPOSTERIORSXY 
dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c("x","y")

## first rescale the coordinates to the original scale 
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,
                                                        coordsHabitatGridCenter = G,
                                                        scaleToGrid = FALSE)$coordsDataScaled

f1 <- f
f1[f1%in%c(1,2)] <-  1
f[] <- as.factor(as.character(f[]))

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = myResultsSXYZ$sims.list$sxy,
  plot.check = TRUE)

## extract density
yearnames <- c("2017", "2018", "2019", "2020", "2021")
leg <- c(FALSE, FALSE, FALSE, FALSE, TRUE)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("density_scale_checkBUFFERS.pdf",12,6)


######  ALL INDIVIDUALS  #####

par(mfrow = c(4,5),
    mar = c(0,1,0,1),
    oma = c(0.5,1.5,2,2),
    bty = "n")

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = myResultsSXYZ$sims.list$z[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = leg[t])
  plot(Xbuf2, add = TRUE)
  plot(corepol, border = "red", add = TRUE)
  plot(nuc, border = "blue", add = TRUE)
  
  
}


######  CUBS  #####

ZZcubs <- myResultsSXYZ$sims.list$z
ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZcubs[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = leg[t])
  plot(Xbuf2, add = TRUE)
  plot(corepol, border = "red", add = TRUE)
  plot(nuc, border = "blue", add = TRUE)
}

######  SUBADULTS  #####

ZZsub <- myResultsSXYZ$sims.list$z
ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZsub[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = leg[t])
  plot(Xbuf2, add = TRUE)
  plot(corepol, border = "red", add = TRUE)
  plot(nuc, border = "blue", add = TRUE)
}


######  ADULTS  #####

ZZad <- myResultsSXYZ$sims.list$z
ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZad[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = leg[t])
  plot(Xbuf2, add = TRUE)
  plot(corepol, border = "red", add = TRUE)
  plot(nuc, border = "blue", add = TRUE)
}

mtext(yearnames, at = c(0.1,0.3,0.5,0.7,0.9), outer = TRUE, line = -1, side = 3)
mtext(c("Adult", "Subadult", "Cub", "All"), at = c(0.17,0.45,0.67,0.92), outer = TRUE, line = 0, side = 2, adj = 1)

dev.off()


## -------------------------------------------------
##      APPENDIX SCRIPT: PLOT SEX PROPORTION    
## ------------------------------------------------- 


sa <- readOGR("D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2_WGS84_31N_all")
eur <- readOGR("D:/MargSalas/Oso/Datos/GIS/Countries", "esp_fr_2")
eur <- spTransform(eur, crs(sa))

par(mfrow = c(1,1))
plot(Xbuf, col = NA, border = NA)
plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2))
plot(nuc2, col = "darkorchid4", border = "darkorchid4", add = TRUE)
plot(Xbuf2_peri2, col = "darkseagreen4", border = "darkseagreen4", add = TRUE)
plot(eur, xlim = c(bbox(sa)[1,1], bbox(sa)[1,2]), ylim = c(bbox(sa)[2,1], bbox(sa)[2,2]), border = adjustcolor("black", alpha.f = 0.2), add = TRUE)

plot(logDistcore, add = TRUE) # Check
mtext("Big core", side = 3, line = -3, cex = 2)




