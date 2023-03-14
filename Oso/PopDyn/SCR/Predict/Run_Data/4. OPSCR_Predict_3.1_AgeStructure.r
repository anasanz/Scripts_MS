
## -------------------------------------------------
##                 Plot predictions
## ------------------------------------------------- 

rm(list = ls())

library(nimbleSCR)
library(nimble)
library(rgdal)

source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")


# Load results from model (to get 5th year)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_param.RData")
sampmat1 <- do.call(rbind, nimOutput)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")
sampmat2 <- do.call(rbind, nimOutputSXY)

sampmat <- cbind(sampmat1, sampmat2)

# Constants:
M.aug <- 300 # Augmented individuals estimation model
Tt <- 5 # Nyears estimation model (2017:2021)

M.new <- 800 # New augmentation limit future prediction
t.new <- 5 # Extra years future prediction

# Load results from predictions

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter")
load("proj_pcr.core1-5000.RData")
z.proj.core <- z.proj.core2
sxy.proj.core <- sxy.proj.core2
age.cat.proj.core <- age.cat.proj.core2


load("proj_pcr.all1-5000.RData")
z.proj.all <- z.proj.all2
sxy.proj.all <- sxy.proj.all2
age.cat.proj.all <- age.cat.proj.all2


# Load buffer core area and habitat grid (to subset in sampling area)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")


## ---- 1. Join all results (past+future) in a single matrix with simlist format ----

# I add another dimension for the number of predictions
ndim <- 2 # 1 with pcr.core, 2 with pcr.all

sxy.allyears <- array(NA, c(dim(z.proj.core)[1], M.new, 2, Tt + t.new, ndim)) # Dataframes all years to store
z.allyears <- age.cat.allyears <- array(NA, c(dim(z.proj.core)[1], M.new, Tt + t.new, ndim))

# Different number of iterations. Predictions only from iter 1 to 5000, but it is in order.
# So the first 5000 values of the matrix of past years are equivalent to the 5000 iterations of the projection.

# Fill in data from past years

for(n in 1:ndim){ # Same for both dimensions of predictions (past years remain the same)
  for(t in 1:Tt){
    for(ite in 1:dim(z.proj.core)[1]){
      # 300 augmented individuals. Store in first five years
      sxy.allyears[ite,1:M.aug,,t,n] <- myResultsSXYZ$sims.list$sxy[ite,,,t] # Only 300 augmented individuals from past
      z.allyears[ite,1:M.aug,t,n] <- myResultsSXYZ$sims.list$z[ite,,t]
      age.cat.allyears[ite,1:M.aug,t,n] <- myResultsSXYZ$sims.list$age.cat[ite,,t]
    }}}

# Fill in data of future years
## PREDICTION 1: PCR ESTIMATED IN THE CORE (Dimension 1)
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.core)[1]){
    # 800 augmented individuals. We take year 5 from results estimation model. Store in 5 last years
    sxy.allyears[ite,,,(t+Tt),1] <- sxy.proj.core[ite,,,(t+1)] 
    z.allyears[ite,,(t+Tt),1] <- z.proj.core[ite,,(t+1)]
    age.cat.allyears[ite,,(t+Tt),1] <- age.cat.proj.core[ite,,(t+1)]
  }}

## PREDICTION 2: PCR ESTIMATED IN ALL STATE SPACE (Dimension 2)
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.all)[1]){
    # 800 augmented individuals. We take year 5 from results estimation model. Store in 5 last years
    sxy.allyears[ite,,,(t+Tt),2] <- sxy.proj.all[ite,,,(t+1)] 
    z.allyears[ite,,(t+Tt),2] <- z.proj.all[ite,,(t+1)]
    age.cat.allyears[ite,,(t+Tt),2] <- age.cat.proj.all[ite,,(t+1)]
  }}

## ---- 2. Age categories ----
# Do this for pcr estimated within core buffer area (Prediction 1)

# Unscale sxy

dimnames(sxy.allyears)[[3]] <- c('x','y')
sxy.allyears.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.allyears,## this are your sxy
                                             coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                             scaleToGrid = FALSE)$coordsDataScaled

## ---- 2.1. CUBS ----
ZZcubs <- z.allyears[,,,1]
ZZcubs[!age.cat.allyears[,,,1] %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

cub_trapBuf <- matrix(NA,nrow=dim(ZZcubs)[1],ncol=dim(ZZcubs)[3]) # nrow = iterations, ncol = year
cub_allss <- matrix(NA,nrow=dim(ZZcubs)[1],ncol=dim(ZZcubs)[3])

for(ite in 1:dim(ZZcubs)[1]){
  for(t in 1:dim(ZZcubs)[3]){
    
    which.alive <- which(ZZcubs[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,1] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    cub_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    cub_allss[ite,t] <- length(which.alive) # The sum of the points alive is the abundance that year and iteration in the state space. Store
    
  }}


## ---- 2.2. SUBADULTS ----

ZZsub <- z.allyears[,,,1]
ZZsub[!age.cat.allyears[,,,1] %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

subad_trapBuf <- matrix(NA,nrow=dim(ZZsub)[1],ncol=dim(ZZsub)[3]) # nrow = iterations, ncol = year
subad_allss <- matrix(NA,nrow=dim(ZZcubs)[1],ncol=dim(ZZcubs)[3])

for(ite in 1:dim(ZZsub)[1]){
  for(t in 1:dim(ZZsub)[3]){
    
    which.alive <- which(ZZsub[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,1] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    subad_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    subad_allss[ite,t] <- length(which.alive)
    }}


## ---- 2.3. ADULTS ----

ZZad <- z.allyears[,,,1]
ZZad[!age.cat.allyears[,,,1] %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

ad_trapBuf <- matrix(NA,nrow=dim(ZZad)[1],ncol=dim(ZZad)[3]) # nrow = iterations, ncol = year
ad_allss <- matrix(NA,nrow=dim(ZZad)[1],ncol=dim(ZZad)[3])

for(ite in 1:dim(ZZad)[1]){
  for(t in 1:dim(ZZad)[3]){
    
    which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,1] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    ad_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    ad_allss[ite,t] <- length(which.alive)
    }}


## ---- 3. Estimate abundance proportion of each age class ----

prop_years <- list()

prop <- matrix(NA,nrow=dim(ZZad)[1],ncol=3) # nrow = iterations, ncol = age categories
prop <- array(NA,c(dim(ZZad)[1],3,2))

ite = 1
t = 1

for(t in 1:dim(ZZad)[3]){
  for(ite in 1:dim(ZZad)[1]){
    prop[ite, 1] <- cub_trapBuf[ite,t]/NIn_trapBuf[ite,t]
    prop[ite, 2] <- subad_trapBuf[ite,t]/NIn_trapBuf[ite,t]
    prop[ite, 3] <- ad_trapBuf[ite,t]/NIn_trapBuf[ite,t]
  }
  prop_years[[t]] <- prop
}

## ---- 3.5. Plot ----

source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
years <- c(2017, 2018, 2019, 2020, 2021)

setwd("D:/MargSalas/Oso/Results/Plots/model3.1")
pdf("age_structure.pdf", 9, 7)

par(mfrow = c(2,3))

for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
  
  plot(1, ylim = c(0.5, ncol(prop_years[[t]])+0.5), 
       xlim = c(0,1), 
       type ="n", yaxt="n", 
       #xaxt="n", 
       xlab = " ", ylab = "", main = years[t],
       cex.axis = 0.8)
  
  axis(2, c(1:ncol(prop_years[[t]])), labels = c("Cub","Subadult","Adult"), las = 2, cex.axis = 1)
  
  for(i in 1:ncol(prop_years[[t]])){
    plot.violins3(list(prop_years[[t]][ ,i]),
                  x = i,
                  at = i,
                  violin.width = 0.2,
                  plot.ci = 0.95,
                  col = c("darksalmon"),
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)}}

dev.off()

# Plot bonito

setwd("D:/MargSalas/Oso/Results/Plots/model3.1")
pdf("age_structure2.pdf", 6, 10)

par(mfrow = c(1,1))

plot(1, ylim = c(0, 19 + 0.5), 
     xlim = c(0,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "",
     cex.axis = 0.8)

axis(2, c(1:20), las = 2, cex.axis = 1)

cubs_at <- c(1,5,9,13,17)
subad_at <- c(2,6,10,14,18)
ad_at <- c(3,7,11,15,19)
ats <- list(cubs_at, subad_at, ad_at)

polygon(x = c(0,0,1,1), y = c(0,4,4,0), col = adjustcolor("grey", alpha.f = 0.5), border = NA)

color_cat1 <- c("lightgoldenrod1", "darkolivegreen4", "darkslategrey")
color_cat2 <- c("lightgoldenrod1", "darkgoldenrod1", "sienna4") # no mucho
color_cat3 <-c("aquamarine1", "aquamarine4","darkslategrey")
color_cat4 <- c("aquamarine1","lightseagreen", "darkslategrey")
color_cat5 <- c("lightgoldenrod1", "darksalmon", "sienna4") # no mucho

for(i in 1:ncol(prop_years[[1]])){ # Look into cubs
  
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    plot.violins3(list(prop_years[[t]][ ,i]),
                  x = i,
                  at = ats[[i]][t],
                  violin.width = 0.4,
                  plot.ci = 0.95,
                  col = color_cat4[i],
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)}}
dev.off()



## ---- 3.6. Compare with age structure of data ----
# This is to see where the model placed the undetected individuals (if is especially in one age cat)

zdatAGE.det <- zdatAGE[1:61,]
zdatAGE.det[is.na(zdatAGE.det)] <- 0
zdatAGE


prop_det <- data.frame(matrix(NA, nrow = 5, ncol = 3))
rownames(prop_det) <- c("2017", "2018", "2019", "2020", "2021")
colnames(prop_det) <- c("Cub", "Subadult", "Adult")
t = 1

for(t in 1:5){
  alive.age <- age.cat.z[,t]*zdatAGE.det[,t]
  prop_det[t,1] <- length(alive.age[which(alive.age == 3 | alive.age == 2)])/length(alive.age[which(alive.age != 0)])
  prop_det[t,2] <- length(alive.age[which(alive.age == 5 | alive.age == 4)])/length(alive.age[which(alive.age != 0)])
  prop_det[t,3] <- length(alive.age[which(alive.age == 6)])/length(alive.age[which(alive.age != 0)])
}
rowSums(prop_det)

