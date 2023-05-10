
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

M.new <- 700 # New augmentation limit future prediction
t.new <- 5 # Extra years future prediction

# Load results from predictions

# ** IMPORTANT
# - The files proj_pcr.core and proj_pcr.all come from the projections of 5000 random iterations (vector "itera")
#   with a PCR estimated for all individuals (per capita recruitment per males and females)
# - The files proj_pcr.core.fem and proj_pcr.all.fem come from the projections of 5000 random iterations (vector "itera")
#   with a PCR estimated for ONLY FEMALES (per capita recruitment per female)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter")
#load("proj_pcr.core.RData")
#load("proj_pcr.all.RData")

load("proj_pcr.core.fem.RData")
load("proj_pcr.all.fem.RData")

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


## ---- 2. Subset abundance in buffer ----

# Unscale sxy

dimnames(sxy.allyears)[[3]] <- c('x','y')
sxy.allyears.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.allyears,## this are your sxy
                                              coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                              scaleToGrid = FALSE)$coordsDataScaled

# Array to store abundance in the buffer each iteration and year

#NIn <- matrix(NA,nrow = dim(z.allyears)[1], ncol=dim(z.allyears)[3]) # nrow = iterations, ncol = year
NIn <- array(NA, c(dim(z.allyears)[1], dim(z.allyears)[3], ndim))

for(n in 1:ndim){
  for(ite in 1:dim(z.allyears)[1]){
    for(t in 1:dim(z.allyears)[3]){
      
      which.alive <- which(z.allyears[ite,,t,n]==1) # Select only the individuals alive (z=1)
      
      which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,n] # Retrieve the activity center for those individuals
      
      sp <- SpatialPoints(which.aliveSXY, proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
      
      which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
      
      NIn[ite,t,n] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      
    }
  }
}

#average number of individuals without the buffer for each year 
colMeans(NIn[,,1])
colMeans(NIn[,,2])


# Sum of individuals alive in total each year (without buffer)
colMeans(apply(z.allyears[,,,1],c(1,3),function(x) sum(x==1, na.rm = TRUE)))
colMeans(apply(z.allyears[,,,2],c(1,3),function(x) sum(x==1, na.rm = TRUE)))


## ---- Plot ----

# Plot names:
# Running the script with proj_pcr.core.RData (or all) -> prediction_abundanceCore_pcr.allind.pdf
# Running the script with proj_pcr.core.fem.RData (or all) -> prediction_abundanceCore_pcr.fem.pdf

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1/Predictions")
pdf("prediction_abundanceCore_pcr.fem.pdf", 7, 5)

plot(1, ylim = c(30,max(NIn)+50), 
     xlim = c(0.5, ncol(NIn)+ncol(NIn)+0.5) , 
     type ="n", 
     #yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

#polygon(x = c(5.5,5.5,22,22), y = c(0,max(NIn)+150,max(NIn)+150,0), col = adjustcolor("grey", alpha.f = 0.2), border = NA)

axis(1, c(1:ncol(NIn)), labels = c(2017:2026), 
     at = c(1,2,3,4,5,7,10,13,16,19),las = 2, cex.axis = 1)

# First five years (present)
for (i in 1:5){
  plot.violins3(list(NIn[ ,i,1]),
                x = i,
                at = i,
                violin.width = 0.3,
                plot.ci = 0.95,
                col = c("yellow4"),
                add = T,
                alpha = 0.8,
                scale.width = FALSE,
                border.col = "yellow4",
                horizontal = FALSE)}

# Scenario 1 (pcr.core)
at.pcr.core = c(1,2,3,4,5,6.5,9.5,12.5,15.5,18.5)
for (i in 6:10){
  plot.violins3(list(NIn[ ,i,1]),
                x = i,
                at = at.pcr.core[i],
                violin.width = 0.3,
                plot.ci = 0.95,
                col = c("yellow3"),
                add = T,
                alpha = 0.8,
                scale.width = FALSE,
                border.col = "yellow3",
                horizontal = FALSE)}

# Scenario 2 (pcr.all)
at.pcr.all = at.pcr.core + 1
for (i in 6:10){
  plot.violins3(list(NIn[ ,i,2]),
                x = i,
                at = at.pcr.all[i],
                violin.width = 0.3,
                plot.ci = 0.95,
                col = c("wheat4"),
                add = T,
                alpha = 0.8,
                scale.width = FALSE,
                border.col = "wheat4",
                horizontal = FALSE)}


dev.off()


