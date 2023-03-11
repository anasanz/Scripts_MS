
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

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y')
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled

# Constants:
M.aug <- 300 # Augmented individuals estimation model
Tt <- 5 # Nyears estimation model (2017:2021)

M.new <- 800 # New augmentation limit future prediction
t.new <- 5 # Extra years future prediction

# Load results from predictions

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter")
#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter")
load("proj_pcr.core.RData")
z.proj.core <- z.proj.core2
sxy.proj.core <- sxy.proj.core2
age.cat.proj.core <- age.cat.proj.core2

dim(z.proj.core2)

# Load buffer core area and habitat grid (to subset in sampling area)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")


## ---- 1. Join all results (past+future) in a single matrix with simlist format ----

sxy.allyears <- array(NA, c(dim(z.proj.core)[1], M.new, 2, Tt + t.new)) # Dataframes all years to store
z.allyears <- age.cat.allyears <- array(NA, c(dim(z.proj.core)[1], M.new, Tt + t.new))

# Different number of iterations. Predictions only from iter 1 to 5000, but it is in order.
# So the first 5000 values of the matrix of past years are equivalent to the 5000 iterations of the projection.

# Fill in data from past years

for(t in 1:Tt){
  for(ite in 1:dim(z.proj.core)[1]){
    # 300 augmented individuals. Store in first five years
    sxy.allyears[ite,1:M.aug,,t] <- myResultsSXYZ$sims.list$sxy[ite,,,t] # Only 300 augmented individuals from past
    z.allyears[ite,1:M.aug,t] <- myResultsSXYZ$sims.list$z[ite,,t]
    age.cat.allyears[ite,1:M.aug,t] <- myResultsSXYZ$sims.list$age.cat[ite,,t]
  }}

# Fill in data of future years

for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.core)[1]){
    # 800 augmented individuals. We take year 5 from results estimation model. Store in 5 last years
    sxy.allyears[ite,,,(t+Tt)] <- sxy.proj.core[ite,,,(t+1)] 
    z.allyears[ite,,(t+Tt)] <- z.proj.core[ite,,(t+1)]
    age.cat.allyears[ite,,(t+Tt)] <- age.cat.proj.core[ite,,(t+1)]
  }}

## ---- 2. Subset abundance in buffer ----

# Unscale sxy

dimnames(sxy.allyears)[[3]] <- c('x','y')
sxy.allyears.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.allyears,## this are your sxy
                                              coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                              scaleToGrid = FALSE)$coordsDataScaled

# Matrix to store abundance in the buffer each iteration and year

NIn <- matrix(NA,nrow = dim(z.allyears)[1], ncol=dim(z.allyears)[3]) # nrow = iterations, ncol = year

for(ite in 1:dim(z.allyears)[1]){
  for(t in 1:dim(z.allyears)[3]){
    
    which.alive <- which(z.allyears[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY, proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    NIn[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    
  }
}

#average number of individuals without the buffer for each year 
colMeans(NIn)

# Sum of individuals alive in total each year (without buffer)
NALL <- apply(z.allyears,c(1,3),function(x) sum(x==1, na.rm = TRUE))
colMeans(NALL)

## ---- Plot ----

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("prediction_pcrCore.pdf", 7, 7)

plot(1, ylim = c(30,max(NIn)+100), 
     xlim = c(0.5, ncol(NIn)+0.5) , 
     type ="n", 
     #yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

polygon(x = c(5.5,5.5,11,11), y = c(0,max(NIn)+150,max(NIn)+150,0), col = adjustcolor("grey", alpha.f = 0.2), border = NA)

axis(1, c(1:ncol(NIn)), labels = c(2017:2026), las = 2, cex.axis = 1)

for (i in 1:ncol(NIn)){
  plot.violins3(list(NIn[ ,i]),
                x = i,
                at = i,
                violin.width = 0.4,
                plot.ci = 0.95,
                col = c("darksalmon"),
                add = T,
                alpha = 0.5,
                scale.width = FALSE,
                border.col = "darkred",
                horizontal = FALSE)}

dev.off()
