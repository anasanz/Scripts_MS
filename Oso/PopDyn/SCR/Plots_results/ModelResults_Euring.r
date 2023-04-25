## -------------------------------------------------
##           Plot Model results PTT EURING
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

# Load buffer core area and habitat grid (to subset in sampling area)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

# Constants:

M.aug <- 300 # Augmented individuals estimation model
Tt <- 5 # Nyears estimation model (2017:2021)

## ---- 1. Abundance ----

# Store:

sxy.estimation <- array(NA, c(dim(sampmat)[1], M.aug, 2, Tt)) # Dataframes all years to store
z.estimation <- age.cat.estimation <- array(NA, c(dim(sampmat)[1], M.aug, Tt))

for(t in 1:Tt){
  for(ite in 1:dim(sampmat)[1]){
    sxy.estimation[ite,1:M.aug,,t] <- myResultsSXYZ$sims.list$sxy[ite,,,t] # Only 300 augmented individuals from past
    z.estimation[ite,1:M.aug,t] <- myResultsSXYZ$sims.list$z[ite,,t]
    age.cat.estimation[ite,1:M.aug,t] <- myResultsSXYZ$sims.list$age.cat[ite,,t]
  }}

## Subset abundance in buffer 

# Unscale sxy

dimnames(sxy.estimation)[[3]] <- c('x','y')
sxy.estimation.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.estimation,## this are your sxy
                                             coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                             scaleToGrid = FALSE)$coordsDataScaled

# Array to store abundance in the buffer each iteration and year

#NIn <- matrix(NA,nrow = dim(z.allyears)[1], ncol=dim(z.allyears)[3]) # nrow = iterations, ncol = year
NIn <- matrix(NA, nrow = dim(z.estimation)[1], ncol = dim(z.estimation)[3])

for(ite in 1:dim(z.estimation)[1]){
  for(t in 1:dim(z.estimation)[3]){
      
      which.alive <- which(z.estimation[ite,,t]==1) # Select only the individuals alive (z=1)
      
      which.aliveSXY <- sxy.estimation.uns[ite,which.alive,,t] # Retrieve the activity center for those individuals
      
      sp <- SpatialPoints(which.aliveSXY, proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
      
      which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
      
      NIn[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      
    }
  }


#average number of individuals without the buffer for each year 
ab <- round(colMeans(NIn),0)

## ---- Plot ----


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("Abundance.pdf", 7, 5)

plot(1, ylim = c(10,max(NIn) + 10), 
     xlim = c(0.5, ncol(NIn)+0.5) , 
     type ="n", 
     #yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

#polygon(x = c(5.5,5.5,22,22), y = c(0,max(NIn)+150,max(NIn)+150,0), col = adjustcolor("grey", alpha.f = 0.2), border = NA)

axis(1, c(1:ncol(NIn)), labels = c(2017:2021), 
     at = c(1,2,3,4,5), las = 1, cex.axis = 1)

# First five years (present)
for (i in 1:Tt){
  plot.violins3(list(NIn[ ,i]),
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

text(x = c(1.13, 2.13, 3.13, 4.13, 5.13), y = ab, labels = ab)
dev.off()

## ---- 2. Survival ----

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("Survival.pdf", 7, 5)
par(mfrow = c(1,1),
    mar = c(5,6,4,2) + 0.1)

params_phi <- c("phi.cub", "phi.sub", 'phi.ad')

plot(1, ylim = c(0.5, length(params_phi)+0.5), 
     xlim = c(0.4,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Survival probability",
     cex.axis = 0.8)

axis(2, c(1:length(params_phi)), labels = c("Cub(1+2)", "Subad(3+4)", "Adult(>4)"), las = 2, cex.axis = 1)

for (i in 1:length(params_phi)){
  plot.violins3(list(myResults$sims.list[names(myResults$sims.list) %in% params_phi[i]][[1]]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("yellow4"),
                add = T,
                alpha = 0.8,
                scale.width = FALSE,
                border.col = "yellow4",
                horizontal = TRUE)}

surv <- c(round(myResults$mean$phi.cub,2), round(myResults$mean$phi.sub,2), round(myResults$mean$phi.ad,2))
text(x = surv + 0.035, y = c(1,2,3), labels = surv)
dev.off()
