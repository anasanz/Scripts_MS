

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)

sourceCpp("GetSpaceUse_PD.cpp")
source("getDensityInput")




#give dim names to your posteriors sxy
dimnames(YOURPOSTERIORSXY)[[3]] <- c("x","y")

## first rescale the coordinates to the original scale 
myResultsSXYZ_M$sims.list$sxy <- scaleCoordsToHabitatGrid(
  coordsData = YOURPOSTERIORSXY,
  coordsHabitatGridCenter = myHabitat$habitat.xy,
  scaleToGrid = FALSE)$coordsDataScaled

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = rrRegions,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = rrRegions,## here put the same than regions argument. 
  s = YOURPOSTERIORSXY,
  plot.check = TRUE)


##extract density
DensityCountriesRegions <- list()
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = myResultsSXYZ_MF$sims.list$z[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:nyears){
ACDens[[t]] <- densityInputRegions$regions.r
ACDens[[t]][] <- NA
ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
plot(ACDens[[t]])
}