
#--FOLLOWING ROYLE BOOK PAGE 136

sigma <- 2 # Transformed (in km): You transform it by multiplying the value from the model * resolution (5km in our case)

q<-qchisq(0.95,2)

radius<-sigma*sqrt(q)

area<-pi*radius^2



#Use the distribution in R

lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)

upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2)+100, nrow = 4, byrow = TRUE) 

s <- c(50, 50) # Currrent activity center location



baseIntensities <- rep(1,4)

habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)

numRows <- nrow(habitatGrid)

numCols <- ncol(habitatGrid)

numWindows <- 4


# in the same unit than nimbleSCR

#for(ite in 1:100){#run through iterations

sd <- 10

# The log probability density of moving from (1,1) to (1.2, 0.8)

dist <- list()


for(rep in 1:500){ # Repetitions to take into account that it is not necessarily always 10 (same as with sigma, we don't always detect individuals at a given distance)
  
  dist[[rep]] <-
    
    rbernppACmovement_normal(
      
      n               = 1,
      
      lowerCoords     = lowerCoords,
      
      upperCoords     = upperCoords,
      
      s               = s,
      
      sd              = sd,
      
      baseIntensities = baseIntensities,
      
      habitatGrid     = habitatGrid,
      
      numGridRows     = numRows,
      
      numGridCols     = numCols,
      
      numWindows      = numWindows)
}


##

predictedAC <- do.call(rbind,dist)

predictedDist <- rgeos::gDistance(SpatialPoints(predictedAC),SpatialPoints(matrix(s,nrow=1,ncol=2)),byid = T)


# this is the distances in scaled to the habitat grid.

# you need to

hist(predictedDist)

mean(predictedDist)

median(predictedDist)
