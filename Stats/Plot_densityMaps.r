
remotes::install_github("mikemeredith/makeJAGSmask")

library(makeJAGSmask)

data(simSCR)
str(simSCR, 2)

# Get the bits we need
sigma <- simSCR$sims.list$sigma
S <- simSCR$sims.list$S    # coordinates of activity centres on the pixel scale
w <- simSCR$sims.list$w    # w = 1 if real, 0 if phantom
S[w == 0] <- NA            # replace coordinates with NA for phantoms

habMat <- simSCR$JAGSmask$habMat
trapMat <- simSCR$JAGSmask$trapMat
upperLimit <- simSCR$JAGSmask$upperLimit

niter <- dim(w)[1]

## 1. PLOTTING ACTIVITY CENTRES

image(x=(1:65)+0.5, y=(1:59)+0.5, z=habMat, asp=1, xlab="", ylab="", col=c("grey", "white"))
for(i in 1:5)
  points(S[,i,], cex=0.1, col=i+1)
for(i in 6:30)
  points(S[,i,], cex=0.1, col=1)
points(trapMat, pch=3, col='red')


## 2. POSTERIOR PROBABILITY DENSITY
## Count the proportion of ACs that fall in each pixel
# To count the number of AC:
totup <- function(index, mat) {
  index <- matrix(index, length(index)/2, 2)
  for(i in 1:nrow(index))
    mat[index[i,1], index[i,2]] <- mat[index[i,1], index[i,2]] + 1
  return(mat)
}

# For one animal:
ACmat <- array(0, dim=dim(habMat))  # initial all-0 matrix
ACmat1 <- totup(S[,2,], ACmat)      # only tot up #2
ACmat <- ACmat1 / niter             # convert to proportion
ACmat[habMat == 0] <- NA            # NA out the non-habitat (it will be white)
ACcol <- heat.colors(12, rev=TRUE)
image(x=(1:65)+0.5, y=(1:59)+0.5, z=ACmat, asp=1, xlab="", ylab="", col=ACcol)
points(trapMat, pch=3, cex=0.5)     # add the traps

## ANIMAL DENSITY: PROBABILITY OF AN ANIMAL AVAILABLE FOR DETECTION IN EACH PIXEL
# Takes into account movement (sigma)
# Array for the posterior predictive animal locations
AL <- array(NA, dim=dim(S))
nind <- dim(S)[2]
# Work through all iterations and all animals
for(i in 1:niter)
  for(j in 1:nind)
    if(w[i,j] == 1) {
      hab <- 0
      while(hab == 0) {
        AL[i,j, ] <- S[i,j, ] + rnorm(2, 0, sigma[i])
        # habitat check
        if(all(AL[i,j,] >= 1) && all(AL[i,j,] < upperLimit))
          hab <- habMat[trunc(AL[i,j,1]),trunc(AL[i,j,2])]
      }
    }
