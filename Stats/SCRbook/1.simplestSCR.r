## -------------------------------------------------
##                 Simplest SCR model
## ------------------------------------------------- 

set.seed(36372) 
Area <- 1                     # area of state-space (unit square)
x <- cbind(rep(seq(.1,.9,.2), each = 5),  # trap locations
           rep(seq(.1,.9,.2), times = 5))

p0 <- 0.3                     # baseline capture probability
sigma <- 0.05                 # Gaussian scale parameter
mu <- 50                      # population density
N <- rpois(1, mu*Area)        # population size
s <- cbind(runif(N, 0, 1),    # activity centers in unit square
           runif(N, 0, 1))
K <- 5 

##  Capture data
# Format: one row per individual, one column per trap 
# (traps are like the occasions in normal cr?)

y <- matrix (NA, N, nrow(x))  

## Calculate probability of captures for each individual
# Per individual (N = 69): 
#       - 1. Vector of the distance from each trap to the activity center
#       - 2. Vector of the capture probability for each trap (p = 25)
#       - 3. Simulate capture history (detection in the 25 traps accordying to p)


for(i in 1:N) { 
  d.ij <- sqrt((x[,1] - s[i,1])^2 +     # 1. distance between x (traps) and s[i] (activity centers)
              (x[,2] - s[i,2])^2)

  p.ij <- p0*exp(-d.ij^2/ (2*sigma^2))  # 2. capture probability as a function from the 

  y[i,] <- rbinom(nrow(x), K, p.ij)     # 3. capture history for animal i
}

#### DOUBT: K = 5? Up to five times each individual detected in each trap?
#           Is k the time occasions??
