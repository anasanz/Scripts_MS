## -------------------------------------------------
##            Single season closed SCR
## ------------------------------------------------- 

# Simulate data (OG) and run model 

library(nimble)
library(dplyr)
library(ggplot2)

## ---- Data simulation ----

# Traps

tr <- seq(15, 85, length = 10)
traps <- data.frame(X = rep(tr, each = length(tr)),
                    Y = rep(tr, times = length(tr))) # 100 coord. traps

viz_traps <- traps %>% 
  ggplot(aes(x = X, y = Y)) +
  geom_point(pch = 3) + 
  xlim(0, 100) +
  ylim(0, 100)
viz_traps

# Generate population

set.seed(10)
xlim <- c(0, 100)
ylim <- c(0, 100) # area 100 * 100 = 1e4
A <- (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])/10000
A

mu <- 50 # density
N <- rpois(1, mu*A) # generate population
N 

#Generate activity centers.

s <- data.frame(s.x = runif(N, xlim[1], xlim[2]), 
                s.y = runif(N, ylim[1], ylim[2]))

viz_traps_ac <- viz_traps + 
  geom_point(data = s, aes(x = s.x, y = s.y), pch = 16, color = "red")
viz_traps_ac


# Generate detections.

sigma <- 5
lambda0 <- 0.4 # p detection is 0.4 when activity center is at the trap
J <- nrow(traps) # nb of traps
K <- 5 # nb capture occasions

yy <- array(NA, c(N, J, K))
dim(yy[,,1]) # Cada "hoja" es una occasion de 50 individuos*100 trampas


for(j in 1:J) {
  # Distance from trap j to the 50 activity centers:
  dist <- sqrt((traps$X[j] - s$s.x)^2 + (traps$Y[j] - s$s.y)^2) 
  # Expected number of detections (lambda) at trap j for each of the activity centers:
  lambda <- lambda0 * exp(-dist^2 / (2 * sigma^2))    # Detection function
  # In this case is a poisson model, so you use lambda as the expected number of detections at the trap location
  # (and lambda0 would be the expected numper of detections if the AC is at the trap location)
  # If you would use a binomial model, it would be the probability of having a detection
  
  # For each occasion k and trap j
  for(k in 1:K) {
    yy[,j,k] <- rpois(N, lambda) # yy is the number of detections (lambda is the same in each occasion, different per trap because of the distance)
  }
}

n <- apply(yy, c(2,3), sum) # Suma de detectiones/occasion en cada trampa

# Plot detections

tot <- apply(n, 1, sum)
dat <- data.frame(traps, tot = tot)

viz_traps_ac +
  geom_point(data = dat, aes(x = X, y = Y, size = tot), alpha = 0.3) +
  scale_size(range = c(0, 20)) +
  labs(x = "",
       y = "",
       size = "# detections")
