## -------------------------------------------------
##      Size estimation PTS using Bayesian
##              Close population model
## ------------------------------------------------- 

## -------------------------------------------------
##                1. M0: Constant p
## ------------------------------------------------- 


library(jagsUI)
library(rjags)

## ---- Simulation ----

data.fn <- function(N = 100, p = 0.5, T = 3){
  yfull <- yobs <- array(NA, dim = c(N, T))
  for (j in 1:T){
    yfull[,j] <- rbinom(n = N, size = 1, prob = p)
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs))
}

data <- data.fn()

str(data)

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model_m0.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)			# Inclusion indicators (probability that exists)
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p		# Can only be detected if z=1
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT <1 min)
out_sim_M0 <- jags(win.data, inits, params, "model_m0.txt", n.chains = nc, 
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(out_sim_M0, dig = 3)
hist(out_sim_M0$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 150))
abline(v = data$C, lwd = 3)
abline(v = mean(out_sim_M0$sims.list$N), col = "blue", lwd = 3)