## -------------------------------------------------
##                Simulation study
## ------------------------------------------------- 

# Provide proof that our p and n lead to unbiased estimates
# Generate data using the estimated population size and detection probabilities and run the model again
# Using our final model: Variation of Mh, with constant p in time, different detection per sex

rm(list=ls())

library(jagsUI)
library(rjags)

## ---- Load our model results to set population parameters ----

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/model_results")
load("out_sim_M0_psex.RData")
sum <- out_sim_M0_psex$summary

Npop <- sum[rownames(sum) == "N", colnames(sum) == "mean"]
Npop_UPCI <- sum[rownames(sum) == "N", colnames(sum) == "97.5%"]
Npop_LOWCI <- sum[rownames(sum) == "N", colnames(sum) == "2.5%"]

p1 <- sum[rownames(sum) == "p[1]", colnames(sum) == "mean"]
p2 <- sum[rownames(sum) == "p[2]", colnames(sum) == "mean"]
psi <- sum[rownames(sum) == "psi", colnames(sum) == "mean"]


## ---- Simulation ----~

# Data simulation function

data.fn <- function(N = 100, p = c(0.1,0.2), T = 5, psi = 0.5){
  
  sex <- rbinom(n = N, size = 1, prob = psi) # 0 females, 1 males
  
  yfull <- yobs <- array(NA, dim = c(N, T))
  
  for (i in 1:N){
    for (j in 1:T){
      yfull[i,j] <- rbinom(n = 1, size = 1, prob = p[sex[i]+1])
    }}
  
  ever.detected <- apply(yfull, 1, max)
  
  C <- sum(ever.detected)
  
  yobs <- yfull[ever.detected == 1,]
  sexobs <- sex[ever.detected == 1]
  
  cat(C, "out of", N, "animals present were detected.\n")
  
  return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs, sex = sex, sexobs = sexobs))
}

# Store iter
niter = 100

Ndf <- list()

for (i in 1:niter){
  
  data <- data.fn(N = Npop, p = c(p1,p2), T = 5, psi = psi)
  
  # Augment data set by 150 potential individuals
  nz <- 150
  yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))
  sexAug <- c(data$sexobs, rep(NA, nz))
  
  # Specify model in BUGS language
  setwd("D:/MargSalas/Ganga/Data/CMR")
  sink("model_m0_pSex.txt")
  cat("
    model {
    
    # Priors
    omega ~ dunif(0, 1)
    psi ~ dunif(0, 1)
    
    for (s in 1:2){
      p[s] ~ dunif(0, 1)
    }
    
    # Likelihood
    for (i in 1:M){
       sex[i] ~ dbern(psi)	
       z[i] ~ dbern(omega)			# Inclusion indicators (probability that exists)
       for (j in 1:T){
          yaug[i,j] ~ dbern(p.eff[i,j])
          p.eff[i,j] <- z[i] * p[sex[i]+1]		# Can only be detected if z=1
          } #j
       } #i
    
    # Derived quantities
    N <- sum(z[])
    }",fill = TRUE)
  sink()
  
  # Bundle data
  win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug), sex = sexAug)
  
  # Initial values
  inits <- function() list(z = rep(1, nrow(yaug)), p = runif(2, 0, 1),
                           sex = c(rep(NA,nrow(data$yobs)), rbinom(nz,1,0.5)))
  
  # Parameters monitored
  params <- c("N", "p", "omega", "psi")
  
  # MCMC settings
  ni <- 2500
  nt <- 2
  nb <- 500
  nc <- 3
  
  # Call WinBUGS from R (BRT <1 min)
  out_sim_M0_psex <- jags(win.data, inits, params, "model_m0_pSex.txt", n.chains = nc, 
                          n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
  Ndf[[i]] <- out_sim_M0_psex

}

setwd("D:/MargSalas/Ganga/Results/CR/optimize_sampling/results")
save(Ndf, file = "ResultsN_Simulation_accuracy_100iter.RData")

# Plot

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/model_results")
load("ResultsN_Simulation_accuracy_100iter.RData")

# Calculate bias

bias_mean <- list()
bias_mode <- list()

for (i in 1:niter){
  # Mean
  mean_ab <- Ndf[[i]]$mean$N
  bias_mean[[i]] <- mean_ab - Npop
  # Mode
  dens_obs1 <- density(Ndf[[i]]$sims.list$N)
  mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
  bias_mode[[i]] <- mode_ab - Npop
}

bias <- data.frame(mean = unlist(bias_mean), mode = unlist(bias_mode))
mean(bias[,2])

setwd("D:/MargSalas/Ganga/Results/Plots/Optimize")
pdf("Simulation_accuracy_correct_boxplot_mean.pdf", 7,5)

boxplot(bias[,1], pch = 19, cex = 0.5, medlwd = 1.3, col = "#35978f", xlab = "Bias", ylab = " ")
points(mean(bias[,1]), pch = 18, col = "white")

dev.off()


setwd("D:/MargSalas/Ganga/Results/Plots/Optimize")
pdf("Simulation_accuracy_correct_boxplot_mode.pdf", 7,5)

boxplot(bias[,2], pch = 19, cex = 0.5, medlwd = 1.3, col = "#35978f", xlab = "Bias", ylab = " ")
points(mean(bias[,2]), pch = 18, col = "white")

dev.off()




