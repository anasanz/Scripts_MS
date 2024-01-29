## -------------------------------------------------
##      Size estimation PTS using Bayesian
##              Close population model
## ------------------------------------------------- 

rm(list=ls())

## -------------------------------------------------
##                1. M0: Constant p in time
##              Different detection per sex
## ------------------------------------------------- 


library(jagsUI)
library(rjags)

## ---- Simulation ----

# Store iter
niter = 100
Ndf <- data.frame(matrix(NA, nrow = niter, ncol = 4))
colnames(Ndf) <- c("mean", "sd", "2.5%", "97.5%")

p1df <- data.frame(matrix(NA, nrow = niter, ncol = 4))
colnames(p1df) <- c("mean", "sd", "2.5%", "97.5%")

p2df <- data.frame(matrix(NA, nrow = niter, ncol = 4))
colnames(p2df) <- c("mean", "sd", "2.5%", "97.5%")

#for (i in 1:niter){
  data.fn <- function(N = 100, p = c(0.2,0.7), T = 3, psi = 0.5){
    
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
  
  data <- data.fn()
  
  str(data)
  
  data$yobs
  sum(data$sex)
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
 
   # Store N
 # Ndf[i,] <- out_sim_M0_psex$summary[rownames(out_sim_M0_psex$summary) %in% "N", colnames(out_sim_M0_psex$summary) %in% c("mean", "sd", "2.5%", "97.5%")]
 # p1df[i,] <- out_sim_M0_psex$summary[rownames(out_sim_M0_psex$summary) %in% "p[1]", colnames(out_sim_M0_psex$summary) %in% c("mean", "sd", "2.5%", "97.5%")]
 # p2df[i,] <- out_sim_M0_psex$summary[rownames(out_sim_M0_psex$summary) %in% "p[2]", colnames(out_sim_M0_psex$summary) %in% c("mean", "sd", "2.5%", "97.5%")]
  
#}

# It is corrected now, I was not sppecifying correctly the sex of observed individuals
# Plot posterior summaries

setwd("D:/MargSalas/Ganga/Results/CR")
pdf("Summary_posteriors_100iter_pSex.pdf")

par(mfrow = c(3,1))
plot(density(Ndf$mean), main = "meanN_100iter", col = "darkcyan", lwd = 1.2, xlab = " ", ylab = " ", bty = "n", ylim = c(0,0.08), xlim = c(60,110))
polygon(c(density(Ndf$mean)$x, 0), c(density(Ndf$mean)$y, 0), col="darkcyan", border = "darkcyan") 

polygon(c(density(Ndf$`2.5%`)$x, 0), c(density(Ndf$`2.5%`)$y, 0), border = "lightblue", lwd = 2) 
polygon(c(density(Ndf$`97.5%%`)$x, 0), c(density(Ndf$`97.5%`)$y, 0), border = "darkblue", lwd = 2) 

abline(v = data$N, lwd = 2)


plot(density(p1df$mean), main = "meaP1_100iter", col = "darkcyan", lwd = 1.2, xlab = " ", ylab = " ", bty = "n", ylim = c(0,7), xlim = c(0,1))
polygon(c(density(p1df$mean)$x, 0), c(density(p1df$mean)$y, 0), col="darkcyan", border = "darkcyan") 

polygon(c(density(p1df$`2.5%`)$x, 0), c(density(p1df$`2.5%`)$y, 0), border = "lightblue", lwd = 2) 
polygon(c(density(p1df$`97.5%%`)$x, 0), c(density(p1df$`97.5%`)$y, 0), border = "darkblue", lwd = 2) 

abline(v = 0.2, lwd = 2)

plot(density(p2df$mean), main = "meaP2_100iter", col = "darkcyan", lwd = 1.2, xlab = " ", ylab = " ", bty = "n", ylim = c(0,8), xlim = c(0,1))
polygon(c(density(p2df$mean)$x, 0), c(density(p2df$mean)$y, 0), col="darkcyan", border = "darkcyan") 

polygon(c(density(p2df$`2.5%`)$x, 0), c(density(p2df$`2.5%`)$y, 0), border = "lightblue", lwd = 2) 
polygon(c(density(p2df$`97.5%%`)$x, 0), c(density(p2df$`97.5%`)$y, 0), border = "darkblue", lwd = 2) 

abline(v = 0.7, lwd = 2)

dev.off()
