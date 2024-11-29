## -------------------------------------------------
##          Data simulation HDS model for PTS
##               Specific HDS model 
##           Simulate repeated surveys
## ------------------------------------------------- 

rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)

set.seed(2013)

## ---- 1. Set information common to all simulations, independently of the number of surveys  ----

### Survey characteristics ###

nSites <- 22			# number of line transect surveys
site <- 1:nSites # Variable site

strip.width <- 400 				# SAME THAN IN NATIONAL CENSUS
dist.breaks <- c(0,100,200,300, 400)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	


### Detection ###

# Parameter values determined by informative priors
# Same model than the national PTS monitoring program (informative priors in all)

mu.alpha <- 3.87567
sig.alpha <- 0.0875
alpha.sig <- rnorm(1, mu.alpha, sig.alpha)

mu.hSquare <- -0.00781
sig.hSquare <- 0.0154
b.hSquare <- rnorm(1, mu.hSquare, sig.hSquare)

mu.jd <- -0.05391
sig.jd <- 0.0240
b.jd <- rnorm(1, mu.jd, sig.jd)

mu.jdSquare <- 0.01176
sig.jdSquare <- 0.0100
b.jdSquare <- rnorm(1, mu.jdSquare, sig.jdSquare)

mu.Vebro <- 0.26297
sig.Vebro <- 0.0690
Vebro <- rnorm(1, mu.Vebro, sig.Vebro) 

# Beta
mu.b <- 0.997
sig.b <- 0.0438 

beta <- rnorm(1, mu.b, sig.b)

# Define HR detection function
g <- function(x, sig, b) 1 - exp(-(x/sig)^-b)

### Abundance ###

mu.bHQ <- 0.125796 # informative prior on bHQ from the farmdindis model with years 2010-2022
sig.bHQ<- 0.2112865 

bHQ <- rnorm(1, mu.bHQ, sig.bHQ) 

hq <- runif(nSites, 0, 3)
hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd

alpha <- 1.5 # Intercept

lambda <- exp(alpha + bHQ*hq_sc)
N <- rpois(nSites,lambda) # Abundance per site
N.tot <- sum(N) # Total number of individuals

## Store data generating parameters ###

sim.data <- list( N.tot = N.tot,
                  bHQ = bHQ,
                  alpha = alpha,
                  b.hSquare = b.hSquare,
                  b.jd = b.jd, 
                  b.jdSquare = b.jdSquare, 
                  Vebro = Vebro,
                  beta = beta)

## ---- 2. Simulation with diferent number of surveys ----

surv <- c(1:8) # Number of repeated surveys

info_surveys <- list()

for (sss in 1:length(surv)){
  
  nSurveys <- surv[[sss]]
    
  ## ---- Detection ----
  # Covariates: Generate scaled covariates from NS model
  
  # HOUR OF DAY
  h_mean <- 0.46054 # Mean and sd from NS model
  h_sd <- 0.2002978
  
  h <- matrix(rnorm(nSites * nSurveys, h_mean, h_sd), nrow = nSites) # 1. Simulate variable
  h_sc <- (h - h_mean) / h_sd # 2. Centered on the variable used in NS (squared mean and sd) -> It should be the same of course
  hSquare_sc <- h_sc^2 # 3. Squared
  
  # JULIAN DAY
  jd_mean <- 140.2123
  jd_sd <- 23.35039
  
  jd <- matrix(rnorm(nSites * nSurveys, jd_mean, jd_sd), nrow = nSites) # 1. Simulate variable
  jd_sc <- (jd - jd_mean) / jd_sd # 2. Centered on the variable used in NS
  jdSquare_sc <- jd_sc^2 #3. Squared
  
  
  # Simulate sigma and detection probabilities for each survey
  sigma <- array(NA, dim = c(nSites, nSurveys))
  for (s in 1:nSurveys) {
    sigma[, s] <- exp(alpha.sig + b.hSquare * hSquare_sc[, s] +
                        b.jd * jd_sc[, s] + b.jdSquare * jdSquare_sc[, s] +
                        Vebro)}
  
  # ---- Simulate continuous distance data ----
  # Add a survey dimension: 
  # - Simulate detection and distance categories for each survey
  # - Detection probabilities (derived from sigma) are survey-specific
  
  y <- array(0, dim = c(nSites, nG, nSurveys))
  
  for(j in 1:nSites) {
    
    if(N[j] == 0) next
    
    for(s in 1:nSurveys) {
      
    # Distance from observer to the individual
    d <- runif(N[j], 0, strip.width) 		# Uniform distribution of animals
    
    # Simulates one distance for each individual in the site (N[j])
    p <- g(x=d, sig=sigma[j,s], b = beta)   		# Detection probability. Sigma is site-time specific
    
    seen <- rbinom(N[j], 1, p)
    
    if(all(seen == 0))
      next
    
    d1 <- d[seen==1] 				# The distance data for seen individuals
    counts <- table(cut(d1, dist.breaks, include.lowest=TRUE))
    y[j,,s] <- counts
    } }
  
  y.sum.sites <- apply(y, c(1, 3), sum)
  
  # ---- Convert data to JAGS format ----
  
  nind <- sum(y.sum.sites)
  
  # Vector with counts
  yLong <- y.sum.sites
  
  # Distance category
  site.dclass <- dclass <- survey.dclass <- NULL
  
  for(j in 1:nSites){ ## HERE see with the model if I map first into sites and then into surveys or viceversa
    for (s in 1:nSurveys) {
      if (yLong[j,s] == 0) # Refers for the ditance classes to the list with years and bins
      next 
      site.dclass <- c(site.dclass, rep(j, yLong[j,s]))
      survey.dclass <- c(survey.dclass, rep(s, yLong[j, s]))
      for (k in 1:nG){
        dclass <- c(dclass, rep(k, y[j, k, s]))	# Distance category index
      } } }
  
  
  # ---- Compile data for JAGS model ----
  
  data1 <- list(nSites = nSites, nSurveys = nSurveys,
                mu.bHQ = mu.bHQ, sig.bHQ = sig.bHQ, # Informative prior for habitat quality
                mu.alpha = mu.alpha, sig.alpha = sig.alpha, # Informative priors for sigma
                mu.hSquare = mu.hSquare, sig.hSquare = sig.hSquare, 
                mu.jd = mu.jd, sig.jd = sig.jd, 
                mu.jdSquare = mu.jdSquare, sig.jdSquare = sig.jdSquare, 
                mu.Vebro = mu.Vebro, sig.Vebro = sig.Vebro, 
                mu.b = mu.b, sig.b = sig.b, # Informative prior for shape parameter
                nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
                y = yLong, nind=nind, dclass=dclass, site.dclass = site.dclass, survey.dclass = survey.dclass,
                hq = hq_sc, # Variables
                hSquare = hSquare_sc,
                jd = jd_sc,
                jdSquare = jdSquare_sc
  )
  
  ## ---- Inits ----
  
  Nst <- apply(yLong,1,max) + 1
  inits <- function(){list(alpha = runif(1),
                           N = Nst)}
  
  ## ---- Params ----
  
  params <- c("Ntotal",
              "alpha", "bHQ", 
              "log.b",
              "b",
              "alpha.sig", "b.hSquare", "b.jd", "b.jdSquare", "log.Vebro")
  
  ## ---- MCMC settings ----
  nc <- 3 ; ni <- 20000 ; nb <- 1000 ; nt <- 2
  
  ## ---- Run model ----
  setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
  source("4.HDS_2Survey_sig[HR_inf_fullModel]_lam[hq[inf]].r")
  
  # With jagsUI 
  out <- jags(data1, inits, params, "4.HDS_2Survey_sig[HR_inf_fullModel]_lam[hq[inf]].txt", n.chain = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
  
  sum <- out$summary
  info_surveys[[sss]] <- sum
}


sim.data
info_surveys

# Save for Rahel
setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Results")
save(sim.data, info_surveys, file = "Data_Results_Sim_occasions.RData")


