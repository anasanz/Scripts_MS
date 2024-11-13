
## -------------------------------------------------
##          Data simulation HDS model for PTS
##           Abundance estimation in 2022
##          Informative prior in detection
## ------------------------------------------------- 

rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)

set.seed(2013)

## ---- Survey characteristics ----

# Half normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

# Number of transects per year (unbalanced)
nSites <- 57			# number of line transect surveys 

strip.width <- 500 				# strip half-width
dist.breaks <- c(0,25,50,100,200,500)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	

# ---- Detection component: Sig = mu.sig (obs) + b*temp ----

# We will use an informative prior in sigma from the estimates obtained in model 2010-2020
# Average sigma accross all observers and years during 2010-2020
#     Intercept of sigma (Mean and SD of observer random effect) + 
#     Effect of temperature  (model beta coefficient)

#setwd("D:/Otros/Ganga/Trend_HDS_model_ch2/6. 10_20_AF/10_20") 
#load("PTALC_HDS_HQ_10_20.RData") # RESULTS MODEL 2010-2020
#summary <- as.data.frame(as.matrix(out$summary))

# Intercept of sigma (Mean and SD of observer random effect)

mu.sig <- 3.672591 # summary$mean[rownames(summary) %in% "mu.sig"]
sig.sig <- 0.2941518 #summary$mean[rownames(summary) %in% "sig.sig"]

log.sig <- rnorm(1, mu.sig, sig.sig)

# Effect of temperature  (model beta coefficient)
mu.bTemp <- 0.08144117 # summary$mean[rownames(summary) %in% "bTemp.sig"]
sig.bTemp <- 0.02336605 # summary$sd[rownames(summary) %in% "bTemp.sig"]

bTemp <- rnorm(1, mu.bTemp, sig.bTemp) 

temp <- rnorm(nSites, 17, 5) # Temp variable
temp_mean <- mean(temp)
temp_sd <- sd(temp)
temp_sc <- (temp - temp_mean) / temp_sd

# Average sigma per site accross all observers and years during 2010-2020
sigma <- exp(log.sig + bTemp * temp_sc ) # Real scale? (e.g., log of 45 meters)

# ---- Abundance component: lam = alpha + b*HQ ----

alpha <- 0.3 # Intercept

bHab <- 0.5 # Beta coef of habitat quality
hq <- runif(nSites, 0, 3)
hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd

lambda <- exp(alpha + bHab*hq)

# Abundance per site
N <- rpois(nSites,lambda)
N.tot <- sum(N) # Total number of individuals

# ---- Simulate continuous distance data ----

y <- matrix(0, nrow = nSites, ncol = length(dist.breaks)-1)

for(j in 1:nSites) {
  if(N[j] == 0) next
  # Distance from observer to the individual
  d <- runif(N[j], 0, strip.width) 		# Uniform distribution of animals
  # Simulates one distance for each individual in the site (N[j])
  p <- g(x=d, sig=sigma[j])   		# Detection probability. Sigma is site-time specific
  seen <- rbinom(N[j], 1, p)
  if(all(seen == 0))
    next
  d1 <- d[seen==1] 				# The distance data for seen individuals
  counts <- table(cut(d1, dist.breaks, include.lowest=TRUE))
  y[j,] <- counts
}

y.sum.sites <- rowSums(y)

# ---- Convert data to JAGS format ----

nind <- sum(y.sum.sites)

# Vector with counts
yLong <- y.sum.sites

# Distance category

site.dclass <- dclass <- NULL

for(j in 1:nSites){
  if (yLong[j] == 0) # Refers for the ditance classes to the list with years and bins
    next 
  site.dclass <- c(site.dclass, rep(j, yLong[j]))
  for (k in 1:nG){
    dclass <- c(dclass, rep(k, y[j, k]))	# Distance category index
    } }



# ---- Compile data for JAGS model ----

data1 <- list(nSites = nSites, mu.sig = mu.sig, sig.sig = sig.sig, mu.bTemp = mu.bTemp, sig.bTemp = sig.bTemp, temp = temp_sc, 
              nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
              y = yLong, nind=nind, dclass=dclass, site.dclass = site.dclass, hq = hq_sc)

## ---- Inits ----

Nst <- yLong + 1
inits <- function(){list(alpha = runif(1), bHQ = runif(1),
                         N = Nst
)}

## ---- Params ----

params <- c("Ntotal",
            "alpha", "bHQ")

## ---- MCMC settings ----
nc <- 3 ; ni <- 1500 ; nb <- 200 ; nt <- 2

## ---- Run model ----
setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
source("1.HDS_sig[inf]_lam[hq].r")

# With jagsUI 
out <- jags(data1, inits, params, "1.HDS_sig[inf]_lam[hq].txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

sum <- out$summary
