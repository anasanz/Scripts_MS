
## -------------------------------------------------
##          Data simulation HDS model for PTS
##           Abundance estimation in 2022
##              Site-specific sigma
##               Informative prior in hq
## ------------------------------------------------- 

### THIS ONE MAYBE BETTER FROM THE FARMDINDIS MODEL (ALREADY DONE SITE-SPEC)

rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)

set.seed(2013)

## ---- Survey characteristics ----

# Half normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

# Number of transects per year (unbalanced)
nSites <- 22			# number of line transect surveys 

strip.width <- 400 				# SAME THAN IN NATIONAL CENSUS
dist.breaks <- c(0,100,200,300, 400)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	

# ---- Detection component: Sig ----

mu.sig <- 3.55 
sig.sig <- 0.2941518 

sigma <- exp(rnorm(nSites, mu.sig, sig.sig)) # En params, no data

# Variable site
site <- 1:nSites

# ---- Abundance component: lam = alpha + b*HQ ----

# We have very scarce observations that will not allow estimating the habitat quality
# so we build an informative prior on bHQ from the farmdindis model with years 2010-2022

# setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
# load("2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData")
# summary <- as.data.frame(as.matrix(out$summary))

mu.bHQ <- 0.125796 # summary$mean[rownames(summary) %in% "bHQ"]
sig.bHQ<- 0.2112865 # summary$sd[rownames(summary) %in% "bHQ"]

bHQ <- rnorm(1, mu.bHQ, sig.bHQ)

hq <- runif(nSites, 0, 3)
hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd

alpha <- 0.3 # Intercept

lambda <- exp(alpha + bHQ*hq)

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
  p <- g(x=d, sig = sigma[j])   		# Detection probability. There is only one sigma (informed by national plan)
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

data1 <- list(nsites = nSites, site = site, mu.bHQ = mu.bHQ, sig.bHQ = sig.bHQ, site.dclass = site.dclass, 
              #mu.sig = mu.sig, sig.sig = sig.sig,
              nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
              y = yLong, nind=nind, dclass=dclass, hq = hq_sc)

## ---- Inits ----

Nst <- yLong + 1
inits <- function(){list(alpha = runif(1),
                         N = Nst
)}

## ---- Params ----

params <- c("Ntotal",
            "alpha", "bHQ", "mu.sig", "sig.sig"
            )

## ---- MCMC settings ----
nc <- 3 ; ni <- 1500 ; nb <- 200 ; nt <- 2

## ---- Run model ----
setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
source("3.1.HDS_sig[j]_lam[hq[inf]].r")

# With jagsUI 
out <- jags(data1, inits, params, "3.HDS_sig[j]_lam[hq[inf]].txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

sum <- out$summary
