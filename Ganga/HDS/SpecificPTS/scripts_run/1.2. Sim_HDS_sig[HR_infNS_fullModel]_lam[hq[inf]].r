## -------------------------------------------------
##          Data simulation HDS model for PTS
##           Abundance estimation in 2022
##          Informative prior in detection (1 sigma)
##               Informative prior in hq
## ------------------------------------------------- 

# Check for fixed effect script in: "D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 2/Bayes/Simulation"
# 1.2. TRIM_PoisGam. S_x[dcat_area_integralBin]_sigma[obs(j,t)_covZone(j)]_lambda[site.alpha.random(j)_year(t)_beta.year(j)].r
rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)

set.seed(2013)

## ---- Survey characteristics ----

# HAZARD RATE DETECTION FUNCTION
g <- function(x, sig, b) 1 - exp(-(x/sig)^-b)

# Number of transects per year (unbalanced)
nSites <- 22			# number of line transect surveys 

strip.width <- 400 				# SAME THAN IN NATIONAL CENSUS
dist.breaks <- c(0,100,200,300, 400)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	

# ---- Detection component: Sig ----

# Covariates: Generate scaled covariates from NS model

# HOUR OF DAY
h_mean <- 0.46054 # Mean and sd from NS model
h_sd <- 0.2002978
  
h <- rnorm(nSites, h_mean , h_sd) # 1. Simulate variable
h_sc <- (h - h_mean) / h_sd # 2. Centered on the variable used in NS (squared mean and sd) -> It should be the same of course
hSquare_sc <- h_sc^2 # 3. Squared


# JULIAN DAY
jd_mean <- 140.2123
jd_sd <- 23.35039

jd <- rnorm(nSites, jd_mean , jd_sd) # 1. Simulate variable
jd_sc <- (jd - jd_mean) / jd_sd # 2. Centered on the variable used in NS
jdSquare_sc <- jd_sc^2 #3. Squared


## DOUBT: The squared variables are centered on the squared mean or the normal mean? 
##        RS: First you have to center the variable and then square it if you want

# Same model than the national PTS monitoring program (informative priors in all)

mu.alpha <- 3.87567
sig.alpha <- 0.0875
alpha <- rnorm(1, mu.alpha, sig.alpha)

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
Vebro <- rnorm(1, mu.Vebro, sig.Vebro) # Fixed effect but informed, becomes a random effect? I think normally it needs to be centered in 0

sigma <- exp(alpha + b.hSquare * hSquare_sc +
                b.jd * jd_sc + b.jdSquare * jdSquare_sc +
                Vebro)

# Variable site
site <- 1:nSites

# Beta
mu.b <- 0.997
sig.b <- 0.0438 

beta <- rnorm(1, mu.b, sig.b)


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
  p <- g(x=d, sig=sigma[j], b = beta)   		# Detection probability. Sigma is site-time specific
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

data1 <- list(nSites = nSites, 
              mu.bHQ = mu.bHQ, sig.bHQ = sig.bHQ, # Informative prior for habitat quality
              mu.alpha = mu.alpha, sig.alpha = sig.alpha, # Informative priors for sigma
              mu.hSquare = mu.hSquare, sig.hSquare = sig.hSquare, 
              mu.jd = mu.jd, sig.jd = sig.jd, 
              mu.jdSquare = mu.jdSquare, sig.jdSquare = sig.jdSquare, 
              mu.Vebro = mu.Vebro, sig.Vebro = sig.Vebro, 
              mu.b = mu.b, sig.b = sig.b, # Informative prior for shape parameter
              nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
              y = yLong, nind=nind, dclass=dclass, site.dclass = site.dclass,
              hq = hq_sc, # Variables
              hSquare = hSquare_sc,
              jd = jd_sc,
              jdSquare = jdSquare_sc
              )

## ---- Inits ----

Nst <- yLong + 1
inits <- function(){list(alpha = runif(1),
                         N = Nst
)}

## ---- Params ----

params <- c("Ntotal",
            "alpha", "bHQ", 
            #"sigma", 
            "log.b",
            "alpha.sig", "b.hSquare", "b.jd", "b.jdSquare", "log.Vebro")

## ---- MCMC settings ----
nc <- 3 ; ni <- 20000 ; nb <- 1000 ; nt <- 2

## ---- Run model ----
setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
source("1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]].r")

# With jagsUI 
out <- jags(data1, inits, params, "1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]].txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)


sum <- out$summary
# Super wide CIs


sim.data <- list( N.tot = N.tot,
                  alpha = alpha,
                  mu.bHQ = mu.bHQ,
                  mu.alpha = mu.alpha,
                  mu.hSquare = mu.hSquare,
                  mu.jd = mu.jd, 
                  mu.jdSquare = mu.jdSquare, 
                  mu.Vebro = mu.Vebro,
                  mu.b = mu.b)
