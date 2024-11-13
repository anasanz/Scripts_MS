rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)



### Same than original trend model, but adding bQ
# mu.lam.site slightly over estimated

#set.seed(2013)

# ---- Data simulation ----

# y[jt] ~ bin(p(sigma),N[jt])
# Sigma[jt] <- obs[jt] + temp[jt] + year

# N[jt] ~ Pois(lambda[jt])
# log(lambda[jt]) <- site.alpha.ran[j] + year.ran[t] + beta*yr[t-1] + beta*hq + W

#####
# ---- Distance sampling data ----

# HAZARD RATE DETECTION FUNCTION
g <- function(x, sig, b) 1 - exp(-(x/sig)^-b)

# Number of transects per year (unbalanced)
nSites <- rep(100,9)			# Same number of transects by year
max.sites <- max(nSites)            # Maximun number of sites is the last year

strip.width <- 200 				# strip half-width, w (in this example only one side of the line transect is surveyed)
dist.breaks <- c(0,25,50,100,200)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	

# Year effect 
yrs <- 1:9 
nyrs <- length(yrs)
year_number <- 0:8 # (RS: start from 0)

###
# ---- Detection component ----

# RANDOM EFFECT IN OBSERVER
obs <- 1:9
nobs <- length(obs)
mu.sig.obs <- log(50)
sig.sig.obs <- 0.25
# Observer effect in sigma
sig.obs <- rnorm(length(obs), mu.sig.obs, sig.sig.obs) 
# Observer covariate
ob.id <- matrix(sample(1:9, max.sites*nyrs, replace = TRUE), nrow = max.sites, ncol = nyrs) # Matix with IDs
ob <- matrix(sig.obs[ob.id],  nrow = max.sites, ncol = nyrs) # Matrix with intercept for simulating data

# YEAR EFFECT IN SIGMA (RANDOM)
sig.sig.year <- 1		
sig.year <- rnorm(nyrs, 0, sig.sig.year) 

# TEMPERATURE COVARIATE
bTemp.sig <- 0.5
temp <- matrix(rnorm(max.sites*nyrs, 50, 7), nrow = max.sites, ncol = nyrs)
#SCALED
temp_mean <- mean(temp)
temp_sd <- sd(temp)
temp_sc <- (temp - temp_mean) / temp_sd


#SIGMA
sigma <- exp(ob + bTemp.sig * temp_sc + 
               matrix(rep(sig.year, each = max.sites), nrow = max.sites, ncol = nyrs))

#BETA
b <- 2

#####
# ----  Abundance component: Site effect, Year effect and Year trend

# Site effect (RANDOM INTERCEPT)
mu.lam.alpha.site <- log(1.5)				
sig.lam.alpha.site <- 0.5				
lam.alpha.site <- rnorm(max.sites, mu.lam.alpha.site, sig.lam.alpha.site) 


# Year effect (RANDOM)
sig.lam.year <- 0.7			
lam.year <- rnorm(nyrs, 0, sig.lam.year) 


#TIME CO-VARIATE (YEAR)

b.lam.year <- 0.3
year <- matrix(NA,nrow = max.sites, ncol = nyrs)
colnames(year) <- yrs
for (i in 0:nyrs){
  year[ ,yrs[i]] <- rep(year_number[i], max.sites)
}

# HQ COVARIATE

bHQ <- 0.5
hq <- matrix(runif(max.sites*nyrs, 0, 3), nrow = max.sites, ncol = nyrs)
#SCALED
hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd


# AUTOCORRELATION AND OVERDISPERSION TERM

rho <- 0.5 # Autoregressive parameter

sig.lam.eps <- 0.2 
eps <- matrix(NA,nrow = max.sites, ncol = nyrs) # Unstructured random variation for overdispersion
for (j in 1:max.sites){
  for (t in 1:nyrs){
    eps[j,t] <- rnorm(1,0,sig.lam.eps)
  }
}

w <- matrix(NA,nrow = max.sites, ncol = nyrs)
lam <- matrix(NA,nrow = max.sites, ncol = nyrs) 

# First year
for(j in 1:max.sites){
  w[j,1] <- eps[j,1] / sqrt(1 - rho * rho)
  lam[j,1] <- exp(lam.alpha.site[j] + 
                    lam.year[1] + 
                    b.lam.year*year[j,1] +
                    bHQ*hq_sc[j,1] +
                    w[j,1])
}

# Later years
for (j in 1:max.sites){
  for (t in 2:nyrs){
    w[j,t] <- rho * w[j,t-1] + eps[j,t]
    lam[j,t] <- exp(lam.alpha.site[j] + 
                      lam.year[t] + 
                      b.lam.year*year[j,t] +
                      bHQ*hq_sc[j,t] +
                      w[j,t])
  }
}

lam.tot <- colSums(lam)


######
# ---- Generate ABUNDANCE per site and year ----

# Abundance
N <- matrix(NA,nrow = max.sites, ncol = nyrs) 
for (j in 1:max.sites){
  for (t in 1:nyrs){
    N[j,t] <- rpois(1,lam[j,t])
  }}

N.tot <- colSums(N)

# Introduce NA (not sampled)
vec <- seq(1,length(N))
na <- sample(vec, 100)
N[na]<-NA

# Cluster size (to correct) per site and year
clus <- list()
for (t in 1:nyrs){
  clus[[t]] <- rpois(nSites[t], 1.5)
} 
clusLong <- ldply(clus,cbind) # 1 long vector with all abundances per site and year
clus3 <- ldply(clus,rbind)
clus_size <- t(clus3) # CLUSTER SIZE per site and year stored in a matrix with columns (this is not really necessary to do, with inventing an average would be fine)
average_clus <- mean(clus_size, na.rm = TRUE)

# Correct N with average cluster size
Nclus <- N * average_clus
Nclus.tot <- colSums(Nclus,na.rm = TRUE) # Total pop.abundance corrected by cluster size


####
# ---- Simulate continuous distance data ----

# Nc = count of individuals detected in each distance interval
yList <- list()
for (i in 1:nyrs){
  yList[[i]] <- array(0, c(nSites[i], length(dist.breaks)-1))
}

for (t in 1:nyrs){
  for(j in 1:max.sites) {
    if(N[j,t] == 0 | is.na(N[j,t]))
      next
    # Distance from observer to the individual
    d <- runif(N[j,t], 0, strip.width) 		# Uniform distribution of animals
    # Simulates one distance for each individual in the site (N[j])
    p <- g(x=d, sig=sigma[j,t], b = b)   		# Detection probability. Sigma is site-time specific
    seen <- rbinom(N[j,t], 1, p)
    if(all(seen == 0))
      next
    d1 <- d[seen==1] 				# The distance data for seen individuals
    counts <- table(cut(d1, dist.breaks, include.lowest=TRUE))
    yList[[t]][j,] <- counts 				# The number of detections in each distance interval
  }}

y.sum.sites <- lapply(yList, function(x) rowSums(x)) # Total count per site each year
y.sum.sites2 <- ldply(y.sum.sites,rbind)
y.sum <- t(y.sum.sites2) # y per site and year stored in a matrix with columns
y.sum[na] <- NA # Add what are real NA generated from Na in N (not 0)


####
# ---- Convert data to JAGS format ----

nind.year <- lapply(yList,sum)
nind <- sum(unlist(nind.year, use.names = F))

y.sum # Matrix with counts 

# Co-variates

temp_sc
ob.id # Matrix with observers
year_number # Vector with year variable

site <- c(1:max.sites)
year <- c(1:nyrs)

# Get one long vector with years, distance category and site

#site <- dclass <- year <- NULL
dclass <- site.dclass <- year.dclass <- NULL # Fixed index to map dclass onto site and year 

for (t in 1:nyrs){
  for(j in 1:max.sites){
    if (y.sum[j,t] == 0 | is.na(y.sum[j,t])) 
      next
    #site <- c(site, rep(j, y.sum[j,t])) # site index: repeat the site as many times as counts in that site (for multi model??)
    # vector of sites through years (disregarding distance class)
    #year <- c(year, rep(t, y.sum[j,t]))
    
    for (k in 1:nG){
      if (yList[[t]][j,k] == 0) # Refers for the ditance classes to the list with years and bins
        next 
      dclass <- c(dclass, rep(k, yList[[t]][j,k]))	# Distance category index
      site.dclass <- c(site.dclass, rep(j, yList[[t]][j,k]))
      year.dclass <- c(year.dclass, rep(t, yList[[t]][j,k]))
    }}
}

# Create one matrix for indexing year when calculating abundance per year in JAGS

allyears <- NULL 
for (i in 1:nyrs){
  allyears <- c(allyears,rep(yrs[i],nSites[i]))
}
m <- data.frame(allyears = allyears)
m$allyears <- as.factor(m$allyears)
indexYears <- model.matrix(~ allyears-1, data = m)

####
# ---- Compile data for JAGS model ----

data1 <- list(nyears = nyrs, nsites = max.sites, nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks,
              year.dclass = year.dclass, site.dclass = site.dclass, y = y.sum, nind=nind, dclass=dclass,
              hqCov = hq_sc, tempCov = temp_sc, ob = ob.id, nobs = nobs, year1 = year_number, site = site, year_index = yrs)


# Inits
Nst <- y.sum + 1
inits <- function(){list(mu.sig = runif(1, log(30), log(50)), sig.sig = runif(1), b = runif(1),
                         mu.lam.site = runif(1), sig.lam.site = 0.2, sig.lam.year = 0.3, bYear.lam = runif(1), bHQ = runif(1),
                         N = Nst)} 

# Params
params <- c( "mu.sig", "sig.sig", "bTemp.sig", "sig.obs", "log.sigma.year", "b", # Save also observer effect
             "mu.lam.site", "sig.lam.site", "sig.lam.year", "bYear.lam", "log.lambda.year", "bHQ", # Save year effect
             "popindex", "sd", "rho", "lam.tot",'Bp.Obs', 'Bp.N'
)


# MCMC settings
nc <- 3 ; ni <- 70000 ; nb <- 3000 ; nt <- 5

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
#setwd("~/Scripts_MS/Ganga/HDS/Farmdindis/Model")
source("2.2.HDS_trendmodel_lam[hq]_sigHR.r")

# With jagsUI 
out <- jags(data1, inits, params, "2.2.HDS_trendmodel_lam[hq]_sigHR.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

setwd("D:/MargSalas/Ganga/Results/HDS/Model_results")
save(out, file = "2.2.Sim_HDS_trendmodel_lam[hq]_sigHR.RData")

summary <- out$summary

data_comp <- list(lam.tot = lam.tot, 
                  mu.sig.obs = mu.sig.obs, sig.sig.obs = sig.sig.obs,
                  bTemp.sig = bTemp.sig,
                  bHQ = bHQ,
                  mu.lam.alpha.site = mu.lam.alpha.site,
                  sig.lam.alpha.site = sig.lam.alpha.site,
                  sig.lam.year = sig.lam.year,
                  b.lam.year = b.lam.year,
                  rho = rho, sig.lam.eps = sig.lam.eps)


