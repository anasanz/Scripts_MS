
rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)


#set.seed(2013)

# ---- Data simulation ----

# Abundance: random intercept in site and random effect in year
# Detection: random intercept in observer, random year effect, temperature covariate and wind co-variate

# 9 years of data
# Balanced number of transects per year (NA in not sampled ones)
# OVERDISPERSION AND SERIAL AUTOCORRELATION: W
# Model:

# y[jt] ~ bin(p(sigma),N[jt])
# Sigma[jt] <- obs[jt] + temp[jt] + wind[jt]

# N[jt] ~ Pois(lambda[jt])
# log(lambda[jt]) <- site.alpha.ran[j] + year.ran[t] + beta*yr[t-1] + W

#####
# ---- Distance sampling data ----

# Half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

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
                    w[j,1])
}

# Later years
for (j in 1:max.sites){
  for (t in 2:nyrs){
    w[j,t] <- rho * w[j,t-1] + eps[j,t]
    lam[j,t] <- exp(lam.alpha.site[j] + 
                      lam.year[t] + 
                      b.lam.year*year[j,t] +
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
    p <- g(x=d, sig=sigma[j,t])   		# Detection probability. Sigma is site-time specific
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
              tempCov = temp_sc, ob = ob.id, nobs = nobs, year1 = year_number, site = site, year_index = yrs)

# ---- JAGS model ----

setwd("S:/PhD/Second chapter/Data/Model")
cat("model{
    
    # PRIORS
    
    # PRIORS FOR LAMBDA
    rho ~ dunif(-1,1) # Autorregresive parameter (serial AC)
    tau <- pow(sd, -2) # Prior for overdispersion in eps
    sd ~ dunif(0, 3)
    
    bYear.lam ~ dnorm(0, 0.001) # Prior for the trend
    
    # Random effects for lambda per site
    mu.lam.site ~ dunif(-10, 10) 
    sig.lam.site ~ dunif(0, 10)
    tau.lam.site <- 1/(sig.lam.site*sig.lam.site)
    
    for (j in 1:nsites){
    log.lambda.site[j] ~ dnorm(mu.lam.site, tau.lam.site)
    }
    
    # Random effects for lambda per year
    sig.lam.year ~ dunif(0, 10) 
    tau.lam.year <- 1/(sig.lam.year*sig.lam.year)
    
    log.lambda.year[1] <- 0
    for (t in 2:nyears){
    log.lambda.year[t] ~ dnorm(0, tau.lam.year)
    }
    
    
    # PRIORS FOR SIGMA
    bTemp.sig ~ dnorm(0, 0.001)
    
    mu.sig ~ dunif(-10, 10) # Random effects for sigma per observer
    sig.sig ~ dunif(0, 10)
    tau.sig <- 1/(sig.sig*sig.sig)
    
    # Random observer effect for sigma
    for (o in 1:nobs){
    sig.obs[o] ~ dnorm(mu.sig, tau.sig)
    }
    
    # Random effects for sigma per year
    sig.sig.year ~ dunif(0, 10) 
    tau.sig.year <- 1/(sig.sig.year*sig.sig.year)
    
    for (t in 1:nyears){
    log.sigma.year[t] ~ dnorm(0, tau.sig.year)
    }
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[site.dclass[i], year.dclass[i], 1:nG]) 
    
    # Bayesian p-value for detection component (Bp.Obs)
    
    dclassnew[i] ~ dcat(fct[site.dclass[i], year.dclass[i], 1:nG]) # generate new observations
    Tobsp[i]<- pow(1- sqrt(fct[site.dclass[i], year.dclass[i],dclass[i]]),2) # Test for observed data
    Tobspnew[i]<- pow(1- sqrt(fct[site.dclass[i], year.dclass[i],dclassnew[i]]),2) # Test for new data
    }
    
    Bp.Obs <- sum(Tobspnew[1:nind]) > sum(Tobsp[1:nind])
    
    # LIKELIHOOD
    
    # FIRST YEAR
    for(j in 1:nsites){ 
    
    sigma[j,1] <- exp(sig.obs[ob[j,1]] + bTemp.sig*tempCov[j,1] + log.sigma.year[year_index[1]])
    
    # Construct cell probabilities for nG multinomial cells (distance categories) PER SITE
    
    for(k in 1:nG){ 
    
    up[j,1,k]<-pnorm(db[k+1], 0, 1/sigma[j,1]^2) ##db are distance bin limits
    low[j,1,k]<-pnorm(db[k], 0, 1/sigma[j,1]^2) 
    p[j,1,k]<- 2 * (up[j,1,k] - low[j,1,k])
    pi[j,1,k] <- int.w[k] / strip.width 
    f[j,1,k]<- p[j,1,k]/f.0[j,1]/int.w[k]                   ## detection prob. in distance category k                      
    fc[j,1,k]<- f[j,1,k] * pi[j,1,k]                 ## pi=percent area of k; drops out if constant
    fct[j,1,k]<-fc[j,1,k]/pcap[j,1] 
    }
    
    pcap[j,1] <- sum(fc[j,1, 1:nG]) # Different per site and year (sum over all bins)
    
    f.0[j,1] <- 2 * dnorm(0,0, 1/sigma[j,1]^2) # Prob density at 0
    
    
    y[j,1] ~ dbin(pcap[j,1], N[j,1]) 
    N[j,1] ~ dpois(lambda[j,1]) 
    
    lambda[j,1] <- exp(log.lambda.site[site[j]] + log.lambda.year[year_index[1]] + bYear.lam*year1[1] + w[j,1]) # year1 is t-1; year_index is t (to index properly the random effect)
    w[j,1] <- eps[j,1] / sqrt(1 - rho * rho)
    eps[j,1] ~ dnorm(0, tau)
    
    # Bayesian p-value on abundance component 
    
    Nnew[j,1]~dpois(lambda[j,1]) ##Create replicate abundances for year 1
    
    FT1[j,1]<-pow(sqrt(N[j,1])-sqrt(lambda[j,1]),2) ### residuals for 'observed' and new abundances in year 1
    FT1new[j,1]<-pow(sqrt(Nnew[j,1])-sqrt(lambda[j,1]),2)
    }
    
    #############
    # LATER YEARS
    for(j in 1:nsites){ 
    for (t in 2:nyears){
    
    sigma[j,t] <- exp(sig.obs[ob[j,t]] + bTemp.sig*tempCov[j,t] + log.sigma.year[year_index[t]])
    
    # Construct cell probabilities for nG multinomial cells (distance categories) PER SITE
    
    for(k in 1:nG){ 
    
    up[j,t,k]<-pnorm(db[k+1], 0, 1/sigma[j,t]^2) ##db are distance bin limits
    low[j,t,k]<-pnorm(db[k], 0, 1/sigma[j,t]^2) 
    p[j,t,k]<- 2 * (up[j,t,k] - low[j,t,k])
    pi[j,t,k] <- int.w[k] / strip.width 
    f[j,t,k]<- p[j,t,k]/f.0[j,t]/int.w[k]                   ## detection prob. in distance category k                      
    fc[j,t,k]<- f[j,t,k] * pi[j,t,k]                 ## pi=percent area of k; drops out if constant
    fct[j,t,k]<-fc[j,t,k]/pcap[j,t] 
    }
    
    pcap[j,t] <- sum(fc[j,t, 1:nG]) # Different per site and year (sum over all bins)
    
    f.0[j,t] <- 2 * dnorm(0,0, 1/sigma[j,t]^2) # Prob density at 0
    
    
    y[j,t] ~ dbin(pcap[j,t], N[j,t]) 
    N[j,t] ~ dpois(lambda[j,t]) 
    
    lambda[j,t] <- exp(log.lambda.site[site[j]] + log.lambda.year[year_index[t]] + bYear.lam*year1[t] + w[j,t])
    w[j,t] <- rho * w[j,t-1] + eps[j,t]
    eps[j,t] ~ dnorm(0, tau)
    
    # Bayesian p-value on abundance component (rest of years)
    
    Nnew[j,t]~dpois(lambda[j,t]) # create replicate abundances for rest of the years
    
    FT1[j,t]<-pow(sqrt(N[j,t])-sqrt(lambda[j,t]),2) # residuals for 'observed' and new abundances for the rest of the years
    FT1new[j,t]<-pow(sqrt(Nnew[j,t])-sqrt(lambda[j,t]),2)
    }
    }
    
    T1p <- sum(FT1[1:nsites,1:nyears]) #Sum of squared residuals for actual data set (RSS test)
    T1newp <- sum(FT1new[1:nsites,1:nyears]) # Sum of squared residuals for new data set (RSS test)
    
    # Bayesian p-value
    Bp.N <- T1newp > T1p
    
    # Derived parameters
    
    for(t in 1:nyears){
    popindex[t] <- sum(lambda[,t])
    }
    
    # Expected abundance per year inside model
    
    lam.tot[1] <- popindex[1] # Expected abundance in year 1
    for (i in 2:nyears){
    lam.tot[i] <- lam.tot[i-1] * # Here I add the starting population size as a baseline for the trend 
    exp(bYear.lam)}
    
    
    }",fill=TRUE, file = "s_sigma(integral)[obs(o,j,t)_covTemp(j,t)_year.random(t)]_lambda[alpha.site.random(j)_year.random(t)_beta.year(j)_w]_BayesP.txt")


# Inits
Nst <- y.sum + 1
inits <- function(){list(mu.sig = runif(1, log(30), log(50)), sig.sig = runif(1),
                         mu.lam.site = runif(1), sig.lam.site = 0.2, sig.lam.year = 0.3, bYear.lam = runif(1),
                         N = Nst)} 

# Params
params <- c( "mu.sig", "sig.sig", "bTemp.sig", "sig.obs", "log.sigma.year", # Save also observer effect
             "mu.lam.site", "sig.lam.site", "sig.lam.year", "bYear.lam", "log.lambda.year", # Save year effect
             "popindex", "sd", "rho", "lam.tot",'Bp.Obs', 'Bp.N'
)


# MCMC settings
nc <- 3 ; ni <- 70000 ; nb <- 3000 ; nt <- 5

# With jagsUI 
# With jagsUI 
out <- jags(data1, inits, params, "s_sigma(integral)[obs(o,j,t)_covTemp(j,t)_year.random(t)]_lambda[alpha.site.random(j)_year.random(t)_beta.year(j)_w]_BayesP.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
summary <- out$summary
print(out)


setwd("S:/PhD/Second chapter/Data/Results/TRIM")
save(out, file = "6.TRIM.RData") # 60000 iter, 4 thining

load("6.TRIM.RData")

out$summary
data_comp <- list(lam.tot = lam.tot, 
                  mu.sig.obs = mu.sig.obs, sig.sig.obs = sig.sig.obs,
                  bTemp.sig = bTemp.sig, 
                  mu.lam.alpha.site = mu.lam.alpha.site,
                  sig.lam.alpha.site = sig.lam.alpha.site,
                  sig.lam.year = sig.lam.year,
                  b.lam.year = b.lam.year,
                  rho = rho, sig.lam.eps = sig.lam.eps)

