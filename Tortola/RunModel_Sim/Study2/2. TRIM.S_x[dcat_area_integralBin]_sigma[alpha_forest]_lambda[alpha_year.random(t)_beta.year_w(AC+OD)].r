rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)


set.seed(2013)

# ---- Data simulation ----

# Model for only one transect transect, to get the independent trend at each transect: 15 years of data

# Abundance: random effect in year, trend and overdispersion
#     Every year the baseline abundance is different, but the same baseline for every section
# Detection: forest

# NO CONVERGE! COMO PUEDE SER? SI EL DE YEAR_RAN EN SIGMA SI!

# y[t] ~ bin(p(sigma),N[t])
# Sigma[t] <- alpha + forest[t] 

# N[t] ~ Pois(lambda[t])
# log(lambda[t]) <- year.ran[t] + beta*yr[t-1] + W

#####
# ---- Distance sampling data ----

# Half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

# Year effect 
yrs <- 1:15 
nyrs <- length(yrs)
year_number <- 0:14 # (RS: start from 0)

# 1 transect and 6 sections per year 
nSection <- 6
sec <- c(1,2,3,4,5,6)

strip.width <- 200 				# strip half-width, w (in this example only one side of the line transect is surveyed)
dist.breaks <- c(0,25,50,100,200)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	


###
# ---- Detection component ----

# FIXED INTERCEPT
alpha <- 0.3 

# FOREST COVARIATE
bforest.sig <- 0.5
forest <- matrix(rnorm(nSection*nyrs, 50, 10), nrow = nSection, ncol = nyrs)
#SCALED
forest_mean <- mean(forest)
forest_sd <- sd(forest)
forest_sc <- (forest - forest_mean) / forest_sd

#SIGMA
sigma <- exp(alpha + bforest.sig * forest_sc)


#####
# ----  Abundance component: Year effect and Year trend        

# YEAR EFFECT IN LAMBDA (RANDOM INTERCEPT)
sig.lam.year <- 0.1		
mu.lam.year <- 0.2
# Year effect in sigma
lam.year <- rnorm(nyrs, mu.lam.year, sig.lam.year) 

#TIME CO-VARIATE (YEAR)
b.lam.year <- 0.2
year <- matrix(NA,nrow = nSection, ncol = nyrs)
colnames(year) <- yrs
for (i in 0:nyrs){
  year[ ,yrs[i]] <- rep(year_number[i], nSection)
}

# AUTOCORRELATION AND OVERDISPERSION TERM
rho <- 0.5 # Autoregressive parameter

sig.lam.eps <- 0.2 
eps <- matrix(NA,nrow = nSection, ncol = nyrs) # Unstructured random variation for overdispersion
for (j in 1:nSection){
  for (t in 1:nyrs){
    eps[j,t] <- rnorm(1,0,sig.lam.eps)
  }
}

w <- matrix(NA,nrow = nSection, ncol = nyrs)
lam <- matrix(NA,nrow = nSection, ncol = nyrs) 

# First year
for(j in 1:nSection){
  w[j,1] <- eps[j,1] / sqrt(1 - rho * rho)
  lam[j,1] <- exp(lam.year[1] + 
                    b.lam.year*year[j,1] +
                    w[j,1])
}

# Later years
for (j in 1:nSection){
  for (t in 2:nyrs){
    w[j,t] <- rho * w[j,t-1] + eps[j,t]
    lam[j,t] <- exp(lam.year[t] + 
                      b.lam.year*year[j,t] +
                      w[j,t])
  }
}

######
# ---- Generate ABUNDANCE per site and year ----

# Abundance
N <- matrix(NA,nrow = nSection, ncol = nyrs) 
for (j in 1:nSection){
  for (t in 1:nyrs){
    N[j,t] <- rpois(1,lam[j,t])
  }}

N.tot <- colSums(N)

# Introduce NA (not sampled a given year)
vec <- seq(1, nyrs)
na <- sample(vec, 2)
N[,na]<-NA



####
# ---- Simulate continuous distance data ----

# Nc = count of individuals detected in each distance interval
yList <- list()
for (i in 1:nyrs){
  yList[[i]] <- array(0, c(nSection, length(dist.breaks)-1))
}

for (t in 1:nyrs){
  for(j in 1:nSection) {
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
y.sum[,na] <- NA # Add what are real NA generated from Na in N (not 0)

# ---- Convert data to JAGS format ----

nind.year <- lapply(yList,sum)
nind <- sum(unlist(nind.year, use.names = F))

y.sum # Matrix with counts 

# Co-variates

forest_sc
year_number # Vector with year variable

section <- c(1:nSection)
year <- c(1:nyrs)

# Get one long vector with years, distance category and site

#site <- dclass <- year <- NULL
dclass <- section.dclass <- year.dclass <- NULL # Fixed index to map dclass onto site and year 

for (t in 1:nyrs){
  for(j in 1:nSection){
    if (y.sum[j,t] == 0 | is.na(y.sum[j,t])) 
      next
    #site <- c(site, rep(j, y.sum[j,t])) # site index: repeat the site as many times as counts in that site (for multi model??)
    # vector of sites through years (disregarding distance class)
    #year <- c(year, rep(t, y.sum[j,t]))
    
    for (k in 1:nG){
      if (yList[[t]][j,k] == 0) # Refers for the ditance classes to the list with years and bins
        next 
      dclass <- c(dclass, rep(k, yList[[t]][j,k]))	# Distance category index
      section.dclass <- c(section.dclass, rep(j, yList[[t]][j,k]))
      year.dclass <- c(year.dclass, rep(t, yList[[t]][j,k]))
    }}
}

# Create one matrix for indexing year when calculating abundance per year in JAGS

allyears <- NULL 
for (i in 1:nyrs){
  allyears <- c(allyears,rep(yrs[i],nSection))
}
m <- data.frame(allyears = allyears)
m$allyears <- as.factor(m$allyears)
indexYears <- model.matrix(~ allyears-1, data = m)

####
# ---- Compile data for JAGS model ----

data1 <- list(nyears = nyrs, nsection = nSection, nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks,
              year.dclass = year.dclass, section.dclass = section.dclass, y = y.sum, nind=nind, dclass=dclass,
              forestCov = forest_sc, year1 = year_number, section = section, year_index = yrs)

# ---- JAGS model ----

setwd("D:/Otros/Tórtola/Model")
cat("model{
    
    # PRIORS
    
    # PRIORS FOR LAMBDA
    rho ~ dunif(-1,1) # Autorregresive parameter (serial AC)
    tau <- pow(sd, -2) # Prior for overdispersion in eps
    sd ~ dunif(0, 3)
    
    bYear.lam ~ dnorm(0, 0.001) # Prior for the trend
    
    # Random effects for lambda per year
    mu.lam.year ~ dunif(-10, 10) 
    sig.lam.year ~ dunif(0, 10)
    tau.lam.year <- 1/(sig.lam.year*sig.lam.year)
    
    log.lambda.year[1] <- 0
    for (t in 2:nyears){
    log.lambda.year[t] ~ dnorm(mu.lam.year, tau.lam.year)
    }
    
    
    # PRIORS FOR SIGMA
    bforest.sig ~ dnorm(0, 0.001)
    alpha.sig ~ dunif(-10, 10) 
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[section.dclass[i], year.dclass[i], 1:nG]) 
    
    # Bayesian p-value for detection component (Bp.Obs)
    
    dclassnew[i] ~ dcat(fct[section.dclass[i], year.dclass[i], 1:nG]) # generate new observations
    Tobsp[i]<- pow(1- sqrt(fct[section.dclass[i], year.dclass[i],dclass[i]]),2) # Test for observed data
    Tobspnew[i]<- pow(1- sqrt(fct[section.dclass[i], year.dclass[i],dclassnew[i]]),2) # Test for new data
    }
    
    Bp.Obs <- sum(Tobspnew[1:nind]) > sum(Tobsp[1:nind])
    
    # LIKELIHOOD
    
    # FIRST YEAR
    for(j in 1:nsection){ 
    
    sigma[j,1] <- exp(alpha.sig + bforest.sig*forestCov[j,1])
    
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
    
    lambda[j,1] <- exp(log.lambda.year[year_index[1]] + bYear.lam*year1[1] + w[j,1]) # year1 is t-1; year_index is t (to index properly the random effect)
    w[j,1] <- eps[j,1] / sqrt(1 - rho * rho)
    eps[j,1] ~ dnorm(0, tau)
    
    # Bayesian p-value on abundance component 
    
    Nnew[j,1]~dpois(lambda[j,1]) ##Create replicate abundances for year 1
    
    FT1[j,1]<-pow(sqrt(N[j,1])-sqrt(lambda[j,1]),2) ### residuals for 'observed' and new abundances in year 1
    FT1new[j,1]<-pow(sqrt(Nnew[j,1])-sqrt(lambda[j,1]),2)
    }
    
    #############
    # LATER YEARS
    for(j in 1:nsection){ 
    for (t in 2:nyears){
    
    sigma[j,t] <- exp(alpha.sig + bforest.sig*forestCov[j,t])
    
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
    
    lambda[j,t] <- exp(log.lambda.year[year_index[t]] + bYear.lam*year1[t] + w[j,t])
    w[j,t] <- rho * w[j,t-1] + eps[j,t]
    eps[j,t] ~ dnorm(0, tau)
    
    # Bayesian p-value on abundance component (rest of years)
    
    Nnew[j,t]~dpois(lambda[j,t]) # create replicate abundances for rest of the years
    
    FT1[j,t]<-pow(sqrt(N[j,t])-sqrt(lambda[j,t]),2) # residuals for 'observed' and new abundances for the rest of the years
    FT1new[j,t]<-pow(sqrt(Nnew[j,t])-sqrt(lambda[j,t]),2)
    }
    }
    
    T1p <- sum(FT1[1:nsection,1:nyears]) #Sum of squared residuals for actual data set (RSS test)
    T1newp <- sum(FT1new[1:nsection,1:nyears]) # Sum of squared residuals for new data set (RSS test)
    
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
    
    
    }",fill=TRUE, file = "Model2_study2.txt")

# Inits
Nst <- y.sum + 1
#inits <- function(){list(mu.lam.year = runif(1), sig.lam.year = runif(1), bYear.lam = runif(1),
#                         alpha.sig = runif(1), bforest.sig = runif(1),
#                        N = Nst)}
inits <- function(){
  
  list(
    
    #mu.A_mean_hyperprior= rnorm(1), mu.D_mean_hyperprior= rnorm(1)
    
  )}

# Params
params <- c( "mu.lam.year", "sig.lam.year", "bYear.lam", "alpha.sig", "bforest.sig", 
             "popindex", "sd", "rho", "lam.tot",'Bp.Obs', 'Bp.N'
)


# MCMC settings
nc <- 3 ; ni <- 100000 ; nb <- 3000 ; nt <- 5

# With jagsUI 
out <- jags(data1, inits, params, "Model1.1_study2.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
summary <- out$summary
print(out)



out$summary
data_comp <- list(N.tot = N.tot, 
                  alpha = alpha,
                  bforest.sig = bforest.sig,
                  mu.lam.year = mu.lam.year,
                  sig.lam.year = sig.lam.year,
                  b.lam.year = b.lam.year,
                  rho = rho, sig.lam.eps = sig.lam.eps)
