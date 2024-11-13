# Load packages

rm(list=ls())

library(rjags)
library(jagsUI)
library(dplyr)
library(rtrim)


setwd("D:/PhD/Otros/Tórtola/Data")

tor <- read.csv("tortola_ds_ready_reobs.csv", sep = ",")
tor[,1] <- "STTUR"
###################################################################
##                       HDS ANALYSIS                           ###
###################################################################

# ---- Information: bins, years, sites ----

strip.width <- 500 				
dist.breaks <- c(0,25,100,500)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-4]
nG <- length(dist.breaks)-1

unique(tor$Year)
yrs <- unique(tor$Year) 
nyrs <- length(yrs)

# ---- Distance observations ----

# Format
all.sites <- unique(tor$site_sec)
max.sites <- length(all.sites)


m <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(m) <- all.sites
colnames(m) <- yrs

# Add counts > 0
tor_yes <- tor[which(tor$count == 1), ]
count <- aggregate(X ~ Year + site_sec, FUN = length, data = tor_yes)

for (i in 1:nrow(count)){
  m[which(rownames(m) %in% count$site_sec[i]), which(colnames(m) %in% count$Year[i])] <- count$X[i]
}

# Add absences (0)
tor_no <- tor[which(tor$count == 0), ]
for (i in 1:nrow(tor_no)){
  m[which(rownames(m) %in% tor_no$site_sec[i]), which(colnames(m) %in% tor_no$Year[i])] <- tor_no$count[i]
}

# Seleccionar aÃ±os del 2005 al 2019
m <- m[,c(4:18)]
# Remove transects only made 4 years or less
m <- m[-which(rowSums(is.na(m)) > 10), ]
sum(rowSums(m[which(rowSums(is.na(m)) > 10), ], na.rm = TRUE)) # Se pierden 553 detecciones de 2906
sum(rowSums(m, na.rm = TRUE))


# Change number of sites 
all.sites <- row.names(m)
max.sites <- length(all.sites)
# Change years
yrs <- colnames(m)
nyrs <- length(yrs)


# Only to check: Count of individuals per year
count.year <- colSums(m,na.rm = TRUE)

# Year
yrs2 <- c(0:(nyrs-1)) # To make it as a continuous variable, otherwise it doesnt work
year <- matrix(NA,nrow = max.sites, ncol = nyrs)
colnames(year) <- yrs
for (i in 1:nyrs){
  year[ ,which(colnames(year) %in% yrs[i])] <- rep(yrs2[i], max.sites)
}

# Observer 
# Format
obs <- matrix(NA, nrow = max.sites, ncol = nyrs)
rownames(obs) <- all.sites
colnames(obs) <- yrs

# Add observers for fields with counts > 0
for (i in 1:nrow(tor)){
  obs[which(rownames(obs) %in% tor$site_sec[i]), which(colnames(obs) %in% tor$Year[i])] <- tor$Observer[i]
}


# Forest
setwd("D:/PhD/Otros/Tórtola/Data")

load("D:/PhD/Otros/Tórtola/Data/data_buff500.rdata")
colnames(data_buff500)[1] <- "Site"
colnames(data_buff500)[2] <- "Section"
data_buff500$site_sec <- paste(data_buff500$Site, data_buff500$Section, sep = "_")
var <- data_buff500[ ,c("site_sec", "Hm_mean", "Hm_max", "FCC_mean", "FCC100%", "Forest%", "richness", "ForestMargin", "ForestMargin%")]
var2 <- var[,which(colnames(var) %in% c("site_sec", "Forest%"))]
# Format
forest <- matrix(NA, nrow = max.sites, ncol = nyrs)
rownames(forest) <- all.sites
colnames(forest) <- yrs

# Add observers for fields with counts > 0
for (i in 1:nrow(tor)){
  forest[which(rownames(forest) %in% var2$site_sec[i]), ] <- var2$`Forest%`[i] 
}

# ---- Specify data in JAGS format ----

# Update dataset "tor" with new sites (deleted sites with less than 5 years sampled)
tor <- tor[which(tor$site_sec %in% rownames(m) & tor$Year %in% colnames(m)), ]
tor_no <- tor_no[which(tor_no$site_sec %in% rownames(m) & tor_no$Year %in% colnames(m)), ]
tor_yes <- tor_yes[which(tor_yes$site_sec %in% rownames(m) & tor_yes$Year %in% colnames(m)), ]


# Distance class and ind
nind <- nrow(tor_yes)
dclass <- tor_yes$Bin

m  # Counts per year and site

# Co-variates

yrs <- 1:nyrs
year_number <- yrs2

# Matrix with observers
ob <- matrix(as.numeric(factor(obs)), nrow = max.sites, ncol = nyrs) # JAGS doesn't accept categorical variables
unique(factor(ob))
obs_id <- unique(factor(ob))[-2]
ob[which(is.na(ob))] <- sample(obs_id, length(which(is.na(ob))), replace = TRUE) # No NA in covariate

nobs <- length(unique(factor(ob)))

# Matrix with forest (put random values where NA)
unique(factor(forest))
for_id <- unique(factor(forest))[-22]
forest[which(is.na(forest))] <- sample(for_id, length(which(is.na(forest))), replace = TRUE) # No NA in covariate

forest_mean <- mean(forest)
forest_sd <- sd(forest)
forest_sc <- (forest - forest_mean) / forest_sd

# Indexes
nsites <- length(unique(tor$Site))
sec <- 1:6
nSection <- length(sec)

# site <- c(1:nsites) # I think I don't need this
# year <- c(1:nyrs)

# Create 2 vectors to index the sitesection random effect

sections <- rep(sec,nsites)
sites <- 1:nsites
sitessections <- rep(sites, each = nSection)

# Fixed index to map dclass onto site and year 
# For the index, create a matrix m where NA are 0 (because I need the same length)

m_index <- m
m_index[which(is.na(m_index))] <- 0

site.dclass <- year.dclass <- NULL

for (t in 1:nyrs){ # sites has to be nested on years because dclass first indexes the sites on the same year
  for (j in 1:max.sites){
    site.dclass <- c(site.dclass, rep(j, m_index[j,t]))
    year.dclass <- c(year.dclass, rep(t, m_index[j,t]))
  } }



# ---- Compile data for JAGS model ----

data1 <- list(nyears = nyrs, SiteSection = max.sites,  nsites = nsites, nSection = nSection, sections = sections, sitessections = sitessections, nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks,
              year.dclass = year.dclass, site.dclass = site.dclass, y = m, nind=nind, dclass=dclass,
              forestCov = forest_sc, ob = ob, nobs = nobs, year1 = year_number, 
              #site = site, 
              year_index = yrs)

# ---- JAGS model ----

setwd("D:/Ana/Model/Tórtola")
cat("model{
    
    # PRIORS
    
    # PRIORS FOR LAMBDA
    rho ~ dunif(-1,1) # Autorregresive parameter (serial AC)
    tau <- pow(sd, -2) # Prior for overdispersion in eps
    sd ~ dunif(0, 3)
    
    bYear.lam ~ dnorm(0, 0.001) # Prior for the trend
    
    # Random effects for lambda per section nested in site
    # Piors random effect
    mu.site ~ dunif(-10, 10) 
    sig.site ~ dunif(0, 10)
    tau.site <- 1/(sig.site*sig.site)
    
    tau.section <- 1/(sig.section*sig.section) # Hyperparameters for variation WITHIN transect per section: Level 2 of ranfom effect (NESTED)
    sig.section ~ dunif(0,500)  
    
    for (j in 1:nsites){
    mu.lambda.site[j] ~ dnorm(mu.site, tau.site)      # Level 1: One mean per transect (site)
    for (s in 1:nSection){
    lam.section[j,s] ~ dnorm(mu.lambda.site[j],tau.section)  # Level 2: Section variation within transect (nested)
    }    
    }
    
    # Random effects for lambda per year
    sig.lam.year ~ dunif(0, 10) 
    tau.lam.year <- 1/(sig.lam.year*sig.lam.year)
    
    log.lambda.year[1] <- 0
    for (t in 2:nyears){
    log.lambda.year[t] ~ dnorm(0, tau.lam.year)
    }
    
    
    # PRIORS FOR SIGMA
    bforest.sig ~ dnorm(0, 0.001)
    
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
    for(j in 1:SiteSection){ 
    
    sigma[j,1] <- exp(sig.obs[ob[j,1]] + bforest.sig*forestCov[j,1] + log.sigma.year[year_index[1]])
    
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
    
    lambda[j,1] <- exp(lam.section[sitessections[j], sections[j]] + log.lambda.year[year_index[1]] + bYear.lam*year1[1] + w[j,1]) # year1 is t-1; year_index is t (to index properly the random effect)
    w[j,1] <- eps[j,1] / sqrt(1 - rho * rho)
    eps[j,1] ~ dnorm(0, tau)
    
    # Bayesian p-value on abundance component 
    
    Nnew[j,1]~dpois(lambda[j,1]) ##Create replicate abundances for year 1
    
    FT1[j,1]<-pow(sqrt(N[j,1])-sqrt(lambda[j,1]),2) ### residuals for 'observed' and new abundances in year 1
    FT1new[j,1]<-pow(sqrt(Nnew[j,1])-sqrt(lambda[j,1]),2)
    }
    
    #############
    # LATER YEARS
    for(j in 1:SiteSection){ 
    for (t in 2:nyears){
    
    sigma[j,t] <- exp(sig.obs[ob[j,t]] + bforest.sig*forestCov[j,t] + log.sigma.year[year_index[t]])
    
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
    
    lambda[j,t] <- exp(lam.section[sitessections[j], sections[j]] + log.lambda.year[year_index[t]] + bYear.lam*year1[t] + w[j,t])
    w[j,t] <- rho * w[j,t-1] + eps[j,t]
    eps[j,t] ~ dnorm(0, tau)
    
    # Bayesian p-value on abundance component (rest of years)
    
    Nnew[j,t]~dpois(lambda[j,t]) # create replicate abundances for rest of the years
    
    FT1[j,t]<-pow(sqrt(N[j,t])-sqrt(lambda[j,t]),2) # residuals for 'observed' and new abundances for the rest of the years
    FT1new[j,t]<-pow(sqrt(Nnew[j,t])-sqrt(lambda[j,t]),2)
    }
    }
    
    T1p <- sum(FT1[1:SiteSection,1:nyears]) #Sum of squared residuals for actual data set (RSS test)
    T1newp <- sum(FT1new[1:SiteSection,1:nyears]) # Sum of squared residuals for new data set (RSS test)
    
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
    
    
    }",fill=TRUE, file = "model7.txt")


# Inits
Nst <- m_index + 1
inits <- function(){list(mu.sig = runif(1, log(30), log(50)), sig.sig = runif(1),
                         mu.site = runif(1), sig.site = 0.2, sig.section = runif(1), bforest.sig = runif(1),
                         sig.lam.year = 0.3, bYear.lam = runif(1), 
                         N = Nst)} 

# Params
params <- c( "mu.sig", "sig.sig", "bforest.sig", "sig.obs", "log.sigma.year", # Save also observer effect
             "mu.site", "sig.site", "sig.section", "sig.lam.year", "bYear.lam", "log.lambda.year", # Save year effect
             "popindex", "sd", "rho", "lam.tot",'Bp.Obs', 'Bp.N'
)


# MCMC settings
nc <- 3 ; ni <- 250000 ; nb <- 3000 ; nt <- 5

# With jagsUI 
# With jagsUI 
out <- jags(data1, inits, params, "model7.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

summary <- out$summary
print(out)

# Error in node lam.section[95,3]
setwd("D:/Ana/Results/Otros/Tortola")
save(out, file = "7tor_1.2.RData")

