rm(list=ls())

library(rjags)
library(jagsUI)
library(dplyr)
library(stringr)


setwd("D:/Otros/Tórtola/Data/Columbid")

tor <- read.csv("stdec_ds_ready_02_21.csv", sep = ",")
tor[,1] <- "STDEC"

###################################################################
##                       HDS ANALYSIS                           ###
###################################################################

# Model 1.1: Same as model 1 but with roughness covariate instead of forest

# ---- Information: bins, years, sites ----

strip.width <- 500 				
dist.breaks <- c(0,25,100,500)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-4]
nG <- length(dist.breaks)-1

yrs <- unique(tor$Year) 
nyrs <- length(yrs)

# 1 transect and 6 sections per year 
nSection <- 6
section <- c(1,2,3,4,5,6)

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

# Only to check: Count of individuals per year
count.year <- colSums(m,na.rm = TRUE)

## ---- SUBSET the data by transect ----

# 1. Keep transects with > 7 time points
m <- as.data.frame(m)
m_subset <- data.frame(matrix(nrow = 0, ncol = ncol(m)))
colnames(m_subset) <- colnames(m)

for (i in 1:nrow(m)) {
  if (sum(!is.na(m[i, ])) > 6) {
    m_subset[nrow(m_subset)+1,] <-  m[i, ] # Add row if it has 7 values or more 
  }}

# 2. Keep transects with at least 5 time points with counts (To be able to calculate trend)

m_subset$site <- str_sub(rownames(m_subset), end=-3) # Add column with site (removing 2 last digits section) to group per site
x <- t(sapply(split(m_subset, m_subset$site), function(x) colSums(x[, c(1:18)]))) # Counts grouped per site
x <- as.data.frame(ifelse(x>0,1,0))

m_subset2 <- data.frame(matrix(nrow = 0, ncol = ncol(m_subset)))
colnames(m_subset2) <- colnames(m_subset)

for (i in 1:nrow(m_subset)) {
  site <- m_subset[i,which(colnames(m_subset) %in% "site")] # Identify the site for the if condition
  if (sum(x[rownames(x) %in% site, ], na.rm = TRUE) > 4) {
    m_subset2[nrow(m_subset2)+1,] <-  m_subset[i, ] # Add row if it has 5 time points with counts or more
  }}

# Transects that will be analyzed
length(unique(m_subset2$site)) 
transect_stdec <- unique(m_subset2$site)

# 3. Check if they fit with transects Tórtola

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
load("1.1TortoData_transects_0221.RData") # Load analyzed transects
transect_stdec <- transect_stdec[transect_stdec %in% transect] 
# Only 20 of the transects that can be analyzed with STDEC were analyzed with STTUR!! :(
# Save them
# setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
# save(transect, file = "1.1CopalData_transects_0221.RData")

# Subset the m df so that we only analyze the transects of sttur
m_subset2 <- m_subset2[which(m_subset2$site %in% transect_stdec), ]


# 4. Load variables

# ROUGHNESS CO-VARIATE
setwd("D:/Otros/Tórtola/Data")
rough_var <- read.csv("TRI_buff_500.csv", sep = ",")
rough_var$site_sec <- paste(rough_var[,1], rough_var[,2], sep = "_")
rough_var <- rough_var[,c(16,15)]

# Format
rough <- matrix(NA, nrow = nrow(m_subset2), ncol = ncol(m_subset2))
rownames(rough) <- rownames(m_subset2)
colnames(rough) <- colnames(m_subset2)

# Add data into the same structure than m_subset2
for (i in 1:nrow(rough_var)){
  rough[which(rownames(rough) %in% rough_var$site_sec[i]), ] <- rough_var$X_mean[i] 
}
rough <- as.data.frame(rough) 
rough$site <- m_subset2$site # For subset

# 3. Make the subset and start loop

transect <- transect_stdec


for (xxx in 1:length(transect)){
  
  data_transect <- m_subset2[which(m_subset2$site %in% transect[xxx]), ]
  
  # Year: Min and max year (if there are years with NA in between, those are considered within the time series)
  
  data_transect_nona <- data_transect[, colSums(is.na(data_transect)) != nrow(data_transect)] # Delete columns with NA
  
  min_year <- min(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) == "site")]))
  max_year <- max(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) == "site")]))
  
  nyrs <- max_year - min_year + 1 # For JAGS
  year_number <- 0:(nyrs-1)
  yrs <- 1:nyrs
  
  # Counts per year and site
  
  m <- data_transect[ ,-which(colnames(data_transect) %in% "site")] # For JAGS
  m <- m[,which(colnames(m) %in% c(min_year:max_year))] # Select data among threshold years (no NA in the limits of the time series)
  
  # Distance class and ind
  
  dtransect_distance <- tor_yes[which(tor_yes$Site %in% transect[xxx]), ]
  
  nind <- nrow(dtransect_distance) # For JAGS
  dclass <- dtransect_distance$Bin
  
  # Fixed index to map dclass onto site and year 
  # For the index, create a matrix m where NA are 0 (because I need the same length)
  
  m_index <- as.matrix(m)
  m_index[which(is.na(m_index))] <- 0
  
  year.dclass <- section.dclass <- NULL
  
  for (t in 1:nyrs){ 
    for (j in 1:nSection){
      section.dclass <- c(section.dclass, rep(j, m_index[j,t]))
      year.dclass <- c(year.dclass, rep(t, m_index[j,t]))
    }}
  
  # rough variable
  rough_transect <- as.matrix(rough[which(rough$site %in% transect[xxx]), which(colnames(rough) %in% colnames(m))])
  
  rough_mean <- mean(rough_transect)
  rough_sd <- sd(rough_transect)
  rough_sc <- (rough_transect - rough_mean) / rough_sd
  
  
  # ---- Compile data for JAGS model ----
  
  data1 <- list(nyears = nyrs, nsection = nSection, nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks,
                year.dclass = year.dclass, section.dclass = section.dclass, y = m, nind=nind, dclass=dclass,
                roughCov = rough_sc, year1 = year_number, section = section, year_index = yrs)
  
  
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
    brough.sig ~ dnorm(0, 0.001)
    
    mu.sig.year ~ dunif(-10, 10) # Random effects for sigma per observer
    sig.sig.year ~ dunif(0, 10)
    tau.sig.year <- 1/(sig.sig.year*sig.sig.year)
    
    # Random year effect for sigma
    for (t in 1:nyears){
    log.sigma.year[t] ~ dnorm(mu.sig.year, tau.sig.year)
    }
    
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
    
    sigma[j,1] <- exp(log.sigma.year[year_index[1]] + brough.sig*roughCov[j,1])
    
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
    
    sigma[j,t] <- exp(log.sigma.year[year_index[t]] + brough.sig*roughCov[j,t])
    
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
    
    
    }",fill=TRUE, file = "Model1.1_study2.txt")
  
  # Inits
  Nst <- m_index + 1
  inits <- function(){list(mu.lam.year = runif(1), sig.lam.year = runif(1), bYear.lam = runif(1),
                           mu.sig.year = runif(1, log(30), log(50)), sig.sig.year = runif(1), brough.sig = runif(1),
                           N = Nst)} 
  
  # Params
  params <- c( "mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", 
               "popindex", "sd", "rho", "lam.tot",'Bp.Obs', 'Bp.N'
  )
  
  
  # MCMC settings
  nc <- 3 ; ni <- 500000 ; nb <- 3000 ; nt <- 5
  
  # With jagsUI 
  out <- jags(data1, inits, params, "Model1.1_study2.txt", n.chain = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
  summary <- out$summary
  print(out)
  
  
  setwd("D:/Otros/Tórtola/Results/Study2/Model1.1")
  save(out, file = paste("1.1StdecData_0221_", transect[xxx], ".RData"), sep = "")
  
}