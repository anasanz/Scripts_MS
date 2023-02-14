# ---- JAGS model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
cat("model{
    
    # PRIORS
    
    # PRIORS FOR LAMBDA
    rho ~ dunif(-1,1) # Autorregresive parameter (serial AC)
    tau <- pow(sd, -2) # Prior for overdispersion in eps
    sd ~ dunif(0, 3)
    
    bYear.lam ~ dnorm(0, 0.001) # Prior for the trend
    
    bHQ ~ dnorm(0, 0.001)

    
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
    
    # PRIOR FOR BETA
    b ~ dunif(0, 100)
    
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
    
    p[j,1,k]<-1-exp(-(midpt[k]/sigma[j,1])^-b)
    pi[j,1,k] <- int.w[k] / strip.width 
    fc[j,1,k]<- p[j,1,k] * pi[j,1,k]                 ## pi=percent area of k; drops out if constant
    fct[j,1,k]<-fc[j,1,k]/pcap[j,1] 
    }
    
    pcap[j,1] <- sum(fc[j,1, 1:nG]) # Different per site and year (sum over all bins)
    
    f.0[j,1] <- 2 * dnorm(0,0, 1/sigma[j,1]^2) # Prob density at 0
    
    
    y[j,1] ~ dbin(pcap[j,1], N[j,1]) 
    N[j,1] ~ dpois(lambda[j,1]) 
    
    lambda[j,1] <- exp(log.lambda.site[site[j]] + log.lambda.year[year_index[1]] + bYear.lam*year1[1] + bHQ*hqCov[j,1] + w[j,1]) # year1 is t-1; year_index is t (to index properly the random effect)
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
    
    p[j,t,k]<-1-exp(-(midpt[k]/sigma[j,t])^-b)
    pi[j,t,k] <- int.w[k] / strip.width 
    fc[j,t,k]<- p[j,t,k] * pi[j,t,k]                 ## pi=percent area of k; drops out if constant
    fct[j,t,k]<-fc[j,t,k]/pcap[j,t] 
    }
    
    pcap[j,t] <- sum(fc[j,t, 1:nG]) # Different per site and year (sum over all bins)
    
    f.0[j,t] <- 2 * dnorm(0,0, 1/sigma[j,t]^2) # Prob density at 0
    
    
    y[j,t] ~ dbin(pcap[j,t], N[j,t]) 
    N[j,t] ~ dpois(lambda[j,t]) 
    
    lambda[j,t] <- exp(log.lambda.site[site[j]] + log.lambda.year[year_index[t]] + bYear.lam*year1[t] + bHQ*hqCov[j,t] + w[j,t])
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
    
    
    }",fill=TRUE, file = "2.2.HDS_trendmodel_lam[hq]_sigHR.txt")
