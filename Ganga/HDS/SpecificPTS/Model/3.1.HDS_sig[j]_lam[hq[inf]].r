## ---- JAGS model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
cat("model{
    
    # PRIORS
    
    # Priors for lambda
    alpha ~ dnorm(0, 0.001)
    bHQ ~ dnorm(0, 0.001)
    
    # Informative priors for sigma
    
    mu.sig ~ dunif(0, 5) # Random effects for sigma per site (IMPORTANT that is positive...I was generating negative sigmas with dunif(-10,10))
    sig.sig ~ dunif(0, 3)
    tau.sig <- 1/(sig.sig*sig.sig)
    
    ## Random site effect for sigma
    for (s in 1:nsites){
    sig[s] ~ dnorm(mu.sig, tau.sig)
    }
    
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[site.dclass[i], 1:nG])  
    }
    
    for(j in 1:nsites){ 
    
    sigma[j] <- exp(sig[j])
    
    # Construct cell probabilities for nG multinomial cells (distance categories) PER SITE
    
    for(k in 1:nG){ 

    up[j,k]<-pnorm(db[k+1], 0, 1/sigma[j]^2) ##db are distance bin limits
    low[j,k]<-pnorm(db[k], 0, 1/sigma[j]^2) 
    p[j,k]<- 2 * (up[j,k] - low[j,k])
    pi[j,k] <- int.w[k] / strip.width 
    f[j,k]<- p[j,k]/f.0[j]/int.w[k]                   ## detection prob. in distance category k                      
    fc[j,k]<- f[j,k] * pi[j,k]                 ## pi=percent area of k; drops out if constant
    fct[j,k]<-fc[j,k]/pcap[j] 
    }
    
    pcap[j] <- sum(fc[j, 1:nG]) # Different per site (sum over all bins)

    f.0[j] <- 2 * dnorm(0,0, 1/sigma[j]^2) 

    y[j] ~ dbin(pcap[j], N[j]) 
    N[j] ~ dpois(lambda[j]) 
    
    lambda[j] <- exp(alpha + bHQ*hq[j])
    }
    
    # Derived parameters
    Ntotal <- sum(N) 
    
    }",fill=TRUE, file = "3.HDS_sig[j]_lam[hq[inf]].txt")
