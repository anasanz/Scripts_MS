## ---- JAGS model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
cat("model{
    
    # PRIORS
    
    # Priors for lambda
    alpha ~ dnorm(0, 0.001)
    bHQ ~ dnorm(mu.bHQ, sig.bHQ) # Informative prior on habitat quality
    
    
    # Priors for sigma
    log.sig ~ dnorm(mu.sig, sig.sig) # Informative priors for sigma
    
    # SIGMA is not specific to transect, only one value
    sigma <- exp(log.sig)
    
    # Prior for beta
    log.b ~ dnorm(mu.b, sig.b)
    b <- exp(log.b)
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nG])  
    }
  
    for(k in 1:nG){ 
    
    p[k]<-1-exp(-(midpt[k]/sigma)^-b)
    pi[k] <- int.w[k] / strip.width 
    fc[k]<- p[k] * pi[k]                 ## pi=percent area of k; drops out if constant
    fct[k]<-fc[k]/pcap
    }
    
    pcap <- sum(fc[1:nG]) # Same in all sites
    
    f.0 <- 2 * dnorm(0,0, 1/sigma^2) 
    
    for(j in 1:length(y)){ 

    y[j] ~ dbin(pcap, N[j]) 
    N[j] ~ dpois(lambda[j]) 
    
    lambda[j] <- exp(alpha + bHQ*hq[j])
    }
    
    # Derived parameters
    Ntotal <- sum(N) 
    
    }",fill=TRUE, file = "1.1.HDS_sig[HR_inf]_lam[hq[inf]].txt")