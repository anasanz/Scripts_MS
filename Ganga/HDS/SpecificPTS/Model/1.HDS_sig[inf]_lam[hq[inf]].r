
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
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nG])  
    }
    
  for(k in 1:nG){ 
    up[k]<-pnorm(db[k+1], 0, 1/sigma^2) ##db are distance bin limits
    low[k]<-pnorm(db[k], 0, 1/sigma^2) 
    p[k]<- 2 * (up[k] - low[k])
    pi[k] <- int.w[k] / strip.width 
    f[k]<- p[k]/f.0/int.w[k]                   ## detection prob. in distance category k                      
    fc[k]<- f[k] * pi[k]                 ## pi=percent area of k; drops out if constant
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
    
    }",fill=TRUE, file = "1.HDS_sig[inf]_lam[hq[inf]].txt")