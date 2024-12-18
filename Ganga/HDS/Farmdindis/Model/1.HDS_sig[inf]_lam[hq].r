
## ---- JAGS model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
cat("model{
    
    # PRIORS
    
    # Priors for lambda
    alpha ~ dnorm(0, 0.001)
    bHQ ~ dnorm(0, 0.001)
    
    # Informative priors for sigma
    
    ###########################################################################
    # RS: either use detection estimates from previous model as fixed parameters
    # meaning, provide them both as data (log.sig=mu.sig, bTemp=mu.bTemp)
    # or give both of them an informative prior, like you do here for bTemp
    
    log.sig ~ dnorm(mu.sig, sig.sig) #: Deterministic coefficient (mu and sig come from model 2010-2020)
    
    bTemp ~ dnorm(mu.bTemp, sig.bTemp)
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[site.dclass[i], 1:nG])  
    }
    
    for(j in 1:length(y)){ 
    
    sigma[j] <- exp(log.sig + bTemp*temp[j])
    
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
    
    }",fill=TRUE, file = "1.HDS_sig[inf]_lam[hq].txt")