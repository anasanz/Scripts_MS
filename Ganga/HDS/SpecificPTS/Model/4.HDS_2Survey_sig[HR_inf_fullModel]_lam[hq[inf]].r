## ---- JAGS model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
cat("model{
    
    # PRIORS
    
    # Priors for lambda
    alpha ~ dnorm(0, 0.001)
    bHQ ~ dnorm(mu.bHQ, sig.bHQ) # Informative prior on habitat quality
    
    
    # Priors for sigma (all informative)
    alpha.sig ~ dnorm(mu.alpha, sig.alpha)
    b.hSquare ~ dnorm(mu.hSquare, sig.hSquare)
    b.jd ~ dnorm(mu.jd, sig.jd)
    b.jdSquare ~ dnorm(mu.jdSquare, sig.jdSquare)
    log.Vebro ~ dnorm(mu.Vebro, sig.Vebro)
    
    # Prior for beta
    log.b ~ dnorm(mu.b, sig.b)
    b <- exp(log.b)
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[site.dclass[i],1:nG, survey.dclass[i]])  
    }
    
    for(j in 1:nSites){ 
    for(s in 1:nSurveys){
    
    sigma[j,s] <- exp(alpha.sig + b.hSquare * hSquare[j,s] +
                b.jd * jd[j,s] + b.jdSquare * jdSquare[j,s] +
                log.Vebro)

  
    for(k in 1:nG){ 
    
    p[j,k,s]<-1-exp(-(midpt[k]/sigma[j,s])^-b)
    pi[j,k,s] <- int.w[k] / strip.width 
    fc[j,k,s]<- p[j,k,s] * pi[j,k,s]                 ## pi=percent area of k; drops out if constant
    fct[j,k,s]<-fc[j,k,s]/pcap[j,s]
    }
    
    pcap[j,s] <- sum(fc[j, 1:nG, s]) 
    y[j,s] ~ dbin(pcap[j,s], N[j])
    }
    
    N[j] ~ dpois(lambda[j]) 
    lambda[j] <- exp(alpha + bHQ*hq[j])
    }
    
    # Derived parameters
    Ntotal <- sum(N) 
    
    }",fill=TRUE, file = "4.HDS_2Survey_sig[HR_inf_fullModel]_lam[hq[inf]].txt")