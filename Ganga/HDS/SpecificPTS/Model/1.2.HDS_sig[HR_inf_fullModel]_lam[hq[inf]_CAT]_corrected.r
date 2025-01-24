## ---- JAGS model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
cat("model{
    
    # PRIORS
    
    # Priors for lambda
    alpha ~ dnorm(0, 0.001)
    
    bHQ1 ~ dnorm(mu.bHQ1, sig.bHQ1) # Informative prior on habitat quality
    bHQ2 ~ dnorm(mu.bHQ2, sig.bHQ2) # Informative prior on habitat quality
    bHQ3 ~ dnorm(mu.bHQ3, sig.bHQ3) # Informative prior on habitat quality


    # Priors for sigma (all informative)
    alpha.sig ~ dnorm(mu.alpha, sig.alpha)
    b.hSquare ~ dnorm(mu.hSquare, sig.hSquare)
    b.jd ~ dnorm(mu.jd, sig.jd)
    b.jdSquare ~ dnorm(mu.jdSquare, sig.jdSquare)

    # Prior for beta
    log.b ~ dnorm(mu.b, sig.b)
    b <- exp(log.b)
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[site.dclass[i],1:nG])  
    }
    
    for(j in 1:length(y)){ 
    
    sigma[j] <- exp(alpha.sig + b.hSquare * hSquare[j] +
                b.jd * jd[j] + b.jdSquare * jdSquare[j])

  
    for(k in 1:nG){ 
    
    p[j,k]<-1-exp(-(midpt[k]/sigma[j])^-b)
    pi[j,k] <- int.w[k] / strip.width 
    fc[j,k]<- p[j,k] * pi[j,k]                 ## pi=percent area of k; drops out if constant
    fct[j,k]<-fc[j,k]/pcap[j]
    }
    
    pcap[j] <- sum(fc[j, 1:nG]) # Same in all sites

    y[j] ~ dbin(pcap[j], N[j]) 
    N[j] ~ dpois(lambda[j]) 
    
    lambda[j] <- exp(alpha +  bHQ1*hqCov1[j] + bHQ2*hqCov2[j] + bHQ3*hqCov3[j])
    }
    
    # Derived parameters
    Ntotal <- sum(N) 
    
    }",fill=TRUE, file = "1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]_CAT]_corrected.txt")