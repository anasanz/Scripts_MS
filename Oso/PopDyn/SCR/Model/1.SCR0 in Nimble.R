
SCR0<-nimbleCode({
  
    ###data augmentation parameter
    psi ~ dunif(0,1)
    ##detection parameters - sigma (movement)
    sigma ~ dunif(0,5) #adjust to units of trap array and space use of species
    ##detection parameters - lam0 (baseline detection probability)
    p0 ~ dunif(0,1)
    
  ####observation model
      for (i in 1:M){
        #activity centers (min and max coordinates for state space)
        S[i,1] ~ dunif(xmin, xmax)
        S[i,2] ~ dunif(ymin, ymax)
        #alive state
        z[i] ~ dbern(psi)
        
        for (j in 1:J){
          #distance, detection function, effective detection probability
          D[i,j] <- sqrt((X[j,1]-S[i,1])^2 + (X[j,2]-S[i,2])^2)
          p[i,j] <- p0 * exp(-D[i,j]^2/(2*sigma^2))
          p.eff[i,j] <- z[i]*p[i,j]
          #detection model
          y[i,j] ~ dbinom(p.eff[i,j], K)# K = number of occasions
        }}
    #total abundance in state-space
    N <- sum(z[1:M]) # Vectorized right?? Entire range. You need to specify in nimble
})


