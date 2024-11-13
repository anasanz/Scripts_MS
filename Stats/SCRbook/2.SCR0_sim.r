
## -------------------------------------------------
##                 SCR0 MODEL
## ------------------------------------------------- 

rm(list=ls())

library(SCRbayes)
library(rjags)
library(jagsUI)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions/scrbook")
source("simSCR0.r")

## ---- Simulate data ----

set.seed(2013)
# Create 5 x 5 grid of trap locations with unit spacing 
traplocs <- cbind(sort(rep(1:5,5)),rep(1:5,5)) 
ntraps <- nrow(traplocs) 
# Compute distance matrix: Distance from one trap to others
Dmat <- e2dist(traplocs,traplocs)

# Define state-space of point process. (i.e., where animals live). 
# "buffer" just adds a fixed buffer to the outer extent of the traps. 
# 
buffer <- 2 
xlim <- c(min(traplocs[,1] - buffer),max(traplocs[,1] + buffer)) 
ylim <- c(min(traplocs[,2] - buffer),max(traplocs[,2] + buffer))

N <- 100 # population size 
K <- 20 # number nights of effort

# Simulate activity centers
sx <- runif(N,xlim[1],xlim[2])
sy <- runif(N,ylim[1],ylim[2]) 
S <- cbind(sx,sy)

# Compute distance matrix: 
D <- e2dist(S,traplocs) # distance of each individual from each trap

# Define parameters of encounter probability
alpha0 <- -2.5
sigma <- 0.5  # scale parameter of half-normal
alpha1 <- 1/(2*sigma*sigma) # convert to coefficient on distance

# Compute Probability of encounter: p(x,s) = p0 * exp(-1/sigma^2*D^2)
probcap <- plogis(-2.5)*exp(-alpha1*D*D) # One per individual and trap

# Generate the encounters of every individual in every trap 
Y <- matrix(NA, nrow = N, ncol = ntraps) 
for(i in 1:nrow(Y)){ 
  Y[i,] <- rbinom(ntraps,K,probcap[i,])
  } # Each element of y is the frequency of encounter (out of K) 
    # of individuals in traps


# Generate the encounters of every individual in every trap and occasion
Y <- array(NA,dim=c(N,ntraps,K))
for(i in 1:nrow(Y)){ 
  for(j in 1:ntraps){ 
    Y[i,j,1:K] <- rbinom(K,1,probcap[i,j])
} }


## ---- Simulating function ----

# To generate this data there is the function simSCR0:
#### 2D ####

data <- simSCR0(N = 100, K = 20, alpha0 = -2.5, sigma = 0.5, 
                discard0 = TRUE, array3d = FALSE, rnd = 2013)
# Discard 0 = TRUE removes all non-detected individuals

Y <- data$Y # ENcounter histories
traplocs <- data$traplocs # Trap locations

#### 3D ####
# array 3d = TRUE: individual encounter histories
#     -> One matrix of nind*ntraps per occasion

data1 <- simSCR0(N = 100, K = 20, alpha0 = -2.5, sigma = 0.5, 
                discard0 = FALSE, array3d = TRUE, rnd = 2013)

# To recover the 2-d matrix from the 3-d array, 
# and subset the 3-d array to individuals that were captured, we do this:

# Sum over the “sample occasions” dimension (3rd margin of the array) 
Y1 <- data1$Y
Y2d <- apply(Y1,c(1,2),sum)
# Compute how many times each individual was captured 
ncaps <- apply(Y2d,1,sum)
# Keep those individuals that were captured 
Y2 <- Y1[ncaps > 0, ,]

## ---- Fit simSCR0 in jags ----

#### MODEL WITH KNOWN N ####

# Use simulator to grab a data set and then harvest the elements of the resulting object for further analysis #

data <- simSCR0(discard0 = FALSE, rnd = 2013) 
y <- data$Y 
traplocs <- data$traplocs
# In this case nind=N because we’re doing the known-N problem # 
nind <- nrow(y) 
X <- data$traplocs 
J <- nrow(X) # number of traps
K <- data$K 
xlim <- data$xlim 
ylim <- data$ylim

# Model 
setwd("D:/MargSalas/Oso/SCR/Model")

cat(" model{ 
alpha0 ~ dnorm(0,.1) 
logit(p0) <- alpha0 
alpha1 ~ dnorm(0,.1) 
sigma <- sqrt(1/(2*alpha1)) 

for(i in 1:N){ # note N here -- N is KNOWN in this example (but we don't know where they are!-> s) 
  s[i,1] ~ dunif(xlim[1],xlim[2]) 
  s[i,2] ~ dunif(ylim[1],ylim[2]) 
    for(j in 1:J){ 
      d[i,j] <- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5) 
      y[i,j] ~ dbin(p[i,j],K) 
      p[i,j] <- p0*exp(- alpha1*d[i,j]*d[i,j]) }
  } 
    }", fill = TRUE, file = "SCR0a.txt")

## Starting values for activity centers, s
sst <- cbind(runif(nind,xlim[1],xlim[2]),runif(nind,ylim[1],ylim[2]))

for(i in 1:nind){ 
  if(sum(y[i,])==0) next 
  # for the observed individuals, we replace those values by each individual’s mean trap coordinate of its encounters
  sst[i,1] <- mean( X[y[i,]>0,1] ) 
  sst[i,2] <- mean( X[y[i,]>0,2] )
}

# Data
data <- list (y=y, X=X, K=K, N=nind, J=J, xlim=xlim, ylim=ylim) 

# Initial values
inits <- function(){ list (alpha0=rnorm(1,-4,.4), alpha1=runif(1,1,2), s=sst)}

# Parameters
parameters <- c("alpha0","alpha1","sigma") 

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Run model

out <- jags(data, inits, parameters, "SCR0a.txt", n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
print(out)

#### MODEL WITH UNKNOWN N -> DATA AUGMENTATION ####

## Simulate the data and extract the required objects ## 
data <- simSCR0(discard0=TRUE,rnd=2013) 
y <- data$Y 
nind <- nrow(y) 
X <- data$traplocs 
K <- data$K 
J <- nrow(X) 
xlim <- data$xlim 
ylim <- data$ylim

## Data augmentation 
M <- 200 
y <- rbind(y, matrix(0, nrow = M-nind, ncol=ncol(y))) 
z <- c(rep(1, nind),rep(0, M-nind)) # Augmented values

## Starting values for s 
sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2])) 

# Helpful to start the si for each observed individual at or near the trap(s) 
# it was captured. So for the first rows of M (nind, captured individuals) I get the mean

for(i in 1:nind){ 
  sst[i,1] <- mean( X[y[i,]>0,1] ) 
  sst[i,2] <- mean( X[y[i,]>0,2] )
}

# Model

setwd("D:/MargSalas/Oso/SCR/Model")

cat(" model{ 
alpha0 ~ dnorm(0,.1) 
logit(p0) <- alpha0 
alpha1 ~ dnorm(0,.1) 
sigma <- sqrt(1/(2*alpha1))
psi ~ dunif(0,1)

for(i in 1:M){
  z[i] ~ dbern(psi)
  s[i,1] ~ dunif(xlim[1],xlim[2]) 
  s[i,2] ~ dunif(ylim[1],ylim[2]) 
    for(j in 1:J){ 
      d[i,j] <- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5) 
      y[i,j] ~ dbin(p[i,j],K) 
      p[i,j] <- z[i]*p0*exp(- alpha1*d[i,j]*d[i,j]) }
}
  N <- sum(z[])
  D <- N/64
    }", fill = TRUE, file = "SCR0b_DA.txt")

# Data
data <- list (y=y, X=X, K=K, M=M, J=J, xlim=xlim, ylim=ylim) 
# Initial values: sst previously calculated
inits <- function(){ list (alpha0=rnorm(1,-4,.4), alpha1=runif(1,1,2), s=sst, z=z) }
# Parameters
parameters <- c("alpha0","alpha1","sigma","N","D")
# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

out_da <- jags(data, inits, parameters, "SCR0b_DA.txt", n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
print(out_da)
