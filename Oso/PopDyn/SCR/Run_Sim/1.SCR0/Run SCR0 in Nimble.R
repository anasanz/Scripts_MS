
library(nimble)
library(MCMCvis)

### generate some basic SCR data

##trap array
X<-as.matrix(expand.grid(seq(-3,3,1), seq(-3,3,1)))
J<-dim(X)[1]

##detection parameters
#sigma=movement
sigma<-1
#p0=baseline detection
p0<-0.2

##state space coordinates
xmin<-min(X[,1]-3*sigma) # Does it multiply by sigma to kind of cover all home ranges? 
ymin<-min(X[,2]-3*sigma)
xmax<-max(X[,1]+3*sigma)
ymax<-max(X[,2]+3*sigma)

##true abundance within state space
N<-50

##random activity centers
Sx<-runif(N, xmin, xmax)
Sy<-runif(N, ymin, ymax)

##distances from activity centers to traps
e2dist <- function (x, y){ # x = activity centers, y = traps
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }

D <- e2dist(cbind(Sx, Sy), X)

##generate detection data
K <- 5 #number of sampling occasions
obs <- matrix(NA, N, J)

for (i in 1:N){ 
  p.eff<-p0*exp(-D[i,]^2/(2*sigma^2))
  obs[i,]<-rbinom(J, K, p.eff) # Simulate one observation per trap (J = 50 values) with specific p and 5 trials
}

##how many individuals detected?
n <- sum(apply(obs, 1, sum)>0) # Sum rows (individuals)

##subset obs to include only individuals detected at least once
seen<-which(apply(obs, 1, sum)>0) # Which rows (individuals)
Y<-obs[seen,]

##augment observed data to size M
M <- 80 
y.in <- rbind(Y, 
            matrix(0, nrow=M-n, ncol=J))

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

##alternative: start observed individuals off as alive, others randomly
##note: in complex models it may be necessary to calculate starting values for
##activity centers as well, based on average observed location

inits <- function(){list(psi= runif(1,0.4, 0.6), 
                       sigma = runif(1,0.5, 1.5),
                       p0 = runif(1,0,1),
                       z = z.in)}

##Compile data for Nimble in two objects: data (things that are random variables)
## and constants. This distinction is not so important in simple models but may
## have implications for more complex models

dat <- list(y = y.in)
consts <- list(M = M, J = J, X = X, 
             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
             K = K)

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Sim")
source('SCR in Nimble.R')

##determine which parameters to monitor
params<-c('N', 'psi', 'sigma', 'p0')

#(1) set up model

model <- nimbleModel(SCR0, constants = consts, 
                     data=dat, check = FALSE)

model$initializeInfo() # Because we didn't set initial values of S?
help(modelInitialization)

#(2) Compile model in c++
#     In complex models, this step can take a while (as well as step 5)
#     Much longer than in JAGS, but the model typically runs much faster
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
#     Also, nimble allows two sets of monitors, these can be thinned at different rates
#     all of which is more important in complex models but not to start with
conf.mcmc<-configureMCMC(model, monitors = params, thin=1)
# For example to change the sampler

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)
# One of the reasons why nimble is faster than JAGS is because it compiles
# in C++

# (6) Run (monitor time just for fun)
system.time(
  (samp <- runMCMC(cmcmc, niter = 2500, nburnin = 1000, nchains=3, inits = inits) )
)

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ<-MCMCsummary(samp)
MCMCtrace(samp)

# Additional steps are used to trouble shoot. If you have an error its
# useful to go through the steps
# To be more efficient
# More of the time, the direct line is enough


#####
# Run with nimbleMCMC
mcmc.output <- nimbleMCMC(code = SCR0,
                            data = dat,
                            constants = consts,
                            inits = inits,
                            monitors = params,
                            thin = 1,
                            niter = 2500,
                            nburnin = 1000,
                            nchains = 3
                            )
str(mcmc.output)
head(mcmc.output$chain1)
MCMCsummary(mcmc.output)
par(mfrow = c(1,1))
MCMCtrace(mcmc.output,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE)

# Every time you run the MCMC you get slightly different values
# to get rid of the randomness: see video Olivier to set the seed
# Very useful to get always the same when you are still building the model!!

