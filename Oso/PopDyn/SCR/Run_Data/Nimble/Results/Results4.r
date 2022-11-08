
## -------------------------------------------------
##                   RESULTS 4
##          Check results different models 
## - OPen models with effort without age (scripts 4)
## ------------------------------------------------- 

rm(list = ls())

library(MCMCvis)
library(coda)
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/DoScale.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/ProcessCodaOutput.R")

## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##                  p0 = effort 
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
load("sampOpenSCR_diftraps_effortCov_1721_FINALDATA_3d.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,21] <- 0
  out.list[[i]][1,7] <- 0
}

out2 <- ProcessCodaOutput(out.list)

MCMCtrace(out.list,   # Convergence: beffort2 and sigma do not mix great
          ind = TRUE,
          pdf = FALSE)

summ_Effort <- MCMCsummary(out.list)

# Observations: 

# - b.effort2 (effort in Spain) does not converge

## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##              p0 = effort + trap
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
load("sampOpenSCR_diftraps_effortTrapCov_1721_FINALDATA_3d.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

MCMCtrace(out.list,   # Convergence: beffort2 and sigma do not mix great
          ind = TRUE,
          pdf = FALSE)

summ_EffortTrap <- MCMCsummary(out.list)

# Observations: 
# - Things that have swifted since we included trap:
#    -> Super positive effect of trap (both)
#    -> Distcore is not affecting anymore 
#    -> Negative effect of trap in Spain on p0
#    * This happens because there were many traps of type "both" in Spain, so
#      effort was counfounding
# - b.effort2 (Spain) doesn't converge, but is the best model, maybe because it includes trap


## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##              p0 = effort + trap + m (0-6)
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
load("sampOpenSCR_diftraps_effortTrapSoccCov_1721_FINALDATA_3d.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

MCMCtrace(out.list,   # Convergence: beffort2 and sigma do not mix great
          ind = TRUE,
          pdf = FALSE)

summ_EffortTrapSocc <- MCMCsummary(out.list)

# Observations: 
# - Doesn't converge at all

## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##              p0 = effort + trap + m (1-7)
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
load("sampOpenSCR_diftraps_effortTrapSocc2Cov_1721_FINALDATA_3d.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

MCMCtrace(out.list,   # Convergence: beffort2 and sigma do not mix great
          ind = TRUE,
          pdf = FALSE)

summ_EffortTrapSocc <- MCMCsummary(out.list)

# Observations: 
# - Doesn't converge at all, REMOVE MONTH

## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##                  REPARAMETRIZED
##p0 = effort (p0[1] = France; p0[2] = Spain; beta = 2 visits France)
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
load("sampOpenSCR_diftraps_effortCov2p0_1721_FINALDATA.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

MCMCtrace(out.list,   # Convergence: beffort2 and sigma do not mix great
          ind = TRUE,
          pdf = FALSE)

summ_Effort2p0<- MCMCsummary(out.list)

# Observations: 
# - Convergence is even worst: P0[2] (Spain, as before) and now also sigma dont converge
# - Maybe better with trap?
# - Effects make sense because trap covariate is not included


## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##              p0 = effort + trap + behav. response
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
load("sampOpenSCR_diftraps_effortTrapBhCov_1721_FINALDATA_3d.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")

MCMCtrace(out.list,   # Convergence: beffort2 and sigma do not mix great
          ind = TRUE,
          pdf = TRUE)

summ_EffortTrapBh <- MCMCsummary(out.list)

## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##              p0 = effort + trap + behav. response
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
load("sampOpenSCR_diftraps_effortTrapBhCov_1721_FINALDATA_3d.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")

MCMCtrace(out.list,   # Convergence: beffort2 and sigma do not mix great
          ind = TRUE,
          pdf = TRUE)

summ_EffortTrapBh <- MCMCsummary(out.list)


## -------------------------------------------------
##    openSCRdenscov + different trap arrays
##             2017 - 2021 FINAL DATA 
##              p0 = effort + trap + behav. response
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/5.OPSCR_sigma")
load("Results_Model5-2.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/5.OPSCR_sigma")

MCMCtrace(out.list,   
          ind = TRUE,
          pdf = TRUE)

summ_EffortTrapBh_sigSex <- MCMCsummary(out.list)
