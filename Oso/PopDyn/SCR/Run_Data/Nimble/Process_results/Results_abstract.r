
## -------------------------------------------------
##                 RESULTS ABSTRACT EURING
## ------------------------------------------------- 


# Two last models of OPSCR:
# - Model with age structure
# - Model with effort covariates


rm(list = ls())

library(MCMCvis)
library(coda)
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/DoScale.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/ProcessCodaOutput.R")


## -------------------------------------------------
##        openSCRdenscov + Age structure + different trap arrays 
##          2017 - 2021 FINAL DATA (Cross check Elena)
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
load("sampOpenSCR_diftraps_age_final1721.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

out2 <- ProcessCodaOutput(out.list)
out1$mean
out2$mean
# Convergence

library(MCMCvis)

MCMCtrace(out.list,
          ind = TRUE,
          pdf = TRUE, 
          open_pdf = FALSE, 
          filename = 'OpenSCR_diftraps_age_2021_Traceplot', 
          wd = 'D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021')

out$mean

summOpen1 <- MCMCsummary(out.list)




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

