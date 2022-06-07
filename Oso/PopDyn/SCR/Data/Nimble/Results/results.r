
## -------------------------------------------------
##          Check results different models
## ------------------------------------------------- 

rm(list = ls())

library(MCMCvis)

## ---- Load data from single year SCR with density covariate ----

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/Nimble/Results/SCRdenscov_year")
load("sampSCRdenscov2017.RData")
samp17 <- samp
summ17 <- MCMCsummary(samp17)
MCMCtrace(samp17)

load("sampSCRdenscov2018.RData")
samp18 <- samp
summ18 <- MCMCsummary(samp18)
MCMCtrace(samp18)


setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/Nimble/Results/openSCRdenscov")
load("sampOpenSCR.RData")
sampOpen <- samp2
summOpen <- MCMCsummary(sampOpen)
MCMCtrace(sampOpen)

write.csv(summOpen, file = "summOpen.csv")

summOpen

