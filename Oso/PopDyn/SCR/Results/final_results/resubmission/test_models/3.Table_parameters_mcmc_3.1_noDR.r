
## -------------------------------------------------
##   Table with parameter results + recruitment
## ------------------------------------------------- 

rm(list = ls())

library(MCMCvis)
library(nimbleSCR)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.1_noDR_allparams_FINAL")
load("myResults_RESUB_3-3.1_noDR_param.RData")
sum <- summary(nimOutput)


param <- c("psi", "omega", "phi.cub", "phi.sub", "phi.ad", "beta.dens", "sigD", "sigma[1]", "sigma[2]", 
           "p.cub", "p.sub", "p.ad", "b.bh", "trapBetas[1]", "trapBetas[2]", "trapBetas[3]")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.1_noDR_allparams_FINAL")
MCMCtrace(nimOutput, 
          params = param,
          ind = TRUE,
          pdf = TRUE,
          filename = "mcmc_m1")

res <- as.data.frame(matrix(nrow = length(param), ncol = 6))
colnames(res) <- c("Parameters", "Mean", "SD", "2.5%", "50%", "9.7%")
res$Parameters <- param

for(i in 1:nrow(res)){
  res[i,2] <- sum$statistics[rownames(sum$statistics) %in% param[i], colnames(sum$statistics) %in% "Mean"]
  res[i,3] <- sum$statistics[rownames(sum$statistics) %in% param[i], colnames(sum$statistics) %in% "SD"]
  res[i,4] <- sum$quantiles[rownames(sum$quantiles) %in% param[i], colnames(sum$quantiles) %in% "2.5%"]
  res[i,5] <- sum$quantiles[rownames(sum$quantiles) %in% param[i], colnames(sum$quantiles) %in% "50%"]
  res[i,6] <- sum$quantiles[rownames(sum$quantiles) %in% param[i], colnames(sum$quantiles) %in% "97.5%"] 
}

res[,c(2:6)] <- round(res[,c(2:6)],3)

res$Parameters[res$Parameters %in% "sigma[1]"] <- "sigmaF"
res$Parameters[res$Parameters %in% "sigma[2]"] <- "sigmaM"
res$Parameters[res$Parameters %in% "trapBetas[1]"] <- "b.effort1"
res$Parameters[res$Parameters %in% "trapBetas[2]"] <- "b.effort2"
res$Parameters[res$Parameters %in% "trapBetas[3]"] <- "b.effort3"
res$Parameters[res$Parameters %in% "beta.dens"] <- "b.DistCore"

# Transform sigma and tau to the scale of the state space: Do it by multiplying by the resolution, 5

res[which(res$Parameters %in% c("sigmaF", "sigmaM", "sigD")), c(2:6)] <- res[which(res$Parameters %in% c("sigmaF", "sigmaM", "sigD")), c(2:6)] * 5

res1 <- res

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Tables")
#openxlsx::write.xlsx(res, 'params.xlsx')

## ---- Recruitment ----

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("pcr_all_fem.RData")


mean_pcr <- apply(pcr_all_fem$pcrmat,1,mean) # Mean pcr per ite

lci <- quantile(mean_pcr,probs = 0.025) 
uci <- quantile(mean_pcr,probs = 0.975)


