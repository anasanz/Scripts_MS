## -------------------------------------------------
##               2.1. Optimize sampling
##                  Remove occasions
## ------------------------------------------------- 

# Subset different parts of the dataset to assess the precision of estimates
# under different sampling scenarios

rm(list=ls())

library(tidyverse)
library(tidyr)
library(magrittr)
library(jagsUI)
library(rjags)

## ---- Load data ----

setwd("D:/MargSalas/Ganga/Data/CMR")

ganga <- readxl::read_excel("sandgrouse_Individual ID+sex_FINAL.xlsx") %>%
  select(7,8,16) %>%
  rename("id" = "INDIVIDUAL ID") %>%
  rename( "Occasion"= "Id_Visit")

ganga <- ganga[which(!is.na(ganga$id)), ] # Select only identified individuals and relevant columns for capture history
ganga <- ganga[which(ganga$Year == 2022), -1] # Select year 2022
ganga$detect <- 1 # 122 detections

# Load sex
setwd("D:/MargSalas/Ganga/Data/CMR")
load("id_sex_sandgrouse_2022.RData") # It is in the same order than capt.hist
id_sex$sex[which(id_sex$sex == "F")] <- 0
id_sex$sex[which(id_sex$sex == "M")] <- 1
id_sex$sex[which(id_sex$sex == "X")] <- NA

## ---- Run subsampling removing different number of occasions ----

n_occasions <- 1:max(ganga$Occasion)
oc_remove <- c(1,2,3)

results_oc_remove <- list()
results_combi <- list()

for (r in 1:length(oc_remove)){ # This is equivalent to the percentage when I remove a % of samples
  
  comb_remove <- combn(n_occasions, oc_remove[r]) # Combinations when removing x number of occasions
  
  for (i in 1:dim(comb_remove)[2]){ # This is equivalent to the iterations when I remove a % of samples
    
    rd <- comb_remove[,i] # This is the occasion/ocassions to remove
    
    ## Subsample
    ganga_sub <- ganga[-which(ganga$Occasion %in% rd),]
    
    ## Create capture history (description in Arrange script)
    capt.hist <- ganga_sub %>%
      distinct() %>%
      spread(Occasion, detect, fill = 0) %>% 
      group_by(id) %>%
      unite("ch", 2:tail(names(.),1), sep = "")
    
    capt.hist$ch <- as.character(capt.hist$ch)
    
    ## Arrange for model
    
    # Determine sex of subsamples
    capt.hist <- left_join(capt.hist, id_sex, by = "id")
    sex <- as.numeric(capt.hist$sex)
    
    # Format capture history
    capt.hist <- as.data.frame(capt.hist$ch)
    colnames(capt.hist)[1] <- "ch"
    n_oc <- nchar(capt.hist$ch[1]) ### Number of characters ch is number of occasions
    
    data_ganga <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = n_oc)
    T = ncol(data_ganga)
    
    # Augment data set by 150 potential individuals
    nz <- 150
    yaug <- rbind(data_ganga, array(0, dim = c(nz, T)))
    sexAug <- c(sex, rep(NA, nz))
    
    # Bundle data
    win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug), sex = sexAug)
    
    # Initial values
    inits <- function() list(z = rep(1, nrow(yaug)), p = runif(2, 0, 1),
                             sex = c(rep(NA,nrow(data_ganga)), rbinom(nz,1,0.5)))
    # Parameters monitored
    params <- c("N", "p", "omega", "psi", "z", "sex")
    
    # MCMC settings
    ni <- 2500
    nt <- 1
    nb <- 500
    nc <- 3
    
    ## Call WinBUGS from R
    out_sim_M0_psex_sub <- jags(win.data, inits, params, "model_m0_pSex.txt", n.chains = nc, 
                                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    results_combi[[i]] <- out_sim_M0_psex_sub
  }
  results_oc_remove[[r]] <- results_combi
}

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/optimize_sampling/results")
save(results_oc_remove, file = "results_oc_remove.RData")


