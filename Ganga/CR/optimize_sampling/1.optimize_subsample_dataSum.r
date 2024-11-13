## -------------------------------------------------
##               1.1. Optimize sampling
##                      Remove %
##            Save data summaries while running
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

# Load model because it doesn't work
setwd("D:/MargSalas/Ganga/Data/CMR")
sink("model_m0_pSex.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
psi ~ dunif(0, 1)

for (s in 1:2){
  p[s] ~ dunif(0, 1)
}

# Likelihood
for (i in 1:M){
   sex[i] ~ dbern(psi)	
   z[i] ~ dbern(omega)			# Inclusion indicators (probability that exists)
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[sex[i]+1]		# Can only be detected if z=1
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}",fill = TRUE)
sink()


## ---- Run subsampling different % ----

per <- c(0.5, 0.75, 0.80, 0.90, 0.95)
per_table <- c("50%", "75%", "80%", "90%", "95%", "100%")
  

iter <- 100

results_per <- list()
results_iter <- list()

# To store data summaries:
explore <- data.frame(per = NA, iter = NA, Ndetect = NA, Nrecapt = NA, perc_recapt = NA)


for(p in 1:length(per)){
  
  for (i in 1:iter){
    
    ## Subsample
    rd <- sample(1:122, size = 122*per[p])
    ganga_sub <- ganga[rd,]
    
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
    n_oc <- 5
    
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
    results_iter[[i]] <- out_sim_M0_psex_sub
    
    # Store data summary
    capt.hist.2 <- capt.hist %>%
      mutate(oc1 = str_sub(ch, 1,1),
             oc2 = str_sub(ch, 2,2),
             oc3 = str_sub(ch, 3,3),
             oc4 = str_sub(ch, 4,4),
             oc5 = str_sub(ch, 5,5)) %>%
      mutate_at(c('oc1', 'oc2','oc3', 'oc4', 'oc5'), as.numeric) %>%
      select(-ch)
    n_recap <- length(which(apply(capt.hist.2[,c(1:5)], 1, sum) > 1))
    perc_recapt <- round(length(which(apply(capt.hist.2[,c(1:5)], 1, sum) > 1))/nrow(capt.hist.2), 2)
    
    explore[nrow(explore) + 1, ] <- c(per_table[p], i, nrow(capt.hist), n_recap, perc_recapt)
    
  }
  
  results_per[[p]] <- results_iter
  
}

## ---- Run also a 100 times the original model with all samples ----

for (i in 1:iter){
  
  ## Create capture history (description in Arrange script)
  capt.hist <- ganga %>%
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
  n_oc <- 5
  
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
  results_iter[[i]] <- out_sim_M0_psex_sub
  
  # Store data summary
  capt.hist.2 <- capt.hist %>%
    mutate(oc1 = str_sub(ch, 1,1),
           oc2 = str_sub(ch, 2,2),
           oc3 = str_sub(ch, 3,3),
           oc4 = str_sub(ch, 4,4),
           oc5 = str_sub(ch, 5,5)) %>%
    mutate_at(c('oc1', 'oc2','oc3', 'oc4', 'oc5'), as.numeric) %>%
    select(-ch)
  
  n_recap <- length(which(apply(capt.hist.2[,c(1:5)], 1, sum) > 1))
  perc_recapt <- round(length(which(apply(capt.hist.2[,c(1:5)], 1, sum) > 1))/nrow(capt.hist.2), 2)
  
  explore[nrow(explore) + 1, ] <- c("100%", i, nrow(capt.hist), n_recap, perc_recapt)
  
}

# Store as 6th element of the list
results_per[[6]] <- results_iter

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/optimize_sampling/results")
save(results_per, file = "results_per2.RData")
save(explore, file = "data_sum_subsample.RData")

