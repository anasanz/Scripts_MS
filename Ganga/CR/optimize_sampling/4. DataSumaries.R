
## -------------------------------------------------
##  Data summaries: Removing occasions and subsample
## ------------------------------------------------- 

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

## ---- 1. Data summary removing different number of occasions ----

n_occasions <- 1:max(ganga$Occasion)
oc_remove <- c(1,2,3)

explore <- data.frame(N_ocRem = NA, ocID = NA, Ndetect = NA, Nrecapt = NA, perc_recapt = NA)

# Removing 1 occasion: 
  r = 1
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
    
    capt.hist.2 <- capt.hist %>%
      mutate(oc1 = str_sub(ch, 1,1),
             oc2 = str_sub(ch, 2,2),
             oc3 = str_sub(ch, 3,3),
             oc4 = str_sub(ch, 4,4)) %>%
      mutate_at(c('oc1', 'oc2','oc3', 'oc4'), as.numeric) %>%
      select(-ch)
    
    n_recap <- length(which(apply(capt.hist.2[,c(2:5)], 1, sum) > 1))
    perc_recapt <- round(length(which(apply(capt.hist.2[,c(2:5)], 1, sum) > 1))/nrow(capt.hist.2), 2)

    explore[nrow(explore) + 1, ] <- c(oc_remove[r], as.character(rd), nrow(capt.hist), n_recap, perc_recapt)
  }
  
  # Removing 2 occasions: 
  r = 2
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
    
    capt.hist.2 <- capt.hist %>%
      mutate(oc1 = str_sub(ch, 1,1),
             oc2 = str_sub(ch, 2,2),
             oc3 = str_sub(ch, 3,3)) %>%
      mutate_at(c('oc1', 'oc2','oc3'), as.numeric) %>%
      select(-ch)
    
    n_recap <- length(which(apply(capt.hist.2[,c(2:4)], 1, sum) > 1))
    perc_recapt <- round(length(which(apply(capt.hist.2[,c(2:4)], 1, sum) > 1))/nrow(capt.hist.2), 2)
    rd_name <- paste(rd[1], "&", rd[2], sep = "")
    explore[nrow(explore) + 1, ] <- c(oc_remove[r], rd_name, nrow(capt.hist), n_recap, perc_recapt)
  }
  
  # Removing 3 occasions: 
  r = 3
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
    
    capt.hist.2 <- capt.hist %>%
      mutate(oc1 = str_sub(ch, 1,1),
             oc2 = str_sub(ch, 2,2)) %>%
      mutate_at(c('oc1', 'oc2'), as.numeric) %>%
      select(-ch)
    
    n_recap <- length(which(apply(capt.hist.2[,c(2:3)], 1, sum) > 1))
    perc_recapt <- round(length(which(apply(capt.hist.2[,c(2:3)], 1, sum) > 1))/nrow(capt.hist.2), 2)
    rd_name <- paste(rd[1], "&", rd[2], "&", rd[3], sep = "")
    explore[nrow(explore) + 1, ] <- c(oc_remove[r], rd_name, nrow(capt.hist), n_recap, perc_recapt)
  }
  
  explore_remove_oc <- explore[-1,]
  explore_remove_oc$Nrecapt <- as.numeric(explore_remove_oc$Nrecapt)
  explore_remove_oc$perc_recapt <- as.numeric(explore_remove_oc$perc_recapt)
  
  sum_remove_oc <- explore_remove_oc %>%
    group_by(N_ocRem) %>%
    summarise(mean_recapt = mean(Nrecapt),
              mean_percent = mean(perc_recapt),
              min_samples = min(Ndetect),
              max_samples = max(Ndetect))

## ---- 2. Data summary subsampling ----
# I had to run the summaries at the same time than the models, as it was a random draw of the samples

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/optimize_sampling/results")
load("data_sum_subsample.RData")

explore_sub <- explore[-1,]
explore_sub$Nrecapt <- as.numeric(explore_sub$Nrecapt)
explore_sub$perc_recapt <- as.numeric(explore_sub$perc_recapt)

sum_sub <- explore_sub %>%
  group_by(per) %>%
  summarise(mean_recapt = mean(Nrecapt),
            mean_percent = mean(perc_recapt))

sum <- data.frame(group = c(sum_sub$per, sum_remove_oc$N_ocRem), 
           Nrecaptures = c(sum_sub$mean_recapt, sum_remove_oc$mean_recapt),
           id_recaptured = c(sum_sub$mean_percent, sum_remove_oc$mean_percent))

setwd("D:/MargSalas/Ganga/Results/Plots/Optimize")
write.csv(sum, file = "data_sumaries.csv")
