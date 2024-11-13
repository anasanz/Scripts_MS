
## -------------------------------------------------
##        Run open CMR on bear data 2010 - 2021
##    Derive survival per sex (Females and Males)
## ------------------------------------------------- 

rm(list=ls())

library(RMark)
library(dplyr)
library(tidyr)
library(magrittr)


# --- Load data and info to join sex ----
setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Data_os_96_20.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_cr <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified
os_cr$Method[which(os_cr$Method == "sampling_station")] <- "Sampling_station"
os_cr$Obs_type[which(os_cr$Obs_type == "hair")] <- "Hair"

# Only data from sampling stations and hair
os_cr <- os_cr[which(os_cr$Method %in% "Sampling_station" & os_cr$Obs_type %in% "Hair"), ]

# Load info database
setwd("D:/MargSalas/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";")
info <- info[,c(4,5)]


#os_cr <- os_cr_all[which(os_cr_all$Region == "Catalunya"), ]

# ---- Format CR MARK ----

os_cr$Occasion <- os_cr$Year

os_cr <- os_cr[ ,c(1,23)]
os_cr$detect <- 1


capt.hist <- os_cr %>%
  # remove duplicates, which may occur when individuals are caught multiple times in an event
  # For example, your event may be a year and an individual may be caught multiple times in a year.
  distinct() %>%
  # spread out data. The fill = 0 adds rows for combinations of id and event where individuals were not observerd
  spread(Occasion, detect, fill = 0) %>% 
  # For every individual....
  group_by(Confirmed_Individual) %>%
  # Paste together 0's and 1's
  # Unite is similar to paste. Here we are pasting the strings together from the second column (first capture event)
  # to the last capture event ("tail(names(.),1)").
  # we don't want any characters separating 0's and 1's, so we use: sep = ""
  unite("ch", 2:tail(names(.),1), sep = "")

colnames(capt.hist)[1] <- "ID"
capt.hist.sex <- left_join(capt.hist, info, by = "ID")
capt.hist.sex <- capt.hist.sex[-46,-1]
capt.hist.sex$Sex <- as.factor(capt.hist.sex$Sex)

#####################
##tell RMark that we want to use "sex" as a grouping variable (stratification)
dipproc = process.data(capt.hist.sex, groups="Sex")
##going forward, we will use the object "dipproc" in the models

## set up model structure for both parameters survival (phi) and detection (p)
#intercept only
null <- list(formula = ~1)

#separate parameters for males and females
#note that the word after ~ needs to match the name you used to build groups in 
#the process.data() function
sex <- list(formula = ~Sex)

#changing parameters over time (ie, capture occasions)
time <- list(formula = ~time)

## Run possible models with p and phi varying by sex and time
phi0.p0 <- mark(dipproc, model.parameters = 
                  list(Phi = null, p = null), brief = T)
phi.sex.p0 <- mark(dipproc, model.parameters = 
                     list(Phi = sex, p = null), brief = T)
phi.time.p0 <- mark(dipproc, model.parameters = 
                      list(Phi = time, p = null), brief = T)
phi0.p.sex <- mark(dipproc, model.parameters = 
                     list(Phi = null, p = sex), brief = T)
phi.sex.p.sex <- mark(dipproc, model.parameters = 
                        list(Phi = sex, p = sex), brief = T)
phi.time.p.sex <- mark(dipproc, model.parameters = 
                         list(Phi = time, p = sex), brief = T)
phi0.p.time <- mark(dipproc, model.parameters = 
                      list(Phi = null, p = time), brief = T)
phi.sex.p.time <- mark(dipproc, model.parameters = 
                         list(Phi = sex, p = time), brief = T)

#################
##build model selection table and look at results of top models
modlist<-collect.models() 
modlist

## Top model results on the link and real scale
phi0.p.sex$results$beta
phi0.p.sex$results$real

phi.sex.p.sex$results$beta
phi.sex.p.sex$results$real

setwd("D:/MargSalas/Oso/SCR/Exp_analysis")
write.csv(phi0.p.sex$results$real, file = "Mark_phi0psex.csv")
write.csv(phi.sex.p.sex$results$real, file = "Mark_phiSexpSex.csv")
