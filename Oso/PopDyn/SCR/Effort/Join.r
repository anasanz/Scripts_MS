
## -------------------------------------------------
##                 Join data
## ------------------------------------------------- 

rm(list = ls())


library(tidyverse)
library(sf)

##○ !!!!!
# One thing to take into account are the observations that come from both hair and camera.
# Repeated??? How do we handle them? Is it important??

# Other thing. What about the independent observation of cubs added in Spain?
# What do I do with the French data? :(:(



## ---- Load data ----

##  2017-2019: All (Spanish and France) from Maelis

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data")
load("data.RData")

# Encounters 
det_2017 <- data[["det"]][["2017"]] %>%
  filter(suivi == "systematic")

det_2018 <- data[["det"]][["2018"]] %>%
  filter(suivi == "systematic")

det_2019 <- data[["det"]][["2019"]] %>%
  filter(suivi == "systematic")

# Traps positions
traps_2017 <- data[["traps"]][["2017"]] %>%
  filter(suivi == "systematic")

traps_2018 <- data[["traps"]][["2018"]] %>%
  filter(suivi == "systematic")

traps_2019 <- data[["traps"]][["2019"]] %>%
  filter(suivi == "systematic")

## 2020
# Spain:
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/Spain/Data")
load("dataSpain20.RData")

det_2020_SP <- dataSpain20[["det"]]
traps_2020_SP <- dataSpain20[["traps"]]

load("dataSpain21.RData")

det_2021_SP <- dataSpain21[["det"]]
traps_2021_SP <- dataSpain21[["traps"]]

## ---- Encounter data frame ----

edf <- data.frame(session = c(rep(1, nrow(det_2017)),
                              rep(2, nrow(det_2018)),
                              rep(3, nrow(det_2019)),
                              rep(4, nrow(det_2020_SP))), # add the session column
                  ind = c(det_2017$id,
                          det_2018$id,
                          det_2019$id,
                          det_2020_SP$id), # row names are individual IDs
                  occ = c(det_2017$month,
                          det_2018$month,
                          det_2019$month,
                          det_2020_SP$month) %>% as.factor() %>% as.numeric, # occasions are month may = 1
                  trap = c(det_2017$trap,
                           det_2018$trap,
                           det_2019$trap,
                           det_2020_SP$trap), # col names are trap IDs
                  sex = fct_c(det_2017$sex, # note the use of fct_c !!
                              det_2018$sex,
                              det_2019$sex,
                              as.factor(det_2020_SP$sex))) 

# Can join it together, although maybe one file per year at the end

# Unify names, put the spanish version
unique(edf$ind[edf$session %in% c(1:3)])
unique(edf$ind[edf$session %in% c(4)])

edf$ind <- iconv(edf$ind,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents

edf$ind[which(edf$ind == "Pepite")] <- "Pepito"
edf$ind[which(edf$ind == "Cannellito")] <- "Canellito"
edf$ind[which(edf$ind == "Callisto")] <- "Callista"
edf$ind[which(edf$ind == "New18_17")] <- "Negre"
edf$ind[which(edf$ind == "New18_04")] <- "Plumita"
edf$ind[which(edf$ind == "New18_10")] <- "New 18-10"
edf$ind[which(edf$ind == "New18_14")] <- "Nheueto"
edf$ind[which(edf$ind == "Gribouille")] <- "Griboulle"
edf$ind[which(edf$ind == "New18_16")] <- "New 18-16"
edf$ind[which(edf$ind == "New18_06")] <- "Salada"
edf$ind[which(edf$ind == "New18_18")] <- "Deserta"
edf$ind[which(edf$ind == "New19_01")] <- "Cirera"

unique(edf$ind)

