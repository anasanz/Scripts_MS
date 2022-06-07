## -------------------------------------------------
##                 SCR0 - Nimble 
##                  French 2017
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)
library(sf)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France")
load("data.RData")

## ---- Trap data ----

traps_2017 <- data[["traps"]][["2017"]] %>%
  filter(suivi == "systematic")

tdf2017 <- data.frame(trap = as.numeric(traps_2017$trap), # trap id number
                      X = st_coordinates(traps_2017)[,1], # longitude
                      Y = st_coordinates(traps_2017)[,2],# latitude
                      sep = rep("/",length(traps_2017$trap)), # separator
                      site = as.factor(traps_2017$site), # type of trap
                      pays = as.factor(traps_2017$pays), # country
                      suivi = as.factor(traps_2017$suivi)) # type of monitoring







##Compile data for Nimble in two objects: data (things that are random variables
## left side of a ~, can be recalculated within model) and constants

dat <- list(y = y.in)
consts <- list(M = M, J = J, X = X, 
               xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
               K = K)