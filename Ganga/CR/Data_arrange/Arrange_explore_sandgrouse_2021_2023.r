## -------------------------------------------------
##            Arrange CR dataset sandgrouse
##                     2021-2023
## ------------------------------------------------- 

rm(list=ls())

library(tidyverse)
library(tidyr)
library(magrittr)

setwd("D:/MargSalas/Ganga/Data/CMR")

ganga <- readxl::read_excel("sandgrouse_Individual ID+sex-16Aug2023.xlsx") %>%
  select(7,8,16) %>%
  rename("id" = "INDIVIDUAL ID") %>%
  rename( "Occasion"= "Id_Visit")

ganga$Occasion[ganga$Year == 2023] <- ganga$Occasion[ganga$Year == 2023] + 1 # Change occasion id 2023

sum <- ganga %>% 
  group_by(Year, Occasion) %>%
  summarise(count = n_distinct(id))

# Capture history
# Select only indentified individuals and relevant columns for capture history
ganga <- ganga[which(!is.na(ganga$id)), ]

## ---- Check recaptures for survival estimation 2021-2023 ----

g <- ganga[,-2] # For this check I remove the occasion column (Interested only in primary occasions = years)
g$detect <- 1 # 122 detections

capt.hist <- g %>%
  # remove duplicates, which may occur when individuals are caught multiple times in an event
  # For example, your event may be a year and an individual may be caught multiple times in a year.
  distinct() %>%
  # spread out data. The fill = 0 adds rows for combinations of id and event where individuals were not observerd
  spread(Year, detect, fill = 0) %>% 
  # For every individual....
  group_by(id) %>%
  # Paste together 0's and 1's
  # Unite is similar to paste. Here we are pasting the strings together from the second column (first capture event)
  # to the last capture event ("tail(names(.),1)").
  # we don't want any characters separating 0's and 1's, so we use: sep = ""
  unite("ch", 2:tail(names(.),1), sep = "")

# Exploration recaptures

capt.hist.2 <- capt.hist %>%
  mutate(year1 = str_sub(ch, 1,1),
         year2 = str_sub(ch, 2,2),
         year3 = str_sub(ch, 3,3)) %>%
  mutate_at(c('year1', 'year2','year3'), as.numeric) %>%
  select(-ch)

length(which(apply(capt.hist.2[,c(2:4)], 1, sum) == 1))/nrow(capt.hist.2) # 70% of individuals haven't been recaptured
length(which(apply(capt.hist.2[,c(2:4)], 1, sum) == 3)) # Only 1 individual captured all times

# 1 0 1: Gives information to differenciate between survival and detection probability
which(capt.hist.2$year1 == 1 & capt.hist.2$year2 == 0 & capt.hist.2$year3 == 1) # None has been: 1 0 1

# Yearly information
length(which(capt.hist.2$year1 == 1)) # 38 captured in year 1 (2021)
length(which(capt.hist.2$year1 == 1 & capt.hist.2$year2 == 1)) # From these, only 16 recaptured in year 2

length(which(capt.hist.2$year2 == 1)) # 66 captured in year 2 (2022)
length(which(capt.hist.2$year1 == 0 & capt.hist.2$year2 == 1)) # 50 are new! weren't captured before
length(which(capt.hist.2$year2 == 1 & capt.hist.2$year3 == 1)) # From these, only 18 recaptured in year 3

length(which(capt.hist.2$year3 == 1)) # 50 captured in year 3 (2023)
length(which(capt.hist.2$year1 == 0 & capt.hist.2$year2 == 0 & capt.hist.2$year3 == 1)) # 32 are new! :O whaaat?


# COnclusiones: 
#   - En el año 2023 el esfuerzo de muestreo se relajó, pero aun asi se identificaron mas individuos que en 2021
#   - El % de recapturas es bajo, porque cada año se identifican nuevos individuos


