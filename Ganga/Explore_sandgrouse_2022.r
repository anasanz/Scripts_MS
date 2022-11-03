
## -------------------------------------------------
##            Arrange CR dataset sandgrouse
##                      2022
## ------------------------------------------------- 

rm(list=ls())

library(tidyverse)
library(tidyr)
library(magrittr)

setwd("D:/MargSalas/Ganga/Data")

ganga <- readxl::read_excel("sandgrouse_Individual ID+sex_29092022.xlsx") %>%
  select(7,8,16) %>%
  rename("id" = "INDIVIDUAL ID") %>%
  rename( "Occasion"= "Id_Visit")

sum <- ganga %>% 
  group_by(Year, Occasion) %>%
  summarise(count = n_distinct(id))

# Capture history
# Select only indentified individuals and relevant columns for capture history
ganga <- ganga[which(!is.na(ganga$id)), ] 

## ---- Check capture histories for abundance 2022 ----

# Select year 2022
g <- ganga[which(ganga$Year == 2022), ]
g <- g[,-1]
g$detect <- 1


capt.hist <- g %>%
  # remove duplicates, which may occur when individuals are caught multiple times in an event
  # For example, your event may be a year and an individual may be caught multiple times in a year.
  distinct() %>%
  # spread out data. The fill = 0 adds rows for combinations of id and event where individuals were not observerd
  spread(Occasion, detect, fill = 0) %>% 
  # For every individual....
  group_by(id) %>%
  # Paste together 0's and 1's
  # Unite is similar to paste. Here we are pasting the strings together from the second column (first capture event)
  # to the last capture event ("tail(names(.),1)").
  # we don't want any characters separating 0's and 1's, so we use: sep = ""
  unite("ch", 2:tail(names(.),1), sep = "")

capt.hist$ch <- as.character(capt.hist$ch)

#☺ Check how many are caotured twice

data_ganga <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = length(unique(g$Occasion)))
T = unique(g$Occasion)

capt.hist[which(apply(data_ganga,1,sum) > 1),]

save(capt.hist, file = "cr_sandgrouse_2022.RData")

## ---- Check capture histories for survival 2021-2023 ----

ganga 

