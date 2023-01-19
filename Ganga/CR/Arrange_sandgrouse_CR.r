
## -------------------------------------------------
##            Arrange CR dataset sandgrouse
## ------------------------------------------------- 

rm(list=ls())

library(dplyr)
library(tidyr)
library(magrittr)

setwd("D:/MargSalas/Ganga/Data")

g <- read.csv("sandgrouse_rawdata.csv", sep = ";")

# Name sampling occassions
g$Occasion <- g$Date
g$Occasion <- as.factor(g$Occasion)
levels(g$Occasion)[1] <- 2
levels(g$Occasion)[2] <- 3
levels(g$Occasion)[3] <- 1
g$Occasion <- as.numeric(as.character(g$Occasion))
# Select only indentified individuals and relevant columns for capture history
g <- g[which(g$INDIVIDUAL.ID != ""), ] 

# Save dataset coordinates
coord_g <- g[,c(4,5,8,13,16)]

library(stringr)
coord_g$X.Coord <- str_sub(coord_g$X.Coord, 5)
#write.csv(coord_g, "coord_g.csv")


# Save detection history
g <- g[,c(13,16)]
g$detect <- 1


capt.hist <- g %>%
  # remove duplicates, which may occur when individuals are caught multiple times in an event
  # For example, your event may be a year and an individual may be caught multiple times in a year.
  distinct() %>%
  # spread out data. The fill = 0 adds rows for combinations of id and event where individuals were not observerd
  spread(Occasion, detect, fill = 0) %>% 
  # For every individual....
  group_by(INDIVIDUAL.ID) %>%
  # Paste together 0's and 1's
  # Unite is similar to paste. Here we are pasting the strings together from the second column (first capture event)
  # to the last capture event ("tail(names(.),1)").
  # we don't want any characters separating 0's and 1's, so we use: sep = ""
  unite("ch", 2:tail(names(.),1), sep = "")

capt.hist$ch <- as.character(capt.hist$ch)

#save(capt.hist, file = "cr_sandgrouse.RData")

#☺ Check how many are caotured twice

data_ganga <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = length(unique(g$Occasion)))
T = unique(g$Occasion)

capt.hist[which(apply(data_ganga,1,sum) > 1),]
