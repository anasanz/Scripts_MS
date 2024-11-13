
## -------------------------------------------------
##            Arrange CR dataset sandgrouse
##                      2022
## ------------------------------------------------- 

rm(list=ls())

library(tidyverse)
library(tidyr)
library(magrittr)

setwd("D:/MargSalas/Ganga/Data/CMR")

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

data_ganga_old <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = length(unique(g$Occasion)))
T = unique(g$Occasion)

capt.hist[which(apply(data_ganga,1,sum) > 1),]

#save(capt.hist, file = "cr_sandgrouse_2022.RData")
capt.hist.old <- capt.hist
recaptures_old <- as.data.frame(capt.hist.old[which(apply(data_ganga_old,1,sum) > 1),])


## -------------------------------------------------
##            Arrange CR dataset sandgrouse
##                 2022 ANALYSIS NEW SAMPLES
## ------------------------------------------------- 

rm(list=ls())

library(tidyverse)
library(tidyr)
library(magrittr)

setwd("D:/MargSalas/Ganga/Data/CMR")

ganga <- readxl::read_excel("sandgrouse_Individual ID+sex_FINAL.xlsx") %>%
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
g$detect <- 1 # 122 detections


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

recaptures_new <- as.data.frame(capt.hist[which(apply(data_ganga,1,sum) > 1),])

#save(capt.hist, file = "cr_sandgrouse_2022_new.RData")

## ---- Create table with sex of each individual ----

ganga_sex <- readxl::read_excel("sandgrouse_Individual ID+sex_FINAL.xlsx") %>%
  select(7,8,15,16) %>%
  rename("id" = "INDIVIDUAL ID") %>%
  rename( "Occasion"= "Id_Visit")

ganga_sex <- ganga_sex[which(!is.na(ganga_sex$id)), ] 
ganga_sex <- ganga_sex[ganga_sex$Year == 2022, ]
ganga_sex$id <- as.factor(ganga_sex$id)
ganga_sex <- ganga_sex %>% group_by(id) %>% arrange(id)

# Check if there are inconsistencies (i.e., one individual assigned to 2 sexes)
# And store sex in table with same order as capture histories
id_sex <- data.frame(id = capt.hist$id, sex = NA)

id <- unique(ganga_sex$id)

for(i in 1:nrow(ganga_sex)){
  sex_id <- ganga_sex[which(ganga_sex$id %in% id[i]),]
  if ((length(unique(sex_id$sex))) == 1) {
    id_sex$sex[id_sex$id %in% id[i]] <- sex_id$sex[1]
  } else {
    id_sex$sex[id_sex$id %in% id[i]] <- "inconsistent"
  }
}

# Fix inconsistencies (few, so manually)
id_check <- id_sex$id[which(id_sex$sex == "inconsistent")]
ganga_sex[which(ganga_sex$id %in% id_check), ] # All females
id_sex$sex[id_sex$id %in% c("SG19", "SG45", "SG57","SG69")] <- "F"

# Check if there are individuals with unknown sex
unique(id_sex$sex)
id_sex$id[which(id_sex$sex %in% "X")] # SG71 unknown sex :(

capt.hist[capt.hist$id %in% "SG71",] # Captured once, maybe I can delete it?

nrow(id_sex[which(id_sex$sex == "M"), ])
nrow(id_sex[which(id_sex$sex == "F"), ])


save(id_sex, file = "id_sex_sandgrouse_2022.RData")

