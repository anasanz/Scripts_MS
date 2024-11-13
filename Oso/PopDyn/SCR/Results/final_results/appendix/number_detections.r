## -------------------------------------------------
##                Summary SAMPLES 
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("edf1721.RData")
edf <- edf[-which(edf$ind %in% c("Nere", "Goiat")), ]

## ---- Nº of samples per age class ----

load("age_detectedInd.RData")
age.cat<-age
age.cat[age==0]<-0
age.cat[age==1]<-1
age.cat[age==2]<-1
age.cat[age==3]<-2
age.cat[age==4]<-2
age.cat[age>=5]<-3

# Add age information into detections
edf$age <- NA
for(i in 1:nrow(edf)){
  edf$age[i] <- age.cat[which(rownames(age.cat) %in% edf$ind[i]), edf$session[i]]
}

sum <- edf %>% group_by(session, age) %>%
  summarise(n())

sum2 <- as.data.frame(spread(sum, key = session, value = 'n()'))
colnames(sum2) <- c("agecat",colnames(age))

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Tables")
#write.csv(sum2, file = "number_detections.csv")
rowMeans(sum2[,2:6])
apply(sum2[,2:6], 1, mean)
apply(sum2[,2:6], 1, sd)

## ---- Nº of samples per year ----

smp <- edf %>% group_by(session) %>% summarise(n())
