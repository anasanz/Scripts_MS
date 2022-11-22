

library(tidyverse)

## -------------------------------------------------
##          Summary numer of available locations       
## ------------------------------------------------- 

# Cargar el excel de Elena (resumen nºobs/individuo/categoría de dispersión)

setwd("D:/MargSalas/Oso/Dispersal")
dat <- openxlsx::read.xlsx("TaulesResum_EPS.xlsx", sheet = 2)
colnames(dat)[4:13] <- dat[1,4:13]
dat <- dat[-1,]

# MALES
m <- dat[which(dat$Sex == "M"), ]
males <- unique(m$Confirmed_Individual_cub) 
length(males) # 52 males

m_natal <- m %>% filter(Age_class_cub2 == "Natal") 
m_subadult <- m %>% filter(Age_class_cub2 == "Subadult") 
m_adultrep <- m %>% filter(Age_class_cub2 == "Adult_reproductor") 
m_adult <- m %>% filter(Age_class_cub2 == "Adult") 

males[which(males %in% c(m_natal$Confirmed_Individual_cub) &
              males %in% (m_subadult$Confirmed_Individual_cub) &
              males %in% (m_adultrep$Confirmed_Individual_cub))]
# 7 males with data at all stages


# FEMALES
f <- dat[which(dat$Sex == "F"), ]
females <- unique(f$Confirmed_Individual_cub) 
length(females) # 59 females

f_natal <- f %>% filter(Age_class_cub2 == "Natal") 
f_subadult <- f %>% filter(Age_class_cub2 == "Subadult") 
f_adultrep <- f %>% filter(Age_class_cub2 == "Adult_reproductor") 
f_adult <- f %>% filter(Age_class_cub2 == "Adult") 

females[which(females %in% c(f_natal$Confirmed_Individual_cub) &
              females %in% (f_subadult$Confirmed_Individual_cub) &
              females %in% (f_adultrep$Confirmed_Individual_cub))]

# 12 males with data at all stages

