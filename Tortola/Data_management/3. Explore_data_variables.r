
rm(list=ls())

library(dplyr)

# Load data
setwd("D:/Otros/Tórtola/Data")

load("D:/Otros/Tórtola/Data/data_buff500.rdata")
colnames(data_buff500)[1] <- "Site"
colnames(data_buff500)[2] <- "Section"

# Create site_section
data_buff500$site_sec <- paste(data_buff500$Site, data_buff500$Section, sep = "_")

# Select interesting variables
colnames(data_buff500)
var <- data_buff500[ ,c("Hm_mean", "Hm_max", "FCC_mean", "FCC100%", "Forest%", "richness", "ForestMargin", "ForestMargin%", "site_sec")]

# ROUGHNESS CO-VARIATE
setwd("D:/Otros/Tórtola/Data")
rough <- read.csv("TRI_buff_500.csv", sep = ",")
rough$site_sec <- paste(rough[,1], rough[,2], sep = "_")
rough <- rough[,c(16,15)]

var <- left_join(var, rough, by = "site_sec")

tor <- read.csv("tortola_ds_ready.csv", sep = ",")

# Número de detecctiones por transecto-sección ~ variables

tor_det <- tor %>% 
  group_by(site_sec) %>%
  summarise(detections = n())

tor_det <- tor %>% 
  group_by(site_sec, site_year) %>%
  summarise(detections = n())

hist(tor_det$detections, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))

dat <- left_join(tor_det, var, by = "site_sec")
dat <- as.data.frame(dat)
variables <- colnames(var)[-9]
class(dat)

par(mfrow = c(2,4))
par(mfrow = c(1,1))
for(i in 1:length(variables)){
  plot(dat$detections ~ dat[,colnames(dat) %in% variables[i]], main = variables[i], pch = 18)
}

par(mfrow = c(2,4))
par(mfrow = c(1,1))
for(i in 1:length(variables)){
  hist(dat[,colnames(dat) %in% variables[i]], main = variables[i])
}


# TEMPERATURE CO-VARIATE
setwd("D:/Otros/Tórtola/Data")

tor <- read.csv("tortola_ds_ready.csv", sep = ",")
tor[,1] <- "STTUR"

hist(tor$Temperatura)

# Proportion of observations within different temperature category
nrow(tor[which(tor$Temperatura == 1), ])/nrow(tor)*100 # <0
nrow(tor[which(tor$Temperatura == 2), ])/nrow(tor)*100 # 0-10
nrow(tor[which(tor$Temperatura == 3), ])/nrow(tor)*100 # 10-20
nrow(tor[which(tor$Temperatura == 4), ])/nrow(tor)*100 # 20-30
nrow(tor[which(tor$Temperatura == 5), ])/nrow(tor)*100 # 20-30

library(dplyr)
t2 <- tor %>% 
  group_by(site_sec, Year) %>%
  mutate(count_sitesec = sum(count)) # Sum the number of detections per site_sec and year
  
t3 <- t2[-which(duplicated(t2)), ]  
plot(t3$count_sitesec ~ t3$Temperatura) # Relation between number of observations and temperature


