##################################################################################
####    ARRANGE DATA FOR DS ANALYSIS OF PREDICTOR COLUMBIDAE      ################
##################################################################################

# Load packages


library(dplyr)


## -------------------------------------------------
##   COPAL (Columba palumbus/Paloma Torcaz/Common Wood pigeon)
## ------------------------------------------------- 

rm(list=ls())

setwd("D:/Otros/Tórtola/Data/Columbid")

tor <- read.csv("copal_ds_02_21.csv", sep = ",")
tor <- tor[ ,-1]
tor$count <- 1

setwd("D:/Otros/Tórtola/Data")
all_tr <- read.csv("FS_SOCCs_ampliats_02_21.csv", sep = ";")
colnames(all_tr)[1] <- "Site"
colnames(all_tr)[2] <- "Year"
all_tr <- all_tr[which(all_tr$Periode == 2), ]
all_tr$site_year <- paste(all_tr$Site, all_tr$Year, sep = "_") # All transects_year sampled
all_tr <- arrange(all_tr, Site, Year)

# ---- Modify bin ----
unique(tor$Bin)
for (i in 1:nrow(tor)){
  if (tor$Bin[i] == "banda1") {tor$Bin[i] = 1}
  else if (tor$Bin[i] == "banda2") {tor$Bin[i] = 2}
  else if (tor$Bin[i] == "banda3") {tor$Bin[i] = 3}
}

# ---- Get absences ---- #

allsites <- data.frame(site_year = all_tr$site_year) # Identify all transects IDs sampled

all_tr_sec <- do.call("rbind", replicate(6, all_tr, simplify = FALSE)) # Add sections
all_tr_sec <- arrange(all_tr_sec, site_year)
all_tr_sec$Section <- rep(c(1:6), times = nrow(all_tr))
all_tr_sec$site_sec <- paste(all_tr_sec$Site, all_tr_sec$Section, sep = "_") # All transects_year sampled

d <- left_join(all_tr_sec, tor, by = c("site_year", "site_sec")) # Which count is NA are the transects with 0 detections
absence <- d[which(is.na(d$count)), ]
f <- unique(tor$site_year)

# ---- Join climate data and absences to database ---- #

# Join climate
tor2 <- left_join(tor, all_tr_sec, by = c("site_sec", "site_year")) 

# Sort absence ant tor2 to do rbind
colnames(tor2)[1] <- "Site"
colnames(tor2)[2] <- "Year"
colnames(tor2)[3] <- "Observer"
colnames(tor2)[4] <- "Section"

#tor2 <- tor2[ ,-c(9:13, 20)]
tor2 <- tor2[ ,-c(10:14, 20)]


#absence <- absence[ ,-c(3,4, 15:18)]
absence <- absence[ ,-c(3,4, 14:17)]
colnames(absence)[1] <- "Site"
colnames(absence)[2] <- "Year"
colnames(absence)[3] <- "Observer"
colnames(absence)[10] <- "Section"
absence$count <- 0

colnames_order <- c("Site", "Year", "Section", "site_year", "site_sec", "Observer", "count", "Bin", "distance", 
                    "Temperatura", "Vent", "Cel", "Pluja", "Visibilitat")

tor2 <- tor2[,match(colnames_order, colnames(tor2))]
absence <- absence[,match(colnames_order, colnames(absence))]

dat <- rbind(tor2, absence)
dat <- arrange(dat, Site, Year, Section)

setwd("D:/Otros/Tórtola/Data/Columbid")
#write.csv(dat, "copal_ds_ready_02_21.csv")  

## -------------------------------------------------
##   STDEC (Streptopelia decaocto/Tórtola turca/Eurasian collared dove)
## ------------------------------------------------- 

rm(list=ls())


setwd("D:/Otros/Tórtola/Data/Columbid")

tor <- read.csv("stdec_ds_02_21.csv", sep = ",")
tor <- tor[ ,-1]
tor$count <- 1

setwd("D:/Otros/Tórtola/Data")
all_tr <- read.csv("FS_SOCCs_ampliats_02_21.csv", sep = ";")
colnames(all_tr)[1] <- "Site"
colnames(all_tr)[2] <- "Year"
all_tr <- all_tr[which(all_tr$Periode == 2), ]
all_tr$site_year <- paste(all_tr$Site, all_tr$Year, sep = "_") # All transects_year sampled
all_tr <- arrange(all_tr, Site, Year)

# ---- Modify bin ----
unique(tor$Bin)
for (i in 1:nrow(tor)){
  if (tor$Bin[i] == "banda1") {tor$Bin[i] = 1}
  else if (tor$Bin[i] == "banda2") {tor$Bin[i] = 2}
  else if (tor$Bin[i] == "banda3") {tor$Bin[i] = 3}
}

# ---- Get absences ---- #

allsites <- data.frame(site_year = all_tr$site_year) # Identify all transects IDs sampled

all_tr_sec <- do.call("rbind", replicate(6, all_tr, simplify = FALSE)) # Add sections
all_tr_sec <- arrange(all_tr_sec, site_year)
all_tr_sec$Section <- rep(c(1:6), times = nrow(all_tr))
all_tr_sec$site_sec <- paste(all_tr_sec$Site, all_tr_sec$Section, sep = "_") # All transects_year sampled

d <- left_join(all_tr_sec, tor, by = c("site_year", "site_sec")) # Which count is NA are the transects with 0 detections
absence <- d[which(is.na(d$count)), ]
f <- unique(tor$site_year)

# ---- Join climate data and absences to database ---- #

# Join climate
tor2 <- left_join(tor, all_tr_sec, by = c("site_sec", "site_year")) 

# Sort absence ant tor2 to do rbind
colnames(tor2)[1] <- "Site"
colnames(tor2)[2] <- "Year"
colnames(tor2)[3] <- "Observer"
colnames(tor2)[4] <- "Section"

#tor2 <- tor2[ ,-c(9:13, 20)]
tor2 <- tor2[ ,-c(10:14, 20)]


#absence <- absence[ ,-c(3,4, 15:18)]
absence <- absence[ ,-c(3,4, 14:17)]
colnames(absence)[1] <- "Site"
colnames(absence)[2] <- "Year"
colnames(absence)[3] <- "Observer"
colnames(absence)[10] <- "Section"
absence$count <- 0

colnames_order <- c("Site", "Year", "Section", "site_year", "site_sec", "Observer", "count", "Bin", "distance", 
                    "Temperatura", "Vent", "Cel", "Pluja", "Visibilitat")

tor2 <- tor2[,match(colnames_order, colnames(tor2))]
absence <- absence[,match(colnames_order, colnames(absence))]

dat <- rbind(tor2, absence)
dat <- arrange(dat, Site, Year, Section)

setwd("D:/Otros/Tórtola/Data/Columbid")
write.csv(dat, "stdec_ds_ready_02_21.csv")  

## -------------------------------------------------
##   COOEN (Columba Oenas/Paloma zurita/Stock dove)
## ------------------------------------------------- 

rm(list=ls())


setwd("D:/Otros/Tórtola/Data/Columbid")

tor <- read.csv("cooen_ds_02_21.csv", sep = ",")
tor <- tor[ ,-1]
tor$count <- 1

setwd("D:/Otros/Tórtola/Data")
all_tr <- read.csv("FS_SOCCs_ampliats_02_21.csv", sep = ";")
colnames(all_tr)[1] <- "Site"
colnames(all_tr)[2] <- "Year"
all_tr <- all_tr[which(all_tr$Periode == 2), ]
all_tr$site_year <- paste(all_tr$Site, all_tr$Year, sep = "_") # All transects_year sampled
all_tr <- arrange(all_tr, Site, Year)

# ---- Modify bin ----
unique(tor$Bin)
for (i in 1:nrow(tor)){
  if (tor$Bin[i] == "banda1") {tor$Bin[i] = 1}
  else if (tor$Bin[i] == "banda2") {tor$Bin[i] = 2}
  else if (tor$Bin[i] == "banda3") {tor$Bin[i] = 3}
}

# ---- Get absences ---- #

allsites <- data.frame(site_year = all_tr$site_year) # Identify all transects IDs sampled

all_tr_sec <- do.call("rbind", replicate(6, all_tr, simplify = FALSE)) # Add sections
all_tr_sec <- arrange(all_tr_sec, site_year)
all_tr_sec$Section <- rep(c(1:6), times = nrow(all_tr))
all_tr_sec$site_sec <- paste(all_tr_sec$Site, all_tr_sec$Section, sep = "_") # All transects_year sampled

d <- left_join(all_tr_sec, tor, by = c("site_year", "site_sec")) # Which count is NA are the transects with 0 detections
absence <- d[which(is.na(d$count)), ]
f <- unique(tor$site_year)

# ---- Join climate data and absences to database ---- #

# Join climate
tor2 <- left_join(tor, all_tr_sec, by = c("site_sec", "site_year")) 

# Sort absence ant tor2 to do rbind
colnames(tor2)[1] <- "Site"
colnames(tor2)[2] <- "Year"
colnames(tor2)[3] <- "Observer"
colnames(tor2)[4] <- "Section"

#tor2 <- tor2[ ,-c(9:13, 20)]
tor2 <- tor2[ ,-c(10:14, 20)]


#absence <- absence[ ,-c(3,4, 15:18)]
absence <- absence[ ,-c(3,4, 14:17)]
colnames(absence)[1] <- "Site"
colnames(absence)[2] <- "Year"
colnames(absence)[3] <- "Observer"
colnames(absence)[10] <- "Section"
absence$count <- 0

colnames_order <- c("Site", "Year", "Section", "site_year", "site_sec", "Observer", "count", "Bin", "distance", 
                    "Temperatura", "Vent", "Cel", "Pluja", "Visibilitat")

tor2 <- tor2[,match(colnames_order, colnames(tor2))]
absence <- absence[,match(colnames_order, colnames(absence))]

dat <- rbind(tor2, absence)
dat <- arrange(dat, Site, Year, Section)

setwd("D:/Otros/Tórtola/Data/Columbid")
write.csv(dat, "cooen_ds_ready_02_21.csv") 



