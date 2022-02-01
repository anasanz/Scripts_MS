

# Load packages

rm(list=ls())

library(dplyr)


setwd("D:/PhD/Otros/Tórtola/Data")

tor <- read.csv("tortola_ds_ready.csv", sep = ",")
tor[,1] <- "STTUR"
###################################################################
##                       HDS ANALYSIS                           ###
###################################################################

# ---- Information: bins, years, sites ----

strip.width <- 500 				
dist.breaks <- c(0,25,100,500)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-4]
nG <- length(dist.breaks)-1

unique(tor$Year)
yrs <- unique(tor$Year) 
nyrs <- length(yrs)

# ---- Distance observations ----

# Format
all.sites <- unique(tor$site_sec)
max.sites <- length(all.sites)

m <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(m) <- all.sites
colnames(m) <- yrs

# Add counts > 0
tor_yes <- tor[which(tor$count == 1), ]
count <- aggregate(X ~ Year + site_sec, FUN = length, data = tor_yes)

for (i in 1:nrow(count)){
  m[which(rownames(m) %in% count$site_sec[i]), which(colnames(m) %in% count$Year[i])] <- count$X[i]
}

# Add absences (0)
tor_no <- tor[which(tor$count == 0), ]
for (i in 1:nrow(tor_no)){
  m[which(rownames(m) %in% tor_no$site_sec[i]), which(colnames(m) %in% tor_no$Year[i])] <- tor_no$count[i]
}

# ---- Check how to do subsets for analysis ----
# Do 2 subsets: 
m_02_10 <- m[,c(1:9)]
m_10_19 <- m[,c(9:18)]

# How many transects are there with only less than 3 years
ncol(m)
m2 <- m[-which(rowSums(is.na(m)) > 15), ]
sum(rowSums(m[which(rowSums(is.na(m)) > 15), ], na.rm = TRUE)) # Se pierden 195 detecciones de 2868
sum(rowSums(m, na.rm = TRUE))

ncol(m_02_10)
m2_02_10 <- m_02_10[-which(rowSums(is.na(m_02_10)) > 6), ]
sum(rowSums(m_02_10[which(rowSums(is.na(m_02_10)) > 6), ], na.rm = TRUE)) # Se pierden 105 detecciones de 1004
sum(rowSums(m_02_10, na.rm = TRUE))

ncol(m_10_19)
m2_10_19 <- m_10_19[-which(rowSums(is.na(m_10_19)) > 7), ]
sum(rowSums(m_10_19[which(rowSums(is.na(m_10_19)) > 7), ], na.rm = TRUE)) # Se pierden 219 detecciones de 2075
sum(rowSums(m_10_19, na.rm = TRUE))



