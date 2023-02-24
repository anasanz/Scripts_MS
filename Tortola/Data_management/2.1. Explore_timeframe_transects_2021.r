rm(list=ls())

library(rjags)
library(jagsUI)
library(dplyr)
library(stringr)


setwd("D:/Otros/Tórtola/Data")

tor <- read.csv("tortola_ds_ready_02_21.csv", sep = ",")
tor[,1] <- "STTUR"

###################################################################
##                       HDS ANALYSIS                           ###
###################################################################

# Model 1.1: Same as model 1 but with roughness covariate instead of forest

# ---- Information: bins, years, sites ----

strip.width <- 500 				
dist.breaks <- c(0,25,100,500)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-4]
nG <- length(dist.breaks)-1

yrs <- unique(tor$Year) 
nyrs <- length(yrs)

# 1 transect and 6 sections per year 
nSection <- 6
section <- c(1,2,3,4,5,6)

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

# Only to check: Count of individuals per year
count.year <- colSums(m,na.rm = TRUE)

## ---- SUBSET the data by transect ----

# 1. Keep transects with > 7 time points
m <- as.data.frame(m)
m_subset <- data.frame(matrix(nrow = 0, ncol = ncol(m)))
colnames(m_subset) <- colnames(m)

for (i in 1:nrow(m)) {
  if (sum(!is.na(m[i, ])) > 6) {
    m_subset[nrow(m_subset)+1,] <-  m[i, ] # Add row if it has 7 values or more 
  }}

# 2. Keep transects with at least 5 time points with counts (To be able to calculate trend)

m_subset$site <- str_sub(rownames(m_subset), end=-3) # Add column with site (removing 2 last digits section) to group per site
x <- t(sapply(split(m_subset, m_subset$site), function(x) colSums(x[, c(1:18)]))) # Counts grouped per site
x <- as.data.frame(ifelse(x>0,1,0))

m_subset2 <- data.frame(matrix(nrow = 0, ncol = ncol(m_subset)))
colnames(m_subset2) <- colnames(m_subset)

for (i in 1:nrow(m_subset)) {
  site <- m_subset[i,which(colnames(m_subset) %in% "site")] # Identify the site for the if condition
  if (sum(x[rownames(x) %in% site, ], na.rm = TRUE) > 4) {
    m_subset2[nrow(m_subset2)+1,] <-  m_subset[i, ] # Add row if it has 5 time points with counts or more
  }}

length(unique(m_subset2$site)) # Number of transects that will be analyzed

# 3. Make the subset and start loop

transect <- unique(m_subset2$site)

#setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
#save(transect, file = "1.1TortoData_transects_0221.RData")

# ---- Create data frame to store year of beginning and year of end ----

temp <- data.frame(matrix(NA, ncol = length(yrs) + 3, nrow = length(transect)))
colnames(temp) <- c("transect_ID", "year_start", "year_end", yrs)
temp$transect_ID <- transect

for (xxx in 1:length(transect)){
  
  data_transect <- m_subset2[which(m_subset2$site %in% transect[xxx]), ]
  
  # Year: Min and max year (if there are years with NA in between, those are considered within the time series)
  
  data_transect_nona <- data_transect[, colSums(is.na(data_transect)) != nrow(data_transect)] # Delete columns with NA
  
  min_year <- min(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) == "site")]))
  max_year <- max(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) == "site")]))
  
  temp[xxx,2:3] <- c(min_year,max_year)
  temp[xxx,4:23] <- ifelse(is.na(data_transect[1,-21]), 0, 1)
  
}
setwd("D:/Otros/Tórtola/Data")
#write.csv(temp, file = "timeframe_transects.csv")

# ---- Check if we can add more transects ----

## ---- SUBSET the data by transect ----

# 1. Keep transects with at least 5 time points
m <- as.data.frame(m)
m_subset <- data.frame(matrix(nrow = 0, ncol = ncol(m)))
colnames(m_subset) <- colnames(m)

for (i in 1:nrow(m)) {
  if (sum(!is.na(m[i, ])) > 4) {
    m_subset[nrow(m_subset)+1,] <-  m[i, ] # Add row if it has 7 values or more 
  }}

# 2. Keep transects with at least 4 time points with counts (To be able to calculate trend)

m_subset$site <- str_sub(rownames(m_subset), end=-3) # Add column with site (removing 2 last digits section) to group per site
x <- t(sapply(split(m_subset, m_subset$site), function(x) colSums(x[, c(1:18)]))) # Counts grouped per site
x <- as.data.frame(ifelse(x>0,1,0))

m_subset2 <- data.frame(matrix(nrow = 0, ncol = ncol(m_subset)))
colnames(m_subset2) <- colnames(m_subset)

for (i in 1:nrow(m_subset)) {
  site <- m_subset[i,which(colnames(m_subset) %in% "site")] # Identify the site for the if condition
  if (sum(x[rownames(x) %in% site, ], na.rm = TRUE) > 3) {
    m_subset2[nrow(m_subset2)+1,] <-  m_subset[i, ] # Add row if it has 5 time points with counts or more
  }}

length(unique(m_subset2$site)) # Number of transects that will be analyzed

# 3. Make the subset and start loop

transect2 <- unique(m_subset2$site)
transect_add <- transect2[-which(transect2 %in% transect)] # Transects to add, less restrictive

#setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
#save(transect_add, file = "1.1TortoData_transects_0221_ADD.RData") 

# Check time frame new transects

# ---- Create data frame to store year of beginning and year of end ----

temp <- data.frame(matrix(NA, ncol = length(yrs) + 3, nrow = length(transect_add)))
colnames(temp) <- c("transect_ID", "year_start", "year_end", yrs)
temp$transect_ID <- transect_add

for (xxx in 1:length(transect_add)){
  
  data_transect <- m_subset2[which(m_subset2$site %in% transect_add[xxx]), ]
  
  # Year: Min and max year (if there are years with NA in between, those are considered within the time series)
  
  data_transect_nona <- data_transect[, colSums(is.na(data_transect)) != nrow(data_transect)] # Delete columns with NA
  
  min_year <- min(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) == "site")]))
  max_year <- max(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) == "site")]))
  
  temp[xxx,2:3] <- c(min_year,max_year)
  temp[xxx,4:23] <- ifelse(is.na(data_transect[1,-21]), 0, 1)
}

#setwd("D:/Otros/Tórtola/Data")
#write.csv(temp, file = "timeframe_transects_add.csv")
