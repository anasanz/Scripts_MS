## -------------------------------------------------
##      Arrange data DS Specific PTS 2022
## -------------------------------------------------

rm(list=ls())

setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
dat <- read.csv('Transectes_Específics_Ganga.csv', sep = ";", header = TRUE)

colnames(dat)[which(colnames(dat) == "ï..ID.Transecte")] <- "transectID"
colnames(dat)[which(colnames(dat) == "Hora.Inici")] <- "Start_time"
colnames(dat)[which(colnames(dat) == "Hora.Final")] <- "End_time"
colnames(dat)[which(colnames(dat) == "Observador")] <- "Observer"
colnames(dat)[which(colnames(dat) == "Vent")] <- "Wind"
colnames(dat)[which(colnames(dat) == "N")] <- "Count"
colnames(dat)[which(colnames(dat) == "DistÃ.ncia")] <- "Distance"

dat$Species <- "PTALC"
dat$Year <- 2022
dat$Effort <- 500

dat <- dat[,c(1,22,2,3,4,7,21,14,16,5,23)]

## ---- Estimate hour and julian day as Pepe did in National Surveys ----

dat$Start_time[14] <- "9:00"

# Hour of day from 0-1
dat$time <- data.table::as.ITime(dat$Start_time) # Convert to time
dat$time <- factor(dat$time)

a <- hms(as.character(dat$time)) # Extract hours and minuts
hour <- hour(a) # get hours
minute <- minute(a) # get minutes

convert_to_decimal <- function(hour, minute) {
  decimal_value <- (hour + minute / 60) / 24
  return(decimal_value)
} # Function to convert hours and minutes to a decimal value in the range 0-1

dat$time01 <- convert_to_decimal(hour = hour, minute = minute)

# Julian date
dat$date <- as.Date(dat$Data, format = "%d/%m/%y")
dat$jd <- as.POSIXlt(dat$date)$yday 

dat <- dat[,c(1:3,15,4:5,13,6:11)]

## ---- Add info of distance bands ----

# Checking the detections I think we should keep the same bins than in Farmdindis 
# If we take the ones from national program, all would fall in bin 1 
# But we do keep 400 as truncation distance

# Banda 1: 0-25
# Banda 2: 25-50
# Banda 3: 50-100
# Banda 4: 100-200
# Banda 5: 200- 400

dat$Bin <- NA # Medium point of each bin

for (i in 1:nrow(dat)){
  if (!is.na(dat$Distance[i])){
    if (dat$Distance[i] <= 25) {dat$Bin[i] = 1}
    else if (dat$Distance[i] <= 50) {dat$Bin[i] = 2}
    else if (dat$Distance[i] <= 100) {dat$Bin[i] = 3}
    else if (dat$Distance[i] <= 200) {dat$Bin[i] = 4}
    else if (dat$Distance[i] <= 400) {dat$Bin[i] = 4}
  }}

setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
write.csv(dat, file = "Data_HDS_Specific.csv")

# Detections

dat_det <- dat[which(dat$Count>0),]

## ---- Explore detection data ----

hist(dat$Distance[which(!is.na(dat$Distance))], breaks = c(0,24,49,99,200,500),
     main = "Detection function", col = "grey", freq = FALSE, xlab = "Distance") 

# WIND: I think it doesn't really matter

xtabs(~Wind, dat_det) # Transects with detections in each category

wind <- c(0:3)
par(mfrow = c(2,2)) # Frequency - Distance with different winds
for (i in 1:length(wind)){
  w <- dat_det[which(dat_det$Wind == wind[i]), ]
  hist(w$Distance, breaks = c(0,25,50,99,200), xlab = "Distance bins (x)", col = "grey", main = paste("Det - Wind ",wind[i]),
       freq = FALSE)
}

# CLOUDS: Is not even here but I checked in original data frame and all detection are in clouds (not in sunny days)
