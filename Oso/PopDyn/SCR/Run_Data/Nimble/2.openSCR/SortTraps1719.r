
## -------------------------------------------------
##         Sort out trap data for open model
## ------------------------------------------------- 
## Create a single file with the combination of traps from all years

rm(list = ls())

library(tidyverse)

setwd("D:/MargSalas/Oso/Datos/Effort/France") # Only systematic

load("tdf2017_fr.RData")
load("tdf2018_fr.RData")
load("tdf2019_fr.RData")

# Check duplicates

check <- list()
tdf <- list(tdf2017, tdf2018, tdf2019)

for (s in 1:length(tdf)) {
  
  dup_year <- tdf[[s]][which(duplicated(tdf[[s]][,c(2:3)])), ]
  check_year <- matrix(NA, nrow = 1, ncol = ncol(tdf[[s]]))
  colnames(check_year) <- colnames(tdf[[s]])
  
  for (i in 1:nrow(dup_year)){
    tdf_year <- tdf[[s]]
    newdf <- tdf[[s]][which(tdf[[s]]$X %in% dup_year$X[i] & tdf[[s]]$Y %in% dup_year$Y[i]), ]
    check_year <- rbind(check_year,newdf)
  }
  check[[s]] <- check_year
}

# There are duplicates, and there are 2 traps that have the same coordinates but are of 
# different type  ("hair" and "both"). Same 2 traps every year is it an error? Keep 
# the one of type "both?

# I remove the one that is "hair" and I keep it as "both"
for (s in 1:length(tdf)) {
tdf[[s]] <- tdf[[s]][which(-tdf[[s]]$X == 317807 & tdf[[s]]$site == "hair" |
                  tdf[[s]]$X == 320469 & tdf[[s]]$site == "hair"), ]
}

# I remove the duplicates from within the session (all to the same type of trap)

tdf2017 <- tdf[[1]][,c(1:3)]
tdf2018 <- tdf[[2]][,c(1:3)]
tdf2019 <- tdf[[3]][,c(1:3)]

tdf2017 <- tdf2017[-which(duplicated(tdf2017[,c(2:3)])), ]
tdf2018 <- tdf2018[-which(duplicated(tdf2018[,c(2:3)])), ]
tdf2019 <- tdf2019[-which(duplicated(tdf2019[,c(2:3)])), ]

tdf2017$session <- 1
tdf2018$session <- 2
tdf2019$session <- 3

## ---- 2017 - 2018 ----

tdf <- rbind(tdf2017, tdf2018)

# The duplicated are from 2018
dup <- tdf[duplicated(tdf[,c(2:3)]), ]  
colnames(dup)[4] <- "Session2018"
colnames(dup)[1] <- "TrapID2018"

# Remove duplicated to join them beside
tdf2 <- tdf[-which(duplicated(tdf[,c(2:3)])), ] 
  # Without duplicates, is the unique traps from 2017 and 2018
tdf3 <- left_join(tdf2, dup)

# Sort out
for (i in 1:nrow(tdf3)){
  if(tdf3$session[i] == 2) {
   tdf3$TrapID2018[i] <- tdf3$trap[i]
   tdf3$Session2018[i] <- 2
   tdf3$trap[i] <- NA
   tdf3$session[i] <- NA
  }}
   
colnames(tdf3)[1] <- "TrapID2017"
colnames(tdf3)[4] <- "Session2017"

tdf3$TrapID2019 <- NA
tdf3$Session2019 <- NA

## ---- 2017+2018 -> 2019 ----

# Sort 2019 to join
colnames(tdf2019)[1] <- "TrapID2019"
colnames(tdf2019)[4] <- "Session2019"
tdf2019$TrapID2017 <- NA
tdf2019$TrapID2018 <- NA
tdf2019$Session2017 <- NA
tdf2019$Session2018 <- NA
tdf2019 <- tdf2019[,c(5,2,3,7,6,8,1,4)]

tdf <- rbind(tdf3, tdf2019)

# The duplicated are from 2019
dup <- tdf[duplicated(tdf[,c(2:3)]), ] 

# Remove duplicated to join them beside
tdf2 <- tdf[-which(duplicated(tdf[,c(2:3)])), ] 
# Without duplicates, is the unique traps from 2017 and 2018
tdf3 <- left_join(tdf2, dup, by = c("X","Y"))

# Sort out
tdf3 <- tdf3[,-c(9:12)]

for (i in 1:nrow(tdf3)){
  if(is.na(tdf3$TrapID2019.x[i])){
    tdf3$TrapID2019.x[i] <- tdf3$TrapID2019.y[i]
    tdf3$Session2019.x[i] <- tdf3$Session2019.y[i]
  }}

tdf3 <- tdf3[,-c(9,10)]

# Order columns
tdf3 <- tdf3[,c(2,3,1,4:8)]
colnames(tdf3)[3:8] <- str_sub(colnames(tdf3)[3:8],1,nchar(colnames(tdf3)[3:8])-2)

# Create new ID column
tdf3$newID <- seq(1:nrow(tdf3))
tdf3 <- tdf3[,c(9,1:8)]


## -------------------------------------------------
##            Assign encounter data to newID
## ------------------------------------------------- 

load("edf2017_2019_fr.RData")

# In this dataset there can be several detections per occasion,
# Because of the effort covariate (some traps were visited twice per occasion)
# And also because some traps contain both a hair snare and a camera trap
# This may not be a problem in models that use a poisson to model the detection,
# But we have a binomial model in the openPop, so I will consider only one detection
# per occasion. I do this by removing the duplicates

edf <- edf[-which(duplicated(edf)), ]
session <- c(1,2,3)
session_name <- c("Session2017", "Session2018", "Session2019")
session_trap <- c("TrapID2017", "TrapID2018", "TrapID2019")

new_edf <- matrix(NA,nrow = 1,ncol = 8)
colnames(new_edf) <- c("session","ind","occ","trap","sex","newID","X","Y")

for(t in 1:length(session)){
  edft <- edf[which(edf$session %in% session[t]),] # Select encounter data session
  tdft <- tdf3[!is.na(tdf3[,colnames(tdf3) %in% session_name[t]]), ] # Select trap data session
  
  tdft <- tdft[,colnames(tdft) %in% c("newID", "X", "Y", session_name[t], session_trap[t])] # Sort out to join
  colnames(tdft)[4] <- "trap"
  colnames(tdft)[5] <- "session"
  
  edft <- left_join(edft,tdft, by = c("trap", "session"))
  new_edf <- rbind(new_edf,edft)
}
new_edf <- new_edf[-1,]

## -------------------------------------------------
##         Format both files ready to run
## ------------------------------------------------- 

tdf <- tdf3[,c(1:3)]
#setwd("D:/MargSalas/Oso/Datos/Effort/France/forOpenPop") 
#save(tdf, file ="tdf_all.RData")

new_edf <- new_edf[,-4] # Remove last trap column
new_edf <- new_edf[,c(1,2,3,5,4)]
colnames(new_edf)[4] <- "trap"

setwd("D:/MargSalas/Oso/Datos/Effort/France/forOpenPop") 
#save(new_edf, file ="edf_all.RData")


## -------------------------------------------------
##         Format active traps per year
## ------------------------------------------------- 

## One matrix with 0 and 1 to get active traps per year

tr <- matrix(NA, nrow = length(tdf3$newID), ncol = 3)
year <- c(2017,2018,2019)

for(t in 1:length(year)){
  tdf3_year <- tdf3[,grepl(year[t], colnames(tdf3))]
  tr[,t] <- ifelse(!is.na(tdf3_year[,1]),1,0)
}

tr <- cbind(tdf3$newID,tr)
colnames(tr) <- c("newID",year)

# Check that is right

setwd("D:/MargSalas/Oso/Datos/Effort/France") # Only systematic

load("tdf2017_fr.RData")
load("tdf2018_fr.RData")
load("tdf2019_fr.RData")

tdf2017$X2 <- round(tdf2017$X, digits = 1) # Esto es porque el fucking r no me deja buscar coordenadas si no
tdf2018$X2 <- round(tdf2018$X, digits = 1)
tdf2019$X2 <- round(tdf2019$X, digits = 1)
          
# Examp: La trampa 15 de 2019 (691 newID) solo deberia estar activa este año, por lo que no existe en los otros

tdf2017[which(tdf2017$X2 == 294927.2), ]
tdf2018[which(tdf2018$X2 == 294927.2), ]
tdf2019[which(tdf2019$X2 == 294927.2), ]

