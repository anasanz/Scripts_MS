
## -------------------------------------------------
##      Create independent cub observations
## ------------------------------------------------- 

## Survival of cubs is very high in CJS and JS models, so we wonder if it is 
# because we only detect cubs at a late age, once they have survived the critic 
# period.
## In this script I try to put the observations of the mother with cubs as independent observations of cubs

rm(list = ls())

library(lubridate)
library(stringr)

# Load monitoring data

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_2.xlsx')


# Load table info individuals to fill out sex and ages of cubs
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)
info <- info[,c(4,5,8)]
colnames(info)[1] <- "Probable_Individual"
info <- info[-which(is.na(info)), ]

info[nrow(info)+1, ] <- list("Moonbots/Bonabe", "M", "2006") # This obs will have to be deleted, but to avoid errors now


## ---- Sure cub detection ----

# When the number of cubs observed is the same as the number of cubs estimated for that female that year

os$cub_with_mother <- 0
 
for (i in 1:nrow(os)){
   if (os$Probable_Individual[i] == "Indetermined") next
   if (os$cub_with_mother[i] == 1) next # To not repeat loop on cubs already added
   if (os$N_cubs_observed[i] == os$N_cubs_estimated[i] & !is.na(os$N_cubs_estimated[i]) & !is.na(os$N_cubs_observed[i])) {
     add_row <- os[i,]     # Row to add with the different names of the cubs
     cubs <- t(add_row[1,colnames(add_row) %in% c("ID_cub1", "ID_cub2", "ID_cub3")]) # Cubs to add
     
     for (j in (1:as.numeric(os$N_cubs_observed[i]))) { # Add as many rows as cubs
       
       add_row <- os[i,] # To set a blank row for every cub
       
       # Fill information of the cub in the row
       add_row[,c(2,3)] <- cubs[j]
       
       add_row$Sex <- info$Sex[info$Probable_Individual == cubs[j]]
       add_row$Age <- add_row$Year - as.numeric(info$Year_birth[info$Probable_Individual == cubs[j]]) 
       
       add_row$Date_register_formated <- as.Date(add_row$Date_register, format = "%d/%m/%Y") # For Age class
       add_row$month <- month(ymd(add_row$Date_register_formated)) 
       
       if(is.na(add_row$month)) { # Si el mes es NA, lo pongo como no crítico, que serán pocos casos
         add_row$month <- 8 
       }
         
       if (add_row$Age <= 4 & add_row$Age > 1.5) {
         add_row$Age_class <- "Subadult" } else if (add_row$Age == 1.5 | add_row$Age == 1) {
           add_row$Age_class <- "Cub2" } else if (add_row$Age == 0 & add_row$month > 7)  {
         add_row$Age_class <- "Cub1" } else if (add_row$Age == 0 & add_row$month <= 7) {
           add_row$Age_class <- "Cub0" }
       
       add_row[,c(20:26)] <- NA # Because it is the cub
       
       add_row <- add_row[,-c(36,37)]
       add_row$ID_obs <- paste(add_row$ID_obs,j, sep = ".")
       add_row$cub_with_mother <- 1 
       
       # Add row
       long <- nrow(os)
       os <- rbind(os[1:i,], add_row ,os[(i + 1):(long),])
     }
   }
 }

# Check
check <- os[which(os$cub_with_mother == 1), ]

# How many mother detections have been traslated into cub detections? 544!
length(unique(str_sub(check$ID_obs,1,nchar(check$ID_obs)-2))) # With confidence that all these cubs were detected

## ---- Unsure cub detection ----

unsure <- os[which(os$N_cubs_observed != os$N_cubs_estimated & !is.na(os$N_cubs_observed) & os$N_cubs_observed != ""), ]
unsure <- unsure[which(!is.na(unsure$X)), ]

# For two observations we know that the cub is Noissette thanks to the comments, so we fill it up

for (i in 1:nrow(os)){
  if (os$cub_with_mother[i] == 1) next # To not repeat loop on cubs already added
  if (os$Remarks[i] == "The cub is Noisette" & !is.na(os$Remarks[i])) {
    add_row <- os[i,]     # Row to add with the different names of the cubs
    cubs <- "Noisette" # Cubs to add
    
    for (j in (1:as.numeric(os$N_cubs_observed[i]))) { # Add as many rows as cubs
      
      add_row <- os[i,] # To set a blank row for every cub
      
      # Fill information of the cub in the row
      add_row[,c(2,3)] <- cubs[j]
      
      add_row$Sex <- info$Sex[info$Probable_Individual == cubs[j]]
      add_row$Age <- add_row$Year - as.numeric(info$Year_birth[info$Probable_Individual == cubs[j]]) 
      
      add_row$Date_register_formated <- as.Date(add_row$Date_register, format = "%d/%m/%Y") # For Age class
      add_row$month <- month(ymd(add_row$Date_register_formated)) 
      
      if(is.na(add_row$month)) { # Si el mes es NA, lo pongo como no crítico, que serán pocos casos
        add_row$month <- 8 
      }
      
      if (add_row$Age <= 4 & add_row$Age > 1.5) {
        add_row$Age_class <- "Subadult" } else if (add_row$Age == 1.5 | add_row$Age == 1) {
          add_row$Age_class <- "Cub2" } else if (add_row$Age == 0 & add_row$month > 7)  {
            add_row$Age_class <- "Cub1" } else if (add_row$Age == 0 & add_row$month <= 7) {
              add_row$Age_class <- "Cub0" }
      
      add_row[,c(20:26)] <- NA # Because it is the cub
      
      add_row <- add_row[,-c(36,37)]
      add_row$ID_obs <- paste(add_row$ID_obs,j, sep = ".")
      add_row$cub_with_mother <- 1 
      
      # Add row
      long <- nrow(os)
      os <- rbind(os[1:i,], add_row ,os[(i + 1):(long),])
    }
  }
}

# There are only 18 mother observations where we can't identify the cub, so we leave it as it is.

## Overall, how many have we added??

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
os1 <- read.csv("Seguiment_Ossos_Pirineus_1996_2020_taula_final.csv", header = TRUE, row.names = NULL, sep = ";")
os1 <- os1[,-1] 

length(os1$Age_class[which(os1$Age_class %in% c("Cub0","Cub1","Cub2"))]) # 370
length(os$Age_class[which(os$Age_class %in% c("Cub0","Cub1","Cub2"))]) # 1418, so many!

## ---- Save ----

# This might be an extra step to run after all scripts of data management.
# At least, it would be essential to run capture-recapture models
rownames(os) <- seq(1,nrow(os))

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
openxlsx::write.xlsx(os, 'Seguiment_Ossos_Pirineus_1996_2021_taula_final_2_cubLocations.xlsx')

# For normal CR we join it to radiotracking in next script. For SCR there is no need, we only use hair and camera data

