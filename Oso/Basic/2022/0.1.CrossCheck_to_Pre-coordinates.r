
## -------------------------------------------------
##     Adapt cross check of Elena Pi to standard Pre-coordinates
## ------------------------------------------------- 

rm(list = ls())

#library(sp)
#library(rgdal)
#library(stringr)
library(tidyverse)
library(lubridate)

# Better to open it as csv, for format purposes (ex Date)
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
check <- read.csv("CAT_final_crosscheck_cubs.csv", header = TRUE, row.names = NULL, sep = ";")
Encoding(check$Remarks) <- "UTF-8" # Change encoding for weird characters (Â imported from excel)
Encoding(check$X) <- "UTF-8" 
check$Y <- as.character(check$Y)
Encoding(check$Y) <- "UTF-8"
Encoding(check$Municipality) <- "UTF-8"
Encoding(check$Site) <- "UTF-8"
Encoding(check$Database) <- "UTF-8"
Encoding(check$Region) <- "UTF-8"
colnames(check)[1] <- "N"
colnames(check)[32] <- "Obs_type(FR)"
colnames(check)[38] <- "N(FR)"


## ---- 1. Corrected coordinates ----

# For the ones where the coordinates have been corrected with the French data base,
# Assign the new coordinate system and put "LambertII" as coordinate system of origin

mod <- grep("lamb", check$variables, ignore.case = TRUE) # Index of the coordinates modified

for (i in 1:nrow(check)){
  if (!(i %in% mod)) next # Only for the ones that are modified
  check$X[i] <- check$x_LambertII[i]
  check$Y[i] <- check$y_LambertII[i]
  check$Coordinate_system[i] <- "LambertII"
}

## ---- 2. Corrected registers ----

# Reflect if a register has been corrected by adding "+ modified recap_individualisation.xlxs" in database field

for (i in 1:nrow(check)){
  if((check$variables[i] == "") | (check$variables[i] == "ALL")) next # Only modified registers (if its not modified, or fully added I dont add the text)
    check$Database[i] <- paste(check$Database[i], " + modified recap_individualisation.xlxs", sep = "")
  }


## ---- 3. Added registers ----

for (i in 1:nrow(check)){
  if(check$variables[i] == "") next
  if(check$variables[i] == "ALL"){
    date <- as.Date(check$Date[i], format = "%d/%m/%Y")
    check$Year[i] <- year(ymd(date))
    check$X[i] <- check$x_LambertII[i]
    check$Y[i] <- check$y_LambertII[i]
    check$Coordinate_system[i] <- "LambertII"
    check$Category[i] <- 1
    check$Database[i] <- "recap_individualisation.xlxs"
  }}

## ---- 4. Manage new and old columns ----
## ---- 4.1. Confirmed individual ----

# Substitute Confirmed_Individual by new Confirmed_ID
for (i in 1:nrow(check)){
  if(check$Confirmed_ID[i] == ""){
    check$Confirmed_ID[i] <- check$Confirmed_Individual[i]
  }}

# Fill the empty Probable_Individual with the Confirmed_ID (newly added rows)
for (i in 1:nrow(check)){
  if(check$Probable_Individual[i] == ""){
    check$Probable_Individual[i] <- check$Confirmed_ID[i]
  }}

# Select new confirmed and probable

check <- check %>% 
  select(-Confirmed_Individual) %>%
  rename("Confirmed_Individual" = "Confirmed_ID")

# Modify names as they were before
# Confirmed individual
check$Confirmed_Individual <- str_to_sentence(check$Confirmed_Individual)
check$Confirmed_Individual <- gsub("_","-",check$Confirmed_Individual)
check$Confirmed_Individual <- gsub("slo","SLO",check$Confirmed_Individual)
check$Confirmed_Individual[check$Confirmed_Individual == "Camille"] <- "Camille_AspeOuest"
check$Confirmed_Individual[check$Confirmed_Individual == "Camille-aspeouest"] <- "Camille_AspeOuest"
check$Confirmed_Individual[check$Confirmed_Individual == "Canellito"] <- "Cannellito"
check$Confirmed_Individual[check$Confirmed_Individual == "Callisto"] <- "Callista"
check$Confirmed_Individual[check$Confirmed_Individual == "Douillous"] <- "Douilloux"
check$Confirmed_Individual[check$Confirmed_Individual == "Nheu "] <- "Nheu"


unique(check$Confirmed_Individual)

# Probable individual
check$Probable_Individual <- str_to_title(check$Probable_Individual)
check$Probable_Individual <- gsub("_","-",check$Probable_Individual)
check$Probable_Individual <- gsub("slo","SLO",check$Probable_Individual)
check$Probable_Individual <- gsub("Slo","SLO",check$Probable_Individual)
check$Probable_Individual[check$Probable_Individual == "Camille"] <- "Camille_AspeOuest"
check$Probable_Individual[check$Probable_Individual == "Camille-aspeouest"] <- "Camille_AspeOuest"
check$Probable_Individual[check$Probable_Individual == "Canellito"] <- "Cannellito"
check$Probable_Individual[check$Probable_Individual == "Callisto"] <- "Callista"
check$Probable_Individual[check$Probable_Individual == "Douillous"] <- "Douilloux"
check$Probable_Individual[check$Probable_Individual == "Doiulloux"] <- "Douilloux"
check$Probable_Individual[check$Probable_Individual == "Griboulle"] <- "Gribouille"


# Cubs
check$ID_cub1[check$ID_cub1 == "Canellito"] <- "Cannellito"
check$ID_cub1[check$ID_cub1 == "Doiulloux"] <- "Douilloux"
check$ID_cub2[check$ID_cub2 == "Griboulle"] <- "Gribouille"

# Check that names match with info_individuals

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)
info <- info[,c(4,5,8)]
info[nrow(info)+1,] <- c("Indetermined")
info[nrow(info)+1,] <- c("-")
info[nrow(info)+1,] <- c("")

check$Confirmed_Individual[!(check$Confirmed_Individual %in% info$ID)]
check$Probable_Individual[!(check$Probable_Individual %in% info$ID)]
check$ID_cub1[!(check$ID_cub1 %in% info$ID)]
check$ID_cub2[!(check$ID_cub2 %in% info$ID)]
check$ID_cub3[!(check$ID_cub3 %in% info$ID)]

## ---- 4.2. Sort out date ----

# New way of sorting out dates: 
# - Column Date: Date contact if available. If not, date register
# - Column Date.precision: "Contact", "Register", or 
                        # "Most_accurate" if it was filled by Elena (because we don't know which it is, but it's contact if available)

#check$Date.precision <- ""

## ---- 4.2.Sort out coordinates ----

# There are some of the new column that are NA. Filled them with the old column
for (i in 1:nrow(check)){
  if(is.na(check$Coordinates.precision_LambertII[i])){
    check$Coordinates.precision_LambertII[i] <- check$Coordinates.precision[i]}
}
     
check$Coordinates.precision <- check$Coordinates.precision_LambertII # Keep the new


## ---- 4.3. Sort out rest of the columns ----

check$Method <- str_to_sentence(check$Method)
check$Obs_type <- str_to_sentence(check$Obs_type)
check <- check %>% select(-c("N","x_LambertII","y_LambertII","Coordinates.precision_LambertII","variables","N(FR)","status_det"))
check <- check[,c(1:27,32,28:31)]

## ---- Fix coordinates ED50 :( ----
# Some of the registers of the crosscheck have a wrong coordinate format
# There is only commas in ED50, so I can fix directly there

for (i in grep(",", check$Y)){
  check$Y[i] <- as.numeric(gsub(",","",check$Y[i]))
  check$Y[i] <- sub("(.{7})(.*)", "\\1.\\2", check$Y[i])
  }


setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
openxlsx::write.xlsx(check, 'Seguiment_Ossos_Pirineus_1996_2021_Pre-coordinates.xlsx')



