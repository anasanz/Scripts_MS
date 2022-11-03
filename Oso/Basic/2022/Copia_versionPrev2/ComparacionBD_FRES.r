
## -------------------------------------------------
##    Comparación base de datos Francia - España 
##                   2017-2019
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)
library(sf)
library(rgdal)
library(mapview)
library(lubridate)

# En ambas bases de datos:
# Choose column that is a priority to check
# Priorizar los individuos identificados (priorid 2), y 
# Priorizar más (priorid 1):
#  - los estrictamente sistemáticos ("Itinéraire", "Suivi appareil photo") 
#  - en la temporada de estudio (priorid 1)
#  - De indicios usados en mi estudio (fotos,pelo): 
#     "Poils,appât smola", "vidéo automatique", "Poils,appât térébenthine", "poils (spontanés)", "photographie automatique", "photographie / vidéo"

# Base de datos FR raw (Maelis):

dat <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/Prospind_2017-2019_Maelisv2.xlsx") %>%
  janitor::clean_names() %>%
  mutate(libelle_type_operation_terrain = as_factor(libelle_type_operation_terrain),
         date = as_date(date_temoin),
         an = year(date),
         month = month(date)) %>%
  rename(x = x_lamb_ii_etendu,
         y = y_lamb_ii_etendu) 
dat <- dat %>%
  dplyr::select(an, libelle_nom_individu, date_temoin, date_ours_estime, libelle_commune, lieu_dit, x, y, libelle_type_operation_terrain, libelle_type_indice, libelle_validation_indice,  libelle_itineraire, libelle_ap_photo, date, month, remarques)

dat$priori <- 0
dat$priori[which(dat$libelle_nom_individu != "Indéterminé")] <- 2

for (i in 1:nrow(dat)) {
  if (dat$month[i] < 12 & dat$month[i] > 4 &  # En temporada de estudio
      dat$priori[i] == 2 &  # Individuos identificados
      dat$libelle_type_operation_terrain[i] %in% c("Itinéraire", "Suivi appareil photo") & # Seguimiento sistemático
      dat$libelle_type_indice[i] %in% c("Poils,appât smola", "vidéo automatique", "Poils,appât térébenthine", "poils (spontanés)", "photographie automatique", "photographie / vidéo")) # Indicios de interés 
    { 
  dat$priori[i] <- 1}
}

length(which(dat$priori == 1))

dat <- dat[,-c(14,15)]

# Add relevant columns:
dat$Absent_in_SP <- 0

setwd("D:/MargSalas/Oso/Datos/Data_raw/Check")
openxlsx::write.xlsx(dat, 'Data_France(Maelis)_1719_check.xlsx')


# Base de datos ESP precoordinates (Victor y Ana):

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2021_Pre-coordinates.csv", header = TRUE, row.names = NULL, sep = ";") 

# Formating
Encoding(os$Remarks) <- "UTF-8" # Change encoding for weird characters (Â imported from excel)
Encoding(os$X) <- "UTF-8" 
os$Y <- as.character(os$Y)
Encoding(os$Y) <- "UTF-8"
Encoding(os$Municipality) <- "UTF-8"
Encoding(os$Site) <- "UTF-8"
Encoding(os$Database) <- "UTF-8"
Encoding(os$Region) <- "UTF-8"
colnames(os)[1] <- "Confirmed_Individual"

os$X <- str_trim(os$X, side = c("both")) # Remove white spaces
os$Y <- str_trim(os$Y, side = c("both"))

os$Probable_Individual <- str_trim(os$Probable_Individual, side = c("both")) # Remove white spaces
os$Probable_Individual[which(os$Probable_Individual == "")] <- "Indetermined" # Need that all white spaces are "Indetermined"

# Create month to subset priority 1
os$Date_register_formated <- os$Date_register 
os$Date_register_formated <- as.Date(os$Date_register_formated, format = "%d/%m/%Y")
os$month <- month(ymd(os$Date_register_formated)) 

# Priority
os$priori <- 0
os$priori[which(os$Probable_Individual != "Indetermined" & 
                  os$Year %in% c(2017, 2018, 2019) &
            os$Database %in% c("Seguiment França Ossos 2010-2017.xls", "Seguiment França Ossos 2016-2020.xls"))] <- 2

for (i in 1:nrow(os)) {
  if (os$month[i] < 12 & os$month[i] > 4 &  # En temporada de estudio
      os$priori[i] == 2 &  # Individuos identificados
      os$Method[i] %in% c("Transect", "Sampling_station") & # Seguimiento sistemático
      os$Obs_type[i] %in% c("Photo","Photo/Video", "Hair", "Video")) {  # Indicios de interés
    os$priori[i] <- 1}
    }

length(which(os$priori == 1))

os <- os[,-c(30,31)]

# Add relevant columns:
os$Absent_in_FR <- 0
os$Hair_trap <- 0


setwd("D:/MargSalas/Oso/Datos/Data_raw/Check")
openxlsx::write.xlsx(os, "Seguiment_Ossos_Pirineus_1996_2021_Pre-coordinates_CHECK17-19.xlsx")


