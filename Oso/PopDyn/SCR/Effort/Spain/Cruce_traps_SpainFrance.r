
## -------------------------------------------------
##           Join trap information Spain
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)

setwd("D:/MargSalas/Oso/Datos/Effort_raw/Spain")

## ---- Load traps that Maelis used for 2017-2019 (all active the three years) ----

trap1 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/UTM trampas pelos y camaras_1719_Maelis.xlsx") %>%
  select(-...9)
trap1$Xutm[trap1$`Codi dispositiu` == "PS525TP130"] <- 345667 # Error

## ---- Join with 2017 ----
  
trap2017 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Activitat-Diversitat Mamífers DEFINITIVA.xlsx", sheet = 8) %>%
  rename("TrapID" = "Nombre Estación") %>%
  group_by(TrapID) %>%
  summarise(
    Xutm = unique(Xutm, na.rm = TRUE),
    Yutm = unique(Yutm, na.rm = TRUE)) %>%
  filter(!is.na(Xutm))
trap2017[6,2] <- 340840
trap2017[11,3] <- 4736068
trap2017 <- trap2017[-c(10,15,24,30),] # Remove errors
trap2017 <- as.data.frame(trap2017)
trap2017 <- rbind(trap2017,c("Raspamala TP070", 341184, 4736498)) # Los que están solo en resumen que coinciden con general
trap2017 <- rbind(trap2017,c("Selves arbre TP088", 360359, 4724861))
trap2017$Xutm <- as.numeric(trap2017$Xutm)
trap2017$Yutm <- as.numeric(trap2017$Yutm)

trap2 <- left_join(trap1, trap2017, by = c("Xutm", "Yutm"))
sum(!is.na(trap2$TrapID)) 

# De las 31 trampas de Activitat Mamifers que me mandó Santi, hay 7 que no aparecen en la tabla
# de trampas que le dieron a Maelis.

missing <- as.numeric(trap2017$TrapID %in% trap2$TrapID[complete.cases(trap2$TrapID)])
trap2017$TrapID[which(missing == 0)]
colnames(trap2)[10] <- "Activitat2017"

# Traps missing: 
##Gueron, there are 2 but not the same one
## Con Trotxa Estallo TP072 no se bien lo que pasa

## ---- Join with 2018 ----
trap2018 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Activitat-Diversitat Mamífers DEFINITIVA.xlsx", sheet = 9) %>%
  rename("TrapID" = "Localización Cámara") %>%
  group_by(TrapID) %>%
  summarise(
    Xutm = unique(Xutm, na.rm = TRUE),
    Yutm = unique(Yutm, na.rm = TRUE)) %>%
  filter(!is.na(Xutm))
trap2018[8,3] <- 4736068
trap2018 <- as.data.frame(trap2018)
trap2018 <- rbind(trap2018,c("Fontallada TP167", 366693, 4720270))
trap2018 <- rbind(trap2018,c("Bosc de Marimanha TP170", 340840, 4735635))
#trap2018 <- rbind(trap2018,c("Port de Tavascan TP142", 357078, 4727913))
trap2018$Xutm <- as.numeric(trap2018$Xutm)
trap2018$Yutm <- as.numeric(trap2018$Yutm)

trap2 <- left_join(trap2, trap2018, by = c("Xutm", "Yutm"))
sum(!is.na(trap2$TrapID)) 


# De las 32 trampas de Activitat Mamifers que me mandó Santi, hay 17 que no aparecen en la tabla
# de trampas que le dieron a Maelis. Muchisimas!!!!

missing <- as.numeric(trap2018$TrapID %in% trap2$TrapID[complete.cases(trap2$TrapID)])
trap2018$TrapID[which(missing == 0)]
colnames(trap2)[11] <- "Activitat2018"

## Missing:
## Baiasca TP342, there are others in Baiasca in the general (Maelis), but not this one
## Bessiberris no existe (TP290 en resumen pero tampoco aparece en general)
## Fontallada TP111 es raro, existe en los franceses pero con un numero cambiado y el Aixeus

trap2019 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Activitat-Diversitat Mamífers DEFINITIVA.xlsx", sheet = 10) %>%
  rename("TrapID" = "Localización Cámara") %>%
  group_by(TrapID) %>%
  summarise(
    Xutm = unique(Xutm, na.rm = TRUE),
    Yutm = unique(Yutm, na.rm = TRUE)) %>%
  filter(!is.na(Xutm))
trap2019 <- as.data.frame(trap2019)
trap2019 <- trap2019[-c(6),]
#trap2019 <- rbind(trap2019,c("Bassa Isavarre TP105", 343878, 4723408)) # Los que están solo en resumen que coinciden con general
trap2019 <- rbind(trap2019,c("Bosc de Marimanha TP170", 340840, 4735635))
trap2019 <- rbind(trap2019,c("Pina TP027", 345612, 4730195))
trap2019 <- rbind(trap2019,c("Port de Tavascan TP142", 357078, 4727913))
trap2019 <- rbind(trap2019,c("Raspamala TP070", 341184, 4736498))
# Rovinets TP086 has 2 different coordinates in general and trap2019, dont know which one is good
trap2019 <- rbind(trap2019,c("Salau TP022", 344742, 4734952))
trap2019 <- rbind(trap2019,c("Selves de Lladorre TP090", 361383, 4724994)) # This one is called TP091 en general
trap2019 <- rbind(trap2019,c("Solana de Montgós TP069", 340145, 4736256)) # This one is called TP068 en general
trap2019 <- rbind(trap2019,c("Solana de Sorpe TP135", 338366, 4723003)) 
#trap2019 <- rbind(trap2019,c("Bosc de Bonabé TP098", 341982, 4735489)) 
trap2019$Xutm <- as.numeric(trap2019$Xutm)
trap2019$Yutm <- as.numeric(trap2019$Yutm)

trap2 <- left_join(trap2, trap2019, by = c("Xutm", "Yutm"))
sum(!is.na(trap2$TrapID)) 

duplicated(trap2$`Codi dispositiu`)

# De las 32 trampas de Activitat Mamifers que me mandó Santi, hay 7 que no aparecen en la tabla
# de trampas que le dieron a Maelis. 
trap2019$TrapID %in% trap2$TrapID[complete.cases(trap2$TrapID)]

missing <- as.numeric(trap2019$TrapID %in% trap2$TrapID[complete.cases(trap2$TrapID)])
trap2019$TrapID[which(missing == 0)]
colnames(trap2)[12] <- "Activitat2019"

trap2$Activitat2017 <- ifelse(!is.na(trap2$Activitat2017),"ACTIVA","n.a.")
trap2$Activitat2018 <- ifelse(!is.na(trap2$Activitat2018),"ACTIVA","n.a.")
trap2$Activitat2019 <- ifelse(!is.na(trap2$Activitat2019),"ACTIVA","n.a.")

setwd("D:/MargSalas/Oso/Datos/Effort_raw/Spain")
  
write.csv(trap2, file = "UTM trampas pelos y camaras_CRUCE_Activitat Mamifers.csv")

