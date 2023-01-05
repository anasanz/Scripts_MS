

rm(list = ls())


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")


load("edf1721.RData")
edf <- edf[-which(duplicated(edf)), ]


load("tdf2017.RData")
load("tdf2018.RData")
load("tdf2019.RData")
load("tdf2020.RData")
load("tdf2021.RData")

nocc <- 7 # number of occasion (month) per session (year) : may to november


## ---- 2017 ----

# effort covariate : factor which account for the country and the number of visits per occasion in France
# 1 = one visit per month in France 
# 2 = two visits per month in France
# 3 = trap is in Spain

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2017$trap))
  
  for (i in 1:length(tdf2017$trap)){
    if (tdf2017$pays[i] == "Espagne"){ # if the trap is in Spain
        eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
            eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
  name <- paste("effort", j, sep = ".")
  tdf2017[name] <- as.factor(eff)
 } }


## ---- 2018 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2018$trap))
  
  for (i in 1:length(tdf2018$trap)){
    if (tdf2018$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2018[name] <- as.factor(eff)
  } }

## ---- 2019 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2019$trap))
  
  for (i in 1:length(tdf2019$trap)){
    if (tdf2019$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2019[name] <- as.factor(eff)
  } }

## ---- 2020 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2020$trap))
  
  for (i in 1:length(tdf2020$trap)){
    if (tdf2020$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2020[name] <- as.factor(eff)
  } }

## ---- 2021 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2021$trap))
  
  for (i in 1:length(tdf2021$trap)){
    if (tdf2021$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2021[name] <- as.factor(eff)
  } }

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

save(tdf2017, file = "tdf2017_effort.RData")
save(tdf2018, file = "tdf2018_effort.RData")
save(tdf2019, file = "tdf2019_effort.RData")
save(tdf2020, file = "tdf2020_effort.RData")
save(tdf2021, file = "tdf2021_effort.RData")

## -------------------------------------------------
##         Re-run but with data without the 
##    cubs that went with mother as independent obs
## ------------------------------------------------- 

rm(list = ls())


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")


load("edf1721_nocubsmother.RData")
edf <- edf[-which(duplicated(edf)), ]


load("tdf2017_nocubsmother.RData")
load("tdf2018_nocubsmother.RData")
load("tdf2019_nocubsmother.RData")
load("tdf2020_nocubsmother.RData")
load("tdf2021_nocubsmother.RData")

nocc <- 7 # number of occasion (month) per session (year) : may to november


## ---- 2017 ----

# effort covariate : factor which account for the country and the number of visits per occasion in France
# 1 = one visit per month in France 
# 2 = two visits per month in France
# 3 = trap is in Spain

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2017$trap))
  
  for (i in 1:length(tdf2017$trap)){
    if (tdf2017$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2017[name] <- as.factor(eff)
  } }


## ---- 2018 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2018$trap))
  
  for (i in 1:length(tdf2018$trap)){
    if (tdf2018$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2018[name] <- as.factor(eff)
  } }

## ---- 2019 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2019$trap))
  
  for (i in 1:length(tdf2019$trap)){
    if (tdf2019$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2019[name] <- as.factor(eff)
  } }

## ---- 2020 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2020$trap))
  
  for (i in 1:length(tdf2020$trap)){
    if (tdf2020$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2020[name] <- as.factor(eff)
  } }

## ---- 2021 ----

for (j in 1 : nocc) {
  eff <- rep(0, length(tdf2021$trap))
  
  for (i in 1:length(tdf2021$trap)){
    if (tdf2021$pays[i] == "Espagne"){ # if the trap is in Spain
      eff[i] <- 3}  
    
    else if (j == 1 | j == 2 | j == 5){ # if is in france and if we are in may, june or september
      eff[i] <- 2}  # the site is visited 2 times at this occasion
    else {eff[i] <- 1} # the site is visited once
    
    name <- paste("effort", j, sep = ".")
    tdf2021[name] <- as.factor(eff)
  } }

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

save(tdf2017, file = "tdf2017_nocubsmother_effort.RData")
save(tdf2018, file = "tdf2018_nocubsmother_effort.RData")
save(tdf2019, file = "tdf2019_nocubsmother_effort.RData")
save(tdf2020, file = "tdf2020_nocubsmother_effort.RData")
save(tdf2021, file = "tdf2021_nocubsmother_effort.RData")



