## -------------------------------------------------
##               SCR0 MODEL WOLVERINE DATA  
## ------------------------------------------------- 

rm(list=ls())

library(SCRbayes)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions/scrbook")
source("SCR23darray.r")

setwd("D:/Stats/SCRbook/Files")
load("wolverine.rda")

## ---- Data format ----

# Encounter Data File:
# At which traps and when each individual encounter occured
edf <- wolverine$wcaps
# Each row representing a unique encounter event including 
# the trap identity, the individual identity.


# Trap Deployment File:
# Coordinates of each trap 
# + which sample occasions each trap was operating (binary indicator)
# Occasions = Days
tdf <- wolverine$wtraps 
head(tdf)

# Trap locations 
traps <- wolverine$wtraps 
traplocs <- traps[,2:3] 

# Vector of sample sizes (number of days each trap was active)
K <- apply(traps[,4:ncol(traps)],1,sum)

# Format to include the effort
# Array of nind x ntraps x noc
y3d <- SCR23darray(wolverine$wcaps,wolverine$wtraps) 
# Binomial encounter frequencies nxJ
y <- apply(y3d,c(1,2),sum)







