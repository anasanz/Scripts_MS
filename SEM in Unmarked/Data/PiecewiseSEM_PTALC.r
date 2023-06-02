#####################################
#        Data processing             #
#        Unmarked model              #
#        PiecewiseSEM                #
#   Example: Alauda Arvensis         #
#####################################

rm(list=ls())


#Set path of your file
#path <- "~/Desktop/LBOE2240 Stage de recherche CTFC/FARMDINDIS weather project/spatial_model/"
path <- "D:/Otros/Adrien/SEM/"
## Module load
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(unmarked)
library(MuMIn)
library(tidyverse)
library(lme4)
library(lmerTest)

setwd("D:/MargSalas/Scripts_MS/SEM in Unmarked/Sim")
source('Standardization R function unmarked GLMM.R')
source('GOF function improved.R')

# Upload datasets

records_birds<-fread(paste0(path, "dataset_birds_corrected.csv"))
covariates <- fread(paste0(path, "dataset_covariates_corrected_complete.csv"))

# Make sure that all factors ar well coded (very import to define the unmarked object)

# Birds dataset
records_birds$ID <- as.factor(records_birds$ID)
records_birds$Any <- as.factor(records_birds$Any)
records_birds$Especie <- as.factor(records_birds$Especie)
records_birds$Us <- as.factor(records_birds$Us)


# Covariates dataset
covariates$ID <- as.factor(covariates$ID)
covariates$Ambit <- as.factor(covariates$Ambit)
covariates$year <- as.factor(covariates$year)
covariates$Vent <- as.factor(covariates$Vent)
covariates$Observador <- as.factor(covariates$Observador)
covariates$Nuvolositat <- as.factor(covariates$Nuvolositat)
covariates$Temperatura <- as.numeric(covariates$Temperatura)

#-------------- Refine dataset for your species of interest ----------

# Here, we will choose the skylark (Code = ALARV) as a model species

# Convert the tag in a factor format, to store all of them (need to done only after removing
# the NA observations for cereals/fallows)
records_birds$Dades_ocells_Id_transecte_detall <- as.factor(records_birds$Dades_ocells_Id_transecte_detall)
covariates$Dades_ocells_Id_transecte_detall <- as.factor(covariates$Dades_ocells_Id_transecte_detall)

unique(records_birds$Especie)
# Now, we need to extract data in birds observations
#records_species <- subset(records_birds, Especie == "TERAX" & Sexe == 1) # Males LB
records_species <- subset(records_birds, Especie == "PTALC") # Males LB


# Skylark is a cereal specialist, so we need to refine dataset only to observations associated
# in cereals (not possible to include site covariates in HDS models)
records_species_light <- subset(records_species, Us == "G")

# I can't difference the observation type, so I am going to keep only the seen ones
records_species_light <- subset(records_species_light, Tipus_observacio == "V")

#------------ Change format of data to define HDS model --------------

##Convert tag for distance band in a "true" distance mesure (easier to convert in a umf object)

records_species_light$distance <- "0"
records_species_light$distance[records_species_light$Banda=="1"]<-"12"
records_species_light$distance[records_species_light$Banda=="2"]<-"35"
records_species_light$distance[records_species_light$Banda=="3"]<-"75"
records_species_light$distance[records_species_light$Banda=="4"]<-"150"
records_species_light$distance[records_species_light$Banda=="5"]<-"350"

# and make sure that this variable is numeric
records_species_light$distance <- as.numeric(records_species_light$distance)

# You can't use obserations covariates in the analysis. So, i keep only the columns that I need.
# Because we have multi-year dataset,"Dades_ocells_Id_transecte_detall" is a unique code 
# for the site-year combination
records_species_light <- records_species_light[,.(Dades_ocells_Id_transecte_detall, distance)]

# Change the format of the dataset (rows for observations and columns for distance bands)
# We use the same distance breaks that bands previously defined
yDat <- formatDistData(records_species_light, distCol="distance",
                       transectNameCol="Dades_ocells_Id_transecte_detall", dist.breaks=c(0, 25, 50, 100,200, 500))

# Because rownames of YDat dataset correspond to transect tag, we need to make the same with
# covariates (to make sure that we have good covariates for the good transect-year combination)

covariates <- arrange(covariates, Dades_ocells_Id_transecte_detall)
covariates <- covariates %>% remove_rownames %>% column_to_rownames(var="Dades_ocells_Id_transecte_detall")

for (i in 1:ncol(covariates)){ # ASP: There were some factors here
  if(class(covariates[,i]) == "character"){
    covariates[,i] <- as.factor(covariates[,i])
  }}

# Need to scale the covariates in advance, because names in the model need to be the same than in siteCovs dataframe
# in piecewiseSEM functions

covariates$fallow_height_sc <- scale(covariates$fallow_height)
covariates$guarets_sc <- scale(covariates$guarets)
covariates$spring_precipitation_sc <- scale(covariates$spring_precipitation)

covariates <- covariates[,colnames(covariates) %in% c("ID", "fallow_height_sc", "guarets_sc", "spring_precipitation_sc")]

# Now, we can build the "umf object", the basic format for analysis in Unmarked 
umf <- unmarkedFrameDS(y=as.matrix(yDat), siteCovs = covariates,  survey="line",
                       dist.breaks=c(0,25,50, 100, 200, 500), tlength=rep(500, 1694),
                       unitsIn="m")

# I remove NA because otherwise piecewise functions don't work
umf@siteCovs <- umf@siteCovs[which(complete.cases(umf@siteCovs)),]
umf@tlength <- umf@tlength[which(complete.cases(umf@siteCovs))]
umf@y <- umf@y[which(complete.cases(umf@siteCovs)),]

# Based on this umf, we can see the summary of transects information and make a histogram
summary(umf)
hist(umf,xlab="distance (m)", main = "Histogram of observations restricted to cereals")

## ---- Build the piecewise sem model ----

##model: Spring precipitation -> Cereal Height -> Abundance ALARV <- nº of surrounding cereal crops

#first component: abundance (DS) model
ds.mlu <- distsamp(~ fallow_height_sc  ~ fallow_height_sc + guarets_sc + (1|ID), 
                   data = umf, keyfun='halfnorm', output = "abun")
summary(ds.mlu) # RESULT: No effect of cereal height on detection (p = 0.43)
#                         No effect of cereal height en abundance (p = 0.9)
#                         No effect of n cereal crops en abundance (p = 0.5)

#second component: model on fallow height
siteCovs<-data.frame(fallow_height_sc = umf@siteCovs$fallow_height_sc, 
                     spring_precipitation_sc = umf@siteCovs$spring_precipitation_sc,
                     guarets_sc = umf@siteCovs$guarets_sc,
                     ID = umf@siteCovs$ID)

lm.ml <- lmer( fallow_height_sc ~ spring_precipitation_sc + (1|ID), siteCovs)
summary(lm.ml) # RESULT: Positive effect of spring precipitation in cereal height

##bundle both paths into full piecewise SEM
##names not necessary but helpful for some output
M.focal<-list(lmm=lm.ml, DS=ds.mlu)

##get standardized coefficients, either for each component model
stand_beta(lm.ml)
stand_beta(ds.mlu)

##or for list of models making up full SEM
stand_beta_wrap(M.focal)

# RESULT: Estimates don't really change with standardize vS independent model coefficients

## ---- Identify missing paths ----

## Here I assume missing paths are:
## Spring precipitation -> Abundance (in DS model)

##Build fully saturated model with all missing paths
lm.full <- lmer( fallow_height_sc ~ spring_precipitation_sc + (1|ID), siteCovs) # In this case it is the same

ds.mlu.full <- distsamp(~ fallow_height_sc  ~ fallow_height_sc + guarets_sc + spring_precipitation_sc + (1|ID), 
                        data = umf, keyfun='halfnorm', output = "abun")

##Bundle into saturated piecewise SEM
M.sat<-list(lm.full, ds.mlu.full)

#includes missing path coefficient estimates
GOF_LL(M.focal, M.sat) # No missing path

#returns:
#  AIC of focal model
#  p.overall: p-value for test focal against saturated model
#  df: differences in number of paramters between the two models
#  X2: Chi-square value used to determine p
#  coefficients of paths that were added in the saturated model

##if p-value>0.05, there is support for the focal model over the saturated model

#get coefficients, p-values of missing paths only
missing.path(M.focal, M.sat) 

# There is no missing path