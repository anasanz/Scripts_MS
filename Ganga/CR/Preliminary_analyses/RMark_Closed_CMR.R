## -------------------------------------------------
##        Run closed CMR on sandgrouse data 2021
## ------------------------------------------------- 

rm(list=ls())

library(RMark)

setwd("D:/MargSalas/Ganga/Data")

load("cr_sandgrouse.RData")
capt.hist <- as.data.frame(capt.hist$ch)
colnames(capt.hist)[1] <- "ch"

head(capt.hist)

## -------------------------------------------------
##                 2021
## ------------------------------------------------- 

## ---- Default closed capture-recapture model (model Mb) ----

## model codes for type of CMR model, not covariate/parameter structure
Mb<-mark(data=capt.hist, model="Closed", brief=T)

##estimates on the link scale
Mb$results$beta
##detection probabilities (p and c) are estimated on the logit scale
##f0 (number of individuals never detected) is estimated on the log scale

##estimates on the real scale
Mb$results$real ##Remember: f0 is animals never detected, NOT total abundance

##derived estimate of total abundance (f0 + n observed)
Mb$results$derived


## ---- Mb2: submodels for parameters p and c, both are constant but differ from each other ----

##has to be a list, with element "formula", which has to begin with a ~
p.dot<-list(formula=~1)
c.dot<-list(formula=~1)

##run model Mb again, using the specified structure for p and c
##both are combined into a list with elements p and c under model.parameters
Mb.2<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.dot, c=c.dot), 
           brief=T)

##check that estimates are the same as before
Mb.2$results$real

Mb$results$derived
Mb.2$results$derived

## ---- M0, p=c, constant over time ----

p.dot.shared<-list(formula=~1,share=TRUE)
M0<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.dot.shared), 
         brief=T)
M0$results$real
M0$results$derived

## ---- Mt, p=c, varying over time ----

p.time.shared<-list(formula=~time,share=TRUE)
Mt<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.time.shared), 
         brief=T)
Mt$results$real
Mt$results$derived

## ---- Perform model selection ----

modlist<-collect.models()
modlist ##this collects all models in the workspace, including both (identical)
        ## versions of model Mb

modlist2<-collect.models(lx=c("M0", "Mb", "Mt"))
modlist2 ##this only collects the models you specify under lx

##notation in the table: if c(), that means p=c
##                       if c(~1), that means p different from c

##abundance under best model
Mt$results$derived

## ---- Mh: mixture model with individual heterogeneity in detection ----

#p=c, two groups of individuals
pmixture.shared=list(formula=~mixture, share=T)

#these kinds of models are considered a different model type, "FullHet"
Mh<-mark(data=capt.hist, model="FullHet", 
         model.parameters=list(p=pmixture.shared), brief=T)

Mh$results$real
Mh$results$derived

##You cannot compare a model of type "FullHet" with a model of type "Closed"
##using AIC
##In order to compare all 4 models (M0, Mb, Mt and Mh) ,you'd have to re-run 
## the first 3 using "FullHet" 

## ---- Compare all models with FullHet ----

Mb_FH<-mark(data=capt.hist, model="FullHet", brief=T)
Mb.2_FH<-mark(data=capt.hist, model="FullHet", model.parameters=list(p=p.dot, c=c.dot), 
           brief=T)
M0_FH<-mark(data=capt.hist, model="FullHet", model.parameters=list(p=p.dot.shared), 
         brief=T)
Mt_FH<-mark(data=capt.hist, model="FullHet", model.parameters=list(p=p.time.shared), 
         brief=T)

modlist3 <- collect.models(lx=c("Mh", "Mb_FH", "M0_FH", "Mt_FH"))
modlist3

# Best is still the model with different p over time (occassions)
Mt_FH$results$real
Mt_FH$results$derived

Mt$results$real
Mt$results$derived

## -------------------------------------------------
##                 2022
## ------------------------------------------------- 

## -------------------------------------------------
##        Run closed CMR on sandgrouse data 2021
## ------------------------------------------------- 

setwd("D:/MargSalas/Ganga/Data")

load("cr_sandgrouse_2022.RData")
capt.hist <- as.data.frame(capt.hist$ch)
colnames(capt.hist)[1] <- "ch"

head(capt.hist)

## ---- Default closed capture-recapture model (model Mb) ----

## model codes for type of CMR model, not covariate/parameter structure
Mb<-mark(data=capt.hist, model="Closed", brief=T)

##estimates on the link scale
Mb$results$beta
##detection probabilities (p and c) are estimated on the logit scale
##f0 (number of individuals never detected) is estimated on the log scale

##estimates on the real scale
Mb$results$real ##Remember: f0 is animals never detected, NOT total abundance

##derived estimate of total abundance (f0 + n observed)
Mb$results$derived


## ---- Mb2: submodels for parameters p and c, both are constant but differ from each other ----

##has to be a list, with element "formula", which has to begin with a ~
p.dot<-list(formula=~1)
c.dot<-list(formula=~1)

##run model Mb again, using the specified structure for p and c
##both are combined into a list with elements p and c under model.parameters
Mb.2<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.dot, c=c.dot), 
           brief=T)

##check that estimates are the same as before
Mb.2$results$real

Mb$results$derived
Mb.2$results$derived

## ---- M0, p=c, constant over time ----

p.dot.shared<-list(formula=~1,share=TRUE)
M0<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.dot.shared), 
         brief=T)
M0$results$real
M0$results$derived

## ---- Mt, p=c, varying over time ----

p.time.shared<-list(formula=~time,share=TRUE)
Mt<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.time.shared), 
         brief=T)
Mt$results$real
Mt$results$derived

## ---- Perform model selection ----

modlist<-collect.models()
modlist ##this collects all models in the workspace, including both (identical)
## versions of model Mb

modlist2<-collect.models(lx=c("M0", "Mb", "Mt"))
modlist2 ##this only collects the models you specify under lx

##notation in the table: if c(), that means p=c
##                       if c(~1), that means p different from c

## ALL THE SAME AIC

##abundance under best model
Mt$results$derived

## ---- Mh: mixture model with individual heterogeneity in detection ----

#p=c, two groups of individuals
pmixture.shared=list(formula=~mixture, share=T)

#these kinds of models are considered a different model type, "FullHet"
Mh<-mark(data=capt.hist, model="FullHet", 
         model.parameters=list(p=pmixture.shared), brief=T)

Mh$results$real
Mh$results$derived

##You cannot compare a model of type "FullHet" with a model of type "Closed"
##using AIC
##In order to compare all 4 models (M0, Mb, Mt and Mh) ,you'd have to re-run 
## the first 3 using "FullHet" 

## ---- Compare all models with FullHet ----

Mb_FH<-mark(data=capt.hist, model="FullHet", brief=T)
Mb.2_FH<-mark(data=capt.hist, model="FullHet", model.parameters=list(p=p.dot, c=c.dot), 
              brief=T)
M0_FH<-mark(data=capt.hist, model="FullHet", model.parameters=list(p=p.dot.shared), 
            brief=T)
Mt_FH<-mark(data=capt.hist, model="FullHet", model.parameters=list(p=p.time.shared), 
            brief=T)

modlist3 <- collect.models(lx=c("Mh", "Mb_FH", "M0_FH", "Mt_FH"))
modlist3

# Best is still the model with different p over time (occassions)
Mt_FH$results$real
Mt_FH$results$derived

Mt$results$real
Mt$results$derived



