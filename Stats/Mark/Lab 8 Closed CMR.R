############## WFC 124 - Lab 8 ######################

##### Closed population capture-mark-recapture

##install package once
#install.packages("RMark")
#requires MARK to be installed on your computer

rm(list=ls())

library(RMark)

#### Load data from CMR study in an enclosure (Edwards and Eberhardt 1967)

dat<-read.csv("Rabbit.data.csv")


###CHALLENGE###
#Generate the following CMR summary statistics:
#1. How many trapping occasions were there?
#2. How many individuals were captured?
#3. How many individuals were captured more than once? (Hint: add up how many times each individual was captured. Then see how many have a sum of >1)
#4. How many times was an individual captured, on average?


##1. number of capture occasions
ncol(dat)
  
##2. number of individuals captured
nrow(dat)

##3. number of individuals captured more than once
ncap<-apply(dat,1,sum) ##number of captures per individual
length(which(ncap>1))  ##number of individuals captured more than 1 time

##4. avg number of captures per individual
mean(ncap)


##load rabbit data in form of a data frame with character string encounter histories
##data set comes with package RMark
## this is the format for analysis in RMark
data("edwards.eberhardt")

##give data a shorter name
capt.hist<-edwards.eberhardt
head(capt.hist)


## run default closed capture-recapture model (model Mb)
## model codes for type of CMR model, not covariate/parameter structure
Mb<-mark(data=capt.hist, model="Closed", brief=T)
# ignore all the automatic output

##### look at summary output

summary(Mb)
str(Mb)

#######there's a lot here that we don't really care about
####### easier: just access pieces of output you're interested in

##estimates on the link scale
Mb$results$beta
##detection probabilities (p and c) are estimated on the logit scale
##f0 (number of individuals never detected) is estimated on the log scale

##estimates on the real scale
Mb$results$real
##Remember: f0 is animals never detected, NOT total abundance

##note that p (probability of first capture) and c (probability of recapture)
##are almost identical; that suggests we may not need this additional parameter
##and should fit a model with p=c instead - coming up


##derived estimate of total abundance (f0 + n observed)
Mb$results$derived


###################### fit alternative models

##create submodels for parameters p and c, both are constant but differ from each other
##has to be a list, with element "formula", which has to begin with a ~
p.dot<-list(formula=~1)
c.dot<-list(formula=~1)

##run model Mb again, using the specified structure for p and c
##both are combined into a list with elements p and c under model.parameters
Mb.2<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.dot, c=c.dot), 
          brief=T)

##check that estimates are the same as before
Mb$results$derived
Mb.2$results$derived




##model M0, p=c, constant over time
p.dot.shared<-list(formula=~1,share=TRUE)
M0<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.dot.shared), 
           brief=T)
M0$results$real
M0$results$derived

##model Mt, p=c, varying over time
p.time.shared<-list(formula=~time,share=TRUE)
Mt<-mark(data=capt.hist, model="Closed", model.parameters=list(p=p.time.shared), 
         brief=T)
Mt$results$real
Mt$results$derived


#####  Perform model selection

modlist<-collect.models()
modlist
##this collects all models in the workspace, including both (identical)
## versions of model Mb


modlist2<-collect.models(lx=c("M0", "Mb", "Mt"))
modlist2
##this only collects the models you specify under lx
##notation in the table: if c(), that means p=c
##                       if c(~1), that means p different from c

##abundance under best model
Mt$results$derived


##build mixture model with individual heterogeneity in detection

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


