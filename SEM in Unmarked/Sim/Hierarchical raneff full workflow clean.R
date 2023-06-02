################################################################################
#### hierarchical model w/  random effects #####################################

# this is a DS example
##model: x3 -> x2 -> Abundance <-x1
rm(list=ls())

library(unmarked)
library(lme4)
library(lmerTest)

setwd("D:/MargSalas/Scripts_MS/SEM in Unmarked/Sim")

source('Standardization R function unmarked GLMM.R')
source('GOF function improved.R')

##one random effect: points nested within regions
J<-100 #number of sampling locations
K<-1 #no temporal repeats
region<-as.factor(rep(1:10,10)) #random effect

#random effect region on counts
eps.r<-rnorm(length(levels(region)), 0, 0.5)
#random effect region on x2
eps.rx<-rnorm(length(levels(region)), 0, 0.5)

##detection params
dbreaks<-seq(0,10, 1)
sigma<-3
#mid-point approximation of detection function
p.d<-exp(-(0.5:9.5)^2/(2*sigma^2))#/sum(exp(-(0:10)^2/(2*sigma^2)))
p.overall<-sum(p.d)/10


set.seed(1234)
x1<-rbinom(J*K,1,0.5) #binary predictor

##x3 affects x2, but not (directly) the response; include random effect
x3<-rnorm(J*K, 0, 2) #
mu.x2<-0.5*x3 + eps.rx[as.numeric(region)]
x2<-rnorm(J*K, mu.x2, 1)


# coefficients  
b0<--1 #intercept
b1<-1 #coefficient beta for x1 
b2<--1 #coefficient beta for x2 
b3<-0.5  #coefficient beta for x3

#generate expected abundance, link scale 
p.link<-b0+b1*x1+b2*x2+eps.r[as.numeric(region)]

#convert to exp abundance/occurrence, natural scale
lam<-exp(p.link)

#compile covariates in a data frame; include all exogenous vars
siteCovs<-data.frame(x1=x1, x2=x2, x3=x3, region=region)#, rept=rept)

#generate abundances and occurrences
N<-rpois(J*K, lam)

#generate detection data 
obsDS<-matrix(0, J*K, length(p.d)) # ncol = distance bins
  
for (j in 1:(J*K)){
  if(N[j]==0)next
  #assign distance category
  dists<-table(sample(1:length(p.d), N[j], replace=TRUE))
  
  #generate detections
  obs<-rbinom(length(dists), dists, p.d[as.numeric(names(dists))])
  if(sum(obs)==0)next
  #skip to next site if no detections

  obsDS[j,][as.numeric(names(dists))]<-obs
}


##compile into unmarked frames
umfDS<-unmarkedFrameDS(obsDS, siteCovs = siteCovs, dist.breaks = dbreaks,
                       survey='line', unitsIn = 'm', tlength=rep(100, J*K))


################################################################################
### fit piecewise models #######################################################

#first component: abundance (DS) model
ds.mlu<-distsamp(~1~x1+x2+(1|region), data=umfDS, keyfun='halfnorm',
                output='abun', unitsOut = 'ha')

#second component: model on x2
lm.ml<-lmer(x2~x3+(1|region), siteCovs)

##bundle both paths into full piecewise SEM
##names not necessary but helpful for some output
M.focal<-list(lmm=lm.ml, DS=ds.mlu)


##get standardized coefficients, either for each component model
stand_beta(lm.ml)
stand_beta(ds.mlu)
##or for list of models making up full SEM
stand_beta_wrap(M.focal)

################################################################################
### next, get model with all missing paths #####################################

##Note: identifying missing paths is up to user
##Here I assume missing paths are:
## x1->x2 (in linear mixed model)
## x3 -> Abundance (in DS model)
## I assume x2 (an endogenous variable) cannot affect x1 (an exogenous variable)
## And I assume we don't care about relationships between x1 and x3 (both 
## exogenous variables)
## But these things may be different and functions may need adjusting in a real
## life analysis

##Build fully saturated model with all missing paths
lm.full<-lmer(x2~x3+x1+(1|region), siteCovs)
ds.mlu.full<-distsamp(~1~x1+x2+x3+(1|region), data=umfDS, keyfun='halfnorm',
                      output='abun', unitsOut = 'ha')

##Bundle into saturated piecewise SEM
M.sat<-list(lm.full, ds.mlu.full)

#test fit of focal model (by comparing to saturated model)
#using likelihood test described here: 
#https://jslefche.github.io/sem_book/local-estimation.html#extensions-to-generalized-mixed-effects-models

#includes missing path coefficient estimates
GOF_LL(M.focal, M.sat)
#returns:
#  AIC of focal model
#  p.overall: p-value for test focal against saturated model
#  df: differences in number of paramters between the two models
#  X2: Chi-square value used to determine p
#  coefficients of paths that were added in the saturated model

##if p-value>0.05, there is support for the focal model over the saturated model

#get coefficients, p-values of missing paths only
missing.path(M.focal, M.sat)

##missing path coefficients are very small, have very high p-values

##Of course, in this example this is to be expected as the data were generated
## under the focal model

