########### WFC124 - Lab week 9 - Open population CMR

rm(list=ls())

library(RMark)

########## 
### read in dipper data, this is without formatting for Rmark, just to look at it
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Stats/Mark")
dip<-read.csv("Dipper.data.csv")

head(dip)

###CHALLENGE###
#Data summaries
#1. How many years (primary occasions) are there?
#2. How many birds were captured?
#3. How many birds were captured more than once?
#4. How many times were birds captured on average?

##number primary occasions (ie, years)
ncol(dip)

##number birds captured across all years
nrow(dip)

##how many birds captured > once
ncap<-apply(dip,1,sum)
sum(ncap>1)

##how many times captured, on average
mean(ncap)

#####################
## load and look at data already formatted for RMark
data(dipper)
head(dipper)

##use table() function to see how many males and females were captured
table(dipper$sex)

#####################
##tell RMark that we want to use "sex" as a grouping variable (stratification)
dipproc <- process.data(dipper, groups="sex")
##going forward, we will use the object "dipproc" in the models

## set up model structure for both parameters survival (phi) and detection (p)
#intercept only
null <- list(formula = ~1)

#separate parameters for males and females
#note that the word after ~ needs to match the name you used to build groups in 
#the process.data() function
sex< - list(formula = ~sex)

#changing parameters over time (ie, capture occasions)
time <- list(formula = ~time)


###############
###CHALLENGE###
#Generate the following 8 models:
#phi0.p0 - null model
#phi.sex.p0 - sex as a covariate on survival, null detection
#phi.time.p0 - time as a covariate on survival, null detection
#phi0.p.sex - null survival, sex as a covariate on detection
#phi.sex.p.sex - sex as covariate on survival AND detection
#phi.time.p.sex - time as covariate on survival, sex on detection
#phi0.p.time - null survival, time on detection
#phi.sex.p.time - sex on survival, time on detection


## Run possible models with p and phi varying by sex and time
## We do not run p~time phi~time because some parameters in that model are not identifiable
## We do not need to specify "model" because CJS is the default.  Modeling p and phi without covariates
##is also the default, and we don't have a c parameter to worry about (as we did in closed models),
##so you can specify all parameters explicitly like this: 

phi0.p0 <- mark(dipproc, model.parameters = 
                  list(Phi = null, p = null), brief = T)
phi.sex.p0 <- mark(dipproc, model.parameters = 
                     list(Phi = sex, p = null), brief = T)
phi.time.p0 <- mark(dipproc, model.parameters = 
                      list(Phi = time, p = null), brief = T)
phi0.p.sex <- mark(dipproc, model.parameters = 
                     list(Phi = null, p = sex), brief = T)
phi.sex.p.sex <- mark(dipproc, model.parameters = 
                        list(Phi = sex, p = sex), brief = T)
phi.time.p.sex <- mark(dipproc, model.parameters = 
                         list(Phi = time, p = sex), brief = T)
phi0.p.time <- mark(dipproc, model.parameters = 
                      list(Phi = null, p = time), brief = T)
phi.sex.p.time <- mark(dipproc, model.parameters = 
                         list(Phi = sex, p = time), brief = T)

#################
##build model selection table and look at results of top models
modlist<-collect.models() 
modlist

############### 
##estimates on the link scale (logit)
phi0.p0$results$beta

##compare estimates on the real (or natural) scale -ie, already backtransformed
phi0.p0$results$real

phi0.p.sex$results$real
## in the null model output, ignore the labeling of gFemale and t1/t2
## all parameters are constant and thus refer to males and females
## and to all time steps in the data

##for comparison, look at a model that has a sex effect
phi.sex.p0$results$real
##here, you see an estimate of Phi for males and one for females
##which makes sense given the model structure (phi~sex)

###############
##some formatting to add in a temporal covariate, indicating whether each of the seven years
##was a flood year (Flood = 1; true for year 2 and 3)

##this function makes a design matrix so we can add covariates that varied over the years of the study
dipper.ddl <- make.design.data(dipproc)

##look at the desing matrix
class(dipper.ddl)
head(dipper.ddl$Phi)
head(dipper.ddl$p)

nrow(dipper)
dim(dipper.ddl$Phi)[1]

# HERE I don't understand, if there are initially nrow(dipper) = 294 individuals, 
# why the phi and p matrix have only 42 rows, how do I add an individual covariate?




##We add the covariate Flood (we could call it whatever we want) and set all of its values to 0
##We do that for both phi and p
dipper.ddl$Phi$Flood = 0
dipper.ddl$p$Flood = 0

##look at matrix to see that it now has a column called Flood
head(dipper.ddl$Phi)

##The we change the Flood values to 1 for time = 2 and time = 3 --> flood years in the study
dipper.ddl$Phi$Flood[dipper.ddl$Phi$time == 2 | dipper.ddl$Phi$time == 3] = 1
dipper.ddl$p$Flood[dipper.ddl$p$time == 2 | dipper.ddl$p$time == 3] = 1

##look at matrix to see that the Flood column now has a value of 1 for years
## 2 and 3
head(dipper.ddl$p)


############# 
##run additional models, looking at effect of flood on survival 
flood <- list(formula = ~Flood) #model structure, effect of flood
sex.flood <- list(formula = ~Flood + sex) #effect of both sex and flood

###CHALLENGE###
#Create 3 additional models incorporating the effect of flood on survival
#phi.flood.p0 - flood on survival, null on detection
#phi.flood.p.sex - flood on survival, sex on detection
#phi.sex.flood.p0 - flood+sex on survival, null on detection
#NOTE: include BOTH dipproc and dipper.ddl - eg mark(dipproc, dipper.ddl, model.parameters.......)

phi.flood.p0 <- mark(dipproc, dipper.ddl, model.parameters = list(Phi = flood, p = null), brief = T)

phi.flood.p.sex <- mark(dipproc, dipper.ddl, model.parameters = list(Phi = flood, p = sex), brief = T)

phi.sex.flood.p0 <- mark(dipproc, dipper.ddl, model.parameters = list(Phi = sex.flood, p = null), brief = T)

##compare all models run so far
collect.models()

topresults <- phi.flood.p0$results$real
topresults
#in topresults, first row is survival in non-flood years, second row is survival in flood years 
#unfortunately, "t1" and "t2" notation is not very helpful
#t2 = parameter where covariate "Flood" ==1 (flood year)
#t1 = parameter where covariate "Flood" ==0 (no flood)

#############
### make a plot of survival probabilities for flood and non-flood years
##pull out survival estimates from model results
phi<-topresults[1:2, ]
phi

##make barplot with axis labels, category labels
##we set ylim for next step
##we assign the plot to an object for the next step, also
bp<-barplot(phi$estimate, xlab = "Year", 
            ylab = "Survival probability", 
            names.arg = c("Non-flood", "Flood"), 
            ylim = c(0, 0.7))


###########
##adding error bars to the barplot using the arrows() function

##define lower and upper limit of error bar as the estimate minus and plus the SE
lower <- phi$estimate - phi$se 
upper <- phi$estimate + phi$se

##NOTE: if you're working on a Mac, your output may not contain the "se" column
## in that case, use (lower and upper confidence interval levels)
## for CIs to show, you may have to expand your y axis in the plot command
## to range from 0 to 1
#lower<-phi$lcl
#upper<-phi$ucl

## now, use arrows function to add error bars
arrows(x0 = bp, x1 = bp, y0 = lower, 
       y1 = upper, code = 3, angle = 90 )









