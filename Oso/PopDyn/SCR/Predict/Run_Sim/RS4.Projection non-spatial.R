
## -------------------------------------------------
##      Predict from simulation script RS (run1)
##              NON-SPATIAL PROJECTION
## ------------------------------------------------- 

#rm(list = ls())

samp2<-readRDS('TestOutputAGE.rds')
sampmat<-do.call(rbind, samp2)
dim(sampmat)

niter<-nrow(sampmat)

t.new<-10 #10 yrs into future
max.age<-30
age.cat<-3 # age categories to keep track off
Narray<-array(NA, c(age.cat, t.new, niter)) # ASP: For each iteration, project abundance per each category and year

##get starting values: abundance in each age class and age
##get survival, recruitment

##z, age at T=5 from model
z.which<-grep('z\\[', colnames(sampmat)) # ASP: index columns all z (sampmat matrix)
##ordered: all individuals for yr 1, then all for yr2, etc
z.est<-array(sampmat[,z.which], c(niter,Maug, Tt))
z.nosuper<-apply(z.est,1:2,sum)==0 # ASP: Across all years if that individual in that iteration was never detected (sum of z = 0), it was never part of the super pop

#age.which<-pmatch(paste('age[', 1:Maug, ', ', Tt,']', sep=''), colnames(sampmat))
age.which<-grep('age\\[', colnames(sampmat))
age.est<-array(sampmat[,age.which], c(niter,Maug, Tt))

length(which(age.est==0)) == length(which(z.est==0)) # ASP: # z.est is longer because it includes individuals that are not part of the superpop

age.est[z.est==0]<-0 # ASP: If is not alive it doesn't have an age (how age model is set up)


##get distribution of ages in each iteration at t=5
ttt<-t(apply(age.est[,,5], 1, table)) # ASP: At each iteration (rows, dim = 1), individuals per age (table of age est)
N.age<-matrix(NA, niter, 3) # ASP: Organize into age categories
for (i in 1:niter){
N.age[i,]<-c(sum(ttt[i,2:3]), sum(ttt[i,4:5]), sum(ttt[i,6:ncol(ttt)]))
}

Narray[,1,]<-t(N.age) # ASP: Fill N at each age class the first year for each iteration

##calculate per capita recruitment
R<-matrix(NA, nrow(sampmat), Tt-1)
for (t in 2:Tt){
  zwt<-paste('z[', 1:Maug,', ', t, ']', sep='' )
  zwtm<-paste('z[', 1:Maug,', ', t-1, ']', sep='' )
  
  for (nn in 1:nrow(sampmat)){
    R[nn,t-1]<-sum(ifelse((sampmat[nn,zwt]-sampmat[nn,zwtm])==1, 1, 0))
  }
}
N.which<-grep('N', colnames(sampmat))[1:(Tt-1)]
pcrmat<-matrix(NA, nrow(sampmat), 4)
for (ite in 1:nrow(sampmat)){
  pcrmat[ite,]<-R[ite,]/sampmat[ite,N.which[1:4]]
}
pcr<-apply(pcrmat,1,mean) # ASP: Mean per capita recruitment per iteration

# ASP: Projection, for each iteration

 for(ite in 1:niter){
ageN<-matrix(NA, t.new, max.age) # ASP: Number of individuals of each age per year
###make vectors of survival
phi.vec<-c(rep(sampmat[ite, 'phi.cub'], 2),
           rep(sampmat[ite, 'phi.sub'], 2),
           rep(sampmat[ite, 'phi.ad'], max.age-4)) # ASP: Survival corresponding to each age

##starting population, ages 
ageN[1,]<-c(ttt[ite,2:ncol(ttt)], rep(0, max.age-ncol(ttt)+1)) # ASP: The first row, year 1

##move forward
for (t in 2:t.new){
  ##survivors
  S<-rbinom(length(ageN[t-1,]), ageN[t-1,], phi.vec) #ASP: Simulate 30 values (ages). For the number of individuals of each age, how many survived according to phi vec (their age specific survival)
  if(S[max.age]>0) stop('Up max.age!!')
  
  ##recruits (N at t-1 times pcr) - SHOULD BE ADULTS ONLY!!!
 # B<-rpois(1, pcr[ite]*sum(Narray[1:3,t-1,ite]))
  #binomial, assuming there are always 500 available
  B<-rbinom(1, 500, (pcr[ite]*sum(Narray[1:3,t-1,ite]))/500) # ASP: Number of recruits each year (pcr*nºAdults previous year). 0 age category
  
  ageN[t,]<-c(B, S[1:(max.age-1)]) # ASP: The recruits pass to be the cubs of year 1 at time t -> al rest AGING

  Narray[,t,ite]<-c(sum(ageN[t,1:2]),sum(ageN[t,3:4]), sum(ageN[t,5:max.age]) ) # ASP: Sum the N of the different ages (into age categories)
  
}
 }

Nall<-apply(Narray, 2:3, sum)

plot(1:5, apply(Nall[1:5,], 1,mean), type='l', ylim=c(0,200))
for (ite in 1:niter){
  points(1:5, Nall[1:5,ite], type='l', col='lightgrey')
}
points(1:10,apply(Nall, 1,mean), type='l')
