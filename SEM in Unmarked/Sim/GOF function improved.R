###############################################################################
### goodness of fit ###########################################################

##compare saturated model (with all paths) against nested model 
##only applies if focal model is not saturated already

##Two options:
## (a) Based on p-value of coefficients representing missing paths - skip!
## (b) Based on the LL of the two models
##    LL of a SEM is the sum of the LL of the submodels

GOF_LL<-function(M.focal,M.saturated){ #M1 is focal model, M2 is saturated model
  
  ##Q: include detection parameters in parm count? SEM structure only on 
  ##   state component; but shouldn't matter
  ## Include check for identical structure in detection model!!
  
  ##get LL, k (#params) of M1
  ##get model type to determine how to obtain LL
  mod.type1<-lapply(M.focal, class)
  mod.type2<-lapply(M.saturated, class)
  #remove names from lists as they will cause problems if not identical
  names(mod.type1)<-names(mod.type2)<-NULL
  
  ##check that lmerTest is used
  if(any(unlist(c(mod.type1, mod.type2)) == 'LmerMod')) stop('Function requires linear mixed models to be run
                                     with lmerTest, an extension package to lme4. Please
                                     install lmerTest and re-run your linear mixed model.')
  
 if(!identical(mod.type1, mod.type2)){
   stop('Order/types of submodels must be identical in M.focal 
        and M.saturated')
 }
  
  ll1<-k1<-aicf<-NULL
  for (i in 1:length(M.focal)){
    if(grepl('unmarked', mod.type1[[i]][1])){
    ll1[i]<--1*M.focal[[i]]@negLogLike
    k1[i]<-attr(logLik(M.focal[[i]]), "df")
    aicf[i]<-M.focal[[i]]@AIC
    }
    if(mod.type1[[i]][1] %in% c('lm', 'glm')){
      ll1[i]<-logLik(M.focal[[i]])
      k1[i]<-attr(logLik(M.focal[[i]]), "df")
      aicf[i]<-AIC(M.focal[[i]])
    }
    if(mod.type1[[i]][1] %in% c('glmerMod', 'lmerModLmerTest')){
      ll1[i]<-logLik(M.focal[[i]])
      k1[i]<-attr(logLik(M.focal[[i]]), "df")
      aicf[i]<-AIC(M.focal[[i]])
    }
  }
  
  ##get LL, k (#params) of M2
  ll2<-k2<-NULL
  for (i in 1:length(M.saturated)){
    if(grepl('unmarked', mod.type2[[i]][1])){
      ll2[i]<--1*M.saturated[[i]]@negLogLike
      k2[i]<-attr(logLik(M.saturated[[i]]), "df")
    }
    if(mod.type2[[i]][1] %in% c('lm', 'glm')){
      ll2[i]<-logLik(M.saturated[[i]])
      k2[i]<-attr(logLik(M.saturated[[i]]), "df")
    }
    if(mod.type2[[i]][1] %in% c('glmerMod', 'lmerModLmerTest')){
      ll2[i]<-logLik(M.saturated[[i]])
      k2[i]<-attr(logLik(M.saturated[[i]]), "df")
    }
  }
  
  #k= sum of differences in nparm across submodels
  k<-sum(abs(k1-k2))
  chiLL<--2*(sum(ll1)-sum(ll2))
  p.out<-1-pchisq(chiLL, k)
  
  #get missing path coefficients
  mpc<-missing.path(M.focal, M.saturated)
  out.list=list(AIC.focal=sum(aicf), p.overall=p.out, df=k, X2=chiLL, 
                missing.paths = mpc )
  return(out.list)
}


###extract missing paths from saturated model
missing.path<-function(M.focal, M.saturated){ 
  
  mod.type1<-lapply(M.focal, class)
  mod.type2<-lapply(M.saturated, class)
  names(mod.type1)<-names(mod.type2)<-NULL
  
  ##check that lmerTest is used
  if(any(unlist(c(mod.type1, mod.type2)) == 'LmerMod')) stop('Function requires linear mixed models to be run
                                     with lmerTest, an extension package to lme4. Please
                                     install lmerTest and re-run your linear mixed model.')
  
  if(!identical(mod.type1, mod.type2)){
    stop('Order/types of submodels must be identical in M.focal 
        and M.saturated')
  }
  
  #get beta coefficients
  parms1<-parms2<-list()
  for (i in 1:length(M.focal)){
  if(grepl('unmarked', mod.type1[[i]][1])){
    #function to suppress print
    log <- capture.output({
      res.f <- summary(M.focal[[i]]@estimates)$state;#[,4];
      res.s <-summary(M.saturated[[i]]@estimates)$state;
    })
    parms1[[i]]<-res.f
    parms2[[i]]<-res.s
  }
    if(mod.type1[[i]][1] %in% c('lm', 'glm','glmerMod', 'lmerModLmerTest')){
      parms1[[i]]<-as.data.frame(summary(M.focal[[i]])$coefficients)
      parms2[[i]]<-as.data.frame(summary(M.saturated[[i]])$coefficients)
    }
  }
  

  ##compare to extract new paths in saturated model
  new.paths<-pcf<-list()
  for (i in 1:length(M.focal)){
  new.paths[[i]]<-rownames(parms2[[i]])[is.na(pmatch(rownames(parms2[[i]]), 
                                                     rownames(parms1[[i]])))]
  pcf[[i]]<-parms2[[i]][new.paths[[i]],]
  }
  
  names(pcf)<-paste('Submodel ', 1:2,sep='')
  return(pcf)
}



