################################################################################
##### standardized coefficients and R2 for unmarked and ubms models #############

##wrapper function so stand_beta can be applied to model lists (ie, full paths)

stand_beta_wrap<-function(models){
  if(!is.list(models)){
    outlist<-stand_beta(models)
  } else {
    outlist<-lapply(models, stand_beta)
  }
  return(outlist)
}

##Based only on Latent Theoretic approach (Nakagawa and Schieltze 2012)

stand_beta<-function(model){  #if FALSE, uses coefs() function from piecewiseSEM

  #mark whether this is detection corrected model
  det.cor=FALSE
  if(try(attr(model@class, 'package'), silent=TRUE) == 'unmarked') det.cor=TRUE

  #check that lmerTest was used
  if(class(model)[1] == 'lmerMod') stop('Function requires linear mixed models to be run
                                     with lmerTest, an extension package to lme4. Please
                                     install lmerTest and re-run your linear mixed model.')
  
  if (det.cor){ #if unmarked model
  
  if(!class(model) %in% c("unmarkedFitDS", "unmarkedFitOccu","unmarkedFitPCount",
                          "ubmsFitOccu","ubmsFitPcount","ubmsFitDistsamp")){
    stop('In unmarked model: Currently, function only works for single season occupancy,
         N-mixture or distance sampling models.')
  }
  
  if (class(model) == "unmarkedFitOccu"){
    warning('In unmarked model: Function only valid for logit link; 
            do not use if using linkPsi=cloglog in occu()')
  }
  
  if (class(model) == "unmarkedFitPCount"){
    warning('In unmarked model: Function only valid for Poisson mixture; 
            do not use if using mixture=NB or ZIP in pcount()')
  }
    
    ##determine family of state variable
    fam<-ifelse(grepl('Occu', class(model)), 'binomial', 'poisson')
    
  ##check if state model has random effect, get raneff variance 
    raneff<-names(model@estimates@estimates$state@randomVarInfo$cnms)
    var.ran<-0
   if (length(raneff)>0){
     
    if(model@estimates@estimates$state@randomVarInfo$cnms != "(Intercept)") {stop(
    'In unmarked model: Model has random slope. Function does not support random slopes at the moment.')}
     
     #sum of variances of random intercepts
    var.ran<-sum(exp(model@estimates@estimates$state@randomVarInfo$estimates)^2)

   }

  ##model predictions, on link scale
  pred<-unmarked::predict(model, 'state', backTransform=FALSE)[,1]
  
  ##calculate variation due to covariates
  var.sys<-var(pred)
  
  ##calculate family specific variation
  intercept<-coef(model)[1]
  var.fam<-ifelse(fam == 'binomial', pi^2/3, log(1/exp(intercept) + 1))

  ##get beta coefficients, make data frame akin to GLMM output
  #function to suppress print
  log <- capture.output({
    res <- summary(model@estimates)$state;#[,4];
  })
  
  beta.raw<-res
  
  #get sd x for non-categorical predictors
  #get predictors used in state model
  keep<-rownames(beta.raw)[-1]
  
  pred.var<-as.matrix(model@data@siteCovs[,keep], ncol=length(keep))
  pred.factors<-apply(pred.var, 2, is.factor)
  sd_x<-apply(pred.var, 2, sd)
  sd_x[pred.factors]<-NA #set to NA for categorical variables which have no variance
  
  } else { # end detection corrected part
    
    ############################################################################
    ##determine family of state variable
    if (class(model)[1] %in% c("lmerModLmerTest", 'lm')) {fam<-'gaussian'}else{

      if ('glm' %in% class(model)){
           fam<-model$family$family
        }
         
         if (class(model)[1] == 'glmerMod'){
          fam<-model@call$family
        }
      
    }
  
    #check that model uses valid family
    if(!fam %in% c('gaussian', 'poisson', 'binomial')) stop('In (G)L(M)M: Currently function is
                                                            only valid for Normal, Poisson
                                                            or Binomial response variables.')
    
    #check for random slopes (error) and random intercepts
    raneff<-NULL #default to no random effect
    var.ran<-0
    
    if (class(model)[1] %in% c("glmerMod","lmerModLmerTest")){
      #first check for random slopes
      r.e<-as.data.frame(VarCorr(model))
      #remove residuals for lmer
      r.e<-r.e[r.e$grp != 'Residual',]
      if(any( r.e$var1 != '(Intercept)')) stop ('In (G)L(M)M: Model contains random slopes; function 
                                                does not support random slopes at the
                                                moment.')
      raneff<-r.e$grp
      var.ran<-sum(r.e$sdcor^2) # sum all random intercept variances
    }
    
    #get beta coefficients - works for all GLMM type models
      beta.raw<-as.data.frame(summary(model)$coefficients)
      
    # for non-gaussian, get other variance components
      if (fam !='gaussian'){
      
        if(class(model)[1] == 'glmerMod' ){
          pred<-predict(model, type='link',re.form =NA)}
        if('glm' %in% class(model) ){
          pred<-predict(model, type='link')
        }
      ##calculate variation due to covariates
      var.sys<-var(pred)
      
      ##calculate family specific variation
      intercept<-beta.raw['(Intercept)', 'Estimate']
      var.fam<-ifelse(fam == 'binomial', pi^2/3, log(1/exp(intercept) + 1))
      } # end non-gaussian variance components
      
      
      #get sd x for non-categorical predictors
      pred.var<-as.matrix(model.frame(model)[,rownames(beta.raw)[-1]], 
                          ncol=length(rownames(beta.raw)[-1]))
      pred.factors<-apply(pred.var, 2, is.factor)
      sd_x<-apply(pred.var, 2, sd)
      sd_x[pred.factors]<-NA #set to NA for categorical variables which have no variance
      
  } #end GLMM part
  ##############################################################################
  
  #SD in Y on link scale
  if (fam != 'gaussian'){
  sd_y<-sqrt(var.sys+var.fam+var.ran)} else {
    sd_y<-sd(model.frame(model)[,1])#get data from model
  }

  #standardized betas
  beta.star<-data.frame(beta.raw[-1,], Std.Estimate =c(beta.raw[-1,1]*sd_x/sd_y),
                        check.names = FALSE)
  
  return(beta.star)

  } 




