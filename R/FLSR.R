# FLSR.R - DESC
# /FLSR.R

# Copyright Henning Winker (JRC) & Iago MOSQUEIRA (WMR), 2021
# Authors:  Henning Winker (JRC) <henning.winker@ec.europa.eu>
#           Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# srrTMB {{{

#' Fits Stock Recruitment Relationships (SRR) in TBM
#'
#' @param object Input FLSR = as.FLSR(stock,model) object with current model options
#' \itemize{
  #'   \item bevholtSV   
  #'   \item rickerSV
  #'   \item segreg
  #'   \item geomean
  #' }    
#' @param s steepness parameter of SRR (fixed or prior mean)    
#' @param spr0 unfished spawning biomass per recruit from FLCore::spr0(FLStock) 
#' @param s.est option to estimate steepness
#' @param s.logitsd prior sd for logit(s), default is 1.4 (flat) if s.est = TRUE 
#' @param lsrp lower bound of plausible spawning ratio potential SRP
#' @param usrp upper bound of plausible spawning ratio potential SRP
#' @param plim depreciated plim = usrp
#' @param pmax depreciated pmax = lsrp
#' @param nyears yearMeans from the tail used to compute a,b from the reference spr0 (default all years)
#' @param report.sR0 option to report s and R0 instead of a,b
#' @param inits option to specify initial values of log(r0), log(SigR) and logit(s)
#' @param lower option to specify lower bounds of log(r0), log(SigR) and logit(s) 
#' @param upper option to specify upper bounds of log(r0), log(SigR) and logit(s)
#' @param upper option to specify upper bounds of log(r0), log(SigR) and logit(s)
#' @param SDreport option to converge hessian and get vcov
#' @param verbose if TRUE, it shows tracing
#' @return A list containing elements 'FLSR', of class *FLSR*
#' @export
#' @examples
#' data(ple4)
#' gm <- srrTMB(as.FLSR(ple4,model=geomean),spr0=mean(spr0y(ple4)))
#' bh <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
#' ri <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=spr0y(ple4))
#' hs <- srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),lsrp=0.05,usrp=0.2)
#' srs = FLSRs(gm=gm,bh=bh,ri=ri,hs=hs) # combine
#' plotsrs(srs) 
#' plotsrts(srs)  # relative
#' plotsrs(srs[2:4],b0=TRUE) # through to B0
#' plotsrs(srs[2:4],b0=TRUE,rel=TRUE)  # relative
#' gm@SV # estimates
#' do.call(rbind,lapply(srs,AIC))

srrTMB <- function(object, spr0="missing", s=NULL, s.est=TRUE,s.logitsd=10,lsrp=0.01,usrp=0.35,plim=lsrp,pmax=usrp,nyears=NULL,report.sR0=FALSE,inits=NULL, lower=NULL, upper=NULL,
  SDreport=TRUE,verbose=FALSE) {
  
  silent = ifelse(verbose,1,0)
  
  if(is.null(nyears)) nyears = dim(ssb(object))[2]
  if(is.null(plim)){
    if(is.null(s)& s.est){ plim=0.01} else {plim=0.1}
  }
  
  # IDENTIFY model
  model <- SRModelName(model(object))
  
  if(model=="mean"){
    gmB0 = TRUE
    if(missing(spr0)){
      spr0 = 1
      gmB0 = FALSE}
    if(verbose) cat("spr0 missing for computing B0","\n")
  }
  if(!model=="mean"){
    if(missing(spr0)){
      stop(paste("Required to specify spr0 for model",model))}}
  
  # check 
  if(model%in%c("bevholt","ricker"))
    stop(paste0("Please use ",model,"SV instead"))
  
  if(!model%in%c("mean","segreg","bevholtSV","rickerSV"))
     stop(paste("S-R model:",model,"is not (yet) defined in FLSRTMB"))
  
  if(length(spr0)>1){
    #if(length(ssb(object))+1!=length(spr0)) stop("The spr0 vector must correct to non NA values of ssb(stock)")
    spr0.yr =  trim(spr0,year=dims(ssb(object))[["minyear"]]:dims(ssb(object))[["maxyear"]])
  } else {
    spr0.yr= ssb(object)
    spr0.yr[] = spr0  
  }
  
  spr0.yr = c(spr0.yr)
  spr0ref = mean(spr0.yr[(length(spr0.yr)-nyears+1):(length(spr0.yr))])
  
  
  
  if(model=="segreg"){ # Adjust dynamically
  ll =plim/pmax
  ul = 1
  #s = an(quantile(c(ll,ul),0.75))
  srp = an(quantile(an((object@ssb/object@rec)/spr0ref),c(0.5)))
  srp = max(min(srp,0.9*pmax,srp),plim*1.1)
  s = 1/(srp/plim)
  }
  if(model=="rickerSV"){
   ll = 0.2
   ul = 20
   s = mean(c(ll,ul))
  } 
  if(model=="bevholtSV"){
  ll=0.2  
  ul = 1
  if(is.null(s)& s.est){s=mean(c(ll,ul))} # central value for s = 0.2-1.0
  if(is.null(s)& !s.est){s=0.7}
  }
  
 
  # GET rec, ssb
  rec <- c(rec(object))
  ssb <- c(ssb(object))
  
  
  # Simple geomean option with time-varying spr0y
  if(model=="mean"){
  if(length(unique(spr0.yr))>1){ 
     gmfit = lm(log(rec)~log(spr0.yr))} else {
       gmfit = lm(log(rec)~1)
     } 
  
  fitted(object) <- an(exp(predict(gmfit,data.frame(spr0.yr))))
  residuals(object) <- log(rec(object)) - log(fitted(object))
  params(object) <- FLPar(a= an(exp(predict(gmfit,data.frame(spr0.yr= spr0ref)))))
 
  # Add loglik manually (to learn)
  p = length(gmfit$coefficients)
  N = length(fitted(object))
  sigma <- summary(gmfit)$sigma*sqrt((N-p)/N) # correct for p
  attr(object@logLik,"df") = p+1
  attr(object@logLik,"nobs") = length(fitted(object)) 
  object@logLik[] =  sum(dnorm(0, mean=residuals(gmfit), sd=sigma, log=TRUE))
  #check logLik(gmfit)
  object@vcov = matrix(0,nrow=1,dimnames = list(c("a"),c("a")))
  
  # AR1 rho
  rho = stats::cor(residuals(object)[,-N],residuals(object)[,-1])
  R0 = an(exp(predict(gmfit,data.frame(spr0.yr= spr0ref))))
  B0 = ifelse(gmB0,R0*spr0ref,NA)
  attr(object,"SV") = data.frame(s=NA,sigmaR=summary(gmfit)$sigma,R0=R0,rho=rho,B0=B0)
  }
  
  # TMB models
  if(!model=="mean"){
    
  # Set r0 init
  r0init= data.frame(rec=rec,ssb=ssb)       
  r0init=median(r0init[quantile(r0init$ssb,0.5)>r0init$ssb,]$rec)
  
  
  # SET init and bounds
  if(model=="segreg"){
    if(is.null(inits)) inits <- c(an(quantile(log(rec),0.4)), log(0.3),to_logits(min(ll*5,0.9),ll=ll))
    #if(is.null(inits)) inits <- c(log(r0init), log(0.4),to_logits(s,lim=lim))
  } 
  
  if(is.null(inits)) inits <- c(log(r0init), log(0.3),to_logits(s,ll=ll,ul=ul))
  
  
  
  if(is.null(lower))
    lower <- c(min(log(rec)), log(0.05),-100)
  if(is.null(upper))
    upper <- c(max(log(rec * 20)), log(1.5),100)

  # SET TMB input
  inp <- list(
    # data
    Data = list(ssb = ssb, rec = rec,prior_s = c(to_logits(s,ll,ul),s.logitsd),
    spr0y = spr0.yr,spr0=spr0ref,plim=plim, nyears=length(ssb),slim=ll,smax=ul,
    # model
    Rmodel = which(model==c("bevholtSV","rickerSV","segreg"))-1),
    # inits
    Params = list(log_r0 = inits[1], log_sigR = inits[2],logit_s=inits[3]),
    # bounds
    lower=lower, upper=upper,
    #
    ReportSD = SDreport
  )
  
  # Compile TMB inputs 
  Map = list()
  # Turn off steepness estimation
  if(!s.est) Map[["logit_s"]] = factor( NA ) 

  # CREATE TMB object
  Obj <- TMB::MakeADFun(data = inp$Data, parameters = inp$Params,map=Map,
    DLL = "FLSRTMB", silent = TRUE)

 
  
  Opt <- stats::nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
    control=list("trace"=silent, "eval.max"=1e4, "iter.max"=1e4),
    lower=inp$lower, upper=inp$upper)

  Opt[["diagnostics"]] = data.frame(Est=Opt$par,
    final_gradient=Obj$gr(Opt$par))
  
  Report <- Obj$report()
  
  if(SDreport) {
    SD <- try(TMB::sdreport(Obj))
  }

  # LOAD output in FLSR

  # DEBUG HACK
  model(object) <- switch(model, bevholtSV=bevholt, rickerSV=ricker,segreg=segreg)

  fitted(object) <- c(Report$rec_hat)
  residuals(object) <- log(rec(object)) - log(fitted(object))
  
  params(object) <- FLPar(a=Report$a, b=Report$b)
  if(report.sR0) params(object) <- FLPar(s=Report$s, R0=Report$r0)
  
  
  # Add full loglik
  p = 2
  N = length(residuals(object))
  sigma <- Report$sigR*sqrt((N-p)/N) # correct for p
  attr(object@logLik,"df") = p+1
  attr(object@logLik,"nobs") = length(fitted(object)) 
  object@logLik[] =  sum(dnorm(0, mean=residuals(object), sd=sigma, log=TRUE))
  
  # AR1 rho
  rho = stats::cor(residuals(object)[,-N],residuals(object)[,-1])
  
  
  if(model!="segreg"){
  attr(object,"SV") = data.frame(s=Report$s,sigmaR=Report$sigR,R0=Report$r0,rho=rho,B0=Report$r0*spr0ref)
  } else{
    attr(object,"SV") = data.frame(s=NA,sigmaR=Report$sigR,R0=Report$r0,rho=rho,B0=Report$r0*spr0ref,BlimB0=round(plim*1/Report$s,4))
    
  }
 
  
  if(SDreport & report.sR0==FALSE){
    object@vcov = matrix(SD$cov,nrow=2,dimnames = list(c("a","b"),c("a","b")))
  }
  } # End TMB model

  return(object)
}
