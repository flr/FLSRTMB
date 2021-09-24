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
#' @param object Input FLSR object.
#' @param s steepness parameter of SRR (fixed or prior mean)    
#' @param spr0 unfished spawning biomass per recruit from FLCore::spr0(FLStock) 
#' @param s.est option to estimate steepness
#' @param s.logitsd prior sd for logit(s), default is 1.3 (flat) if s.est = TRUE 
#' @param plim determines the minimum break point of the hockey-stick as ratio blim/b0
#' @param nyears yearMeans from the tail used to compute a,b from the reference spr0 (default all years)
#' @param report.sR0 option to report s and R0 instead of a,b
#' @param inits option to specify initial values of log(r0), log(SigR) and logit(s)
#' @param lower option to specify lower bounds of log(r0), log(SigR) and logit(s) 
#' @param upper option to specify upper bounds of log(r0), log(SigR) and logit(s)
#' @param upper option to specify upper bounds of log(r0), log(SigR) and logit(s)
#' @param SDreport option to converge hessian and get vcov
#'
#' @return A list containing elements 'FLSR', of class *FLSR*
#' @export

srrTMB <- function(object, spr0, s=NULL, s.est=TRUE,s.logitsd=1.3,plim=NULL,nyears=NULL,report.sR0=FALSE,inits=NULL, lower=NULL, upper=NULL,
  SDreport=TRUE) {
  
  if(is.null(nyears)) nyears = dim(ssb(object))[2]
  if(is.null(plim)){
    if(is.null(s)& s.est){ plim=0.01} else {plim=0.2}
  }
  
  if(is.null(s)& s.est){s=0.6} # central value
  if(is.null(s)& !s.est){s=0.8}
  
  # IDENTIFY model
  model <- SRModelName(model(object))

  if(length(spr0)>1){
  #if(length(ssb(object))+1!=length(spr0)) stop("The spr0 vector must correct to non NA values of ssb(stock)")
  spr0.yr =  trim(spr0,year=dims(ssb(object))[["minyear"]]:dims(ssb(object))[["maxyear"]])
  } else {
  spr0.yr= ssb(object)
  spr0.yr[] = spr0  
  }
  
  spr0.yr = c(spr0.yr)
  spr0ref = mean(spr0.yr[(length(spr0.yr)-nyears+1):(length(spr0.yr))])
    
  # GET rec, ssb
  rec <- c(rec(object))
  ssb <- c(ssb(object))

  # SET init and bounds
  if(is.null(inits))
    inits <- c(mean(log(rec)), log(0.4),to_logits(s))
  if(is.null(lower))
    lower <- c(min(log(rec)), log(0.05),-20)
  if(is.null(upper))
    upper <- c(max(log(rec * 20)), log(1.5),20)

  # SET TMB input
  inp <- list(
    # data
    Data = list(ssb = ssb, rec = rec,prior_s = c(to_logits(s),s.logitsd),
    spr0y = spr0.yr,spr0=spr0ref,plim=plim, nyears=length(ssb),
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
    control=list("trace"=1, "eval.max"=1e4, "iter.max"=1e4),
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
  
  if(model!="segreg"){
  attr(object,"SV") = data.frame(s=Report$s,sigmaR=Report$sigR,R0=Report$r0)
  } else{
    attr(object,"SV") = data.frame(sigmaR=Report$sigR,R0=Report$r0)
    
  }
  
  if(SDreport & report.sR0==FALSE){
    object@vcov = matrix(SD$cov,nrow=2,dimnames = list(c("a","b"),c("a","b")))
  }

  return(object)
}
