# FLSR.R - DESC
# /FLSR.R

# Copyright Henning Winker (JRC) & Iago MOSQUEIRA (WMR), 2021
# Authors:  Henning Winker (JRC) <henning.winker@ec.europa.eu>
#           Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# srrTMB {{{

#' Fits Stock Recruitment Relationships (SRR) in TMB
#'
#' @param object Input FLSR = as.FLSR(stock,model) object with current model options
#' \itemize{
#'   \item bevholtSV   
#'   \item bevholtDa
#'   \item rickerSV
#'   \item segreg
#'   \item geomean
#' }    
#' @param s steepness parameter of SRR (fixed or prior mean)    
#' @param spr0 unfished spawning biomass per recruit from FLCore::spr0(FLStock) 
#' @param s.est option to estimate steepness
#' @param s.logitsd prior sd for logit(s), default is 1.4 (flat) if s.est = TRUE
#' @param r0.pr option to condition models on r0 priors (NULL = geomean)
#' @param lplim lower bound of spawning ratio potential SRP, default 0.0001
#' @param uplim upper bound of plausible spawning ratio potential SRP , default 0.3
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
#' hs <- srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),lplim=0.05,uplim=0.2)
#' srs = FLSRs(gm=gm,bh=bh,ri=ri,hs=hs) # combine
#' plotsrs(srs) 
#' plotsrts(srs)  # relative
#' plotsrs(srs[2:4],b0=TRUE) # through to B0
#' plotsrs(srs[2:4],b0=TRUE,rel=TRUE)  # relative
#' gm@SV # estimates
#' do.call(rbind,lapply(srs,AIC))

setGeneric("srrTMB", function(object, ...) standardGeneric("srrTMB"))

#' @rdname srrTMB

setMethod("srrTMB", signature(object="FLSRs"),
  function(object, ...) {
    lapply(object, srrTMB, ...)
  }
)

#' @rdname srrTMB

setMethod("srrTMB", signature(object="FLSR"),
  function(object, spr0="missing",
  s=NULL, s.est=TRUE, s.logitsd=20, r0.pr="missing",
  lplim=0.001, uplim=0.3, plim=lplim, pmax=uplim,
  nyears=NULL, report.sR0=FALSE, inits=NULL,
  lower=NULL, upper=NULL, SDreport=TRUE,verbose=FALSE) {
  
  d.type = c("None")
    
  silent = ifelse(verbose,1,0)
  
  s.inp = s

    
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
    if(verbose) warning("spr0 missing for computing B0","\n")
  }
  if(!model=="mean"){
    if(missing(spr0)){
      stop(paste("Required to specify spr0 for model", model))}}
  
  # check 
  if(model%in%c("bevholt","ricker"))
    model <- paste0(model, "SV")
 
  if(!model%in%c("mean","segreg","bevholtSV","rickerSV","bevholtDa"))
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
    mu = mean(c(plim,pmax))
    s = 1/(mu/plim)
  }
  if(model=="rickerSV"){
    ll = 0.2
    ul = 20
    s = mean(c(ll,ul))
  } 
  if(model%in%c("bevholtSV","bevholtDa")){
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
    #if(length(unique(spr0.yr))>1){ 
    #gmfit = lm(log(rec)~log(spr0.yr))} else {
    gmfit = lm(log(rec)~1)
    #} 
    
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
    r0init=median(r0init[quantile(r0init$ssb,0.6)>r0init$ssb,]$rec)
    
    
    # SET init and bounds
    if(model=="segreg"){
      srp = an(quantile(an((object@ssb/object@rec)/spr0ref),c(0.6)))
      srp = max(min(srp,0.9*pmax,srp),plim*1.1)
      if(is.null(s.inp)){
        sinit = 0.99 #1/1.1#1/(srp/plim)
      } else {
        sinit =  s.inp
      }
      if(is.null(inits)) inits <- c(an(quantile(log(rec),0.4)), log(0.3),to_logits(sinit,ll=ll))
      #if(is.null(inits)) inits <- c(log(r0init), log(0.4),to_logits(s,lim=lim))
    } 
    # ><> depensation
    if(model=="bevholtDa") d.type = "A"
    
    if(is.null(inits)) inits <- c(log(r0init), log(0.3),to_logits(s,ll=ll,ul=ul))
    
    if(missing(r0.pr)){
      prior_r0 = c(1,100,0)
    } else{
      if(is.null(r0.pr)){
        prior_r0= c(exp(mean(log(rec),na.rm=TRUE)),0.2)
      } else {
    prior_r0=r0.pr}
    prior_r0[1] = log(prior_r0[1]) 
    prior_r0[3] = 1  
    }
    
  
    
    if(is.null(lower))
      lower <- c(min(log(rec)), log(0.05),-100,-10)
    if(is.null(upper))
      upper <- c(max(log(rec * 20)), log(1.5),100,10)
    # preliminary HACK for depensation model  
     Rmod = ifelse(model=="bevholtDa","bevholtSV",model)
     # SET TMB input
    inp <- list(
      # data
      Data = list(ssb = ssb, rec = rec,prior_s = c(to_logits(s,ll,ul),s.logitsd), prior_r0 = prior_r0,
                  spr0y = spr0.yr,spr0=spr0ref,plim=plim, nyears=length(ssb),slim=ll,smax=ul,
                  
                  # model
                  Rmodel = which(Rmod==c("bevholtSV","rickerSV","segreg"))-1,
                  depensationModel = which(d.type == c("None","A")) - 1),
      # inits
      Params = list(log_r0 = inits[1], log_sigR = inits[2],logit_s=inits[3], log_d = 0),
      # bounds
      lower=lower, upper=upper,
      #
      ReportSD = SDreport
    )
    
    # Compile TMB inputs 
      Map = list()
      if(d.type == "None")
          Map$log_d = factor(NA)
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
        attr(object,"SD") <- SD
    }
    
    # LOAD output in FLSR
    
    # DEBUG HACK
    model(object) <- switch(model, bevholtSV=bevholt, rickerSV=ricker,segreg=segreg,bevholtDa=bevholtDa)
    
    fitted(object) <- c(Report$rec_hat)
    residuals(object) <- log(rec(object)) - log(fitted(object))
    
    params(object) <- FLPar(a=Report$a, b=Report$b,d=Report$d)
    
    if(d.type=="none") params(object) = params(object)[1:2] 
    
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
      attr(object,"SV") = data.frame(s=Report$s,sigmaR=Report$sigR,R0=Report$r0,rho=rho,B0=Report$r0*spr0ref, d = ifelse(d.type=="None",NA,Report$d))
    } else{
      attr(object,"SV") = data.frame(s=NA,sigmaR=Report$sigR,R0=Report$r0,rho=rho,B0=Report$r0*spr0ref, d = ifelse(d.type=="None",NA,Report$d))
      
    }
    
    
    if(SDreport & report.sR0==FALSE){
      object@vcov = matrix(SD$cov,nrow=2,dimnames = list(c("a","b"),c("a","b")))
    }
  } # End TMB model
  
  attr(object,"settings") = list(s=s.inp,s.est=s.est,s.logitsd=s.logitsd,spr0=spr0,lplim=lplim,uplim=uplim,nyears=nyears,
                                 inits=inits,lower=lower,upper=upper,d.type=d.type)
  
  
  
  return(object)
}
)
