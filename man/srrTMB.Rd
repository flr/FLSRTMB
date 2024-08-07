% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLSR.R
\name{srrTMB}
\alias{srrTMB}
\alias{srrTMB,FLSRs-method}
\alias{srrTMB,FLSR-method}
\title{Fits Stock Recruitment Relationships (SRR) in TMB}
\usage{
srrTMB(object, ...)

\S4method{srrTMB}{FLSRs}(object, ...)

\S4method{srrTMB}{FLSR}(
  object,
  spr0 = "missing",
  s = NULL,
  s.est = TRUE,
  s.logitsd = 20,
  r0.pr = "missing",
  lplim = 0.001,
  uplim = 0.3,
  Blim = "missing",
  plim = lplim,
  pmax = uplim,
  nyears = NULL,
  report.sR0 = FALSE,
  inits = NULL,
  lower = NULL,
  upper = NULL,
  SDreport = TRUE,
  verbose = FALSE,
  rm.yrs = "missing"
)
}
\arguments{
\item{object}{Input FLSR = as.FLSR(stock,model) object with current model options
\itemize{
\item bevholtSV
\item bevholtDa
\item rickerSV
\item segreg
\item geomean
}}

\item{spr0}{unfished spawning biomass per recruit from FLCore::spr0(FLStock)}

\item{s}{steepness parameter of SRR (fixed or prior mean)}

\item{s.est}{option to estimate steepness}

\item{s.logitsd}{prior sd for logit(s), default is 1.4 (flat) if s.est = TRUE}

\item{r0.pr}{option to condition models on r0 priors (NULL = geomean)}

\item{lplim}{lower bound of spawning ratio potential SRP, default 0.0001}

\item{uplim}{upper bound of plausible spawning ratio potential SRP , default 0.3}

\item{Blim}{fixing Blim, only works with segreg}

\item{plim}{depreciated plim = usrp}

\item{pmax}{depreciated pmax = lsrp}

\item{nyears}{yearMeans from the tail used to compute a,b from the reference spr0 (default all years)}

\item{report.sR0}{option to report s and R0 instead of a,b}

\item{inits}{option to specify initial values of log(r0), log(SigR) and logit(s)}

\item{lower}{option to specify lower bounds of log(r0), log(SigR) and logit(s)}

\item{upper}{option to specify upper bounds of log(r0), log(SigR) and logit(s)}

\item{SDreport}{option to converge hessian and get vcov}

\item{verbose}{if TRUE, it shows tracing}

\item{rm.yrs}{remove recruitment years from model fit}
}
\value{
A list containing elements 'FLSR', of class \emph{FLSR}
}
\description{
Fits Stock Recruitment Relationships (SRR) in TMB
}
\examples{
data(ple4)
gm <- srrTMB(as.FLSR(ple4,model=geomean),spr0=mean(spr0y(ple4)))
bh <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
ri <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=spr0y(ple4))
hs <- srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),lplim=0.05,uplim=0.2)
srs = FLSRs(gm=gm,bh=bh,ri=ri,hs=hs) # combine
plotsrs(srs) 
plotsrts(srs)  # relative
plotsrs(srs[2:4],b0=TRUE) # through to B0
plotsrs(srs[2:4],b0=TRUE,rel=TRUE)  # relative
gm@SV # estimates
do.call(rbind,lapply(srs,AIC))
}
