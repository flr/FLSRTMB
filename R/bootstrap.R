# bootstrap.R - DESC
# /home/mosqu003/Projects/FLR/pkgs/mine/FLSRTMB/R/bootstrap.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# bootstrapSR {{{

#' Bootstrap fits of mutliple stock-recruits relationships
#'
#' Definition ...
#'
#' The returned 'FLPar' object contains
#'
#' @param x An object of class 'FLStock'.
#' @param iter Number of bootstrap iterations, 'numeric'.
#' @param models Name(s) of model(s) to fit, 'character'. See Details.
#' @param verbose Should progress be reported, 'logical'.
#'
#' @return A list with elements An object or class 'FLPar' containing the estimated paramaters.
#'
#' @name bootstrapSR
#'
#' @author Iago Mosqueira (WMR)
#' @seealso \link{FLSR}, link{srrTMB}
#' @keywords models
#' @examples
#' data(ple4)
#' bootstrapSR(ple4, iter=50, model=c("bevholt", "segreg"))

bootstrapSR <- function(x, iters=100, method=c("best", "loglik", "aic"),
  models=c("bevholt", "ricker", "segreg"), verbose=TRUE, ...) {

  # COMPUTE average of spr0 by year
  spr0x <- yearMeans(spr0y(x))

  # BUILD FLSR
  x <- as.FLSR(x)
  
  # SAMPLES, year * iters
  id <- matrix(sample(seq(dim(x)[2]), dim(x)[2] * iters, replace=TRUE),
    nrow=dim(x)[2], ncol=iters)

  # SELECT models
  models <- match.arg(models, several.ok=TRUE)

  # MATCH with those in FLSRTMB
  mod <- list(bevholt=bevholtSV, ricker=rickerSV, segreg=segreg)[models]

  method <- match.arg(method)

  # BOOTSTRAP

  p <- progressor(along=seq(iters), offset=1)

  res <- foreach(i=seq(iters)) %dopar% {

    y <- x

    rec(y) <- rec(y)[, id[, i]]
    ssb(y) <- ssb(y)[, id[, i]]
  
    fits <- lapply(mod, function(m) {
      model(y) <- m
      srrTMB(y, spr0=spr0x, verbose=FALSE, ...)
    })

    if(verbose)
      p(message = sprintf(paste0("[", i, "]")))

    # COMPUTE LogLik and AIC across fits

    llkhds <- unlist(lapply(fits, 'logLik'))
    aics <- unlist(lapply(fits, 'AIC'))
    
    # - FIND BEST model

    # best: largest logLik
    if (method == "best") {
      best <- which.max(llkhds)
    # logLik
    } else if(method == "loglik") {
      probs <- llkhds / sum(llkhds)
      u <- runif(1, 0, 1)
      best <- which(u <= cumsum(probs))[1]
    # aic: relative likelihood from AIC
    } else if(method == "aic") {
      daic <- aics - min(aics)
      relkhd <- exp(-0.5 * daic)
      aicw <- relkhd / sum(relkhd)
      u <- runif(1, 0, 1)
      best <- which(u <= cumsum(aicw))[1]
    }
    
    # MATCH models: bevholt=1, ricker=2, segreg=3
    m <- match(models[best], c("bevholt", "ricker", "segreg"))

    fit <- fits[[best]]

    # RETURN: params, model, spr0, logLik, SR fit params
    rbind(params(fit), FLPar(m=m, spr0=spr0x, logLik=llkhds[best]),
      FLPar(attr(fit, 'SV')))
  }

  # COMBINE along iters
  params <- Reduce(combine, res)

  return(params)
}
# }}}

# plot_bootstrapSR {{{

plot_bootstrapSR <- function(fits, params) {

  it <- dims(params)$iter

  # CREATE ssb vector for range
  ssbs <- FLQuant(seq(1, max(ssb(fits[[1]])), length=100),
    dim=c(1, 100, 1, 1, 1, it))

  # PREDICT rec at ssbs 
  # recs <- FLCore:::mixed(params$a, params$b, params$m, ssb=ssbs)

  recs <- Reduce(combine, lapply(1:500, function(i)
    FLCore:::mixed(iter(params, i)$a, iter(params, i)$b, iter(params, i)$m,
      ssb=iter(ssbs,i))
    ))

  # ADD error
  recs <- exp(log(recs) + rnorm(it, recs %=% 0, sd=(c(params$sigmaR))) -
    (0.5 * c(params$sigmaR) ^ 2))

  # CREATE df for plotting
  dat <- model.frame(FLQuants(rec=recs, ssb=ssbs), drop=TRUE)

  # CREATE quantiles df
  qdat <- subset(dat, iter == 1)

  qdat$q05 <- tapply(dat$rec, dat$ssb, FUN = function(x)
    quantile(x, prob = 0.05), simplify = TRUE)
  qdat$q50 <- tapply(dat$rec, dat$ssb, FUN = function(x)
    quantile(x, prob = 0.50), simplify = TRUE)
  qdat$q95 <- tapply(dat$rec, dat$ssb, FUN = function(x)
    quantile(x, prob = 0.95), simplify = TRUE)

  # DROP extreme values
  dat <- subset(dat, rec <= max(rec(run)) * 1.25)

  # PLOT fits
  plotsrs(fits) +
    annotate("text", x=-Inf, y=Inf, hjust = -0.2, vjust = 1.5,
      label=.table_srmodels(params)) +
    # QUANTILE smoothers
    geom_smooth(data=qdat, aes(y=q50),
      colour="black", fill=NA, linewidth=0.5,
      method='loess', formula=y~x, se=FALSE) +
    geom_smooth(data=qdat, aes(y=q05),
      colour="black", fill=NA, linewidth=0.5, linetype=2,
      method='loess', formula=y~x, se=FALSE) +
    geom_smooth(data=qdat, aes(y=q95),
      colour="black", fill=NA, linewidth=0.5, linetype=2,
      method='loess', formula=y~x, se=FALSE) -
    # POINTS
    geom_point(data=dat, aes(x=jitter(ssb, 2), rec), 
      colour="gray", fill="white", alpha=0.1, size=0.5) +
    theme(legend.position="none") +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = scientific)
}
# }}}

# .table_srmodels {{{
.table_srmodels <- function(x) {

  tab <- table(x$m)
  mod <- c("Bevholt:", "Ricker:", "Segreg:")[as.numeric(names(tab))]
  val <- round(c(tab) / sum(tab), 2)

  paste(mod, val, collapse="\n")
}
# }}}
