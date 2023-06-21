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

bootstrapSR <- function(x, iters=100, method=c("best", "rejection"),
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

    llkhds <- unlist(lapply(fits, 'logLik'))

    # FIND BEST model
    if(method == "rejection") {
      probs <- llkhds / sum(llkhds)
      u <- runif(1, 0, 1)
      best <- fits[[which(u <= cumsum(probs))[1]]]
      }
    else if (method == "best") {
      best <- fits[[which.max(llkhds)]]
    }

    # MATCH models: bevholt=1, ricker=2, segreg=3
    m <- match(models[which.max(llkhds)], c("bevholt", "ricker", "segreg"))

    rbind(params(best), FLPar(m=m, spr0=spr0x), FLPar(attr(best, 'SV')))
  }

  # COMBINE along iters
  params <- Reduce(combine, res)

  return(params)
}
# }}}
