#' from_logits()
#'
#' convert steepness from logit
#' @param logit_h logit(steepness)
#' @return steepness h 
#' @export
from_logits <- function(logit_h){
  out=  0.2001 + 0.7998*1/(1+exp(-logit_h))
  out}

#' to_logits()
#'
#' convert steepness to logit
#' @param h steepness h
#' @return logit transformed steepness h 
#' @export
to_logits <- function(h){
  -log(0.7998/(h-0.2001)-1) 
}

#' spr0y()
#'
#' Function to compute annual spr0 
#' @param object class FLStock
#' @param byage if TRUE it return spr0_at_age
#' @return FLQuant with annual spr0y  
#' @export
#' @author Laurence Kell
spr0y<-function(object,byage=FALSE){
  survivors=exp(-apply(m(object),2,cumsum))
  survivors[-1]=survivors[-dim(survivors)[1]]
  survivors[1]=1
  expZ=exp(-m(object[dim(m(object))[1]]))
  if (!is.na(range(object)["plusgroup"]))
    survivors[dim(m(object))[1]]=survivors[dim(m(object))[1]]*(-1.0/(expZ-1.0))
  
  fec=mat(object)*stock.wt(object)*exp(-m(object)*m.spwn(object))
  
  if(byage) rtn = fec * survivors
  if(!byage) rtn =  apply(fec * survivors, 2, sum)
  rtn}



#' spr0y()
#'
#' Function to compute annual sprF as function of F_a
#' @param object class FLStock
#' @param byage if TRUE it return sprF_at_age
#' @return FLQuant with annual sprFy  
#' @export
#' @author Henning Winker and Laurence Kell 
sprFy = function(object,byage = FALSE){
  survivors = exp(-apply(m(object)+harvest(object), 2, cumsum))
  survivors[-1] = survivors[-dim(survivors)[1]]
  survivors[1] = 1
  expZ = exp(-(m(object[dim(m(object))[1]])+harvest(object[dim(m(object))[1]])))
  if (!is.na(range(object)["plusgroup"])) 
    survivors[dim(m(object))[1]] = survivors[dim(m(object))[1]] * 
    (-1/(expZ - 1))
  fec = mat(object) * stock.wt(object) * exp(-(m(object)+harvest(object)) * m.spwn(object))
  if(byage) rtn = fec * survivors
  if(!byage) rtn =  apply(fec * survivors, 2, sum)
  rtn}
