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

#' resilience()
#'
#' Function to compute r, generation time (gt) and annual reproductive rate (alpha)
#' @param object class FLStock
#' @param s steepness of the stock recruitment relationship
#' @return FLQuants with FLQuant r, gt and alpha
#' @export
#' @author Henning Winker and Laurence Kell 
resilience <- function(object,s=0.7){ 
  spr0 = spr0y(object)
  spr0_a = spr0y(object,byage=T)
  # Reproductive output Rs for bonyfish
  rs = 4*s/(spr0*(1-s))
  wm = stock.wt(object)*mat(object)
  nage = dim(object)[1]
  nyr = dim(object)[2]
  age = dims(object)$min:dims(object)$max 
  
  r = FLQuant(NA,dimnames=list(year=dims(object)$minyear:dims(object)$maxyear),units="1/year")  
  gt = FLQuant(NA,dimnames=list(year=dims(object)$minyear:dims(object)$maxyear),units="years")  
  
  # Make Leslie matrix 
  for(y in 1:dim(object)[2]){
    if(is.na(spr0[,y])) next() 
    L=mat.or.vec(nage,nage)
    L[1,] =   an(rs[,y])*wm[,y]
    #fill rest of Matrix with Survival
    for(i  in 2:nage) L[i,(i-1)] = exp(-m(object)[i,1]) 
    # Plus group
    L[nage,(i)] = exp(-m(object)[nage,y]) 
    # Net reproductive rate
    r[,y]=log(as.numeric(eigen(L)$values[1]))
    gt[,y] = sum(age*spr0_a[,y])/spr0[,y]
  }
  
  # return intrinsic rate of population increase r and generation GT
  return(FLQuants(r=r,gt=gt))
}
