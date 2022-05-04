
# srrjitter {{{

#' Jitter of S-R fits
#'
#' @param fit  fit of srrTBM() 
#' @param steps number of jitter steps
#' @return list
#' @export
#' @examples
#' data(ple4)
#' hs = srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),lplim=0.07,uplim=0.5)
#' plotsrs(hs)
#' jitter = srrjitter(hs)
#' plotsrs(jitter$groups)
#' plotsrs(jitter$best) # Best
#' # Relax lower bound
#' plotsrs(FLSRs(init=hs,best=jitter$best))

srrjitter <- function(fit,steps=100){
srrfit = fit 
settings = srrfit@settings
model = SRModelName(fit@model)
if(model=="segreg"){
lplim = settings$lplim  
uplim = settings$uplim
ll =lplim/uplim
ul = 1
s.range = from_logits(seq(-7,6,13/(steps-1)),ll=lplim/uplim,ul=0.99)[1:steps]
}
if(model=="ricker"){
    model(srrfit) = "rickerSV"
    s.range = from_logits(seq(-7,6,13/(steps-1)),ll=0.21,ul=10)[1:steps]
    
}    
if(model=="bevholt"){
  model(srrfit) = "bevholtSV"
  s.range = from_logits(seq(-7,6,13/(steps-1)),ll=0.21,ul=0.99)[1:steps]
}

jf = FLSRs(
     lapply(s.range,function(x){
      srrTMB(srrfit,spr0=settings$spr0,
            s=x,s.est=settings$s.est,
            s.logitsd = settings$s.logitsd,
            lplim=settings$lplim,uplim=settings$uplim,
            SDreport = FALSE,verbose=FALSE)}))

names(jf) = ac(1:steps)    

aic = as.vector(do.call(c,lapply(jf, function(x){
  AIC(x)
})))

LL = as.vector(do.call(c,lapply(jf, function(x){
  logLik(x)[[1]]
})))

blim =  as.vector(do.call(c,lapply(jf, function(x){
  params(x)[[2]]
})))

group=match(round(LL,3),unique(round(LL,3)))
run = 1:steps


igroup=NULL 
for(i in 1:length(unique(group))){
  igroup=c(igroup,run[group%in%i][LL[group%in%i]==min(LL[group%in%i])][1])
}

if(model=="segreg"){
out = data.frame(x=lplim/s.range,LL=LL,AIC=aic,run=1:steps,group=factor(group))
} else {
out = data.frame(x=s.range,LL=LL,AIC=aic,run=1:steps,group=factor(group))
}

ibest = which(out$LL==out$LL)[1]

init = fit
best = srrTMB(srrfit,spr0=settings$spr0,
             s=s.range[ibest],s.est=settings$s.est,
             s.logitsd = settings$s.logitsd,
             lplim=settings$lplim,uplim=settings$uplim,
             SDreport = TRUE,verbose=FALSE)
groups = jf[igroup] 
names(groups) = unique(group)

res = list(init=init,best=best,groups=groups,jitter=out)

return(res)
} # }}}


# plotjitter {{{
#' plots Jitter results
#'
#' @param jitter  output from srrjitter() 
#' @return ggplot
#' @export
#' @examples
#' data(ple4)
#' hs = srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),lplim=0.07,uplim=0.2)
#' plotsrs(hs)
#' jitter = srrjitter(hs)
#' plotsrs(jitter$groups)
#' plotjitter(jitter)
#' bh = srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
#' jitbh = srrjitter(bh)
#' plotjitter(jitbh)
#' plotsrs(jitbh$groups)

plotjitter <- function(jitter){
  model = SRModelName(model(jitter$best))
  p = ggplot(jitter$jitter,aes(x,LL,color=group))+theme_bw()+
    geom_point(size=2)+ylim(min(jitter$jitter$LL)-1,max(jitter$jitter$LL)+1)+
    geom_hline(yintercept = logLik(jitter$init)[[1]])
  if(model=="segreg"){
    p = p+ylab("LogLik")+xlab(expression(B[lim]/B[0]))
  } else{
    p = p+ylab("LogLik")+xlab("s")
  }
  return(p)
}
  
  