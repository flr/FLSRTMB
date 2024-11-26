
# plotsrts {{{

#' Plots time series of observed and predicted recruitment
#'
#' @param object Input FLSR of FLSRs object
#' @return ggplot
#' @export
#' @examples
#' data(ple4)
#' bh <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' plotsrts(bh)
#' # Try more models
#' hs <- srrTMB(as.FLSR(ple4,model=segreg),spr0=(spr0y(ple4)),plim=0.05,pmax=0.2)
#' ri <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=mean(spr0y(ple4)))
#' bh.tv <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
#' ri.tv <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=spr0y(ple4))
#' all= FLSRs(hs=hs,bh=bh,ri=ri,bh.tv=bh.tv,ri.tv=ri.tv)
#' plotsrts(all)
#' do.call(c,lapply(all,function(x)AIC(x))) # AIC


plotsrts <- function(object,relative=TRUE){

  
  
  
if(class(object)=="FLSR"){  
df = data.frame(as.data.frame(rec(object)),ssb=as.data.frame(ssb(object))$data)
df2 = data.frame(as.data.frame(fitted(object)),name=object@name)  

p = ggplot(df,aes(year,data))+theme_bw()+geom_point(pch=21,color=1,fill="white")+
  geom_line(data=df2,aes(year,data),size=0.8)+theme(legend.title = element_blank())+
  ylab("Recruitment")+xlab("Year")+theme(legend.position = "none") 
}
  
if(class(object)=="FLSRs"){  
  df = data.frame(as.data.frame(rec(object[[1]])),ssb=as.data.frame(ssb(object[[1]]))$data)
  df2 = do.call(rbind,Map(function(x,y){
    data.frame(as.data.frame(fitted(x)),name=y)
  },x=object,y=as.list(object@names))) 
  
  df2$name = factor(df2$name,levels=unique(df2$name))
  
  p = ggplot(df2,aes(year,data,color=name))+theme_bw()+
    geom_line(size=0.8)+theme(legend.title = element_blank())+
    ylab("Recruitment")+xlab("Year")+geom_point(data=df,aes(year,data) ,pch=21,color=1,fill="white")

  if(length(object)==1) p= p+theme(legend.position = "none")  
}  
return(p)
}  
#}}}


#' plotsrs based on plot(FLSRs)

#' Plots FLSRs, i.e. multiple S-R relationships
#' @param object of class FLSRS 
#' @param path connect points sequentially 
#' @return ggpplot
#' @export 
#' @examples 
#' data(ple4)
#' bh <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' plotsrs(bh)
#' plotsrs(bh,b0=TRUE) # plot through to B0
#' plotsrs(bh,b0=TRUE,rel=TRUE) # plot relative to B0
#' # Try more models
#' hs <- srrTMB(as.FLSR(ple4,model=segreg),spr0=(spr0y(ple4)),plim=0.05,pmax=0.2)
#' ri <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=mean(spr0y(ple4)))
#' bh.tv <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
#' ri.tv <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=spr0y(ple4))
#' all= FLSRs(hs=hs,bh=bh,ri=ri,bh.tv=bh.tv,ri.tv=ri.tv)
#' plotsrs(all)
#' plotsrs(all,path=FALSE)
#' # plot through to b0
#' plotsrs(all,b0=TRUE)
#' # plot all relative relative to B0
#' plotsrs(all,rel=TRUE)
#' plotsrs(all,rel=TRUE,b0=TRUE)

plotsrs <- function(object,path=TRUE,b0=FALSE,rel=FALSE){
            
            if(class(object)=="FLSR"){
              object=FLSRs(object)
            }
               
            x= object
            uns <- units(x[[1]])
            
           
              dat <- Reduce(rbind, Map(function(x, i)
                cbind(sr=i, model.frame(FLQuants(ssb=ssb(x), rec=rec(x)), drop=TRUE)),
                x, names(x)))
            
            # EXTRACT models & pars
            mods <- lapply(x, 'model')
            pars <- lapply(x, 'params')
            inp <- data.frame(ssb=seq(0, max(dat$ssb), length=100), rec=NA)
            
            # RESULTS
            res <- Map( function(x,y) {
              d. <- dat[dat$sr==x,]
              if(b0) sub <- data.frame(ssb=seq(0, max(d.$ssb,y@SV[["B0"]]), length=100), rec=NA)
              if(!b0) sub <- inp
              if(rel){
                data.frame(sr=x, ssb=sub$ssb/y@SV[["B0"]],
                           rec=eval(as.list(mods[[x]])[[3]], c(list(ssb=sub$ssb), as(pars[[x]], 'list')))/y@SV[["R0"]]
                )
              } else {
                  
              data.frame(sr=x, ssb=sub$ssb,
                         rec=eval(as.list(mods[[x]])[[3]], c(list(ssb=sub$ssb), as(pars[[x]], 'list')))
              )
              }
            },x=names(mods),y=object)
            
            res <- Reduce('rbind', res)
            res$sr = factor(res$sr,levels=object@names)
            # GET plot
            p <- ggplot(na.omit(res), aes(x=ssb, y=rec, colour=sr,fill=sr)) +
              geom_line(size=0.8) + theme_bw()
             
              if(length(object)==1){
              p <- p+theme(legend.position="none",legend.title = element_blank()) 
              } else {
                p <- p+theme(legend.position="right",legend.title = element_blank())   
              }
             if(rel){
              p = p+ xlab(expression(SSB/SSB[0])) +
                 ylab(expression(R/R[0])) 
             } else{
               p = p+ xlab(expression(SSB)) +
                 ylab(expression(Recruits)) 
               
             }
            
          if(path){
              rps=aggregate(ssb~sr,dat,length)[,2]
              
                if(rel){
                B0 = an(do.call(c,lapply(object,function(x)x@SV[["B0"]])))
                R0 = an(do.call(c,lapply(object,function(x)x@SV[["R0"]])))
                B0s=rep(B0,rps)
                R0s=rep(R0,rps)
                dat$ssb = dat$ssb/B0s
                dat$rec = dat$rec/R0s
                } 
            
            if(length(unique(dat$ssb))==nrow(dat)/length(unique(dat$sr))){    
              p = p+geom_path(data=dat,cex=0.1,col="black",linetype=1,alpha=0.4)+
                geom_point(data=dat,color=1,cex=c(rep(1.5,nrow(dat)-1),2.),pch=21,fill=c(rep("white",nrow(dat)-1),"black"))
            } else {
              nsr = length(object)
              pchs= cexs= NULL
              for(i in 1:length(rps)){
                pchs = c(pchs,c(rep(21,rps[i]-1),22))
                cexs = c(cexs,c(rep(1.5,rps[i]-1),2))}
              p = p+geom_path(data=dat,cex=0.1,col="black",linetype=1,alpha=0.4)+
                geom_point(data=dat,aes(fill=sr),color=1,cex=cexs,pch=pchs,alpha=0.5)
              
              }
              p <- p+geom_line(data=na.omit(res),size=0.8) 
            } # path
            
    return(p)
} # }}}


#' sprior plot 

#' Plots the logit prior distribution for steepness 
#' @param s steepness, default 0.6 for a approx. uniform prior with s.logistsd = 20
#' @param s.logitsd s steepness, default 20 for a approx. uniform prior with s = 0.6
#' @param ll lower bound of s = 0.2
#' @param ul lower bound of s = 1
#' @return ggplot
#' @export
#' @examples
#' sprior() # approx. uniform with some curving on bounds
#' sprior(s=0.8,s.logitsd=0.5) 
sprior <- function(s=0.6,s.logitsd=20,ll=0.2,ul=1){
  d=0.00001
  x.logit = seq(-10,10,0.01)
  mu.logit = to_logits(s,ll=0.2+d,ul=1-d)
  y = dnorm(x.logit,mu.logit,s.logitsd)
  x= from_logits(x.logit)
  df = data.frame(x=c(0.2,x,1),Density=c(0,y,0))
  ggplot(df, aes(x,Density))+theme_bw()+
    geom_area(aes(y=Density),fill="grey",alpha=1,col=1)+
    geom_vline(xintercept = s,size=0.5,col=2,linetype="dashed")+xlab("Steepness s")
}


#' blimprior plot 

#' Plots the bounds of the hockey-stick break-point  
#' @param lplim steepness, default 0.6 for a approx. uniform prior with s.logistsd = 20
#' @param uplim s steepness, default 20 for a approx. uniform prior with s = 0.6
#' @param s.sdlogit default 20
#' @param par parameter on x-axis default "plim", else "sstar" 
#' @return ggplot
#' @export
#' @examples
#' blimprior() # approx. uniform with some curving on bounds
#' blimprior(lplim=0.001,uplim=0.3,s.logitsd=20) 
#' # Non-bias corrected
#' blimprior(lplim=0.001,uplim=0.3,s.logitsd=20,bias.correct=FALSE)

blimprior <- function(lplim=0.01,uplim=0.3,s.logitsd=50,par="plim",bias.correct=TRUE){
  d=0.00001
  ll = lplim/uplim 
  ul = 1
  mu = mean(c(lplim,uplim))
  s = 1/(mu/lplim)
  # bias.correct
  if(bias.correct){
  plim = lplim
  pmax = seq(0.0001,0.9999,0.0001)
  pmax = pmax[which((lplim/from_logits(-10, ll = plim/pmax, ul = 1)-uplim)^2==min((lplim/from_logits(-10, ll = plim/pmax, ul = 1)-uplim)^2))]
  mu = mean(c(lplim,(uplim+pmax)/2)) #><> fix
  ll =plim/pmax
  ul = 1
  s = 1/(mu/lplim) 
  }
  # checks  
  1/s*lplim  
  lplim/s
  
  x.logit = seq(-10,10,0.01)
  mu.logit = to_logits(s,ll=ll+d,ul=1-d)
  y = dnorm(x.logit,mu.logit,s.logitsd)
  x= from_logits(x.logit,ll=ll,ul=1)
  if(par=="plim"){
  x=  lplim/x
  df = data.frame(x=c(lplim,x,uplim),Density=c(0,y,0))
  p <- ggplot(df, aes(x,Density))+theme_bw()+
    geom_area(aes(y=Density),fill="grey",alpha=1,col=1)+
    xlim(0,min(1.5*uplim,1))+xlab(expression("Bounds"~B[lim]/B[0]))
  } else {
    df = data.frame(x=c(ll,x,ul),Density=c(0,y,0))
    p <- ggplot(df, aes(x,Density))+theme_bw()+
      geom_area(aes(y=Density),fill="grey",alpha=1,col=1)+
      xlim(0,1)+xlab(expression(s^"*"))
  }
  return(p)
}
 

#' dprior plot 

#' Plots the logit prior distribution for depensation 
#' @param d depensation, default 1 for a approx. uniform prior with s.logistsd = 20
#' @param d.logitsd dependation sd, default 20 for a approx. uniform prior with s = 0.6
#' @param ll lower bound of d = 0.25
#' @param ul lower bound of d = 4
#' @return ggplot
#' @export
#' @examples
#' dprior() # approx. uniform with some curving on bounds
#' dprior(d=1,d.logitsd=2) 
dprior <- function(d=1,d.logitsd=100,ll=0.5,ul=3){
  C=0.00001
  x.logit = seq(-10,10,0.01)
  mu.logit = to_logitd(d,ll=ll+C,ul=ul-C)
  y = dnorm(x.logit,mu.logit,d.logitsd)
  x= (from_logitd(x.logit,ll=ll+C,ul=ul-C))
  df = data.frame(x=c(x),Density=c(y))
  ggplot(df, aes(x,Density))+theme_bw()+
    scale_x_continuous(limits = c(0, NA))+
    geom_area(aes(y=Density),fill="grey",alpha=1,col=1)+
    geom_vline(xintercept = d,size=0.5,col=2,linetype="dashed")+
    xlab("Depensation d")
}

