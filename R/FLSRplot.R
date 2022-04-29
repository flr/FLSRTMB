
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
#' plotsrs(all,path=TRUE)

plotsrts <- function(object){

if(class(object)=="FLSR"){  
df = data.frame(as.data.frame(rec(object)),ssb=as.data.frame(ssb(object))$data)
df2 = data.frame(as.data.frame(fitted(object)),name=object@name)  

p = ggplot(df,aes(year,data))+theme_bw()+geom_point(pch=21,color=1,fill="white")+
  geom_line(data=df2,aes(year,data),size=0.8)+theme(legend.title = element_blank())+
  ylab("Recruitment")+xlab("Year")  
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


