

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
#' plotsrs(all,path=TTRUE)

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


#' plotsrs {{{ based on plot(FLSRs)

#' Plots FLSRs, i.e. multiple S-R relationships
#' @param object of class FLSRS 
#' @param path connect points sequentially 
#' @return ggpplot
#' @export 

plotsrs <- function(object,path=TRUE) {
            x= object
            uns <- units(x[[1]])
            
           
              dat <- Reduce(rbind, Map(function(x, i)
                cbind(sr=i, model.frame(FLQuants(ssb=ssb(x), rec=rec(x)), drop=TRUE)),
                x[1], names(x[1])))
            
            # EXTRACT models & pars
            mods <- lapply(x, 'model')
            pars <- lapply(x, 'params')
            inp <- data.frame(ssb=seq(0, max(dat$ssb), length=100), rec=NA)
            
            # RESULTS
            res <- lapply(names(mods), function(x) {
              data.frame(sr=x, ssb=inp$ssb,
                         rec=eval(as.list(mods[[x]])[[3]], c(list(ssb=inp$ssb), as(pars[[x]], 'list')))
              )
            })
            
            res <- Reduce('rbind', res)
            res$sr = factor(res$sr,levels=object@names)
            # GET plot
            p <- ggplot(na.omit(res), aes(x=ssb, y=rec, colour=sr)) +
              geom_line(size=0.8) + theme_bw()+
              xlab(as.expression(paste0("SSB (", sub('\\*', '%.%', uns$ssb), ")"))) +
              ylab(as.expression(paste0("Recruits (", sub('\\*', '%.%', uns$rec), ")"))) +
              theme(legend.position="right",legend.title = element_blank()) 
            
            if(path){
              p = p+geom_path(data=dat,cex=0.1,col="black",linetype=1,alpha=0.4)+
                geom_point(data=dat,color=1,cex=c(rep(1.5,nrow(dat)-1),2.),pch=21,fill=c(rep("white",nrow(dat)-1),"black"))
              
            } else {
              p <- p+geom_point(data=dat,color=1,fill="white",pch=21,cex=1.8)
            }
            
            return(p)
          }
          
# }}}


