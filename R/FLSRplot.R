

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
#' hs <- srrTMB(as.FLSR(ple4,model=segreg),spr0=(spr0y(ple4)),plim=0.02,pmax=0.3)
#' ri <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=mean(spr0y(ple4)))
#' bh.tv <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
#' ri.tv <- srrTMB(as.FLSR(ple4,model=rickerSV),spr0=spr0y(ple4))
#' all= FLSRs(hs=hs,bh=bh,ri=ri,bh.tv=bh.tv,ri.tv=ri.tv)
#' plotsrts(all)
#' do.call(c,lapply(all,function(x)AIC(x))) # AIC
#' plot(all)+theme(legend.position="right")

plotsrts <- function(object){

if(class(object)=="FLSR"){  
df = data.frame(as.data.frame(rec(object)),ssb=as.data.frame(ssb(object))$data)
df2 = data.frame(as.data.frame(fitted(object)),name=object@name)  

p = ggplot(df,aes(year,data))+geom_point(pch=21,color=1,fill="white")+
  geom_line(data=df2,aes(year,data))+theme(legend.title = element_blank())+
  ylab("Recruitment")+xlab("Year")  
}
if(class(object)=="FLSRs"){  
  df = data.frame(as.data.frame(rec(object[[1]])),ssb=as.data.frame(ssb(object[[1]]))$data)
  df2 = do.call(rbind,Map(function(x,y){
    data.frame(as.data.frame(fitted(x)),name=y)
  },x=object,y=as.list(object@names))) 
  
  df2$name = factor(df2$name,levels=unique(df2$name))
  
  p = ggplot(df,aes(year,data))+geom_point(pch=21,color=1,fill="white")+
    geom_line(data=df2,aes(year,data,color=name))+theme(legend.title = element_blank())+
    ylab("Recruitment")+xlab("Year")
}  
return(p)
}  
#}}}