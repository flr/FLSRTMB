# FLSRTMB

*Beta version FLSRTMB: fitting stock-recruitment with TMB in FLR*  

### Authors: Henning Winker (EC-JRC) & Iago Mosqueira (WUR)*

# Features
+ Rapidly fits stock-recruitment models with very convergence properties 
+ Uses FLR classes as input and output 
+ Enables fitting spawner-recruitment model with time-varying spr0(y)  
+ Enables use of steepness priors from fishlife
+ Provides options for conditioned hockey-stick with a break point b > plim (e.g. plim = 0.1B0) 

# Installation
Installing `FLSRTMB` requires the librabry(devtools), which can be install by 'install.packages('devtools')' and a R version >= 3.5. `FLSRTMB` also requires the latest version of `library(FLCore)` and suggests using `library(ggplotFL)` for plotting. All can be installed from github.

`devtools::install_github("flr/FLCore")`

`devtools::install_github("flr/ggplotFL")`

`devtools::install_github("henning-winker/FLSRTMB")`

`library(FLSRTMB)`

Compiling C++ in windows can be troublesome. As an alternative to installing from github, a windows package binary zip file can be downloaded [here](https://github.com/Henning-Winker/FLSRTMB/tree/main/BinaryPackage/win).

# Quick test drive

Example data are loaded for North Sea from `FLCore`

`data(ple4)`

The Beverton and Holt model is parameterized as a function of steepness (s) and unfished spawning potential ratio SPR0.
To create the `FLSR` input object, the `model=bevholtSV` is selected. 

### Beverton-Holt Model

`sr <- as.FLSR(ple4,model=bevholtSV)`

The function spr0y computes annual spr0. A good starting point can be the average  

`spr0 <- yearMeans(spr0y(ple4))`

Fit a bevholt model without constraints 

`bh = srrTMB(sr,s.est=T, spr0=spr0)`

`plot(bh)`

Note the output `params` are a and b by default but can be changed to R0 and s.

`params(srrTMB(sr,s.est=T, spr0=spr0,report.sR0 = TRUE))`


For many stock spr0 time-varying due to variations in weight-at-age, M-at-age or maturity-at-age 
plot(spr0y(ple4))+ylab("spr0")

To account for reduced or increased unfished stock sizes given these variation, `spr0y()` can be directly inputted

`bh.y = srrTMB(sr,s.est=T, spr0=spr0y(ple4))`

Note that by default a,b are computed from the average spr0 across all years. However, if the analyst suspects a non reversable systematic change, the `spr0` reference can also be taken from the most recent years

`bh.y3 = srrTMB(sr,s.est=T, spr0=spr0y(ple4),nyears=3)`

Compare fits 

`plot(FLSRs(spr0=bh,spr0y=bh.y,spr0y3 = bh.y3))+theme(legend.position="right")`

Option to estimate steepness with prior, e.g. from meta-analysis, such as FishLife (Thorson, 2020)
This requires to specify `s` with the  prior mean and `s.logitsd` is the sd on logit scale 

`bh.prior = srrTMB(sr,s.est=T,s=0.7,s.logitsd=0.4, spr0=spr0y(ple4))`

`plot(FLSRs(spr0=bh,spr0y=bh.y,s.prior = bh.prior))+theme(legend.position="right")`

It is also possible to fix steepness

`s = c(0.75,0.8,0.85,0.9,0.95)`

`bhs <- FLSRs(sapply(s, function(x) { return( srrTMB(sr,s=x,s.est=F, spr0=spr0y(ple4)))}))`

`bhs@names = c(paste("s =",round(s,3)))`

`plot(bhs)+theme(legend.position="right")`

### Hockey Stick (segmented regression)

First fit a simple hockey-stick without constraints 

`hs = srrTMB(as.FLSR(ple4,model=segreg),s.est=T, spr0=spr0y(ple4))`

Now add the constraint that plim = Blim/B0 > 0.1, which can ensure that Blim is within plausible risk adverse biological limits relative B0

`hs.b01 = srrTMB(as.FLSR(ple4,model=segreg),s.est=T, spr0=spr0y(ple4),plim=0.1)`

Compare

`plot(FLSRs(hs=hs,hs.b01=hs.b01))`

# Licence

European Commission Joint Research Centre D.02. Released under the EUPL 1.1.

