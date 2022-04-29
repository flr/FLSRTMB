# test-FLSR.R - DESC
# /test-FLSR.R

# Copyright Henning Winker (JRC) & Iago MOSQUEIRA (WMR), 2021
# Authors:  Henning Winker (JRC) <henning.winker@ec.europa.eu>
#           Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# TEST
library(ggplotFL)
library(FLCore)
library(FLSRTMB)

# Example data are loaded for North Sea from FLCore

data(ple4)

#The Beverton and Holt model is parameterized as a function of steepness (s) and unfished spawning potential ratio SPR0. To create the FLSR input object, the model=bevholtSV is selected.

#Beverton-Holt Model

sr <- as.FLSR(ple4,model=bevholtSV)

#The function spr0y computes annual spr0. A good starting point can be the average

spr0 <- yearMeans(spr0y(ple4))

# Fit a bevholt model without constraints

bh = srrTMB(sr,s.est=T, spr0=spr0)

plot(bh)

#Note the output params are a and b by default but can be changed to R0 and s.

params(srrTMB(sr,s.est=T, spr0=spr0,report.sR0 = TRUE))

ri = srrTMB(as.FLSR(ple4,model=rickerSV),s.est=T, spr0=spr0)


#For many stock spr0 time-varying due to variations in weight-at-age, M-at-age or maturity-at-age plot(spr0y(ple4))+ylab("spr0")

#To account for reduced or increased unfished stock sizes given these variation, spr0y() can be directly inputted

bh.y = srrTMB(sr,s.est=T, spr0=spr0y(ple4))

#Note that by default a,b are computed from the average spr0 across all years. However, if the analyst suspects a non reversable systematic change, the spr0 reference can also be taken from the most recent years

bh.y3 = srrTMB(sr,s.est=T, spr0=spr0y(ple4),nyears=10)

# Compare fits

plot(FLSRs(spr0=bh,spr0y=bh.y,spr0y3 = bh.y3))+theme(legend.position="right")

# Option to estimate steepness with prior, e.g. from meta-analysis, such as FishLife (Thorson, 2020) This requires to specify s with the prior mean and s.logitsd is the sd on logit scale

bh.prior = srrTMB(as.FLSR(ple4,model=bevholtSV),s.est=T,s=0.7,s.logitsd=0.4, spr0=spr0y(ple4))

plot(FLSRs(spr0=bh,spr0y=bh.y,s.prior = bh.prior))+theme(legend.position="right")

# It is also possible to fix steepness

s = c(0.75,0.8,0.85,0.9,0.95) 

bhs <- FLSRs(sapply(s, function(x) { return( srrTMB(sr,s=x,s.est=F, spr0=spr0y(ple4)))})) 

bhs@names = c(paste("s =",round(s,3))) 

plot(bhs)+theme(legend.position="right")

#-------------------------------------
# Hockey Stick (segmented regression)
#-------------------------------------

# First fit a simple hockey-stick without constraints

hs = srrTMB(as.FLSR(ple4,model=segreg),s.est=T, spr0=spr0y(ple4))

# Now add the constraint that plim = Blim/B0 > 0.1, which can ensure that Blim is within plausible risk adverse biological limits relative B0

hs.b01 = srrTMB(as.FLSR(ple4,model=segreg),s.est=T, spr0=spr0y(ple4),plim=0.1)

# Compare

plot(FLSRs(hs=hs,hs.b01=hs.b01))+theme(legend.position="right")

# Compute Blim/B0 
ab = (params(hs.b01))
R0 = ab[1]*ab[2]
B0 = R0*yearMeans(spr0y(ple4))
Blim = ab[2]
Blim/B0


