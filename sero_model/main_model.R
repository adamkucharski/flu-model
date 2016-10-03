# Model of serological dynamics - uses PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015)

 setwd("~/Documents/flu-model/sero_model/")
# setwd("~/Dropbox/git/flu-model/sero_model")

library(reshape2)
library(mvtnorm)
library(MASS)
library(coda)
library(RColorBrewer)

library(foreach)
library(doMC)
registerDoMC(4)  #change the 2 to your number of CPU cores
getDoParWorkers()

rm(list=ls(all=TRUE))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data and functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#source("load_data.R") # Reformat HaNam data and save to file
# make_fluscape_rdata(pathfssvn="~/fluscape/trunk/") # Reformat Flu B data and save to file
source("sero_functions.R")
source("posterior_analysis_flu.R")
source("sero_funcs_steven.r") # Load Flu B format

compile.c() # Compile c code


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN INFERENCE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Run cross-sectional inference
foreach(kk1=c(2011:2012)) %dopar% {
  #if(kk==2013){kk1=c(2007:2012)}else{kk1=kk}
  flutype0="B"
  if(flutype0=="H3"){ dy1=c(2007:2012) }else{ dy1=c(2011,2012) }
  data.infer(year_test=kk1,mcmc.iterations=1e6,loadseed=1,flutype=flutype0,fix.param=c("disp_k","wane","muShort"))
}

# - - - - - - - - - - - - - - - - - 
# Run longtudinal inference on H3 or B data

flutype0="H1"
if(flutype0=="H3"){ dy1=c(2007:2012) }
if(flutype0=="B"){ dy1=c(2011,2012) } 
if(flutype0=="H1"){ dy1=c(2009:2011) } 

#load.flu.map.data() # define spline from H3 antigenic map data
load("datasets/spline_fn.RData") # load spline function for map **NEED TO LOAD THIS before inference run**
for(kk in 1:4){
#foreach(kk=1:4) %dopar% {
  #if(kk==2013){kk1=c(2007:2012)}else{kk1=kk}
  # Fits to spline if am.spl is defined [** AK: need to tidy this up **]
  data.infer(year_test=dy1,mcmc.iterations=1e4,loadseed=kk,flutype=flutype0,fix.param=c("disp_k","error"),fit.spline=NULL) #,"map.fit"
  
}

data.infer(year_test=dy1,mcmc.iterations=1e4,loadseed=kk,flutype=flutype0,fix.param=c("disp_k","error","muShort"),fit.spline=NULL) #,"map.fit"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# PLOT POSTERIORS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Plot posteriors for longtudinal data (including attack rates - FIG 3)
for(kk in 1){
  
  plot.posteriors(year_test=dy1,loadseed=kk,flu.type=flutype0,f.lim=T,plotmap = F)
  
}

# - - - - - - - - - - - - - - - - - 
# Plot posteriors for cross-sectional data

for(kk in c(2011:2012)){
  flutype0="H1"
  dy1=kk
  plot.posteriors(year_test=dy1,loadseed=1,flu.type=flutype0,f.lim=)
}

# plot.compare(define.year.vec=c(2007:2012) ) #c(c(2007:2012),"2007_2008_2009_2010_2011_2012"))
# plot.posteriors(simDat=T,loadseed="SIM",year_test=c(2007:2012),plotmap=T)

# - - - - 
# Plot specific titre vs estimates (FIG 1) and antibody kinetics (FIG 2)

plot.posterior.titres.select(loadseed=1,year_test=c(2007:2012),flu.type="H3",simDat=F,btstrap=50,part_pick=c(31,57,25),year_pick=c(2008:2010))
plot.antibody.changes(btstrap=200)

# - - - - - - - - - - - - - - - - - 
# SUPPLEMENTARY FIGURES
# Plot titre vs estimates

load("datasets/spline_fn.RData") # load spline function for map **NEED TO LOAD THIS before next run**
plot.posterior.titres(loadseed=1,flu.type="H3",simDat=F,year_test=c(2007:2012),btstrap=10)

#H1 titres
plot.posterior.titres(loadseed=1,flu.type="H1",simDat=F,year_test=c(2009:2011),btstrap=2)


# Plot convergence for MCMC chains
plot.multi.chain.posteriors(burnCut=0.25)





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation results
foreach(kk=1:4) %dopar% {
  simulation.infer(seed_i=kk,mcmc.iterations=5e5) # Run inference
}

# Plot posteriors and attack rate comparisons for simulation plots
for(kk in 1:4){
  plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),year_test=c(2007:2012),plotmap=F,f.lim=F)
}

dy1=c(2007:2012)
kk=1
plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="",simDat=T,year_test=dy1,btstrap=10)

