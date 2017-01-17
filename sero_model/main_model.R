# Model of serological dynamics - uses extended PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015-)
# Main execution code

setwd("~/Documents/flu-model/sero_model/")
# setwd("~/Dropbox/git/flu-model/sero_model")

library(reshape2)
library(mvtnorm)
library(MASS)
library(coda)
library(RColorBrewer)
library(magrittr)
library(plot3D)
library(colorspace)
 
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

flutype0="H3HN"
if(flutype0=="H3FS"){ dy1=c(2009) }
if(flutype0=="H3HN"){ dy1=c(2007:2012) }
if(flutype0=="B"){ dy1=c(2011,2012) } 
if(flutype0=="H1"){ dy1=c(2009:2011) }
load.flu.map.data()
load("datasets/spline_fn.RData") # load spline function for map **NEED TO LOAD THIS before inference run**


# >>> Run code up to here to set everything up


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN INFERENCE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Run cross-sectional inference on H3 HaNam data
foreach(kk1=c(2007:2012)) %dopar% {
  data.infer(year_test=kk1,mcmc.iterations=20,loadseed=1,flutype=flutype0,fix.param=c("tau1","wane","muShort",fit.spline=am.spl,switch0=20))
}

# >>> TESTING RUNS

# Run longitudinal inference on H3 HaNam data
foreach(kk=1:4) %dopar% {

  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=4e5,loadseed=kk,
             flutype=flutype0,fix.param=c("tau1","vary.init"),
             fit.spline=am.spl,switch0=2,linearFn=T,vp1=0.4) #,"map.fit"

}

# Run cross-sectional inference on H3 FluScape data
flutype0="H3FS"
if(flutype0=="H3FS"){ dy1=c(2009) }
#for(kk in 1:4){
foreach(kk=1:4) %dopar% {
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=3e5,loadseed=kk,flutype=flutype0,
             fix.param=c("tau1","muShort","wane","sigma2","vary.init"),fit.spline=am.spl,switch0=2,linearFn=T) #,"map.fit"
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# PLOT POSTERIORS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Plot posteriors for longtudinal data (including attack rates - FIG 3) for H3 Vietnam
for(kk in 1:3){
  
  plot.posteriors(year_test=c(2007:2012),loadseed=kk,flu.type="H3HN",
                  fr.lim=T,plotmap = F,plot.corr = T,linearFn=T)
  
}

# Plot posteriors for H3 FluScape data
plot.posteriors(year_test=c(2009),loadseed=1,flu.type="H3FS",fr.lim=T,plotmap = F,plot.corr = F,linearFn=T)


# - - - - - - - - - - - - - - - - - 
# Plot posteriors for cross-sectional data

for(kk in c(2007:2012)){
  flutype0="H3HN"
  dy1=kk
  plot.posteriors(year_test=dy1,loadseed=1,flu.type=flutype0,fr.lim=F)
}

# plot.compare(define.year.vec=c(2007:2012) ) #c(c(2007:2012),"2007_2008_2009_2010_2011_2012"))
# plot.posteriors(simDat=T,loadseed="SIM",year_test=c(2007:2012),plotmap=T)

# - - - - 
# Plot specific titre vs estimates (FIG 1) and antibody kinetics (FIG 2) for H3 Vietnam
plot.posterior.titres.select(loadseed=1,year_test=c(2007:2012),flu.type="H3HN",simDat=F,btstrap=200,part_pick=c(57,31,25),year_pick=c(2008:2010),linearFn=T)


# Plot specific titre vs estimates for H3 FluScape
plot.posterior.titres.select(loadseed=1,year_test=c(2009),flu.type="H3FS",simDat=F,btstrap=200,part_pick=c(92,84,111),year_pick=c(2009),linearFn=T)


# Plot specific antibody kinetics for H3 Vietnam -- DEPRECATED
#plot.antibody.changes(btstrap=200)

# Rewind and run historical landscapes (FIG 2) -- NEED TO UPDATE FOR LINEAR/EXP FUNCTIONS
run.historical.landscapes(loadseed=1,ymax=6.05,linearFn=T,d.step=0.25)


# - - - - - - - - - - - - - - - - - 
# SUPPLEMENTARY FIGURES
# Plot titre vs estimates

plot.posterior.titres(loadseed=1,flu.type="H3HN",simDat=F,
                      year_test=c(2007:2012),btstrap=5,plotRes=T,linearFn=T) # Note linear function


#H3 FluScape titres
plot.posterior.titres(loadseed=1,flu.type="H3FS",simDat=F,year_test=c(2009),btstrap=50,plotRes=T,linearFn=T)


#H1 titres -- DEPRECATED
#plot.posterior.titres(loadseed=1,flu.type="H1",simDat=F,year_test=c(2009:2011),btstrap=2)


# >>> IMPORTANT FOR TESTING RUNS
# Plot convergence for MCMC chains for H3 Vietnam

plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3HN",loadpick=c(1:3),fr.lim=F,linearFn=T,runsPOST=3e5)



# Plot convergence for MCMC chains for H3 FluScape
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3FS",year_test=c(2009),loadpick=c(1:4),fr.lim=T,linearFn=T)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation study and inference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Generate simulated data and infer parameters -- simulation parameters are defined in sero_functions.R
# flu.type defines which dataset format (i.e. test strains, test years) the simulated data will produce

foreach(kk=1:4) %dopar% {
  
  simulation.infer(seed_i=kk,mcmc.iterations=1e5, flu.type="H3HN", strain.fix=T,
                   fit.spline=am.spl,vp1=0.4,linearFn=T) # Generate random data and run inference (strain.fix=T -> use Ha Nam strains)

}

# Plot convergence for MCMC chains for H3 Vietnam simulated data
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3HN",simDat=T,year_test=c(2007:2012),
                            linearFn=T,loadpick = c(1:4))


# Plot convergence for MCMC chains for H3 China simulated data
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3FS",simDat=T,year_test=2009)


# Plot simulation study posteriors and attack rate comparisons for simulation plots
for(kk in 1:3){
  
  plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),flu.type="H3HN",year_test=c(2007:2012),plotmap=F,fr.lim=T,linearFn=T) #H3 Vietnam
  #plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),year_test=c(2009),plotmap=F,fr.lim=T) #H3 FluScape
  
}

# Plot simulation study titres against inferred model
dy1=c(2007:2012)
kk=1
plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="H3",simDat=T,year_test=c(2007:2012),btstrap=10,linearFn=T) #H3 Vietnam
#plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="H3FS",simDat=T,year_test=c(2009),btstrap=10) #H3 FluScape

