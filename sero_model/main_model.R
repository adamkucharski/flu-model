# Model of serological dynamics - uses extended PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015-)
# Main execution code

# Set up directories
setwd("~/Documents/flu-model/sero_model/")
# Set up directories
if( !file.exists("posterior_sero_runs") ){ dir.create("posterior_sero_runs") }
if( !file.exists("plot_simulations") ){ dir.create("plot_simulations") }

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
library(doRNG)
registerDoMC(4)  #change to the number of CPU cores
getDoParWorkers()

rm(list=ls(all=TRUE))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data and functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# source("load_data_Vietnam.R") # Reformat Vietnam data and save to file
# source("load_data_China.R") # Reformat FluScape and save to file

source("sero_functions.R")
source("posterior_analysis_flu.R")

compile.c() # Compile c code

flutype0="H3HN"
if(flutype0=="H3FS"){ dy1=c(2009) }
if(flutype0=="H3HN"){ dy1=c(2007:2012) }
if(flutype0=="B"){ dy1=c(2011,2012) } 
if(flutype0=="H1"){ dy1=c(2009:2011) }
load.flu.map.data()
load("datasets/spline_fn.RData") # load spline function for map **NEED TO LOAD THIS before inference run**

set.seed(5) # Set seed


# >>> Run code up to here to set everything up


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN INFERENCE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# MAIN RUNS


# Run longitudinal inference on H3 Vietnam data
flutype0="H3HN"
if(flutype0=="H3HN"){ dy1=c(2007:2012) }
foreach(kk=1:4) %dorng% {
  
  fix.param.in = c("tau1","vary.init") #if(kk<=2){c("tau1") } else{c("tau1","wane") }
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=1e3,loadseed=kk,
             flutype=flutype0,fix.param = fix.param.in, #choose parameters to fix
             fit.spline=am.spl,switch0=2,linearFn=T,vp1=0.4) 
  
}

# Run cross-sectional inference on H3 FluScape Neut data
flutype0="H3FS"
if(flutype0=="H3FS"){ dy1=c(2009) }
#for(kk in 1:4){
foreach(kk=1:4) %dorng% {
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=2e5,loadseed=kk,flutype=flutype0,
             fix.param=c("tau1","muShort","wane","sigma2","vary.init"),fit.spline=am.spl,vp1=0.4,switch0=2,linearFn=T) 
  
}

# Run cross-sectional inference on H3 FluScape HI data
flutype0="H3FS_HI"
if(flutype0=="H3FS_HI"){ dy1=c(2009) }
#for(kk in 1:4){
foreach(kk=1:4) %dorng% {
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=2e5,loadseed=kk,flutype=flutype0,
             fix.param=c("tau1","muShort","wane","sigma2","vary.init"),fit.spline=am.spl,vp1=0.4,switch0=2,linearFn=T) 
  
}

# Sensitivity analysis
# Run longitudinal inference on H3 Vietnam data -- no waning
flutype0="H3HN"
if(flutype0=="H3HN"){ dy1=c(2007:2012) }
foreach(kk=5:8) %dorng% {
  
  fix.param.in = c("tau1","muShort","wane","sigma2","vary.init") #if(kk<=2){c("tau1") } else{c("tau1","wane") }
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=1e3,loadseed=kk,
             flutype=flutype0,fix.param = fix.param.in, #choose parameters to fix
             fit.spline=am.spl,switch0=2,linearFn=T,vp1=0.4) 
  
}

# Run longitudinal inference on H3 Vietnam data -- leave one out analysis
# Start index at 100
flutype0="H3HN"
if(flutype0=="H3HN"){ dy1=c(2007:2012) }
foreach(kk=1:8) %dorng% {
  
  if(kk<=4){fix.param.in = c("tau1","vary.init") } else{fix.param.in = c("tau1","muShort","wane","sigma2","vary.init") }
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=2e5,loadseed=kk,
             flutype=flutype0,fix.param = fix.param.in, #choose parameters to fix
             fit.spline=am.spl,switch0=2,linearFn=T,vp1=0.4,leave_out_10 = T) 
  
}


# Run cross-sectional inference on H3 Vietnam data - DEPRECATED
# foreach(kk1=c(2007:2012)) %dorng% {
#   data.infer(year_test=kk1,mcmc.iterations=20,loadseed=1,flutype=flutype0,fix.param=c("tau1","wane","muShort",fit.spline=am.spl,switch0=20))
# }

# Run cross-sectional model on H3 Vietnam data - DEPRECATED
# data.infer(year_test=dy1,mcmc.iterations=1e3,loadseed=5,flutype=flutype0,fix.param=c("tau1","wane","muShort",fit.spline=am.spl,switch0=20))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# PLOT POSTERIORS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# Plot posteriors for H3 FluScape Neut data
plot.posteriors(year_test=c(2009),loadseed=1,flu.type="H3FS",fr.lim=T,plotmap = F,plot.corr = T,linearFn=T)

# Plot posteriors for H3 FluScape HI data
plot.posteriors(year_test=c(2009),loadseed=2,flu.type="H3FS_HI",fr.lim=T,plotmap = F,plot.corr = T,linearFn=T)

# Plot posteriors for longtudinal data (including attack rates - FIG 3C-D) for H3 Vietnam
for(kk in 1){
  
  plot.posteriors(year_test=c(2007:2012),loadseed=kk,flu.type="H3HN",
                  fr.lim=T,plotmap = F,plot.corr = T,linearFn=T)
  
}

# - - - - - - - - - - - - - - - - - 
# Plot posteriors for cross-sectional data - DEPRECATED

# for(kk in c(2007:2012)){
#   flutype0="H3HN"
#   dy1=kk
#   plot.posteriors(year_test=dy1,loadseed=1,flu.type=flutype0,fr.lim=F)
# }


# - - - - 
# Plot specific titre vs estimates (FIG 1) and antibody kinetics (FIG SUPP) for H3 Vietnam
plot.posterior.titres.select(loadseed=1,year_test=c(2007:2012),flu.type="H3HN",simDat=F,btstrap=200,part_pick=c(57,31,25),year_pick=c(2008:2010),linearFn=T)

# Plot specific titre vs estimates for H3 FluScape Neuts - NOT USED
#plot.posterior.titres.select(loadseed=1,year_test=c(2009),flu.type="H3FS",simDat=F,btstrap=200,part_pick=c(92,84,111),year_pick=c(2009),linearFn=T)

# Plot specific titre vs estimates for H3 FluScape HI
plot.posterior.titres.select(loadseed=1,year_test=c(2009),flu.type="H3FS_HI",simDat=F,btstrap=200,part_pick=c(27,14,111),year_pick=c(2009),linearFn=T)


# Rewind and run historical landscapes (FIG 2) 
run.historical.landscapes(loadseed=1,ymax=6.05,linearFn=T,d.step=0.2)


# - - - - - - - - - - - - - - - - - 
# SUPPLEMENTARY FIGURES

# - - - 
# TITRE DISTRIBUTIONS
# Plot H3 FluScape Neut titres
plot.posterior.titres(loadseed=1,flu.type="H3FS",simDat=F,year_test=c(2009),btstrap=100,plotRes=T,linearFn=T)

# Plot H3 FluScape HI titres
plot.posterior.titres(loadseed=1,flu.type="H3FS_HI",simDat=F,year_test=c(2009),btstrap=100,plotRes=T,linearFn=T)

# Plot Vietnam titre vs estimates
plot.posterior.titres(loadseed=1,flu.type="H3HN",simDat=F,
                      year_test=c(2007:2012),btstrap=100,plotRes=T,linearFn=T) 

# Plot Vietnam titre vs estimates -- No waning
plot.posterior.titres(loadseed=5,flu.type="H3HN",simDat=F,
                      year_test=c(2007:2012),btstrap=50,plotRes=T,linearFn=T) 

# Plot out-of-sample-performance
hold.out.analysis(loadseed=101,year_test=c(2007:2012),flu.type="H3HN",simDat=F,btstrap=1e2,plotRes=T,linearFn=T) # With waning
hold.out.analysis(loadseed=107,year_test=c(2007:2012),flu.type="H3HN",simDat=F,btstrap=1e2,plotRes=T,linearFn=T) # Without waning


# - - - 
# CONVERGENCE CHECKS
# Plot convergence for MCMC chains for H3 Neut FluScape
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3FS",year_test=c(2009),loadpick=c(1:4),fr.lim=T,linearFn=T)

# Plot convergence for MCMC chains for H3 HI FluScape
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3FS_HI",year_test=c(2009),loadpick=c(1:4),fr.lim=T,linearFn=T)

# Plot convergence for MCMC chains for H3 Vietnam
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3HN",loadpick=c(1:4),fr.lim=F,linearFn=T) 

# Plot convergence for MCMC chains for H3 Vietnam -- No waning
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3HN",loadpick=c(5:8),fr.lim=F,linearFn=T,no_wane=T)


# PLot posterior estimate for number of infections to compare H3 HI and Neut for FluScape data (Fig S9)
plot_posterior_infection_number()

# Fitted model schematic (Fig S10)
fitted_model_schematic(bstrap=1e4)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation study and inference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Generate simulated data and infer parameters -- simulation parameters are defined in sero_functions.R
# flu.type defines which dataset format (i.e. test strains, test years) the simulated data will produce

#t.check = Sys.time()

foreach(kk=1:12) %dorng% {
  
  simulation.infer(seed_i=kk,mcmc.iterations=1e4, flu.type="H3HN", strain.fix=T,
                   fit.spline=am.spl,vp1=0.4,linearFn=T) # Generate random data and run inference (strain.fix=T -> use Vietnam strains)

}

#print(Sys.time() - t.check)

# Plot convergence for MCMC chains for H3 Vietnam simulated data
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3HN",simDat=T,year_test=c(2007:2012),
                            linearFn=T,loadpick = c(1:4))


# Plot convergence for MCMC chains for H3 China simulated data - DEPRECATED
#plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3FS",simDat=T,year_test=2009,linearFn=T)


# Plot simulation study posteriors and attack rate comparisons for simulation plots (Fig 3A-B)
for(kk in 1){
  
  plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),flu.type="H3HN",year_test=c(2007:2012),plotmap=F,fr.lim=T,linearFn=T) #H3 Vietnam
  #plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),year_test=c(2009),plotmap=F,fr.lim=T) #H3 FluScape - DEPRECATED
  
}

# Plot simulation study titres against inferred model
dy1=c(2007:2012)
kk=2
plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="H3",simDat=T,year_test=c(2007:2012),btstrap=100,linearFn=T) #H3 Vietnam
#plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="H3FS",simDat=T,year_test=c(2009),btstrap=10) #H3 FluScape - DEPRECATED


# Plot test strains (Fig 3A inset)
plot_hist_strains()

# Plot estimated vs true attack rates from simulated data (Fig 3B inset)
plot.multi.true.vs.estimated(burnCut=0.25,flu.type="H3HN",simDat=T,loadpick=c(1:12))



