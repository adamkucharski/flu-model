# Model of serological dynamics - uses PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015)

# setwd("~/Documents/flu-model/sero_model/")
# setwd("~/Dropbox/git/flu-model/sero_model")

library(reshape2)
library(mvtnorm)
library(MASS)
library(coda)

library(foreach)
library(doMC)
registerDoMC(4)  #change the 2 to your number of CPU cores
getDoParWorkers()

rm(list=ls(all=TRUE))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data and functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#source("load_data.R") # Reformat HaNam data and save to file
source("sero_functions.R")
source("posterior_analysis_flu.R")
source("sero_funcs_steven.r") # Load Flu B format


compile.c() # Compile c code

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define simple function of seed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

simulation.infer <- function(seed_i) {

  loadseed=paste("SIM_",seed_i,sep="")
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # SIMULATION MODEL
  # Generate simulated data 
  #tau1=back-boost  / tau2=suppress / disp_k=dispersion (deprecated) 
  #sigma1=long-term cross-reactivity / sigma 2=short-term CR
  
  thetaSim = c(mu=3,tau1=0,tau2=0.1,wane=1,sigma=0.3,sigma2=0.1,muShort=5,error=0.05,disp_k=1)
  npartM=50
  
  # Generate 2D map
  inf_years.in=seq(1970,2012,1)
  sim.map.in = generate.antigenic.map(inf_years.in)
  
  simulate_data(test_years=c(2009:2012), # this needs to be vector
                inf_years=inf_years.in,strain_years=seq(1970,2012,1),n_part=npartM,
                roundv=T, # Generate integer titre data
                thetastar=thetaSim,
                #antigenic.map.in = sim.map.in,
                #pmask=c("wane","sigma2"), # Specify what is included
                p.inf=0.15,seedi=loadseed)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # INFERENCE MODEL
  # Run MCMC for simulated data set

  loadseed="SIM" # ** Fix for initial testing **
  load(paste("R_datasets/Simulated_data_",loadseed,"_1.RData",sep="")) # Load simulation data for inference step that follows
  
  # Set initial theta
  theta0=c(mu=NA,tau1=NA,tau2=NA,wane=NA,sigma=NA,muShort=NA,error=NA,disp_k=1)
  theta0[["mu"]]=4 # basic boosting
  theta0[["tau1"]]=0.1 # back-boost
  theta0[["tau2"]]=0.1 # suppression via AGS
  theta0[["wane"]]=-log(0.5)/1 # short term waning - half life of /X years
  theta0[["sigma"]]=0.3 # long-term cross-reaction
  theta0[["sigma2"]]=0.1 # short-term cross-reaction
  theta0[["muShort"]]=3 # short term boosting
  theta0[["error"]]=0.1 # measurement error
  theta0[["disp_k"]]=0.1 # overdispersion (deprecated)
  theta=theta0
  vp1=0.02 #probability individual infection history resampled - this is adaptive in model
  
  define.year=c(2009:2012) # years to include in inference
  
  # browser()
  
  sim.map.in0 = generate.antigenic.map(inf_years.in)
  
  # RUN MCMC
  # Note: NEED TO RE-INITIALISE DATAFRAME IF REPEAT RUN (i.e. reload dataset above)
  run_mcmc(
    test.yr=define.year,
    test_years,
    inf_years,
    strain_years,
    n_part,
    test.list=test.listSim, # use simulated data as input
    theta=theta0,
    runs=1e5, # number of MCMC runs
    varpart_prob=vp1,
    hist.true=NULL,
    switch1=10, # ratio of infection history resamples to theta resamples. This is fixed
    #pmask=c("wane","sigma2"), # ,"map.fit" specify parameters to fix
    seedi=loadseed,
    #antigenic.map.in= sim.map.in0, # Define random initial map to fit
    linD=F)
  
}

# Run code on network

system.time(

  for(ii in 1:1){  # Use multiple seeds for simulation code
    
    fnSeedLoop(ii)
    
  } # End loop over seeds
    
)

# Do some of these over the network
library("didewin")
didewin::didewin_config_global(credentials="~/.smbcredentials",
                               home="~/dide/home",
                               temp="~/dide/tmp")

make_trees <- function(n, nspp) {
  lapply(seq_len(n), function(...) ape::rtree(nspp))
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Inference using cross-sectional vs longitudinal data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

data.infer <- function(year_test,mcmc.iterations,loadseed=1,fix.param=NULL,mushort0=5) {
  # INFERENCE MODEL
  # Run MCMC for specific data set

  load("R_datasets/HaNam_data.RData")

  # Plot simulation data vs history
  #source("simulation_plots.R")
  
  # Set initial theta
  theta0=c(mu=NA,tau1=NA,tau2=NA,wane=NA,sigma=NA,muShort=NA,error=NA,disp_k=NA)
  theta0[["mu"]]=3 # basic boosting
  theta0[["tau1"]]=0.05 # back-boost
  theta0[["tau2"]]=0.1 # suppression via AGS
  theta0[["wane"]]=-log(0.5)/0.5 + runif(1,c(-1,1)) # short term waning - half life of /X years
  theta0[["sigma"]]=0.2 # cross-reaction
  theta0[["muShort"]]=mushort0 # short term boosting
  theta0[["error"]]=0.1 # measurement error
  theta0[["disp_k"]]=0.01 # dispersion parameter - NOT CURRENTLY USED
  theta=theta0
  vp1=0.02 #probability individual infection history resampled - this is adaptive in model
  
  define.year=year_test # years to include in inference
  
  # browser()
  
  # RUN MCMC
  # Note: NEED TO RE-INITIALISE DATAFRAME IF REPEAT RUN (i.e. reload dataset above)
  run_mcmc(
    test.yr=define.year,
    test_years,
    inf_years,
    strain_years,
    n_part,
    test.list,
    theta=theta0,
    runs=mcmc.iterations, # number of MCMC runs
    varpart_prob=vp1,
    hist.true=NULL,
    switch1=10, # ratio of infection history resamples to theta resamples. This is fixed
    pmask=fix.param, #c("disp_k"), #c("wane"), #,"muShort"), # specify parameters to fix
    seedi=loadseed,
    linD=F)
}


# - - - - - - - - - - - - - - - - - 
# Run cross-sectional inference
foreach(kk1=c(2007:2012)) %dopar% {
  #if(kk==2013){kk1=c(2007:2012)}else{kk1=kk}
  data.infer(kk1,mcmc.iterations=1e2,loadseed=1,fix.param=c("disp_k","wane","muShort"),mushort0=1e-5)
}

# Multi-year inference
foreach(kk=1:4) %dopar% {
  #if(kk==2013){kk1=c(2007:2012)}else{kk1=kk}
  kk1=c(2007:2012)
  data.infer(kk1,mcmc.iterations=1e5,loadseed=kk,fix.param=c("disp_k","map.fit"))
}



# - - - - - - - - - - - - - - - - - 
# Plot posteriors
for(kk in 1:4){
  plot.posteriors(define.year=c(2007:2012),loadseed=kk)
}

for(kk in c(2007:2012)){
  plot.posteriors(define.year=kk,loadseed=1)
}

plot.compare(define.year.vec=c(2007:2012) ) #c(c(2007:2012),"2007_2008_2009_2010_2011_2012"))



