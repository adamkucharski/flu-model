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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define simulation model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

simulation.infer <- function(seed_i,mcmc.iterations=1e3) {

  loadseed=paste("SIM_",seed_i,sep="")
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # SIMULATION MODEL
  # Generate simulated data 
  #tau1=back-boost  / tau2=suppress / disp_k=dispersion (deprecated) 
  #sigma1=long-term cross-reactivity / sigma 2=short-term CR
  
  thetaSim = c(mu=3,tau1=0.02,tau2=0.1,wane=1,sigma=0.3,muShort=5,error=0.05,disp_k=1,sigma2=0.1)
  npartM=70
  
  # Generate 2D map
  define.year=c(2007:2012) # test years
  inf_years.in=seq(1968,2012,1)
  sim.map.in = generate.antigenic.map(inf_years.in)
  
  attack.yr=read.csv("datasets/sim_attack.csv")[,1]
  
  #attack.yr=rlnorm(inf_years.in,meanlog=log(0.15)-1^2/2,sdlog=0.5)
  #write.csv(attack.yr,"datasets/sim_attack.csv")
  
  simulate_data(test_years=define.year, # this needs to be vector
                inf_years=inf_years.in,strain_years=seq(1968,2012,1),n_part=npartM,
                roundv=T, # Generate integer titre data
                thetastar=thetaSim,
                #antigenic.map.in = sim.map.in,
                #pmask=c("wane","sigma2"), # Specify what is included
                p.inf=attack.yr,seedi=loadseed)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # INFERENCE MODEL
  # Run MCMC for simulated data set

  load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep="")) # Load simulation data for inference step that follows
  
  # Set initial theta
  theta0=c(mu=NA,tau1=NA,tau2=NA,wane=NA,sigma=NA,muShort=NA,error=NA,disp_k=1,sigma2=0.1)
  theta0[["mu"]]=2 # basic boosting
  theta0[["tau1"]]=0.1 # back-boost
  theta0[["tau2"]]=0.1 # suppression via AGS
  theta0[["wane"]]=1  # -log(0.5)/1 # short term waning - half life of /X years
  theta0[["sigma"]]=0.3 # long-term cross-reaction
  theta0[["sigma2"]]=0.1 # short-term cross-reaction
  theta0[["muShort"]]=3 # short term boosting
  theta0[["error"]]=0.1 # measurement error
  theta0[["disp_k"]]=0.1 # overdispersion (deprecated)
  theta=theta0
  vp1=0.02 #probability individual infection history resampled - this is adaptive in model

  # browser()
  sim.map.in0 = 0.3*(cbind(inf_years,inf_years)-min(inf_years)) #generate.antigenic.map(inf_years.in) # Define uniform initial map to fit
  
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
    runs=mcmc.iterations, # number of MCMC runs
    varpart_prob=vp1,
    hist.true=NULL,
    switch1=10, # ratio of infection history resamples to theta resamples. This is fixed
    pmask=c("disp_k"), # ,"map.fit" specify parameters to fix
    seedi=loadseed,
    #antigenic.map.in= sim.map.in0, # Define random initial map to fit
    linD=F)
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run code on network (SR code)
fn.network<-function(){
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
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Inference using cross-sectional vs longitudinal data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

data.infer <- function(year_test,mcmc.iterations=1e3,loadseed=1,flutype="H3",fix.param=NULL , fit.spline=NULL) {
  # INFERENCE MODEL
  # Run MCMC for specific data set
  if(flutype=="H3"){
    load("R_datasets/HaNam_data.RData")
    #am.spl<-load.flu.map.data() # define spline from antigenic map data
  }else{
    load("R_datasets/Fluscape_data_List.RData")
  }

  # Plot simulation data vs history
  #source("simulation_plots.R")
  
  # Set initial theta
  theta0=c(mu=NA,tau1=NA,tau2=NA,wane=NA,sigma=NA,muShort=NA,error=NA,disp_k=NA,sigma2=NA)
  theta0[["mu"]]=3 # basic boosting
  theta0[["tau1"]]=0.05 # back-boost
  theta0[["tau2"]]=0.1 # suppression via AGS
  theta0[["wane"]]=-log(0.5)/0.5 + if(sum(fix.param=="wane")==0){runif(1,c(-1,1))}else{0} # short term waning - half life of /X years -- add noise to IC if fitting
  theta0[["sigma"]]=0.3 + if(sum(fix.param=="wane")==0){0.1*runif(1,c(-1,1))}else{0} # cross-reaction
  theta0[["sigma2"]]=0.1 + if(sum(fix.param=="wane")==0){0.1*runif(1,c(-1,1))}else{0} # short-term cross-reaction
  theta0[["muShort"]]=6 + if(sum(fix.param=="wane")==0){2*runif(1,c(-1,1))}else{0} # short term boosting
  theta0[["error"]]=0.05 # measurement error
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
    seedi=paste(loadseed,"_",flutype,sep=""), # record output
    antigenic.map.in=inf_years, # define specific map structure (or initial structure if fitting)
    am.spline=fit.spline, # decide whether to fit antigenic map along "am.spl" spline function
    linD=F)
  
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN INFERENCE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Run cross-sectional inference
foreach(kk1=c(2011:2012)) %dopar% {
  #if(kk==2013){kk1=c(2007:2012)}else{kk1=kk}
  flutype0="B"
  if(flutype0=="H3"){ dy1=c(2007:2012) }else{ dy1=c(2011,2012) }
  
  data.infer(year_test=kk1,mcmc.iterations=1e5,loadseed=1,flutype=flutype0,fix.param=c("disp_k","wane","muShort"))
}

# - - - - - - - - - - - - - - - - - 
# Run longtudinal inference on H3 or B data

#load.flu.map.data() # define spline from H3 antigenic map data
load("datasets/spline_fn.RData") # load spline function for map

for(kk in 1){
#foreach(kk=1:4) %dopar% {
  #if(kk==2013){kk1=c(2007:2012)}else{kk1=kk}
  flutype0="H3"
  if(flutype0=="H3"){ dy1=c(2007:2012) }else{ dy1=c(2011,2012) }
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=2e5,loadseed=kk,flutype=flutype0,fix.param=c("disp_k"),fit.spline=1) #,"map.fit"
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# PLOT POSTERIORS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot posteriors for longtudinal data (FIG 3)
for(kk in 1:4){
  
  flutype0="H3"
  if(flutype0=="H3"){ dy1=c(2007:2012) }else{ dy1=c(2011,2012) }
  
  plot.posteriors(year_test=dy1,loadseed=kk,flutype=flutype0,f.lim=T,plotmap = T)
  
}

# - - - - - - - - - - - - - - - - - 
# Plot posteriors for cross-sectional data

for(kk in c(2011:2012)){
  flutype0="B"
  dy1=kk
  plot.posteriors(year_test=dy1,loadseed=1,flutype=flutype0,f.lim=T)
}

# plot.compare(define.year.vec=c(2007:2012) ) #c(c(2007:2012),"2007_2008_2009_2010_2011_2012"))
# plot.posteriors(simDat=T,loadseed="SIM",year_test=c(2007:2012),plotmap=T)

# - - - - - - - - - - - - - - - - - 
# Plot titre vs estimates

flutype="H3"
if(flutype=="H3"){ dy1=c(2007:2012) }else{ dy1=c(2011,2012) }
plot.posterior.titres(loadseed=1,flu.type=flutype,simDat=F,year_test=dy1,btstrap=10)

# - - - - 
# Plot specific titre vs estimates (FIG 1) and antibody kinetics (FIG 2)

plot.posterior.titres.select(loadseed=1,year_test=c(2007:2012),flu.type="H3",simDat=F,btstrap=50,part_pick=c(31,57,25),year_pick=c(2008:2010))
plot.antibody.changes()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation results
foreach(kk=1:4) %dopar% {
  simulation.infer(seed_i=kk,mcmc.iterations=5e5) # Run inference
}

for(kk in 1:4){
  plot.posteriors(simDat=T,loadseed=paste("SIM_",kk,sep=""),year_test=c(2007:2012),plotmap=F,f.lim=T)
}

dy1=c(2007:2012)
kk=1
plot.posterior.titres(loadseed=paste("SIM_",kk,sep=""),flu.type="",simDat=T,year_test=dy1,btstrap=10)

