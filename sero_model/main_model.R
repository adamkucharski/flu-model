# Model of serological dynamics - uses PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015)

# setwd("~/Documents/flu-model/sero_model/")
# setwd("./sero_model")

library(reshape2)
library(foreach)
library(doMC)
library(mvtnorm)
registerDoMC(4)  #change the 2 to your number of CPU cores
getDoParWorkers()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data and functions (Fonville et al.)
#source("load_data.R")
source("sero_functions.R")

compile.c() # Compile c code
# logabs2(seq(-5, 5, by=2))


# - - - - - - - - - - - - -
# Generate simulated data
thetaSim=c(mu=4,sigma=0.3,tau1=0.1,tau2=0.1,wane=0.1)
npartM=50
simulate_data(test_years=seq(2010,2010),
              inf_years=seq(1990,2010,1),
              strain_years=seq(1990,2010,2),
              n_part=npartM,thetastar=thetaSim,p.inf=0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# INFERENCE MODEL
# Run MCMC for specific data set

#load("R_datasets/HaNam_data.RData")
load("R_datasets/Simulated_data.RData")

# Plot simulation data vs history
source("simulation_plots.R")

# Set initial theta
theta0=c(mu=NA,sigma=NA,tau1=NA,tau2=NA,wane=NA)
theta0[["mu"]]=4
theta0[["sigma"]]=0.5
theta0[["tau1"]]=0.1
theta0[["tau2"]]=0.1
theta0[["wane"]]=0.1
theta=theta0
vp1=0.1 #probability individual infection history resampled

run_mcmc(test.yr=2010,runs=10000,hist.true=NULL,varpart_prob=vp1,test_years,inf_years,strain_years,n_part,test.list,theta0)
#historytabSim

# Plot posteriors and compare to simulation
source("simulation_diagnostics.R",local=TRUE)

