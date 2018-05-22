# Model of serological dynamics - uses extended PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015-)
# Main execution code

# Set up directories
setwd("~/Documents/flu-model/sero_model/")
# setwd("~/Dropbox/git/flu-model/sero_model")

# Set up directories
if( !file.exists("posterior_sero_runsA") ){ dir.create("posterior_sero_runsA") }
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
registerDoMC(4)  #change the 2 to your number of CPU cores
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

set.seed(5) # Set seed for reproducibility


# >>> Run code up to here to set everything up


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN INFERENCE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


foreach(kk=1:10) %dorng% {
  
  simulation.infer(seed_i=kk,mcmc.iterations=1e2, flu.type="H3HN", strain.fix=T,
                   fit.spline=am.spl,vp1=0.4,linearFn=T) # Generate random data and run inference (strain.fix=T -> use Vietnam strains)
  
}


# Plot convergence for MCMC chains for H3 Vietnam simulated data
plot.multi.chain.posteriors(burnCut=0.25,flu.type="H3HN",simDat=T,year_test=c(2007:2012),
                            linearFn=T,loadpick = c(1:4))
