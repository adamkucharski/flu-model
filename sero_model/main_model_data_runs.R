# Model of serological dynamics - uses extended PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015-)
# Main execution code

# Set up directories
setwd("~/Documents/flu-model/sero_model/")
# setwd("~/Dropbox/git/flu-model/sero_model")

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
registerDoMC(cores = detectCores())  #change the 2 to your number of CPU cores
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
# RUN INFERENCE ON DATA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# Run longitudinal inference on H3 Vietnam data
# Index from 100
foreach(kk=101:112) %dorng% {
  
  #fix.param.in = c("tau1","vary.init") 
  if(kk<=106){fix.param.in = c("tau1","vary.init") } else{fix.param.in = c("tau1","wane","vary.init","sigma2","muShort") }
  # Fits to spline if am.spl is defined
  data.infer(year_test=dy1,mcmc.iterations=2e5,loadseed=kk,
             flutype=flutype0,fix.param = fix.param.in, #choose parameters to fix
             fit.spline=am.spl,switch0=2,linearFn=T,vp1=0.4,leave_out_10 = T) 
  
}



