# Simulate artificial serological data for DENV/ZIKV
# Author: AJ Kucharski (2016)

# setwd("~/Documents/flu-model/sero_model/zika_model/")
# setwd("~/Dropbox/git/flu-model/sero_model/zika_model/")
source("simulation_sera_model.R")

# - - - - - - - - - - - - - - - - - 
# Simulate serology for ZIKV/DENV

# inf.years.sera = range of years for outbreaks
# zikv.attack = attack rate in final year for strain 5
# p.inf.in = mean attack for each strain (drawn from lognormal with sd.val.in=1 by default)
# dmatrix.in = impose cross-reaction structure

# By default cross reaction structure assumes type 1 error (given infection) and type 2 (from cross-reaction)
# error = probability negative given infection (1-sensitivity)
# sigma = probability positive from cross-reaction (1-specificity)
theta.serology=c(error=0.1,sigma=0.3) 
 
simulate_sera_data(strains=5,inf.years.sera=c(1974:2015),time.series.in=NULL,theta=theta.serology,p.inf.in=0.04*c(1,1,1,1,1),sd.val.in=1.5,seedi=1,roundv=F,dmatrix.in=NULL,zikv.attack=0.4)

# Plot results
plot_simulated_sera_data(strains=5,seedi=1)



# Run MCMC inference
source("inference_sera_model.R")