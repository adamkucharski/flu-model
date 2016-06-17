# Simulate artificial serological data for DENV/ZIKV
# Author: AJ Kucharski (2016)

# setwd("~/Documents/flu-model/sero_model/")
# setwd("~/Dropbox/git/flu-model/sero_model")
source("simulation_sera_model.R")

# - - - - - - - - - - - - - - - - - 
# Simulate serology for ZIKV/DENV

# By default cross reaction structure assumes type 1 error (given infection) and type 2 (from cross-reaction)

theta.serology=c(error=0.05,sigma=0.4) 
# error = probability negative given infection (1-sensitivity)
# sigma = probability positive from cross-reaction (1-specificity)

# inf.years.sera = range of years for outbreaks
# zikv.attack = attack rate in final year for strain 5
# p.inf.in = mean attack for each strain (drawn from lognormal with sd.val.in=1 by default)
# dmatrix.in = impose cross-reaction structure


 
simulate_sera_data(strains=5,inf.years.sera=c(1980:2015),time.series.in=NULL,theta=theta.serology,p.inf.in=c(0.05,0.05,0.05,0.1,0.1),sd.val.in=1.5,seedi=1,roundv=F,dmatrix.in=NULL,zikv.attack=0.5)

# Plot results
plot_simulated_sera_data(strains=5,seedi=1)