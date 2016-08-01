# Simulate artificial serological data for DENV/ZIKV
# Author: AJ Kucharski (2016)

library(foreach)
library(doMC)
registerDoMC(4)  #change the 2 to your number of CPU cores
getDoParWorkers()

rm(list=ls(all=TRUE))

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
per_sample0=50
seedRuns=10

for(seedK in 1:seedRuns){

  simulate_sera_data(strains=5,inf.years.sera=c(1985:2016),time.series.in=NULL,theta=theta.serology,
                     p.inf.in=0.05*c(1,1,1,1,1),sd.val.in=1.5,seedi=seedK,roundv=F,dmatrix.in=NULL,zikv.attack=0.5,per_sample=per_sample0)

  # Plot results
  #plot_simulated_sera_data(strains=5,seedi=seedK)

}

# - - - - - - - - - - - - - - - - - 
# Run MCMC inference

for(scenario in 1:4){
  foreach(seedK=c(1:seedRuns)) %dopar% {
    inference_model(seedK,strains=5,runsMCMC=1e4,scenario,per_sample=per_sample0)
  }
}

# - - - - - - - - - - - - - - - - - 
# Plot output

plot.posteriors(per_sample=per_sample0,scenario=3,seedK=1)
plot.performance(per_sample=per_sample0,age_out=15,strains,scenarioN=4,runs=seedRuns)


