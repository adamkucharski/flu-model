# Model of serological dynamics - uses PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015)

library(reshape2)
library(foreach)
library(doMC)
library(mvtnorm)
registerDoMC(4)  #change the 2 to your number of CPU cores
getDoParWorkers()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data (Fonville et al.)

setwd("~/Dropbox/Imperial/Fluscape_2/sero_model/")
source("load_data.R")
source("sero_functions.R")
setwd("~/Dropbox/Imperial/Fluscape_2/sero_model/")

compile.c() # Compile c code

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set up parameters

test_year=2012

theta0=c(mu=NA,sigma=NA,tau1=NA,tau2=NA,wane=NA)
theta0[["mu"]]=2
theta0[["sigma"]]=0.25
theta0[["tau1"]]=0.1
theta0[["tau2"]]=0.1
theta0[["wane"]]=0.1
theta=theta0

historyii=rbinom(inf.n, 1, 0.2) # dummy infection history

# Index variables
jj_year=c(1:test.n)[test.years==test_year]

# Specific MCMC parameters
runs=1000

nparam=length(theta)
npcov=rep(1,nparam)
cov_matrix_theta0 = diag(npcov)

varpart_prob=0.1 # Probability individual infection history resampled
#varind_prob=0.1

thetatab=matrix(NA,nrow=(runs+1),ncol=length(theta))
colnames(thetatab)=names(theta)
thetatab[1,]=theta

historytab=matrix(NA,nrow=n_part,ncol=inf.n)

# Pick plausible initial conditions
for(ii in 1:n_part){
  historytab[ii,]=setuphistIC(ii,jj_year)
}

colnames(historytab)=as.character(inf_years)

# Preallocate matrices
likelihoodtab=matrix(-Inf,nrow=(runs+1),ncol=n_part)
accepttabT=rep(NA,(runs/2))
accepttabH=rep(NA,(runs))

  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run MCMC

for (m in 1:runs){
    
    # Adaptive covariance matrix
    if(m==1){
      epsilon0=0.0001
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      #varpart_prob0=varpart_prob
    }else{
      epsilon0=min(0.1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m))
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      
      #varpart_prob0=min(1,exp(log(varpart_prob0)+(accept_rate-0.234)*0.999^m))
      #varpart_prob=varpart_prob0
    }

    # - - - - - - - - - - - - - - - -
    # Resample parameters
    
    #aTime=Sys.time() #TIMER 1

    if(m %% 2==1){
      theta_star = SampleTheta(thetatab[m,], m,cov_matrix_theta) #resample theta
      history_star = historytab
      pickA=c(1:n_part)
      
    }else{
      pickA=NULL
      while(length(pickA)<2){
        pickA=sample(n_part, round(varpart_prob*n_part)) # check that not length zero
      }
      history_star = SampleHistory(historytab,pickA) #resample history
      theta_star =thetatab[m,]
    }
    
    dmatrix=outputdmatrix(theta_star) # Arrange parameters

    # - - - - - - - - - - - - - - - -
    # LIKELIHOOD function - Can introduce joint fitting here

    lik_val=likelihoodtab[m,]
    for(ii in pickA){
      # Set history to zero after test date
      
      lik_val[ii]=estimatelik(ii,jj_year,as.numeric(history_star[ii,]),dmatrix,theta_star)
      #if(is.na(lik_val[ii])){lik_val[ii]=-Inf}
    }
    

    # - - - - - - - - - - - - - - - -
    # Metropolis Hastings step

    output_prob = ComputeProbability(sum(likelihoodtab[m,]),sum(lik_val)) 

    if(runif(1) < output_prob){
      thetatab[m+1,] = theta_star
      if(m %% 2==0){historytab = history_star} # Only change if resampled
      likelihoodtab[m+1,] = lik_val
      if(m %% 2==1){accepttabT[(m+1)/2]=1}
      
    }else{
      thetatab[m+1,] = thetatab[m,]
      likelihoodtab[m+1,] = likelihoodtab[m,]
      if(m %% 2==1){accepttabT[(m+1)/2]=0}
    }

    if(m<50){
      accept_rateT=0.234
    }else{
      accept_rateT=sum(accepttabT[1:((m+1)/2)])/((m+1)/2)
    }
    
    #Sys.time()-aTime  #TIMER 2

    if(m %% min(runs,1000) ==0){
      #print(c(m,accept_rate,round(pf_likelihoodtab[m])))
      save(likelihoodtab,thetatab,historytab,file=paste("posterior_sero_runs/outputR.RData",sep=""))
    }

}

par(mfrow=c(2,1))
par(mar = c(3,3,2,2))

# Plot profile likelihood
plot(rowSums(likelihoodtab),type="l")

# Plot histogram of boosting
hist(as.data.frame(thetatab)$mu)


