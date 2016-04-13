#Functions

# - - - - - - - - - - - - - - - -
# Set initial condition (for infection history) as infection if titre >=X

setuphistIC<-function(ii,jj,inf.n,test.list,testyear_index, inf_years){ # ii=participant | jj=test year
  
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  spyear=unique(as.numeric(test.jj[3,])) # year of samples taken
  
  hist0=rep(0,inf.n)   
  #hist0[sample(c(1:inf.n),round(0.1*inf.n))]=1
  
  # Check test data available
  if(length(test.jj[,1])>1){
    
    # Set up test strains
    titredat=as.numeric(test.jj[2,]) # Define titre data
    maxt=(titredat==max(titredat))
    
    # set max titre strain to infection=1 in history
    #if(sum(maxt)>1){
    #  hist0[inf_years==spyear[sample(c(1:length(maxt))[maxt],1)]]=1
    #}else{
    #  
    #  hist0[(inf_years==spyear[titredat==max(titredat)])]=1
    #}
    
    # Use simple cutoff for titres
    for(i in 1:length(spyear)){
      if(max(titredat[(as.numeric(test.jj[3,])==spyear[i])])>=4 & runif(1)>0.5 ){
        hist0[(inf_years==spyear[i])]=1
      }
    }
    
    #hist0=sample(c(0,1),inf.n,replace=T,prob=c(0.9,0.1)) # Constrain max number of infections to 10% attack rate?
    
  }
  pos.hist=(hist0>0)
  
  min.range=max(1,testyear_index[1]-20) # Add infection within past 20 years
  if(sum(hist0[1:min(testyear_index)])==0){hist0[sample(min.range:testyear_index[1],1)]=1} # Make sure at least one infection previous to test year
  hist0
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Likelihood given infection history and parameters (done in C)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Functions to set up parameters for model

outputdmatrix<-function(theta,inf_years,locmat=NULL){
  if (is.null(locmat)) {
    (dmatrix=sapply(inf_years,function(x){exp(-theta[["sigma"]]*abs(inf_years-x))})) # note that second entry is actually sample year
  } else {
    stop("Non-null locmat not yet implemented in outputdmatrix")
    # Up to here on the distance matrix
    (dmatrix=sapply(inf_years,function(x){exp(-theta[["sigma"]]*abs(inf_years-x))}))
  }
  dmatrix
}

# Compile c code
compile.c<-function(){
  require("Rcpp")
  setwd("./c_code")
  #system("R CMD SHLIB c_model2.c")
  system("R CMD SHLIB c_model2_sr.c")
  #dyn.load("c_model2.so")
  dyn.load("c_model2_sr.so") # Note edit to remove ./ for cluster runs
  # sourceCpp("./cpp_steven.cpp")
  setwd("..")
}

# - - - - - - - - - - - - - - - -
# Define expected titre function

func1 <- function(x,titredat,dd,theta,testyear_index) {
  if (!is.numeric(x)){stop("argument x must be numeric")}
  out <- .C("c_model2_sr",
            n=as.integer(length(x)),
            itot=as.integer(sum(x)),
            nsample=as.integer(length(titredat)),
            x=as.double(x),
            x1=as.double(rep(0,length(x))),
            titre=as.double(titredat),
            titrepred=as.double(rep(0,length(titredat))),
            dd=as.double(dd),
            ntheta=as.integer(length(theta)),
            theta=as.double(theta),
            inputtestyr=as.integer(testyear_index)
  )
  # browser()
  # out$titrepred - out2$titrepred
  return(out$titrepred)
}


# - - - - - - - - - - - - - - - -
# Calculate likelihood for given participant and test year

estimatelik<-function(ii,jj,historyii,dmatrix,theta_star,test.list,testyearI){ # ii=participant | jj=test year
  
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  # Check test data available
  if(length(test.jj[,1])==1){0}else{
    
    # Set up test strains
    test.part=as.numeric(test.jj[4,]) # index of sample strains data available for
    titredat=test.jj[2,] # Define titre data
    
    d.ij=dmatrix[test.part,] # Define cross-immunity matrix for sample strain
    d_vector=melt(t(d.ij))$value
    
    #if(ii==62){historyii[10]=1}

    expect=func1(historyii,titredat,d_vector,theta_star,testyearI) # Output expectation

    # Calculate likelihood - ** have added summation for k>8 **
    largett=(titredat>=8)

    sum(dpois(as.numeric(titredat[!largett]), expect[!largett], log = TRUE))+
      sum(ppois(8, lambda=expect[largett], lower=FALSE,log=TRUE))
  }
  
}


# - - - - - - - - - - - - - - - -
# Simulation infection history data

simulate_data<-function(test_years,historytabPost=NULL, inf_years,strain_years,n_part=20,thetastar=theta0,p.inf=0.2,seedi=1){ # ii=participant | jj=test year
  
  # Variables needed: test_years,inf_years,strain_years,n_part
  #strain_years=seq(1968,2010,4)
  
  # Set year of birth
  age.yr=sample(1:80,n_part,replace = TRUE)
  
  test.n=length(test_years)
  inf.n=length(inf_years)
  nstrains=length(strain_years)
  sample.index=strain_years-min(strain_years)+1
  
  # Check inputs are correct
  if(sum(max(test_years)==inf_years)==0){
    print("need infection years >= test years")
    return
  }
  
  dmatrix=outputdmatrix(thetastar,inf_years)
  
  #Set per year incidence, to create correlation between participants
  log.sd=1
  attack.yr=rlnorm(inf.n,meanlog=log(p.inf)-log.sd^2/2,sdlog=log.sd)
  
  # Simulate random infection history for each participant
  if(is.null(historytabPost)){
    historytabSim=NULL
    for(ii in 1:n_part){
      hist0=(runif(inf.n)<attack.yr)+0
      alive=((max(test_years)-age.yr[ii])<=inf_years)
      historytabSim=rbind(historytabSim,hist0*alive)
    }
  }else{
    historytabSim=historytabPost
  }
  
  
  # Simulate titres for each participant
  # ** NEED TO INCLUDE TEST YEAR IN FUNCTION IF DECAY ADDED **
  
  test.list=list()
  
  for(ii in 1:n_part){
    
    subjectn=ii
    i.list=list()
    historyii=historytabSim[ii,]
    
    for(jj in 1:test.n){
      
      d.ij=dmatrix[sample.index,] # Define cross-immunity matrix for sample strain
      d_vector=melt(t(d.ij))$value
      testyr=test_years[jj]
      testyearI=c(1:inf.n)[inf_years==testyr]
      
      expect=func1(historyii,sample.index,d_vector,thetastar,testyearI) # Output expectation
      
      #titredat=sapply(expect,function(x){rpois(1,x)}) # Generate titre
      titredat=round(expect)
      titredat=sapply(titredat,function(x){min(x,8)})
      
      
      i.list[[jj]]=rbind(test.year=rep(testyr,nstrains),
                         titredat,
                         strain_years,
                         sample.index
      )
    }
    #i.list[[1]][2,]
    #
    
    test.list[[ii]]=i.list
  }
  # Export data
  #browser()
  if(is.null(historytabPost)){
    save(test_years,inf_years,strain_years,n_part,test.list,age.yr,historytabSim,file=paste("R_datasets/Simulated_data_",seedi,".RData",sep=""))
  }else{
    save(test_years,inf_years,strain_years,n_part,test.list,age.yr,file=paste("R_datasets/Simulated_dataPost_",seedi,".RData",sep=""))
  }
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Resample infection history - included ageA table in case needed later

SampleHistory<-function(historyA,pick,inf.n,ageA,inf_years){
  
  infvector=c(1:inf.n)
  infvector2=rev(infvector)
  
  for(ii in pick){
    #ls_pick=foreach(ii=(1:length(pick))) %dopar% {  # Parallel loop - slower to farm out
    rand1=runif(1)
    x=historyA[ii,]
    
    # Remove infection
    if(rand1<1/3){
      infectID=infvector[(as.numeric(x)>0)]
      if(length(infectID)>0){
        x[sample(c(infectID,infectID),1)]=0
      }
    }
    
    # Add new infection
    if(rand1>1/3 & rand1<2/3){
      ninfecID=infvector[(as.numeric(x)==0)]
      if(length(ninfecID)>0){
        x[sample(c(ninfecID,ninfecID),1)]=1
      }
    }
    
    # Move infection
    if(rand1>2/3){
      infectID=infvector[(as.numeric(x)>0)]
      ninfecID=infvector[(as.numeric(x)==0)]
      if(length(infectID)>0 & length(ninfecID)>0){
        x[sample(c(infectID,infectID),1)]=0
        x[sample(c(ninfecID,ninfecID),1)]=1
      }
    }
    
    # Add prior on birth year - exponentially less likely to update if infections outside
    #if(inf.n>ageA[ii]){
    #  a1=0.01*exp(1)*exp(-sum(x[1:(inf.n-ageA[ii])])) # EDIT infvector2 tweak this parameter to penalise more/less
    #  if( a1 > runif(1) ){
    #    historyA[ii,]=x
    #  }
    #}
    
    historyA[ii,]=x
    
  } # end loop over individuals
  historyA
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Resample age - add 1, 0, -1 with equal probability

SampleAge<-function(pick,ageA){
  
  b1=sapply(ageA[pick],function(x){x+sample(c(-1:1),1)})
  ageA[pick]=b1
  ageA
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Convert infection history to binary - not currently used

convert_binary <- function(x){sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

ComputeProbability<-function(marg_likelihood,marg_likelihood_star){
  
  # uniform priors
  p_theta_star = 1; p_theta = 1
  
  # probability symmetic
  q_theta_given_theta_star = 1; q_theta_star_given_theta = 1
  
  val = exp((marg_likelihood_star-marg_likelihood))*(p_theta_star/p_theta)*(q_theta_given_theta_star/q_theta_star_given_theta) 
  min(val, 1)
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

SampleTheta<-function(theta_in,m,covartheta){
  
  # sample new parameters from nearby: 
  theta_star = as.numeric(exp(rmvnorm(1,log(theta_in), covartheta)))
  names(theta_star)=names(theta_in)
  
  return(thetaS=theta_star)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Metropolis-Hastings algorithm

run_mcmc<-function(
  test.yr,
  test_years,
  inf_years,
  strain_years,
  n_part,
  test.list,
  theta,
  runs,
  varpart_prob,
  hist.true=NULL,
  switch1=2,
  seedi=1
  ){
  
  # DEBUG set params <<<
  # hist.true=NULL; test.yr=c(2008); runs=1; switch1=2; varpart_prob=0.05
  
  test.n=length(test_years); inf.n=length(inf_years); nstrains=length(strain_years)
  sample.index=strain_years-min(strain_years)+1
  historyii=rbinom(inf.n, 1, 0.1) # dummy infection history
  
  # Index variables
  jj_year=match(test.yr,test_years); testyear_index=match(test.yr,inf_years)
  sample.n=length(jj_year)
  
  # Specific MCMC parameters
  # Up to here
  #browser()
  pmask=rep(TRUE,length(theta))
  theta_fitted <- theta[pmask]
  
  nparam=length(theta); npcov=rep(1,nparam)
  cov_matrix_theta0 = diag(npcov)
  
  thetatab=matrix(NA,nrow=(runs+1),ncol=length(theta)); colnames(thetatab)=names(theta)
  thetatab[1,]=theta
  
  historytab=matrix(NA,nrow=n_part,ncol=inf.n)
  historytabCollect=historytab
  age.tab=matrix(NA,nrow=n_part,ncol=1)
  
  # Pick plausible initial conditions
  if(is.null(hist.true)){
    for(ii in 1:n_part){
      historytab[ii,]=setuphistIC(ii,jj_year[1],inf.n,test.list,testyear_index,inf_years) # Pick first test year
    }
  } else { historytab=hist.true }
  
  colnames(historytab)=as.character(inf_years)
  
  # Plausible intial ages - based on earliest strain in history
  #age.tab=sapply(
  #  apply(historytab,1,function(x){min(c(inf.n:1)[x==1])}),
  #  function(y){ sample(y:80, 1, replace=T) })
  
  # Preallocate matrices
  likelihoodtab=matrix(-Inf,nrow=(runs+1),ncol=n_part)
  accepttabT=NULL
  accepttabH=NULL
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Run MCMC
  
  for (m in 1:runs){
    
    # Adaptive covariance matrix
    if(m==1){
      epsilon0=0.01
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      varpart_prob0=varpart_prob
    }else{
      epsilon0=max(0.00001,min(1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m)))
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      varpart_prob0=max(0.01,min(1,exp(log(varpart_prob0)+(accept_rateH-0.234)*0.999^m)))
    }
    
    # - - - - - - - - - - - - - - - -
    # Resample parameters
    
    #aTime=Sys.time() #TIMER 1
    
    if(m %% switch1==0 | m==1){ # m==1 condition as have to calculate all liks on first step
      theta_star = SampleTheta(thetatab[m,], m,cov_matrix_theta) #resample theta
      #age_star = age.tab
      history_star = historytab
      pickA=c(1:n_part)
      
    }else{
      pickA=NULL
      pickA=sample(n_part, ceiling(varpart_prob0*n_part)) # check that not length zero
      #age_star = age.tab #SampleAge(pickA,age.tab) #resample age (not for now)
      history_star = SampleHistory(historytab,pickA,inf.n,age_star,inf_years) #resample history
      theta_star =thetatab[m,]
    }
    
    dmatrix=outputdmatrix(theta_star,inf_years) # Arrange parameters
    
    # - - - - - - - - - - - - - - - -
    # LIKELIHOOD function - Only calculate for updated history
    
    lik_val=likelihoodtab[m,]
    for(ii in pickA){
      # Set history to zero after test date
      lik.ii=rep(NA,sample.n)
      for(kk in 1:sample.n){
        #DEBUG DEBUG set params <<<  ii=1;kk=2;historyii=as.numeric(history_star[ii,])
        lik.ii[kk]=estimatelik(ii,jj_year[kk],as.numeric(history_star[ii,]),dmatrix,theta_star,test.list,testyear_index[kk])
      }
      lik_val[ii]=sum(lik.ii)
      #if(is.na(lik_val[ii])){lik_val[ii]=-Inf}
    }
    
    # - - - - - - - - - - - - - - - -
    # Metropolis Hastings step
    
    #print(c(m,sum(likelihoodtab[m,]),sum(lik_val))) # PRINT LIKELIHOOD
    
    output_prob = ComputeProbability(sum(likelihoodtab[m,]),sum(lik_val)) 
    
    if(is.na(output_prob) & m==1){stop('check initial parameter values')}
    
    if(runif(1) < output_prob){
      thetatab[m+1,] = theta_star
      if(m %% switch1!=0){historytab = history_star} # Only change if resampled
      #if(m %% switch1==0){age.tab = age_star} # Only change if resampled
      likelihoodtab[m+1,] = lik_val
      if(m %% switch1==0){accepttabT=c(accepttabT,1)}
      if(m %% switch1!=0){accepttabH=c(accepttabH,1)}
      
    }else{
      thetatab[m+1,] = thetatab[m,]
      likelihoodtab[m+1,] = likelihoodtab[m,]
      if(m %% switch1==0){accepttabT=c(accepttabT,0)}
      if(m %% switch1!=0){accepttabH=c(accepttabH,0)}
    }
    
    
    if(m<max(100)){
      accept_rateT=0.234
      accept_rateH=0.234
    }else{
      accept_rateT=sum(accepttabT)/length(accepttabT)
      accept_rateH=sum(accepttabH)/length(accepttabH)
    }
    
    #Sys.time()-aTime  #TIMER 2
    if(m %% min(runs,20) ==0){
      historytabCollect=rbind(historytabCollect,historytab)
    }
    
    if(m %% min(runs,100) ==0){
      print(c(m,accept_rateH,varpart_prob0,round(sum(likelihoodtab[m,]))))
      save(likelihoodtab,thetatab,n_part,test.list,historytab,historytabCollect,age.tab,file=paste("posterior_sero_runs/outputR",test.yr[1],"_",seedi,".RData",sep=""))
    }
    
  } #End runs loop
  
}
