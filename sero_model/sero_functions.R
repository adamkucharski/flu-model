#Functions

# - - - - - - - - - - - - - - - -
# Set initial condition (for infection history) as infection if titre >=X

setuphistIC<-function(ii,jj,inf.n,test.list){ # ii=participant | jj=test year
  
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
      if(sum(maxt)>1){
        hist0[inf_years==spyear[sample(c(1:length(maxt))[maxt],1)]]=1
      }else{
        
        hist0[(inf_years==spyear[titredat==max(titredat)])]=1
      }

      # Use simple cutoff for titres
      #for(i in 1:length(spyear)){
      #  if(max(titredat[(as.numeric(test.jj[3,])==spyear[i])])>=4){
      #    hist0[(inf_years==spyear[i])]=1
      #  }
      #}
      
  }
  pos.hist=(hist0>0)
  if(sum(hist0)==0){hist0[sample(1:inf.n,1)]=1} # Make sure at least one infection
  hist0
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Likelihood given infection history and parameters (done in C)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Functions to set up parameters for model

outputdmatrix<-function(theta,inf_years){
  (dmatrix=sapply(inf_years,function(x){exp(-theta[["sigma"]]*abs(inf_years-x))})) # note that second entry is actually sample year
}

# Compile c code
compile.c<-function(){
  setwd("c_code")
  system("R CMD SHLIB c_model2.c")
  system("R CMD SHLIB c_model2_sr.c")
  dyn.load("./c_model2.so")
  dyn.load("./c_model2_sr.so")
  setwd("..")
}

# - - - - - - - - - - - - - - - -
# Define expected titre function

func1 <- function(x,titredat,dd,mu) {
  if (!is.numeric(x)){stop("argument x must be numeric")}
  out <- .C("c_model2",
            n=as.integer(length(x)),
            nsample=as.integer(length(titredat)),
            x=as.double(x),
            x1=as.double(rep(0,length(x))),
            titre=as.double(titredat),
            titrepred=as.double(rep(0,length(titredat))),
            dd=as.double(dd),
            mu=as.double(mu)
  )
  out2 <- .C("c_model2_sr",
            n=as.integer(length(x)),
            nsample=as.integer(length(titredat)),
            x=as.double(x),
            x1=as.double(rep(0,length(x))),
            titre=as.double(titredat),
            titrepred=as.double(rep(0,length(titredat))),
            dd=as.double(dd),
            mu=as.double(mu)
  )
  # browser()
  # out$titrepred - out2$titrepred
  return(out$titrepred)
}


# - - - - - - - - - - - - - - - -
# Calculate likelihood for given participant and test year

estimatelik<-function(ii,jj,historyii,dmatrix,thetastar,test.list){ # ii=participant | jj=test year

  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  # Check test data available
  if(length(test.jj[,1])==1){0}else{
  
  # Set up test strains
  test.part=as.numeric(test.jj[4,]) # index of sample strains data available for
  titredat=test.jj[2,] # Define titre data
  
  d.ij=dmatrix[test.part,] # Define cross-immunity matrix for sample strain
  d_vector=melt(t(d.ij))$value
  
  expect=func1(historyii,titredat,d_vector,thetastar[["mu"]]) # Output expectation
  #plot(as.numeric(test.jj[3,]),expect,ylim=c(0,100))
  
  # Calculate likelihood - ** need to add summation if k>8 **
  largett=(titredat>=8)
  
  sum(dpois(as.numeric(titredat[!largett]), expect[!largett], log = TRUE))+
    sum(ppois(8, lambda=expect[largett], lower=FALSE,log=TRUE))
}

}


# - - - - - - - - - - - - - - - -
# Simulation infection history data

simulate_data<-function(test_years,historytabPost=NULL, inf_years,strain_years,n_part=20,thetastar=theta0,p.inf=0.2){ # ii=participant | jj=test year
  
  # Variables needed: test_years,inf_years,strain_years,n_part
  #strain_years=seq(1968,2010,4)
  
  # Set year of birth
  age.yr=sample(1:80,n_part,replace = TRUE)
  
  test.n=length(test_years)
  inf.n=length(inf_years)
  nstrains=length(strain_years)
  sample.index=strain_years-min(strain_years)+1
  
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
      
      expect=func1(historyii,sample.index,d_vector,thetastar[["mu"]]) # Output expectation
      #titredat=sapply(expect,function(x){rpois(1,x)}) # Generate titre
      titredat=round(expect)
      titredat=sapply(titredat,function(x){min(x,8)})
      
      testyr=test_years[jj]
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
  if(is.null(historytabPost)){
    save(test_years,inf_years,strain_years,n_part,test.list,age.yr,historytabSim,file=paste("R_datasets/Simulated_data.RData",sep=""))
  }else{
    save(test_years,inf_years,strain_years,n_part,test.list,age.yr,file=paste("R_datasets/Simulated_dataPost.RData",sep=""))
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
      if(length(infectID)>1){
        x[sample(infectID,1)]=0
      }
    }
    
    # Add new infection
    if(rand1>1/3 & rand1<2/3){
      ninfecID=infvector[(as.numeric(x)==0)]
      if(length(ninfecID)>0){
        x[sample(ninfecID,1)]=1
      }
    }
    
    # Move infection
    if(rand1>2/3){
      infectID=infvector[(as.numeric(x)>0)]
      ninfecID=infvector[(as.numeric(x)==0)]
      if(length(infectID)>0 & length(ninfecID)>0){
        x[sample(infectID,1)]=0
        x[sample(ninfecID,1)]=1
      }
    }
    
    # Add prior on birth year - exponentially less likely to update if infections outside
    if(inf.n>ageA[ii]){
      a1=exp(-0.1*sum((infvector2*x)[1:(inf.n-ageA[ii])])) # tweak this parameter to penalise more/less
      if( a1 > runif(1) ){
        historyA[ii,]=x
      }
    }
  
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
  p_theta_star = 1
  p_theta = 1
  
  # probability symmetic
  q_theta_given_theta_star = 1
  q_theta_star_given_theta = 1
  
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

run_mcmc<-function(test.yr,test_years,inf_years,strain_years,n_part,test.list,theta0,runs,varpart_prob,hist.true=NULL){
  
  test.n=length(test_years)
  inf.n=length(inf_years)
  nstrains=length(strain_years)
  sample.index=strain_years-min(strain_years)+1
  historyii=rbinom(inf.n, 1, 0.1) # dummy infection history
  
  # Index variables
  jj_year=c(1:test.n)[test_years==test.yr]
  
  # Specific MCMC parameters
  nparam=length(theta)
  npcov=rep(1,nparam)
  cov_matrix_theta0 = diag(npcov)
  
  thetatab=matrix(NA,nrow=(runs+1),ncol=length(theta))
  colnames(thetatab)=names(theta)
  thetatab[1,]=theta
  
  historytab=matrix(NA,nrow=n_part,ncol=inf.n)
  age.tab=matrix(NA,nrow=n_part,ncol=1)
  
  # Pick plausible initial conditions
  if(is.null(hist.true)){
    for(ii in 1:n_part){
      historytab[ii,]=setuphistIC(ii,jj_year,inf.n,test.list)
    }
  }else{
    historytab=hist.true
  }
  
  colnames(historytab)=as.character(inf_years)

  # Plausible intial ages - based on earliest strain in history
  age.tab=sapply(
    apply(historytab,1,function(x){min(c(inf.n:1)[x==1])}),
    function(y){ sample(y:80, 1, replace=T) })
  
  # Preallocate matrices
  likelihoodtab=matrix(-Inf,nrow=(runs+1),ncol=n_part)
  accepttabT=rep(NA,(runs/2))
  accepttabH=rep(NA,(runs))
  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Run MCMC
  
  for (m in 1:runs){
    
    # Adaptive covariance matrix
    if(m==1){
      epsilon0=0.001
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      #varpart_prob0=varpart_prob
    }else{
      epsilon0=min(1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m))
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      
      #varpart_prob0=min(1,exp(log(varpart_prob0)+(accept_rate-0.234)*0.999^m))
      #varpart_prob=varpart_prob0
    }
    
    # - - - - - - - - - - - - - - - -
    # Resample parameters
    
    #aTime=Sys.time() #TIMER 1
    
    if(m %% 2==1 | varpart_prob==0){
      theta_star = SampleTheta(thetatab[m,], m,cov_matrix_theta) #resample theta
      age_star = age.tab
      history_star = historytab
      pickA=c(1:n_part)
      
    }else{
      pickA=NULL
      pickA=sample(n_part, ceiling(varpart_prob*n_part)) # check that not length zero
      age_star = SampleAge(pickA,age.tab) #resample history
      history_star = SampleHistory(historytab,pickA,inf.n,age_star,inf_years) #resample history
      theta_star =thetatab[m,]
    }
    
    dmatrix=outputdmatrix(theta_star,inf_years) # Arrange parameters
    
    # - - - - - - - - - - - - - - - -
    # LIKELIHOOD function - Only calculate for updated history
    
    lik_val=likelihoodtab[m,]
    for(ii in pickA){
      # Set history to zero after test date
      
      lik_val[ii]=estimatelik(ii,jj_year,as.numeric(history_star[ii,]),dmatrix,theta_star,test.list)
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
    
    if(m %% min(runs,200) ==0){
      print(c(m,accept_rateT,epsilon0,round(sum(likelihoodtab[m,]))))
      save(likelihoodtab,thetatab,historytab,file=paste("posterior_sero_runs/outputR.RData",sep=""))
    }
    
  }
  
}
