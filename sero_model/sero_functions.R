#Functions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Functions to set up parameters for model

outputdmatrix<-function(theta){
  (dmatrix=sapply(inf_years,function(x){exp(-theta[["sigma"]]*abs(inf_years-x))})) # note that second entry is actually sample year
}

# - - - - - - - - - - - - - - - -
# Set initial condition as infection if titre >=4

setuphistIC<-function(ii,jj){ # ii=participant | jj=test year
  
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  spyear=unique(as.numeric(test.jj[3,]))
  
  hist0=rep(0,inf.n)   
  hist0[sample(c(1:inf.n),round(0.1*inf.n))]=1
  
  # Check test data available
  if(length(test.jj[,1])>1){
    
      # Set up test strains
      titredat=as.numeric(test.jj[2,]) # Define titre data
      
      for(i in 1:length(spyear)){
        if(max(titredat[(as.numeric(test.jj[3,])==spyear[i])])>=4){
          hist0[(inf_years==spyear[i])]=1
        }
      }
      
  }

  if(sum(hist0)==0){hist0[sample(1:inf.n,1)]=1}
  hist0
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Likelihood given infection history and parameters (done in C)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Compile c code
compile.c<-function(){
  setwd("~/Dropbox/Imperial/Fluscape_2/sero_model/c_code")
  system("R CMD SHLIB c_model2.c")
  dyn.load("./c_model2.so")
}

# - - - - - - - - - - - - - - - -
# Define likelihood function

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
  return(out$titrepred)
}


# - - - - - - - - - - - - - - - -
# Calculate likelihood for given participant and test year

estimatelik<-function(ii,jj,historyii,dmatrix,thetastar){ # ii=participant | jj=test year

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
  sum(dpois(as.numeric(titredat), expect, log = TRUE))

}

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Resample infection history

infvector=c(1:inf.n)

SampleHistory<-function(historyA,pick){

  for(ii in pick){
  #ls_pick=foreach(ii=(1:length(pick))) %dopar% {  # Parallel loop - slower to farm out
    
    x=historyA[ii,]
    
    rand1=runif(1)

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
    
    historyA[ii,]=x
    
  }
  
  historyA
  
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Convert infection history to binary

convert_binary <- function(x){sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

ComputeProbability<-function(pf_marg_likelihood,pf_marg_likelihood_star){

  # uniform priors
  p_theta_star = 1
  p_theta = 1
  
  # probability symmetic
  q_theta_given_theta_star = 1
  q_theta_star_given_theta = 1
  
  val = exp((pf_marg_likelihood_star-pf_marg_likelihood))*(p_theta_star/p_theta)*(q_theta_given_theta_star/q_theta_star_given_theta) 
  min(val, 1)
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

SampleTheta<-function(theta_in,m,covartheta){
  
  # sample new parameters from nearby: 
  mean_vector_theta = theta_in
  mean_vector_theta0=mean_vector_theta

  theta_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta0), covartheta)))
  names(theta_star)=names(theta_in)
  
  return(thetaS=theta_star)
  
}

