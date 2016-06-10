simulate_sera_data<-function(test_years,historytabPost=NULL, inf_years,strain_years,n_part=20,thetastar=theta0,p.inf=0.2,seedi=1,roundv=F,linD=F,dmatrix.in=NULL){ # ii=participant | jj=test year
  
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
  
  if(is.null(dmatrix.in)){
    dmatrix=outputdmatrix(thetastar,inf_years,linD)
  }else{
    dmatrix=dmatrix.in
  }
  
  #Set per year incidence, to create correlation between participant infection histories
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
      if(roundv==T){titredat=round(expect)}else{titredat=expect}
      titredat=sapply(titredat,function(x){min(x,8)})

      i.list[[jj]]=rbind(test.year=rep(testyr,nstrains),
                         titredat,
                         strain_years,
                         sample.index
      )
    }
    #i.list[[1]][2,]

    test.list[[ii]]=i.list
  }
  test.listSim=test.list
  
  # Export data
  #browser()
  if(is.null(historytabPost)){
    save(test_years,inf_years,strain_years,n_part,test.listSim,age.yr,historytabSim,file=paste("R_datasets/Simulated_data_",seedi,".RData",sep=""))
  }else{
    save(test_years,inf_years,strain_years,n_part,test.listSim,age.yr,file=paste("R_datasets/Simulated_dataPost_",seedi,".RData",sep=""))
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