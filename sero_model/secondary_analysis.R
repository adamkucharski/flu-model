# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# More elaroate plotting using posterior data
# Plot posteriors and compare MCMC output to titre data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate epidemics from posterior. Calculate antigenic landscape and selection pressure.

simulate.pressure<-function(loadseed=1,year_test=c(2007:2012),flu.type,simDat=F,btstrap=5){
  
  #simDat=F;loadseed=1;year_test=c(2007:2012);plotmap=F;f.lim=T;flu.type="H3"
  
  # SIMULATION CODE - ADD NEW VECTOR OF TEST STRAINS AFTER INFECTION YEAR?
  # CAN IT HANDLE THIS? YES - JUST MAKE INFECTION_YEARS vector REALLY LONG
  
  if(simDat==F){
    if(flu.type=="H3"){
      load("R_datasets/HaNam_data.RData")
    }else{
      load("R_datasets/Fluscape_data_List.RData")
      year_test=c(2011,2012)
    }
    loadseed=paste(loadseed,"_",flutype,sep="")
  }else{
    load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))
    test.list=test.listSim
    hist.true=historytabSim
  }
  
  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,".RData",sep="")) # Note that this includes test.listPost
  
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  runs1=ceiling(0.2*runsPOST)
  
  # Set up matrices to store -- need btstrap >1
  n.strains=length(strain_years) # this loads from main_model.R
  n.test=length(test.yr)
  n.inf=length(inf_years)
  
  store.mcmc.test.data=array(NA, dim=c(btstrap,n_part,n.strains,n.test,2)) # Store expected titres for each test year
  store.mcmc.hist.data=array(NA, dim=c(btstrap,n_part,n.inf,n.test)) # Store history for each test year
  
  # - - - - - - - - - - - -
  # Sample from MCMC runs to get stored matrices of expected titre and estimated infection years
  
  for(sampk in 1:btstrap){
    pickA=sample(c(runs1:runsPOST),1)
    pickAhist=ceiling(pickA/20)+1 # check which history this specifies (i.e. store every 20)
    hist.sample=historytabCollect[((pickAhist-1)*n_part+1):(pickAhist*n_part),1:n.inf]
    theta.max=as.data.frame(thetatab)[pickA,]
    
    for(pickyr in 1:n.test){
      
      # Note here that inf_years and strain_years are loads from main_model.R
      # Output expected titre - could include Poisson measurement uncertainty here?
      
      # MAKE EDITS HERE
      # DEFINE NEW STRAINS
      
      
      simulate_data(test.yr[pickyr],historytabPost=hist.sample,
                    inf_years,
                    strain_years,
                    n_part,thetastar=theta.max,p.inf=0.1,
                    #pmask=c("sigma2"), # For old fitted data, need to specify that sigma2 wasn't fitted
                    linD=F)
      
      load("R_datasets/Simulated_dataPost_1.RData")
      
      # Mask infections after test year
      for(ii0 in 1:n_part){
        
        sim.titre=test.listSim[[ii0]][[1]] # sort sample years - drawn from simulated data above
        hist.sampleB=hist.sample;  hist.sampleB[,as.numeric(colnames(hist.sample))>test.yr[pickyr]]=0 # don't show infections after test year
        store.mcmc.test.data[sampk,ii0,,pickyr,1]= min(inf_years)-1+sort(sim.titre["sample.index",]) # Sampled strain years
        s.titre=sim.titre["titredat",order(sim.titre["sample.index",])]
        store.mcmc.test.data[sampk,ii0,,pickyr,2]=s.titre # Sampled expected titre
        #store.mcmc.test.data[sampk,ii0,,pickyr,2]=rpois(length(s.titre),lambda=s.titre) # Sampled expected titre - include Poisson noise
        store.mcmc.hist.data[sampk,ii0,,pickyr]=hist.sampleB[ii0,] # Sampled history
        
      }  # end loop over participants
      
    } # end loop over test years
  } # end bootstrap loop
  
  # - - - - - - - - - - - - 
  # Plot figures from MCMC posteriors
  
  for(pickyr in 1:n.test){
    
    par(mfrow=c(2,5)); par(mar = c(5,5,1,1))
    # Mask infections after test year
    
    for(ii0 in 1:n_part){
      simtitreX=store.mcmc.test.data[,ii0,,pickyr,1]
      simtitreY=store.mcmc.test.data[,ii0,,pickyr,2]
      hist.sample=store.mcmc.hist.data[,ii0,,pickyr] # for participant ii0 in year pickyr
      
      plot(inf_years,8*hist.sample[1,],type="l",ylim=c(0,9),col='white',xlab="year",ylab="titre")
      
      # Sample from infection history
      for(ksamp in 1:btstrap){
        for(jj in 1:n.inf){
          lines(min(inf_years)-1+c(jj,jj),c(-1,12*hist.sample[ksamp,jj]-1),col=rgb(0.8,0.8,0.8,0.05),lwd=2) # Plot estimated infections
        }
      }
      
      # Calculate credible interval for expected titres
      medP=apply(simtitreY,2,function(x){median(x)})
      ciP1=apply(simtitreY,2,function(x){quantile(x,0.025)})
      ciP2=apply(simtitreY,2,function(x){quantile(x,0.975)})
      polygon(c(simtitreX[1,],rev(simtitreX[1,])),c(ciP1,rev(ciP2)),lty=0,col=rgb(0,0.3,1,0.2))
      lines(simtitreX[1,],medP,pch=1,col='blue')
      points(simtitreX[1,],medP,pch=19,cex=0.5,col='blue')
      
      # Plot true titres
      points(min(inf_years)-1+test.list[[ii0]][[pickyr]][4,],test.list[[ii0]][[pickyr]][2,],pch=1,col='red')
      
      # Plot true infections if simulation
      if(simDat==T){
        histSim1=hist.true[ii0,]
        histSim1[inf_years>year_test[pickyr]]=0 # don't show infections after test year
        lenhis=rep(0,length(histSim1))
        for(jj in 1:length(lenhis)){
          lines(min(inf_years)-1+c(jj,jj),c(10*histSim1[jj]-2,12*histSim1[jj]-2),lwd=2,col=rgb(0,0.5,0))
        }
      }
      
      if(ii0 %% 10==0){
        dev.copy(pdf,paste("plot_simulations/sim",ii0,"P_",pickyr,".pdf",sep=""),width=12,height=6)
        dev.off()
      }
      
    }  # end loop over participants
    
  } # end loop over test years
  
}