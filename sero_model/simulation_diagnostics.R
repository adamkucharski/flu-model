# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation diagnostics
# Compare MCMC output to simulation data

# Convert vector to median and 95% CrI
c.text<-function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

plot.posteriors<-function(simDat=F,loadseed=1,define.year){
  
  loadseed="1" #"1_w12"
  
  load(paste("posterior_sero_runs/outputR_f",paste(define.year,"_",collapse="",sep=""),"s",loadseed,".RData",sep=""))
  par(mfrow=c(3,3))
  par(mar = c(5,5,1,1))
  colA=rgb(0.8,0.8,0.8)
  
  # Plot profile likelihood
  lik.tot=rowSums(likelihoodtab)
  maxlik=max(lik.tot)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  runs1=ceiling(0.25*runsPOST)
  plot(rowSums(likelihoodtab)[runs1:runsPOST],type="l",ylab="likelihood",ylim=c(maxlik-500,maxlik))
  #plot(as.data.frame(thetatab)$mu[runs1:runsPOST],type="l",ylab="mu")
  #plot(as.data.frame(thetatab)$sigma[runs1:runsPOST],type="l",ylab="sigma")
  
  hist(as.data.frame(thetatab)$error[runs1:runsPOST],main=NULL,col=colA,xlab="error",prob=TRUE,xlim=c(0,0.1))
  
  # Plot histogram of boosting
  hist(as.data.frame(thetatab)$mu[runs1:runsPOST],main=NULL,col=colA,xlab="mu",prob=TRUE,xlim=c(0,5))
  if(simDat==T){abline(v=thetaSim[["mu"]],col="red")}
  
  hist(as.data.frame(thetatab)$sigma[runs1:runsPOST],main=NULL,col=colA,xlab="sigma",xlim=c(0,0.5))
  if(simDat==T){abline(v=thetaSim[["sigma"]],col="red")}
  
  hist(as.data.frame(thetatab)$tau1[runs1:runsPOST],main=NULL,col=colA,xlab="tau1",prob=TRUE,xlim=c(0,0.7))
  if(simDat==T){abline(v=thetaSim[["tau1"]],col="red")}
  
  hist(as.data.frame(thetatab)$tau2[runs1:runsPOST],main=NULL,col=colA,xlab="tau2",prob=TRUE,xlim=c(0,0.7))
  if(simDat==T){abline(v=thetaSim[["tau2"]],col="red")}
  
  hist(as.data.frame(thetatab)$wane[runs1:runsPOST],main=NULL,col=colA,xlab="wane",prob=TRUE,xlim=c(0,3))
  if(simDat==T){abline(v=thetaSim[["wane"]],col="red")}
  
  hist(as.data.frame(thetatab)$muShort[runs1:runsPOST],main=NULL,col=colA,xlab="mu_Short",prob=TRUE,xlim=c(0,20))
  if(simDat==T){abline(v=thetaSim[["muShort"]],col="red")}
  
  hist.sample=length(historytabCollect[,1])/n_part
  ind.infN=rowSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])
  hist(ind.infN,breaks=seq(-0.5,max(ind.infN)+0.5,1),col=colA,xlab="infections",prob=TRUE,main=paste("mean/med=",signif(mean(ind.infN),2),"/",median(ind.infN),sep=""),xlim=c(0,40))
  
  dev.copy(pdf,paste("plot_simulations/posterior",ifelse(simDat==T,paste("mu",thetaSim[["mu"]],"_sigma",thetaSim[["sigma"]],sep=""),""),"_np",n_part,"_yr",paste(define.year,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=12,height=8)
  dev.off()

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compare multiple MCMC outputs for vector of years

plot.compare<-function(simDat=F,loadseed=1,define.year.vec){
  
  loadseed="1" #"1_w12"
  n.test=length(define.year.vec)
  names1=c("test","mu","tau1","tau2","sigma","muShort","error","infections")
  store.val=array(NA,dim=c(3,length(names1),n.test),dimnames=list(NULL,names1,NULL))
  range.p=rbind(c(0,0),c(0,5),c(0,0.5),c(0,0.5),c(0,0.5),c(0,10),c(0,0.1),c(0,30))# define parameter ranges for plots

  for(kk in 1:n.test){
  
    load(paste("posterior_sero_runs/outputR_f",define.year.vec[kk],"_s",loadseed,".RData",sep=""))

    # Store median and 95% CrI
    lik.tot=rowSums(likelihoodtab); maxlik=max(lik.tot); runsPOST=length(lik.tot[lik.tot!=-Inf]); runs1=ceiling(0.25*runsPOST)
    hist.sample=length(historytabCollect[,1])/n_part; ind.infN=rowSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])
    store.val[,,kk]=cbind(rep(kk,3),
                     c.nume(as.data.frame(thetatab)$mu[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$tau1[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$tau2[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$sigma[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$muShort[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$error[runs1:runsPOST]),
                     c.nume(ind.infN)
                     )
  }
  
  # - - - - - - - - - - - - - - - 
  # Plot comparison of parameters
  
  par(mfrow=c(2,4))
  
  for(jj in 2:length(names1)){   # Iterate across parameters
    colA=rgb(0,0,0.8)
    plot(c(1:n.test),c(1:n.test),pch=19,col=rgb(1,1,1),ylim=range.p[jj,],xaxt="n",xlab="test year",ylab="estimate",main=names1[jj])
    axis(1,at=c(1:n.test),labels=define.year.vec)
    grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  
      for(kk in 1:n.test){ # Iterate across test years
        points(kk,store.val[1,jj,kk],pch=19,col=colA)
        lines(c(kk,kk),c(store.val[2,jj,kk],store.val[3,jj,kk]),col=colA)
      }
  
  }
  

  dev.copy(pdf,paste("plot_simulations/posterior_compare",paste(define.year.vec,"_",collapse="",sep=""),".pdf",sep=""),width=12,height=8)
  dev.off()
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot expected titres using sampled posterior estimates against true titres 

plot.posterior.titres<-function(loadseed=1,define.year){
  
  load(paste("posterior_sero_runs/outputR_f",paste(define.year,"_",collapse="",sep=""),"s",loadseed,".RData",sep=""))
  
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  runs1=ceiling(0.2*runsPOST)

  # Set up matrices to store
  btstrap=100
  n.strains=length(strain_years) # this loads from main_model.R
  n.test=length(test.yr)
  n.inf=length(inf_years)
  
  store.mcmc.test.data=array(NA, dim=c(btstrap,n_part,n.strains,n.test,2)) # Store expected titres for each test year
  store.mcmc.hist.data=array(NA, dim=c(btstrap,n_part,n.inf,n.test)) # Store history for each test year
  
  # - - - - - - - - - - - - 
  # Sample from MCMC runs to get stored matrices

  for(sampk in 1:btstrap){
    pickA=sample(c(runs1:runsPOST),1)
    pickAhist=ceiling(pickA/20)+1 # check which history this specifies
    hist.sample=historytabCollect[((pickAhist-1)*n_part+1):(pickAhist*n_part),1:n.inf]
    theta.max=as.data.frame(thetatab)[pickA,]
  
  for(pickyr in 1:n.test){
    
    # Note here that inf_years and strain_years are loads from main_model.R
    # Output expected titre - could include Poisson measurement uncertainty?
    simulate_data(test.yr[pickyr],historytabPost=hist.sample,
                  inf_years,
                  strain_years,
                  n_part,thetastar=theta.max,p.inf=0.1,linD=F)
    
    load("R_datasets/Simulated_dataPost_1.RData")
    
    # Mask infections after test year
    for(ii0 in 1:n_part){
      
      sim.titre=test.listSim[[ii0]][[1]] # sort sample years
      hist.sampleB=hist.sample;  hist.sampleB[,as.numeric(colnames(hist.sample))>test.yr[pickyr]]=0 # don't show infections after test year
      store.mcmc.test.data[sampk,ii0,,pickyr,1]= min(inf_years)-1+sort(sim.titre["sample.index",])
      store.mcmc.test.data[sampk,ii0,,pickyr,2]=sim.titre["titredat",order(sim.titre["sample.index",])]
      store.mcmc.hist.data[sampk,ii0,,pickyr]=hist.sampleB[ii0,]
      
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
      hist.sample=store.mcmc.hist.data[,ii0,,pickyr]
  
      plot(inf_years,8*hist.sample[1,],type="l",ylim=c(0,9),col='white',xlab="year",ylab="titre")
      
      # Sample from infection history
      for(ksamp in 1:btstrap){
        for(jj in 1:n.inf){
          lines(min(inf_years)-1+c(jj,jj),c(-1,10*hist.sample[ksamp,jj]-1),col=rgb(0.8,0.8,0.8,0.01),lwd=2) # Plot estimated infections
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
      
      
      if(ii0 %% 10==0){
        dev.copy(pdf,paste("plot_simulations/sim",ii0,"P_",pickyr,".pdf",sep=""),width=12,height=6)
        dev.off()
      }
    
    }  # end loop over participants
      
  } # end loop over test years

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot simulated data against runs

plot.sim.data<-function(){


  # UPDATE THIS BIT
  #Compare model fits using posterior infection history (historytabPost) and parameters
  simulate_data(define.year,historytabPost=historytab,
                inf_years,
                strain_years,
                n_part=npartM,thetastar=as.data.frame(thetatab)[runsPOST,],p.inf=0.1)
  
  load("R_datasets/Simulated_dataPost.RData")
  par(mfrow=c(2,5))
  par(mar = c(5,5,1,1))
  for(ii0 in 1:n_part){
    lenhis=rep(0,length(historytabSim[ii0,]))
    plot(8*historytabSim[ii0,],type="l",ylim=c(0,9),col='white')
    #for(jj in 1:length(lenhis)){
    #  lines(c(jj,jj),c(0,9*historytabSim[ii0,jj]),col='red')
    #}
    for(jj in 1:length(lenhis)){
      lines(c(jj,jj),c(0,9*historytab[ii0,jj]),col='blue')
    }
    lines(test.list[[ii0]][[1]][4,],test.list[[ii0]][[1]][2,],type="l")
    points(test.list[[ii0]][[1]][4,],test.list[[ii0]][[1]][2,],pch=19)
    if(ii0 %% 10==0){
      dev.copy(pdf,paste("plot_simulations/sim",ii0,"P.pdf",sep=""),width=12,height=6)
      dev.off()
    }
  }

}

