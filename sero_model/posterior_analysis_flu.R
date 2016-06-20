# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation and inference diagnostics
# Plot posteriors and compare MCMC output to titre data

# Convert vector to median and 95% CrI
c.text<-function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

# Rotate antigenic maps to be uniform
scale.map<-function(map.pick){
  f.m=length(map.pick[,1])
  # translate map to finish at (0,0)
  map.pick[,1]=map.pick[,1]-map.pick[f.m,1]
  map.pick[,2]=map.pick[,2]-map.pick[f.m,2]
  
  # rotate map so final coordinate is (0,0) and penultimate is (0,1) 
  r.theta=atan(map.pick[(f.m-1),1]/map.pick[(f.m-1),2]) # angle of rotation
  
  map.pick[,1]=map.pick[,1]*cos(r.theta)-map.pick[,2]*sin(r.theta)
  map.pick[,2]=map.pick[,1]*sin(r.theta)+map.pick[,2]*cos(r.theta)
  
  # flip so 3rd from last is positive
  map.pick[,1]=sign(map.pick[(f.m-2),1])*map.pick[,1]
  
  map.pick
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot MCMC posterior distributions

plot.posteriors<-function(simDat=F,loadseed=1,define.year=c(2007:2012)){
  
  loadseed="SIM"
  load(paste("posterior_sero_runs/outputR_f",paste(define.year,"_",collapse="",sep=""),"s",loadseed,".RData",sep=""))
  par(mfrow=c(3,3))
  par(mar = c(5,5,1,1))
  colA=rgb(0.8,0.8,0.8)
  
  # Define lengths and sizes of inputs
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  maxlik=max(lik.tot[1:runsPOST])

  #plot(as.data.frame(thetatab)$mu[runs1:runsPOST],type="l",ylab="mu")
  
  # - - - - - - - 
  # Calculate ESS by burn-in
  runs1=ceiling(0.25*runsPOST)
  
  calculate.ESS<-function(runs1){
    thetaT=as.data.frame(thetatab)[runs1:runsPOST,]
    ltheta=length(thetaT[["mu"]])
    thin.theta=thetaT[seq(1,ltheta,switch1),]
    
    ESS.calc=effectiveSize(thin.theta)
    ESS.calc
  }

  ESS.calc=calculate.ESS(runs1)
  
  thetaT=as.data.frame(thetatab)[runs1:runsPOST,]
  ltheta=length(thetaT[["mu"]])
  thin.theta=thetaT[seq(1,ltheta,switch1),]
  
  ESS.label<-function(x){paste("ESS=",signif(as.numeric(x),3))}
  
  # - - - - - - - 
  # Plot results
  
  plot(rowSums(likelihoodtab)[runs1:runsPOST],type="l",ylab="likelihood",ylim=c(maxlik-500,maxlik))
  
  # Plot histograms of parameters
  hist(as.data.frame(thetatab)$error[runs1:runsPOST],main=ESS.label(ESS.calc[["error"]]),col=colA,xlab="error",prob=TRUE,xlim=c(0,0.1))
  if(simDat==T){abline(v=thetaSim[["error"]],col="red")}

  hist(thin.theta[["mu"]],main= ESS.label(ESS.calc[["mu"]]),col=colA,xlab="mu",prob=TRUE,xlim=c(0,5))
  if(simDat==T){abline(v=thetaSim[["mu"]],col="red")}
  
  hist(thin.theta[["sigma"]],main=ESS.label(ESS.calc[["sigma"]]),col=colA,xlab="sigma",xlim=c(0,0.5))
  if(simDat==T){abline(v=thetaSim[["sigma"]],col="red")}
  
  hist(thin.theta[["tau1"]],main=ESS.label(ESS.calc[["tau1"]]),col=colA,xlab="tau1",prob=TRUE,xlim=c(0,0.7))
  if(simDat==T){abline(v=thetaSim[["tau1"]],col="red")}
  
  hist(thin.theta[["tau2"]],main=ESS.label(ESS.calc[["tau2"]]),col=colA,xlab="tau2",prob=TRUE,xlim=c(0,0.7))
  if(simDat==T){abline(v=thetaSim[["tau2"]],col="red")}
  
  hist(thin.theta[["wane"]],main=ESS.label(ESS.calc[["wane"]]),col=colA,xlab="wane",prob=TRUE,xlim=c(0,10))
  if(simDat==T){abline(v=thetaSim[["wane"]],col="red")}
  
  hist(thin.theta[["muShort"]],main=ESS.label(ESS.calc[["muShort"]]),col=colA,xlab="mu_Short",prob=TRUE,xlim=c(0,40))
  if(simDat==T){abline(v=thetaSim[["muShort"]],col="red")}
  
  # Plot distribution of infections
  hist.sample=length(historytabCollect[,1])/n_part
  ind.infN=rowSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])
  hist(ind.infN,breaks=seq(0,max(ind.infN)+1,2),col=colA,xlab="infections",prob=TRUE,main=paste("mean/med=",signif(mean(ind.infN),2),"/",median(ind.infN),sep=""),xlim=c(0,40))
  
  dev.copy(pdf,paste("plot_simulations/posterior",ifelse(simDat==T,paste("mu",thetaSim[["mu"]],"_sigma",thetaSim[["sigma"]],sep=""),""),"_np",n_part,"_yr",paste(define.year,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=12,height=8)
  dev.off()
  
  # - - - -
  # Plot antigenic map
  par(mfrow=c(1,1))
  par(mar = c(5,5,1,1))
  map.sample=length(map.tabCollect)

  minMT=min(unlist(map.tabCollect)); maxMT=max(unlist(map.tabCollect))
  MTx=c(-5,15)
  MTy=c(-5,10)
  plot(map.tabCollect[[1]],xlim=MTx,ylim=MTy,type="l",col="white",lwd=2,xlab="strain dimension 1", xaxs="i", yaxs="i",ylab="strain dimension 2")
  #lines(inf_years,inf_years-min(inf_years),col="lightgray",lty=2)
  #plot(map.tabCollect[1,],0*map.tabCollect[1,],type="l",ylim=c(0,1),col="white",yaxt="n",ylab="",yaxs="i",xlab="strain position")
  for(ii in 1:200){
    pick=sample(c(round(0.15*map.sample):map.sample),1)
    
    map.pick = scale.map(map.tabCollect[[pick]])
    sigma.scale=as.numeric(as.data.frame(thetatab)[pick*20,]["sigma"]) # scale to one antigenic unit
    map.pick[,1]=map.pick[,1]/(exp(-sigma.scale))
    map.pick[,2]=map.pick[,2]/(map.pick[2,2]*exp(-sigma.scale))
    
    map.pick = map.tabCollect[[pick]]
    points(map.pick,col=rgb(0.5,0.5,0.9,0.1),pch=19,cex=1)
  }
  
  if(simDat==T){
    map.pick = scale.map(antigenic.map.in)
    map.pick[,1]=map.pick[,1]/(exp(-thetaSim[["sigma"]]))
    map.pick[,2]=map.pick[,2]/(map.pick[2,2]*exp(-thetaSim[["sigma"]]))
    lines(antigenic.map.in,col="black",lwd=2)
  }
  
  ## Alternative way of plotting linear antigenic space
  #plot(map.tabCollect[[1]],0*map.tabCollect[[1]],type="l",ylim=c(0,1),col="white",yaxt="n",ylab="",yaxs="i",xlab="strain position")
  #lines(inf_years,inf_years-min(inf_years),col="lightgray",lty=2)
  #for(ii in 1:500){ 
  #  for(jj in 1:length(map.tabCollect[1,])){
  #    pick=sample(hist.sample-1,1)
  #    lines(c(map.tabCollect[[pick]][jj],map.tabCollect[[pick]][jj]),c(0,1),type="l",col=rgb(0.7,0.7,0.8,0.05))
  #  }
  #}
  
  dev.copy(pdf,paste("plot_simulations/antigenic_map",ifelse(simDat==T,paste("mu",thetaSim[["mu"]],"_sigma",thetaSim[["sigma"]],sep=""),""),"_np",n_part,"_yr",paste(define.year,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=6,height=6)
  dev.off()

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compare multiple MCMC outputs for vector of different years fitted

plot.compare<-function(simDat=F,loadseed=1,define.year.vec=c(2007:2012)){
  
  loadseed="1" #"1_w12"
  n.test=length(define.year.vec)
  names1=c("test","mu","tau1","tau2","disp_k","sigma","muShort","error","infections")
  store.val=array(NA,dim=c(3,length(names1),n.test),dimnames=list(NULL,names1,NULL))
  
  labels.Y=if(length(define.year.vec)>6){c(define.year.vec[1:6],"All")}else{define.year.vec}
  
  for(kk in 1:n.test){
  
    load(paste("posterior_sero_runs/outputR_f",define.year.vec[kk],"_s",loadseed,".RData",sep=""))
    # Store median and 95% CrI
    lik.tot=rowSums(likelihoodtab); maxlik=max(lik.tot); runsPOST=length(lik.tot[lik.tot!=-Inf]); runs1=ceiling(0.25*runsPOST)
    hist.sample=length(historytabCollect[,1])/n_part; ind.infN=rowSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])
    store.val[,,kk]=cbind(rep(kk,3),
                     c.nume(as.data.frame(thetatab)$mu[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$tau1[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$tau2[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$disp_k[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$sigma[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$muShort[runs1:runsPOST]),
                     c.nume(as.data.frame(thetatab)$error[runs1:runsPOST]),
                     c.nume(ind.infN))
  }
  # - - - - - - - - - - - - - - - 
  # Plot comparison of parameters
  par(mfrow=c(2,4))
  
  range.p=rbind(c(0,0),c(0,5),c(0,0.5),c(0,0.5),c(0,5),c(0,0.5),c(0,10),c(0,0.1),c(0,30))# define parameter ranges for plots
  range.p[5,]=c(0,max(store.val[1,5,]))
  for(jj in 2:length(names1)){   # Iterate across parameters
    colA=rgb(0,0,0.8)
    plot(c(1:n.test),c(1:n.test),pch=19,col=rgb(1,1,1),ylim=range.p[jj,],xaxt="n",xlab="test year",ylab="estimate",main=names1[jj])
    axis(1,at=c(1:n.test),labels=labels.Y)
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

plot.posterior.titres<-function(loadseed=1,define.year=c(2007:2012)){
  
  load("R_datasets/HaNam_data.RData")
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
  # Sample from MCMC runs to get stored matrices of expected titre and estimated infection years

  for(sampk in 1:btstrap){
    pickA=sample(c(runs1:runsPOST),1)
    pickAhist=ceiling(pickA/20)+1 # check which history this specifies
    hist.sample=historytabCollect[((pickAhist-1)*n_part+1):(pickAhist*n_part),1:n.inf]
    theta.max=as.data.frame(thetatab)[pickA,]
  
  for(pickyr in 1:n.test){
    
    # Note here that inf_years and strain_years are loads from main_model.R
    # Output expected titre - could include Poisson measurement uncertainty here?
    simulate_data(test.yr[pickyr],historytabPost=hist.sample,
                  inf_years,
                  strain_years,
                  n_part,thetastar=theta.max,p.inf=0.1,linD=F)
    
    load("R_datasets/Simulated_dataPost_1.RData")
    
    # Mask infections after test year
    for(ii0 in 1:n_part){
      
      sim.titre=test.listSim[[ii0]][[1]] # sort sample years
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
# Plot simulated data

plot.sim.data<-function(){

  loadseed="SIM"
  load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))
  
  #Compare model fits using posterior infection history (historytabPost) and parameters
  simulate_data(define.year,historytabPost=historytabSim,
                inf_years,
                strain_years,
                n_part=npartM,thetastar=thetaSim,p.inf=0.1)
  
  load("R_datasets/Simulated_dataPost_1.RData")
  par(mfrow=c(2,5))
  par(mar = c(5,5,1,1))
  for(ii0 in 1:n_part){
    lenhis=rep(0,length(historytabSim[ii0,]))
    plot(8*historytabSim[ii0,],type="l",ylim=c(0,9),col='white')
    #for(jj in 1:length(lenhis)){
    #  lines(c(jj,jj),c(0,9*historytabSim[ii0,jj]),col='red')
    #}
    for(jj in 1:length(lenhis)){
      lines(c(jj,jj),c(0,9*historytabSim[ii0,jj]),col='blue')
    }
    lines(test.listSim[[ii0]][[1]][4,],test.listSim[[ii0]][[1]][2,],col=rgb(0.8,0.8,0.8),lwd=2) # Plot estimated infections
    points(test.listSim[[ii0]][[1]][4,],test.listSim[[ii0]][[1]][2,],pch=19)
    if(ii0 %% 10==0){
      dev.copy(pdf,paste("plot_simulations/sim",ii0,"P.pdf",sep=""),width=12,height=6)
      dev.off()
    }
  }

}

