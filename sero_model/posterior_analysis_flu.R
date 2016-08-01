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

# Transform antigenic maps to be uniform
scale.map<-function(map.pick){
  f.m=1 #length(map.pick[,1])
  # translate map to finish at (0,0)
  map.pick[,1]=map.pick[,1]-map.pick[f.m,1]
  map.pick[,2]=map.pick[,2]-map.pick[f.m,2]
  
  # rotate map so final coordinate is (0,0) and penultimate is (0,1) 
  #r.theta=atan(map.pick[(f.m-1),1]/map.pick[(f.m-1),2]) # angle of rotation
  #map.pick[,1]=map.pick[,1]*cos(r.theta)-map.pick[,2]*sin(r.theta)
  #map.pick[,2]=map.pick[,1]*sin(r.theta)+map.pick[,2]*cos(r.theta)
  
  # flip so 3rd from last is positive
  #map.pick[,1]=sign(map.pick[(f.m-2),1])*map.pick[,1]
  
  map.pick
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot MCMC posterior distributions


plot.posteriors<-function(simDat=F,loadseed=1,flutype="",year_test=c(2007:2012),plotmap=F,f.lim=F){

  #simDat=T;loadseed="SIM";year_test=c(2007:2012);plotmap=F;f.lim=T;flutype="H3"
  
  if(simDat==F){loadseed=paste(loadseed,"_",flutype,sep="")}
  if(simDat==T){load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))}

  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,".RData",sep=""))
  par(mfrow=c(5,2))
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
    thetaT=as.data.frame(thetatab)[runs1:runsPOST,]; ltheta=length(thetaT[["mu"]]); thin.theta=thetaT[seq(1,ltheta,switch1),]
    ESS.calc=effectiveSize(thin.theta); ESS.calc
  }

  ESS.calc=calculate.ESS(runs1)
  thetaT=as.data.frame(thetatab)[runs1:runsPOST,]
  ltheta=length(thetaT[["mu"]])
  thin.theta=thetaT[seq(1,ltheta,switch1),]
  ESS.label<-function(x){paste("ESS=",signif(as.numeric(x),3))}
  
  # - - - - - - - 
  # Plot results
  
  plot(lik.tot[runs1:runsPOST],type="l",ylab="likelihood",ylim=c(maxlik-500,maxlik))
  
  # Plot histograms of parameters
  hist(thin.theta[["error"]],main=ESS.label(ESS.calc[["error"]]),col=colA,xlab="error",prob=TRUE,xlim=c(0,ifelse(f.lim==F,15,1.1*max(thin.theta[["error"]])))); if(simDat==T){abline(v=theta.sim.out[["error"]],col="red")}
  hist(thin.theta[["mu"]],main= ESS.label(ESS.calc[["mu"]]),col=colA,xlab="mu",prob=TRUE,xlim=c(0,ifelse(f.lim==F,15,1.1*max(thin.theta[["mu"]])))); if(simDat==T){abline(v=theta.sim.out[["mu"]],col="red")}
  hist(thin.theta[["sigma"]],main=ESS.label(ESS.calc[["sigma"]]),col=colA,xlab="sigma",xlim=c(0,ifelse(f.lim==F,1,1.1*max(thin.theta[["sigma"]])))); if(simDat==T){abline(v=theta.sim.out[["sigma"]],col="red")}
  hist(thin.theta[["muShort"]],main=ESS.label(ESS.calc[["muShort"]]),col=colA,xlab="mu_Short",prob=TRUE,xlim=c(0,ifelse(f.lim==F,20,1.1*max(thin.theta[["muShort"]])))); if(simDat==T){abline(v=theta.sim.out[["muShort"]],col="red")}
  hist(thin.theta[["sigma2"]],main=ESS.label(ESS.calc[["sigma2"]]),col=colA,xlab="sigma2",xlim=c(0,ifelse(f.lim==F,1,1.1*max(thin.theta[["sigma2"]])))); if(simDat==T){abline(v=theta.sim.out[["sigma2"]],col="red")}
  hist(thin.theta[["tau1"]],main=ESS.label(ESS.calc[["tau1"]]),col=colA,xlab="tau1",prob=TRUE,xlim=c(0,ifelse(f.lim==F,0.2,1.1*max(thin.theta[["tau1"]])))); if(simDat==T){abline(v=theta.sim.out[["tau1"]],col="red")}
  hist(thin.theta[["tau2"]],main=ESS.label(ESS.calc[["tau2"]]),col=colA,xlab="tau2",prob=TRUE,xlim=c(0,ifelse(f.lim==F,0.2,1.1*max(thin.theta[["tau2"]])))); if(simDat==T){abline(v=theta.sim.out[["tau2"]],col="red")}
  hist(thin.theta[["wane"]],main=ESS.label(ESS.calc[["wane"]]),col=colA,xlab="wane",prob=TRUE,xlim=c(0,ifelse(f.lim==F,5,1.1*max(thin.theta[["wane"]])))); if(simDat==T){abline(v=theta.sim.out[["wane"]],col="red")}
  # Plot distribution of infections
  hist.sample=length(historytabCollect[,1])/n_part # need this sample value because table is stacked
  ind.infN=rowSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])
  hist(ind.infN,breaks=seq(0,max(ind.infN)+1,2),col=colA,xlab="infections",prob=TRUE,main=paste("mean/med=",signif(mean(ind.infN),2),"/",median(ind.infN),sep=""),xlim=c(0,40))

  dev.copy(pdf,paste("plot_simulations/posterior",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=8,height=12)
  dev.off()
  
  # - - - -
  # Plot attack rates, scaled by proportion alive

  par(mfrow=c(1,1))
  par(mar = c(5,5,1,1))
  if(flutype=="H3" & simDat==F){
      yob.data=data.frame(read.csv("datasets/HaNam_YOB.csv",header=FALSE)) # Import age distribution
      n.alive=sapply(inf_years,function(x){sum(yob.data<=x)})
  }else{
      yob.data=cbind(rep(1,n_part),rep(1,n_part)) # Import age distribution
      n.alive=n_part+0*inf_years
  }
    
  attack=colSums(historytabCollect[round(0.2*hist.sample*n_part):(hist.sample*n_part),])/(length(ind.infN)*(n.alive/length(yob.data[,1]))) #scale by proportion alive
  attackCI=NULL
  for(jj in 1:length(inf_years)){
    htest <- binom.test(round(n.alive*attack)[jj], n.alive[jj], p = 1,conf.level=0.95)
    meanA=attack[jj]
    conf1=htest$conf.int[1]
    conf2=htest$conf.int[2]
    attackCI=rbind(attackCI,c(meanA,conf1,conf2))
  }
  attackCI=data.frame(attackCI)
  names(attackCI)=c("mean","CI1","CI2")
  colA=rgb(0,0,0.8)
  plot(inf_years,attackCI$mean,pch=19,col=colA,ylim=c(0,1),xlab="year",ylab="attack rate")
  for(kk in 1:length(inf_years)){ # Iterate across test years
     lines(c(inf_years[kk],inf_years[kk]),c(attackCI$CI1[kk],attackCI$CI2[kk]),col=colA)
  }
  if(simDat==T){ #Add simulated attack rates
     #attack.yr=read.csv("datasets/sim_attack.csv")[,1]
     load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep=""))
     attack.yr=colSums(historytabSim)/n_part
     points(inf_years,attack.yr)
  }
    
  dev.copy(pdf,paste("plot_simulations/attack",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=6,height=6)
  dev.off()

  
  # - - - -
  # Plot antigenic map
  if(plotmap==T){
    load("datasets/spline_fn.RData") # load spline function for map
    ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
    
    par(mfrow=c(1,1))
    par(mar = c(5,5,1,1))
    
    vals1=predict(am.spl,scalemap(inf_years))
    map.sample=length(map.tabCollect)

    MTx=c(332,372)
    MTy=c(245,262)
    plot(vals1,type="l",xlim=MTx,ylim=MTy,col="white",lwd=2,xlab="strain dimension 1", xaxs="i", yaxs="i",ylab="strain dimension 2")
    
    #lines(inf_years,inf_years-min(inf_years),col="lightgray",lty=2)
    #plot(map.tabCollect[1,],0*map.tabCollect[1,],type="l",ylim=c(0,1),col="white",yaxt="n",ylab="",yaxs="i",xlab="strain position")
    for(ii in 1:100){
      #pick=sample(c(round(0.15*map.sample):map.sample),1)
      pick=sample(c(1:map.sample),1)
      vals1=predict(am.spl,scalemap(map.tabCollect[[pick]]))
      #sigma.scale=as.numeric(as.data.frame(thetatab)[pick*20,]["sigma"]) # scale to one antigenic unit
      #map.pick = map.tabCollect[[pick]]
      points(vals1,col=rgb(0.5,0.5,0.9,0.1),pch=19,cex=1)
    }
    
    points(ag.coord$AG_y,ag.coord$AG_x)

    dev.copy(pdf,paste("plot_simulations/antigenic_map",ifelse(simDat==T,"SIM",""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=6,height=6)
    dev.off()
  }

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

# plot.posterior.titres(loadseed="SIM",simDat=T,define.year=c(2007:2012))

plot.posterior.titres<-function(loadseed=1,year_test=c(2007:2012),flu.type,simDat=F,btstrap=5){
  
  if(simDat==F){
    if(flu.type=="H3"){
      load("R_datasets/HaNam_data.RData")
    }else{
      load("R_datasets/Fluscape_data_List.RData")
      year_test=c(2011,2012)
    }
    loadseed=1 # DEBUG
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
    pickAhist=ceiling(pickA/20)+1 # check which history this specifies
    hist.sample=historytabCollect[((pickAhist-1)*n_part+1):(pickAhist*n_part),1:n.inf]
    theta.max=as.data.frame(thetatab)[pickA,]
    
  for(pickyr in 1:n.test){
    
    # Note here that inf_years and strain_years are loads from main_model.R
    # Output expected titre - could include Poisson measurement uncertainty here?
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot simulated data

plot.sim.data<-function(){

  loadseed="SIM"
  load(paste("R_datasets/Simulated_data_",loadseed,"_1.RData",sep=""))
  
  #load("R_datasets/Simulated_dataPost_1.RData")
  
  #Compare model fits using posterior infection history (historytabPost) and parameters
  picktest=c(2007:2012)
  n.test=length(picktest)
  
  simulate_data(test_years=picktest,historytabPost=historytabSim,
                inf_years,
                strain_years,
                n_part=npartM,thetastar=thetaSim,p.inf=0.1) #theta.max
  
  load("R_datasets/Simulated_dataPost_1.RData")
  par(mfrow=c(2,5))
  par(mar = c(5,5,1,1))
  for(pickyr in 1:n.test){
    for(ii0 in 1:n_part){
      plot(8*historytabSim[ii0,],type="l",ylim=c(0,9),col='white')
      #for(jj in 1:length(lenhis)){
      #  lines(c(jj,jj),c(0,9*historytabSim[ii0,jj]),col='red')
      #}
      histSim1=historytabSim[ii0,]
      histSim1[inf_years>picktest[pickyr]]=0 # don't show infections after test year
      lenhis=rep(0,length(histSim1))
      
      for(jj in 1:length(lenhis)){
        lines(c(jj,jj),c(0,9*histSim1[jj]),col='blue')
      }
      lines(test.listSim[[ii0]][[pickyr]][4,],test.listSim[[ii0]][[pickyr]][2,],col=rgb(0.8,0.8,0.8),lwd=2) # Plot estimated infections
      points(test.listSim[[ii0]][[pickyr]][4,],test.listSim[[ii0]][[pickyr]][2,],pch=19)
      if(ii0 %% 10==0){
        dev.copy(pdf,paste("plot_simulations/simPlot",ii0,"P_",pickyr,".pdf",sep=""),width=12,height=6)
        dev.off()
      }
    }
  }

}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Rewind history and run with flat incidence

run.titre.time<-function(loadseed=1,year_test=c(2007:2012),flu.type="H3",simDat=F,btstrap=5,n_partSim=2,simTest.year=c(1968:2010)){
  
  if(simDat==F){
    if(flu.type=="H3"){
      load("R_datasets/HaNam_data.RData")
    }else{
      load("R_datasets/Fluscape_data_List.RData")
      year_test=c(2011,2012)
    }
    loadseed=paste(loadseed,"_",flu.type,sep="")
  }else{
    load(paste("R_datasets/Simulated_data_",loadseed,"_1.RData",sep=""))
  }
  
  load(paste("posterior_sero_runs/outputR_f",paste(year_test,"_",collapse="",sep=""),"s",loadseed,".RData",sep="")) # Note that this includes test.listPost
  
  lik.tot=rowSums(likelihoodtab)
  runsPOST=length(lik.tot[lik.tot!=-Inf])
  runs1=ceiling(0.2*runsPOST)
  
  # Set up matrices to store -- need btstrap >1
  strain_years=inf_years # look at strains from every year
  n.strains=length(strain_years) # this loads from main_model.R
  n.inf=length(inf_years)
  hist.sample0=rep(c(1,rep(0,29)),100)[1:n.inf] # CURRENTLY JUST FOR ONE PARTICIPANT
  simTest.year=sort(c(inf_years[hist.sample0==1],inf_years[hist.sample0==1]+1)) # infection year and one year after
  n.test=length(simTest.year)

  store.mcmc.test.data=array(NA, dim=c(btstrap,n_partSim,n.strains,n.test,2)) # Store expected titres for each test year
  store.mcmc.hist.data=array(NA, dim=c(btstrap,n_partSim,n.inf,n.test)) # Store history for each test year
  
  # - - - - - - - - - - - -
  # Sample from MCMC runs to get stored matrices of expected titre and estimated infection years
  
  for(sampk in 1:btstrap){
    pickA=sample(c(runs1:runsPOST),1)
    pickAhist=ceiling(pickA/20)+1 # check which history this specifies
    hist.sample=rbind(hist.sample0,hist.sample0)
    theta.max=as.data.frame(thetatab)[pickA,]
    
    for(pickyr in 1:n.test){ # ITERATE OVER TIME HERE
      
      # Note here that inf_years and strain_years are loads from main_model
      #hist.sample0[ inf_years<simTest.year[pickyr] ] # only take years up to test year -- already included in simulation function!
      simulate_data(simTest.year[pickyr],historytabPost=hist.sample,
                    inf_years,
                    strain_years,
                    n_partSim,thetastar=theta.max,p.inf=0.1,
                    #pmask=c("sigma2"), # For old fitted data, need to specify that sigma2 wasn't fitted
                    linD=F)
      
      load("R_datasets/Simulated_dataPost_1.RData")
      
      # Mask infections after test year
      for(ii0 in 1){
        
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
  # Plot development of titres
  
  par(mfrow=c(1,1)); par(mar = c(5,5,1,1))
  
  plot(inf_years,8*hist.sample[1,],type="l",ylim=c(0,9),col='white',xlab="year",ylab="titre")
  
  for(pickyr in 1:n.test){

    # Mask infections after test year
    
    for(ii0 in 1){
      simtitreX=store.mcmc.test.data[,ii0,,pickyr,1]
      simtitreY=store.mcmc.test.data[,ii0,,pickyr,2]
      hist.sample=store.mcmc.hist.data[,ii0,,pickyr] # for participant ii0 in year pickyr
      
      # Sample from infection history
      for(ksamp in 1:btstrap){
        for(jj in 1:n.inf){
          lines(min(inf_years)-1+c(jj,jj),c(-1,12*hist.sample[ksamp,jj]-1),col=rgb(0.8,0.8,0.8,0.01),lwd=2) # Plot estimated infections
        }
      }
      
      # Calculate credible interval for expected titres
      medP=apply(simtitreY,2,function(x){median(x)})
      ciP1=apply(simtitreY,2,function(x){quantile(x,0.025)})
      ciP2=apply(simtitreY,2,function(x){quantile(x,0.975)})
      polygon(c(simtitreX[1,],rev(simtitreX[1,])),c(ciP1,rev(ciP2)),lty=0,col=rgb(0,0.3,1,0.2))
      lines(simtitreX[1,],medP,pch=1,col='blue')
      points(simtitreX[1,],medP,pch=19,cex=0.5,col='blue')
      
      if(ii0 %% 10==0){
        dev.copy(pdf,paste("plot_simulations/sim",ii0,"P_",pickyr,".pdf",sep=""),width=12,height=6)
        dev.off()
      }
      
    }  # end loop over participants
    
  } # end loop over test years
  
  
}

