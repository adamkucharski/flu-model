# Simulate artificial serological data for DENV/ZIKV
# Use model with cross reaction

library(mvtnorm)

generate_timeseries<-function(strains,inf_years,n_part=20,p.inf.in=0.2,sd.val.in=0.5){

  inf.n=length(inf_years)
  #Set per year incidence, to create correlation between participant infection histories
  time_series=NULL
  for(ii in 1:strains){
    if(length(p.inf.in)>1){p.inf=p.inf.in[ii]}else{p.inf}
    if(length(sd.val.in)>1){sd.val=sd.val.in[ii]}else{sd.val=sd.val.in}
    log.sd=sd.val
    attack.yr=sapply(rlnorm(inf.n,meanlog=log(p.inf)-log.sd^2/2,sdlog=log.sd),function(x){min(x,1)})
    time_series=cbind(time_series,attack.yr)
  }
  time_series
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate serological data

simulate_sera_data<-function(strains,inf.years.sera=c(1980:2015),time.series.in=NULL,theta=c(error=0.05,sigma=0.4),p.inf.in=0.2,sd.val.in=1,seedi=1,roundv=F,dmatrix.in=NULL,zikv.attack=0.5){ # ii=participant | jj=test year
  
  # theta guide:
  # error = probability negative given infection (i.e. 1 - sensitivity)
  # sigma = probability positive from cross-reaction (i.e. 1 - specificity)
  # - - - - - - 
  # Set year of birth - sample uniformly
  #age.yr=sort(sample(1:50,n_part,replace = TRUE))
  inf.n=length(inf.years.sera)
  age.yr=c(0:(inf.n-1))
  n_part=inf.n

  # Generate annual attack rates
  if(is.null(time.series.in)){
    time.series=generate_timeseries(strains,inf.years.sera,n_part,p.inf.in,sd.val.in)
    
    # circulation in last year only for final strain, with attack rate = zikv.attack
    time.series[,strains]=0
    time.series[inf.n,strains]=zikv.attack
    
  }else{
    time.series=time.series.in
  }
  
  # - - - - - - 
  # Simulate random infection history for each participant
  
  historytabSim=NULL
  for(ii in 1:n_part){
    # Identify which years individual was alive for (ensure new borns can be infected)
    mask.alive=c(max(1,inf.n-(age.yr[ii])):inf.n)
    # Calculate probability of infection with each strain
    if(length(mask.alive)>1){
      prob.inf=1-apply(1-time.series[mask.alive,], 2, prod)
    }else{
      prob.inf=1-sapply(1-time.series[mask.alive,], prod)
    }
    historytabSim=rbind(historytabSim,prob.inf)
  }

  # Define cross reaction matrix
  if(is.null(dmatrix.in)){
    dmatrix=diag(strains)*(1-theta[["error"]])+(1-diag(strains))*theta[["sigma"]] # default cross-reaction structure
  }else{
    dmatrix=dmatrix.in
  }
  
  # - - - - - - 
  # Simulate assay results for each participant
  
  test.list.sera=NULL
  
  for(ii in 1:n_part){
    prob.inf=historytabSim[ii,]
  
    # Probability test positive to each strain
    prob.positive = 1-apply(1-t(dmatrix*prob.inf),1,prod)
    test.list.sera = rbind(test.list.sera,prob.positive)
  }
  
  #test.list.sera
  
  # - - - - - - 
  # Export data
  save(test.list.sera,historytabSim,inf.years.sera,strains,n_part,time.series,age.yr,theta,file=paste("R_datasets/Sero_sim_",seedi,"_",strains,".RData",sep=""))

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot simulated serological data

f.y<-function(x){rev(x)} # swap order of ages

plot_simulated_sera_data<-function(seedi,strains){

  load(paste("R_datasets/Sero_sim_",seedi,"_",strains,".RData",sep=""))
  
  par(mfrow=c(3,1))
  
  col.list=list(col1=rgb(0.9,0.6,0),col2=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
  label.age=seq(0,length(age.yr),2)
  lw.1=1.5
  
  # Plot attack rates
  plot(inf.years.sera,time.series[,1],type="l",ylim=c(0,1),col=col.list$col1,xlab="year",ylab="attack rate",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(inf.years.sera,time.series[,2],type="l",col=col.list$col2,lwd=lw.1)
  lines(inf.years.sera,time.series[,3],type="l",col=col.list$col3,lwd=lw.1)
  lines(inf.years.sera,time.series[,4],type="l",col=col.list$col4,lwd=lw.1)
  lines(inf.years.sera,time.series[,5],type="l",col=col.list$col5,lwd=lw.1)
  #grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  title(main=LETTERS[1],adj=0)

  # Plot probability of infection by age
  plot(f.y(historytabSim[,1]),ylim=c(0,1.01),col=col.list$col1,xaxt="n",xlab="age in 2015",ylab="probability ever infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(f.y(historytabSim[,2]),col=col.list$col2,lwd=lw.1)
  lines(f.y(historytabSim[,3]),col=col.list$col3,lwd=lw.1)
  lines(f.y(historytabSim[,4]),col=col.list$col4,lwd=lw.1)
  lines(f.y(historytabSim[,5]),col=col.list$col5,lwd=lw.1)
  axis(1,at=label.age,labels=rev(label.age))
  #grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  title(main=LETTERS[2],adj=0)
  
  # Plot probability of seropositivity by age
  plot(f.y(test.list.sera[,1]),ylim=c(0.2,1.01),col=col.list$col1,xaxt="n",xlab="age in 2015",ylab="probability seropositive",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(f.y(test.list.sera[,2]),col=col.list$col2,lwd=lw.1)
  lines(f.y(test.list.sera[,3]),col=col.list$col3,lwd=lw.1)
  lines(f.y(test.list.sera[,4]),col=col.list$col4,lwd=lw.1)
  lines(f.y(test.list.sera[,5]),col=col.list$col5,lwd=lw.1)
  axis(1,at=label.age,labels=rev(label.age))
  title(main=LETTERS[3],adj=0)
  #grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  
  dev.copy(pdf,paste("plot_simulations/serology_plot.pdf",sep=""),width=7,height=10)
  dev.off()
  
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# MCMC Functions

SampleTheta<-function(theta_initial,m,covartheta,nparam,n_part){
  
  # sample from multivariate normal distribution - no adaptive sampling
  theta_star = as.numeric(exp(mvrnorm(1,log(theta_initial), Sigma=covartheta)))
  names(theta_star)=names(theta_initial)
  
  if(!is.na(match("error",names(theta_star)))){ # Check whether fitting theta or attack vector
    mu1=min(2-theta_star[["sigma"]],theta_star[["sigma"]])
    theta_star[["sigma"]]=ifelse(mu1<0,theta_initial[["sigma"]],mu1)
    
    mu1=min(2-theta_star[["error"]],theta_star[["error"]])
    theta_star[["error"]]=ifelse(mu1<0,theta_initial[["error"]],mu1)
  }else{
    #print("check")
    # Ensure FoI is below 3 (i.e. <95% probability infection per year)
    theta_check=sapply(theta_star,function(x){min(6-x,x)})
    theta_star=theta_check*(as.numeric(theta_check)>0) + theta_initial*(as.numeric(theta_check)<0)
  }
  
  theta_star
}

ComputeProbability<-function(marg_likelihood,marg_likelihood_star){
  # Flat priors on theta => symmetric update probability
  calc.lik = exp(marg_likelihood_star-marg_likelihood)
  min(1, calc.lik)
}

LikelihoodTitres<-function(titre.data,dmatrix,forcetab_star,inf.n,strains,n_sample){
  
  lik=NULL

  for(ii in 1:inf.n){ # Iterate across years
    p.inf=NULL
    for(jj in 1:strains){ # Iterate across strains
      p.inf=c(p.inf,1-exp(-sum(forcetab_star[jj,1:ii]))) # Infected in this period with strain j
    }
    # Calculate serological likelihood
    prob.positive = 1-apply(1-t(dmatrix*p.inf),1,prod)
    lik=rbind(lik,dbinom(round(n_sample[,ii]*titre.data[,ii]),size=n_sample[,ii],prob=prob.positive,log=T))
  }
  
  sum(lik)
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run MCMC function

run_mcmc<-function(
  titre.data, # Note that this starts present and goes back into past (i.e. ordered increasing with age). Data as proportions
  n_sample, # Total samples in each age group and each strain
  inf_years,
  theta,
  strains=5,
  prob_inf=NULL, # initial conditions
  force_constrain, # matrix to constrain circulation years - Note that this is increasing age
  switch1=2,
  runs,
  seedi=1
){

  inf.n=length(inf_years)
  # Preallocate memory
  nparam=length(theta); npcov=rep(1,nparam)
  cov_matrix_theta0 = diag(npcov)
  cov_matrix_force0 = diag(rep(1,inf.n))
  
  thetatab=matrix(NA,nrow=(runs+1),ncol=length(theta)); colnames(thetatab)=names(theta)
  thetatab[1,]=theta
  
  if(!is.null(prob_inf)){
    # Convert attack rate to FoI
    #diff=as.matrix(prob_inf)-as.matrix(cbind(rep(0,strains),prob_inf[,1:(length(prob_inf[1,])-1)]))
    forcetab=-log(1-prob_inf)
  }else{
    forcetab=matrix(0.01,nrow=strains,ncol=inf.n)
  }
  forcetab=force_constrain*forcetab
  forcetabCollect=forcetab

  # Preallocate matrices
  likelihoodtab=matrix(-Inf,nrow=(runs+1))
  accepttabT=NULL
  accepttabH=NULL
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Run MCMC
  
  for (m in 1:runs){
    
    # Adaptive covariance matrix
    if(m==1){
      epsilon0=0.01
      varpart_prob0=0.01
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_force=varpart_prob0*cov_matrix_force0
    }else{
      epsilon0=max(1e-5,min(1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m)))
      varpart_prob0=max(1e-5,min(1,exp(log(varpart_prob0)+(accept_rateH-0.234)*0.999^m))) # force of infection sampling
      
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_force=varpart_prob0*cov_matrix_force0  # ***DEBUGGING***
    }
    
    # - - - - - - - - - - - - - - - -
    # Resample parameters
    
    if(m %% switch1==0 | m==1){ # m==1 condition as have to calculate all liks on first step
      theta_star = SampleTheta(thetatab[m,], m,cov_matrix_theta,nparam=sum(cov_matrix_theta0)) #resample theta
      forcetab_star=forcetab
    }else{
      forcetab_star = t(apply(forcetab,1,function(x){SampleTheta(x, m,cov_matrix_force,nparam=sum(cov_matrix_force))})) #resample theta)
      theta_star = thetatab[m,]
    }
    
    dmatrix = diag(strains)*(1-theta_star[["error"]])+(1-diag(strains))*theta_star[["sigma"]] # default cross-reaction structure
    
    # - - - - - - - - - - - - - - - -
    # LIKELIHOOD function - Only calculate for updated history
    
    lik_val = LikelihoodTitres(titre.data,dmatrix,forcetab_star,inf.n,strains,n_sample)
    # - - - - - - - - - - - - - - - -
    # Metropolis Hastings step
    
    output_prob = ComputeProbability(sum(likelihoodtab[m,]),sum(lik_val)) 
    
    #if(is.na(output_prob)){print(c(m,sum(likelihoodtab[m,]),sum(lik_val),theta_star))} # Print likelihood (For DEBUG)}
    if(is.na(output_prob) & m==1){stop('check initial parameter values')}
    
    if(runif(1) < output_prob){
      thetatab[m+1,] = theta_star
      if(m %% switch1!=0){forcetab = forcetab_star} # Only change if resampled
      likelihoodtab[m+1] = lik_val
      if(m %% switch1==0){accepttabT=c(accepttabT,1)}
      if(m %% switch1!=0){accepttabH=c(accepttabH,1)}
      
    }else{
      thetatab[m+1,] = thetatab[m,]
      likelihoodtab[m+1] = likelihoodtab[m]
      if(m %% switch1==0){accepttabT=c(accepttabT,0)}
      if(m %% switch1!=0){accepttabH=c(accepttabH,0)}
    }
    
    
    if(m<max(100)){
      accept_rateT=0.234 # target acceptance rate for theta
      accept_rateH=0.234 # target acceptance rate for infection history
    }else{
      accept_rateT=sum(accepttabT)/length(accepttabT)
      accept_rateH=sum(accepttabH)/length(accepttabH)
    }
    
    if(m %% min(runs,20) ==0){
      forcetabCollect=rbind(forcetabCollect,forcetab)
    }
    
    if(m %% min(runs,1000) ==0){
      print(c(m,accept_rateT,varpart_prob0,likelihoodtab[m]))
      save(likelihoodtab,thetatab,inf_years,inf.n,strains,titre.data,forcetab,forcetabCollect,switch1,file=paste("posterior_runs/outputR_f",seedi,".RData",sep=""))
    }
    
  } #End runs loop
  
}

