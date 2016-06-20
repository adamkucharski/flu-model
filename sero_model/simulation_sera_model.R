# Simulate artificial serological data for DENV/ZIKV
# Use model with cross reaction

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
  # error = probability negative given infection
  # sigma = probability positive from cross-reaction
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

  #plot(time.series[,1],type="l")

  # - - - - - - 
  # Simulate random infection history for each participant
  
  historytabSim=NULL
  for(ii in 1:n_part){
    # Identify which years individual was alive for (+1 to ensure new borns can be infected)
    mask.alive=c(max(1,inf.n-(age.yr[ii]+1)):inf.n)
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

plot_simulated_sera_data<-function(seedi,strains){

  load(paste("R_datasets/Sero_sim_",seedi,"_",strains,".RData",sep=""))
  
  par(mfrow=c(3,1))
  
  col.list=list(col1=rgb(0.9,0.6,0),col2=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
  
  f.y<-function(x){rev(x)} # swap order of ages
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


