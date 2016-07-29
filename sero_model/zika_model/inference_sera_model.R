# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run inference from simulated serological data
# Author: AJ Kucharski (2016)

# Convert simulated samples into probabilities
load("R_datasets/Sero_sim_1_5.RData")
per_sample=100
n_sample=matrix(per_sample,nrow=n_strains,ncol=inf.n)
titre.prob=t(test.list.sera)

prob_inf0=t(time.series) # Fix initial conditions to correct history
prob_inf0=as.matrix(prob_inf0[,c(length(prob_inf0[1,]):1)])

titre.data=NULL
for(jj in 1:strains){
  tstrain=NULL
  for(kk in 1:length(inf.years.sera)){
    rand=rbinom(1,size=n_sample[jj,kk],prob=titre.prob[jj,kk])/per_sample
    tstrain=cbind(tstrain,rand)
  }
  titre.data=rbind(titre.data,tstrain)
}

# Impose ICs and run MCMC
force_constrain=matrix(1,nrow=strains,ncol=length(inf.years.sera))
force_constrain[strains,2:length(inf.years.sera)]=0 # Constrain ZIKV to final year only

run_mcmc(
  titre.data, # Note that this starts present and goes back into past (i.e. ordered with age). Data as proportions
  n_sample, # Total samples in each age group and each strain
  inf.years.sera,
  theta=c(error=0.05,sigma=0.5),
  strains,
  prob_inf=prob_inf0, # initial conditions
  force_constrain, # matrix to constrain circulation years - Note that this is increasing age
  runs=1e5,
  switch1=5,
  seedi=1
)


load("posterior_runs/outputR_f1.RData")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot posteriors

par(mfrow=c(2,4))
colA=rgb(0.8,0.8,0.8)
col.list=list(col1=rgb(0.9,0.6,0),col2=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
lw.1=1.5
label.age=seq(0,length(age.yr),2)

# Plot posteriors
maxlik=(likelihoodtab==max(likelihoodtab))
run2=length(likelihoodtab)
run1=round(0.2*run2)
plot(likelihoodtab[run1:run2],type="l")
theta_post=data.frame(thetatab[run1:run2,])
hist(1-theta_post[["sigma"]],main="",col=colA,xlab="specificity",prob=TRUE,xlim=c(0.5,1)); abline(v=1-theta.serology[["sigma"]],col="red")
hist(1-theta_post[["error"]],main="",col=colA,xlab="sensitivity",prob=TRUE,xlim=c(0.5,1)); abline(v=1-theta.serology[["error"]],col="red")

# - - - 
# Plot infection curves -  last sample from MCMC against observed data

siml_force=t(time.series); 

nblock=length(forcetabCollect)/(length(inf_years)*strains) # get blocks
post_force=forcetabCollect[((nblock-1)*strains+1):(nblock*strains),]
#post_force=1-exp(-post_force) # Convert FoI to probability
#post_force=post_force[,c(length(inf_years):1)] # Swap arrangement to match

lik=NULL

for(ii in 1:(inf.n)){ # Iterate across years
  p.inf=NULL
  for(jj in 1:strains){ # Iterate across strains
    p.inf=c(p.inf,1-exp(-sum(post_force[jj,1:(ii)]))) # Infected in this period with strain j
  }
  # Convert to probability of infection
  prob.positive = p.inf # 1-apply(1-t(dmatrix*p.inf),1,prod)
  lik=rbind(lik,prob.positive)
}

plot(f.y(historytabSim[,1]),ylim=c(0,1.01),col=col.list$col1,xaxt="n",xlab="age in 2015",ylab="probability ever infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
lines(f.y(lik[,1]),col=col.list$col1,lwd=lw.1,lty=2)
axis(1,at=label.age,labels=rev(label.age))

plot(f.y(historytabSim[,2]),ylim=c(0,1.01),col=col.list$col2,xaxt="n",xlab="age in 2015",ylab="probability ever infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
lines(f.y(lik[,2]),col=col.list$col2,lwd=lw.1,lty=2)
axis(1,at=label.age,labels=rev(label.age))

plot(f.y(historytabSim[,3]),ylim=c(0,1.01),col=col.list$col3,xaxt="n",xlab="age in 2015",ylab="probability ever infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
lines(f.y(lik[,3]),col=col.list$col3,lwd=lw.1,lty=2)
axis(1,at=label.age,labels=rev(label.age))

plot(f.y(historytabSim[,4]),ylim=c(0,1.01),col=col.list$col4,xaxt="n",xlab="age in 2015",ylab="probability ever infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
lines(f.y(lik[,4]),col=col.list$col4,lwd=lw.1,lty=2)
axis(1,at=label.age,labels=rev(label.age))

plot(f.y(historytabSim[,5]),ylim=c(0,1.01),col=col.list$col5,xaxt="n",xlab="age in 2015",ylab="probability ever infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
lines(f.y(lik[,5]),col=col.list$col5,lwd=lw.1,lty=2)
axis(1,at=label.age,labels=rev(label.age))

dev.copy(pdf,paste("plot_simulations/posterior_plot.pdf",sep=""),width=10,height=6)
dev.off()
