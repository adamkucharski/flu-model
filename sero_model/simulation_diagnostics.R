# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation diagnostics
# Compare MCMC output to simulation data

load("posterior_sero_runs/outputR.RData")
par(mfrow=c(3,3))
par(mar = c(5,5,1,1))
colA=rgb(0.8,0.8,0.8)

# Plot profile likelihood
lik.tot=rowSums(likelihoodtab)
maxlik=max(lik.tot)
runsPOST=length(lik.tot[lik.tot!=-Inf])
runs1=ceiling(0.5*runsPOST)
plot(rowSums(likelihoodtab)[runs1:runsPOST],type="l",ylab="likelihood",ylim=c(maxlik-50,maxlik))
plot(as.data.frame(thetatab)$mu[1:runsPOST],type="l",ylab="mu")
plot(as.data.frame(thetatab)$sigma[1:runsPOST],type="l",ylab="sigma")

# Plot histogram of boosting
hist(as.data.frame(thetatab)$mu[runs1:runsPOST],main=NULL,col=colA,xlab="mu",prob=TRUE,xlim=c(1,5))
abline(v=thetaSim[["mu"]],col="red")

hist(as.data.frame(thetatab)$sigma[runs1:runsPOST],main=NULL,col=colA,xlab="sigma",xlim=c(0,0.5))
abline(v=thetaSim[["sigma"]],col="red")

hist(as.data.frame(thetatab)$tau1[runs1:runsPOST],main=NULL,col=colA,xlab="tau1",prob=TRUE,xlim=c(0,0.5))
abline(v=thetaSim[["tau1"]],col="red")

plot(as.data.frame(thetatab)$mu[runs1:runsPOST],as.data.frame(thetatab)$sigma[runs1:runsPOST],pch=18,xlab="mu",ylab="sigma")

plot(as.data.frame(thetatab)$tau1[runs1:runsPOST],as.data.frame(thetatab)$sigma[runs1:runsPOST],pch=18,xlab="tau1",ylab="sigma")

plot(rowSums(historytab-historytabSim),type="l",ylab="history residual")

dev.copy(pdf,paste("plot_simulations/posterior_mu",thetaSim[["mu"]],"_sigma",thetaSim[["sigma"]],"_npart",n_part,".pdf",sep=""),width=12,height=8)
dev.off()


# UPDATE THIS BIT
#Compare model fits using posterior infection history (historytabPost) and parameters
simulate_data(test_years,historytabPost=historytab,
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

