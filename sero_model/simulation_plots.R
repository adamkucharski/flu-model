
# - - -
# Plot true infections against simulated titres
par(mfrow=c(2,5))
par(mar = c(5,5,1,1))
for(ii0 in 1:n_part){
  lenhis=rep(0,length(historytabSim[ii0,]))
  plot(8*historytabSim[ii0,],type="l",ylim=c(0,9),col='white')
  for(jj in 1:length(lenhis)){
    lines(c(jj,jj),c(0,9*historytabSim[ii0,jj]),col='red')
  }
  lines(test.list[[ii0]][[1]][4,],test.list[[ii0]][[1]][2,],type="l")
  points(test.list[[ii0]][[1]][4,],test.list[[ii0]][[1]][2,],pch=19)
  if(ii0 %% 10==0){
    dev.copy(pdf,paste("plot_simulations/sim",ii0,".pdf",sep=""),width=12,height=6)
    dev.off()
  }
}

#historytabSim[ii0,]
#setuphistIC(ii0,jj_year,inf.n,test.list)

#thetaSim[["mu"]]=3
#thetaSim[["sigma"]]=0.25
dmatrix=outputdmatrix(thetaSim,inf_years) # Arrange parameters
#estimatelik(ii0,1,as.numeric(historytabSim[ii0,]),dmatrix,thetaSim,test.list)

#lik_val=NULL
#for(ii0 in 1:n_part){
#  lik_val=c(lik_val,estimatelik(ii0,1,as.numeric(historytabSim[ii0,]),dmatrix,thetaSim,test.list))
#}
#sum(lik_val)

# - - -
