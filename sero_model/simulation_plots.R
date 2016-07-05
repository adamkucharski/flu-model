







plot1<-function(){

  # This code is deprecated
  # - - -
  # Plot true infections against simulated titres for first strain
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


}



map1<-function(){
# - - - -
# Plot antigenic map
if(plotmap==T){
  par(mfrow=c(1,1))
  par(mar = c(5,5,1,1))
  map.sample=length(map.tabCollect)
  minMT=min(unlist(map.tabCollect)); maxMT=max(unlist(map.tabCollect))
  MTx=c(-5,15)
  MTy=c(-5,15)
  plot(map.tabCollect[[1]],xlim=MTx,ylim=MTy,type="l",col="white",lwd=2,xlab="strain dimension 1", xaxs="i", yaxs="i",ylab="strain dimension 2")
  #lines(inf_years,inf_years-min(inf_years),col="lightgray",lty=2)
  #plot(map.tabCollect[1,],0*map.tabCollect[1,],type="l",ylim=c(0,1),col="white",yaxt="n",ylab="",yaxs="i",xlab="strain position")
  for(ii in 1:200){
    #pick=sample(c(round(0.15*map.sample):map.sample),1)
    pick=sample(c(1:map.sample),1)
    map.pick = scale.map(map.tabCollect[[pick]])
    #sigma.scale=as.numeric(as.data.frame(thetatab)[pick*20,]["sigma"]) # scale to one antigenic unit
    #map.pick = map.tabCollect[[pick]]
    lines(map.pick,col=rgb(0.5,0.5,0.9,0.1),pch=19,cex=1)
  }
  
  if(simDat==T){
    map.pick = scale.map(antigenic.map.in)
    lines(map.pick,col="black",lwd=2)
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
  dev.copy(pdf,paste("plot_simulations/antigenic_map",ifelse(simDat==T,paste("mu",thetaSim[["mu"]],"_sigma",thetaSim[["sigma"]],sep=""),""),"_np",n_part,"_yr",paste(year_test,"_",collapse="",sep=""),loadseed,".pdf",sep=""),width=6,height=6)
  dev.off()
}
  
}