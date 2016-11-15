
# This code is deprecated

# Calculate values based on two fold rise or more in that year
sconverttab=NULL

for(kk in 2:length(test.yr) ){ # Only valid for 2008-2011 (no test strains for 2011)
  pyear=0
  nyear=0
  for(ii in 1:n_part){
    t.part1=test.list[[ii]][[kk-1]]
    t.part2=test.list[[ii]][[kk]]
    
    if(length(t.part1[,1])>1 & length(t.part2[,1])>1){ # Check if year to compare
      # Check to match test strains
      matchd1d2 = intersect(names(t.part1),names(t.part2)[t.part2[3,]==test.yr[kk]])
      
      if(length(matchd1d2) > 0){
        
        diffT = t.part2[2,match(matchd1d2,names(t.part2))] - t.part1[2,match(matchd1d2,names(t.part1))] # Compare titres
        nyear = nyear +1
        if(max(diffT) >= 2){pyear = pyear + 1}
      }
      
    }
  }
  sconverttab=rbind(sconverttab, c(pyear/nyear,nyear))
}

attackCIsero=NULL
for(jj in 1:5){
  if(jj == 5){attackCIsero=rbind(attackCIsero,c(-1,-1,-1))}else{ # As NA in final year
    htest <- binom.test(round(sconverttab[jj,1]*sconverttab[jj,2]), sconverttab[jj,2], p = 1,conf.level=0.95)
    meanA=sconverttab[jj,1]
    conf1=htest$conf.int[1]
    conf2=htest$conf.int[2]
    attackCIsero=rbind(attackCIsero,c(meanA,conf1,conf2))
  }
}
attackCIsero=data.frame(attackCIsero)
names(attackCIsero)=c("mean","CI1","CI2")





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



# 2D plots
# Store bootstrap runs
storetitreS = array(NA,dim = c(length(year_x),length(cross_x)))
storetitreL = array(NA,dim = c(length(year_x),length(cross_x)))

for(ii in 1:length(year_x)){
  for(cc in 1:length(cross_x)){
    SS1=NULL
    SL1=NULL
    for(sampk in 1:btstrap){
      # - - - - - - - - - - - -
      # Sample from MCMC posterior to get trajectory
      
      pickA=sample(c(runs1:runsPOST),1)
      theta.max=as.data.frame(thetatab)[pickA,]
      
      SS1 = c(SS1, theta.max$mu*exp(-theta.max$sigma*abs(cross_x[cc])) )
      
      SL1 = c(SL1, min(8, theta.max$mu*exp(-theta.max$sigma*abs(cross_x[cc])) + theta.max$muShort*exp(-theta.max$sigma2*abs(cross_x[cc])) * exp(- year_x[ii] * theta.max$wane) ) )
      
    }
    
    storetitreS[ii,cc] = median(SS1)
    storetitreL[ii,cc] = median(SL1)
    
  }
}

storetitreS2=data.frame(melt(storetitreS))
names(storetitreS2)=c("x","y","z")

storetitreL2=data.frame(melt(storetitreL))
names(storetitreL2)=c("x","y","z")

par(mfrow=c(2,1))
p1 = ggplot(storetitreS2, aes(x, y, z = z)) +
  #stat_contour(aes(colour = ..level..))
  geom_raster(aes(fill = z)) +
  geom_contour(colour = "white")
#geom_tile(aes(fill = z)) #+ stat_contour(aes(colour = ..level..))

p2 = ggplot(storetitreL2, aes(x, y, z = z)) +
  theme(plot.title=element_text(hjust=0,face="bold"),panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(size = 1,rgb(0.95,0.95,0.95)),panel.background=element_blank(),panel.border=element_rect(colour = "black", fill=FALSE))+
  geom_raster(aes(fill = z)) +
  geom_contour(colour = "white")
#geom_tile(aes(fill = z)) #+ stat_contour(aes(colour = ..level..))

theme(plot.title=element_text(hjust=0,face="bold"),panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(size = 1,rgb(0.95,0.95,0.95)),panel.background=element_blank(),panel.border=element_rect(colour = "black", fill=FALSE))+
  xlab("age") +
  ylab("proportion seropositive 2015") +
  ylim(0,1)

grid.arrange(p1,p2,  ncol=1) #,main=paste(ifelse(wvacc==1,"With vaccine","Without vaccine"),sep=""))
