# Model of serological dynamics - uses PLOS Biology model (Kucharski et al. 2015)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data (Fonville et al.)

#setwd("~/Dropbox/Imperial/Fluscape_2/sero_model/")

source("load_data.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Pick sample year and plot

for(ii in 1:n_part){
  
  subjectn=ii
  alt_year=2011
  sample_yr=2012
  
  dataplot1=data1[data1$Subject.number==subjectn,]
  points1=dataplot1[dataplot1$Sample.year==sample_yr,]
  
  dataplot2=data1[data1$Subject.number==subjectn,]
  points2=dataplot2[dataplot2$Sample.year==alt_year,]
  
  sortyr=order(strain_years)
  
  #par(mar = c(5,4,4,2) + 0.1)
  plot(strain_years+runif(nstrains,0,0),as.numeric(points1[,strain_names])+runif(nstrains,0,0),xlab="test strain (year of circulation)",ylab="log titre",ylim=c(0,7),col=rgb(0,0,1),pch=19,cex=1)
  points(strain_years+runif(nstrains,0,0),as.numeric(points2[,strain_names])+runif(nstrains,0,0),col=rgb(1,0,0),pch=1,cex=1,lwd=2)
  #axis(1, at=c(1:nstrains), labels=strain_names[sortyr], las=2)
  grid(nx = NULL, ny = F, col = "lightgray",lty = "dotted")
  text(x=1970,y=7,labels=alt_year,col="red")
  text(x=1970,y=6,labels=sample_yr,col="blue")
  
  dev.copy(pdf,paste("data_figures/Plot_A",subjectn,"year_",alt_year,"_",sample_yr,".pdf",sep=""),width=8,height=6)
  dev.off()

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Animate

library(animation)

draw.curve<-function(sample_yr,subN){
  
  dataplot1=data1[data1$Subject.number==subN,]
  points1=dataplot1[dataplot1$Sample.year==sample_yr,]
  
  plot(strain_years,as.numeric(points1[,strain_names]),xlab="test strain (year of circulation)",ylab="log titre",xlim=c(1968,2012),ylim=c(0,7),col=rgb(0,0,1),pch=19,cex=0.7,main=paste("Participant",subN))
  grid(nx = NULL, ny = F, col = "lightgray",lty = "dotted")
  lines(c(sample_yr,sample_yr),c(-1,8),col=rgb(0,0.5,0),lwd=2)

}

draw.curve(sample_yr=2007,subN=subjectn)

trace.animate <- function() {
  lapply(seq(2007,2012), function(i) {
    draw.curve(i,subjectn)
  })
}

#save all iterations into one GIF
saveGIF(trace.animate(), interval = 1, movie.name=paste("data_figures/titres_",subjectn,".gif",sep=""))

