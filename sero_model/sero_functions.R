# Model of serological dynamics - uses extended PLOS Biology model (Kucharski et al. 2015)
# Author: AJ Kucharski (2015-)
# Functions

# - - - - - - - - - - - - - - - - 
# Load antigenic map data
# Set up antigenic coords

load.flu.map.data<-function(){
  
  ag.coord=read.csv("datasets/antigenic_coords.csv", as.is=T,head=T)
  strain_names=ag.coord$viruses
  
  # Convert into years for each strain
  strain_years=as.numeric(sapply(strain_names,function(x){
    a1=max(which(strsplit(x, "")[[1]]=="/"))
    lstr=nchar(x)
    yr1=substr(x, a1+1, lstr)
    
    if(nchar(yr1)>4){yr1=substr(yr1, 1, 4)}
    year=yr1
    if(nchar(yr1)==2 & as.numeric(yr1)>15){year=paste("19",yr1,sep="")}
    if(nchar(yr1)==2 & as.numeric(yr1)<15){year=paste("20",yr1,sep="")}
    year
  }
  ))
  
  # Convert antigenic coords into cluster centroids
  clust1=c("HK68","EN72","VI75","TX77","BK79","SI87","BE89","BE92","WU95","SY97","FU02","CA04","WI05","PE09")
  clust2=c("BI/16190/68","BI/21793/72","BI/1761/76", "BI/2271/76","NL/233/82", "NL/620/89","NL/823/92","A/BEIJING/32/92","WU/359/95","A/SYDNEY/5/97","FU/411/02","A/CALIFORNIA/7/2004","A/WISCONSIN/67/2005","A/PERTH/16/2009")
  
  clust.names=data.frame(cbind(clust1,clust2),stringsAsFactors=F)
  
  ag.coord1=ag.coord[match(clust.names$clust2,ag.coord$viruses),] # pick representative cluster strains
  
  # Plot strains vs clusters
  ag.coord=ag.coord[ag.coord$viruses!="NL/823/92",] # remove BE 89 dead end (Fonville method)
  am.spl=smooth.spline(ag.coord$AG_y,ag.coord$AG_x)
  xx=c(330:370)
  prd1=predict(am.spl,xx)
  save(am.spl,file="datasets/spline_fn.RData")
}


# - - - - - - - - - - - - - - - -
# Set initial conditions - with small number of ICs

setuphistIC_SA<-function(ii,jj,inf.n,test.list,testyear_index, test_years, inf_years){ # ii=participant | jj=test year
  
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  spyear=unique(as.numeric(test.jj[3,])) # year of samples taken
  
  hist0=rep(0,inf.n)   

  # Check test data available - may be issue if age column added too
  if(length(test.jj[,1])>1){
    
    # Set up test strains
    titredat=as.numeric(test.jj[2,]) # Define titre data
    maxt=(titredat==max(titredat))
    
    # Use simple cutoff for titres -- set high titres = 1 in history
    for(i in 1:length(spyear)){
      if(max(titredat[(as.numeric(test.jj[3,])==spyear[i])])>=4 & runif(1)>0.1 ){
        hist0[(inf_years==spyear[i])]=1
      }
    }
  }

  # - - 

  min.range = max(1,testyear_index[1]-5) # Add infection within past 5 years
  inf_index = inf_years-min(inf_years)+1
  inf_pick = sample(c(1:inf.n)[inf_index>=min.range & inf_index<= testyear_index[1]],1)  # pick strain within plausible region to add
  if(sum(hist0[inf_years < min(test_years)])==0){hist0[inf_pick]=1} # Make sure at least one infection previous to test year
  hist0
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Likelihood given infection history and parameters (done in C)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Functions to set up parameters for model
# Calculate dmatrix by year
outputdmatrix<-function(theta,inf_years,linearD=F,locmat=NULL){
  if (is.null(locmat)) {
    # Exponential decay
    (dmatrix=sapply(inf_years,function(x){exp(-theta[["sigma"]]*abs(inf_years-x))})) # note that second entry refers to sample year in titre calc
    # Linear decay
    if(linearD==T){
      (dmatrix=sapply(inf_years,function(x){y=1-1000*theta[["sigma"]]*abs(inf_years-x)/inf_years; y[y<0]=0; y }) ) # note that second entry refers to sample year in titre calc
    }
  } else {
    stop("Non-null locmat not yet implemented in outputdmatrix")
    # Up to here on the distance matrix
    (dmatrix=sapply(inf_years,function(x){exp(-theta[["sigma"]]*abs(inf_years-x))}))
  }
  dmatrix
}

# Calculate dmatrix from antigenic map data (either fitted or specified) - **This makes above function redundant**

scalemap<-function(xx,inf_years){
  if(max(inf_years)>2012){stop("need infection range to be inside antigenic map")}
  map.range <- c(1968:2012)
  alen <- c(333.83,370.28); alenA <- (alen[2]-alen[1])/(max(xx)-min(xx)); alenB <- alen[1]-alenA*min(xx)
  s1 <- alenA*xx+alenB-alen[1]; s1*length(inf_years)/length(map.range) +alen[1]
}

outputdmatrix.fromcoord <- function(thetasigma,inf_years,anti.map.in,linearD=F){ #anti.map.in can be vector or matrix - rows give inf_years, columns give location

  # Check if map is 1D or 2D
  #if(length(anti.map.in)==length(inf_years)){
  if(linearD==F){
    # Exponential decay function
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
      (dmatrix=sapply(anti.map.in,function(x){exp(-thetasigma*abs(anti.map.in-x))}))
    }else{ # If spline function defined, calculate directly from input
      (dmatrix=apply(anti.map.in,1,function(x){exp(-thetasigma*sqrt(
        colSums(apply(anti.map.in,1,function(y){(y-x)^2}))
        ))}))
    }
  }else{
    # Linear decay function
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
      (dmatrix=sapply(anti.map.in,function(x){y=abs(anti.map.in-x); y   })) 
      
    }else{ # If spline function defined, calculate directly from input
      (dmatrix=apply(anti.map.in,1,function(x){y=sqrt(
        
        colSums(apply(anti.map.in,1,function(y){(y-x)^2}))
        
      ); y # 1-1*thetasigma* //  y[y<0]=0; HAVE REMOVED BASE
      }))
    }
  }
}

# Generate random walk 2D antigenic map 
generate.antigenic.map <- function(inf_years){
  
  map=matrix(0,ncol=2,nrow=length(inf_years))
  for(ii in 2:length(inf_years)){
    map[ii,]=map[ii-1,]+ (1 - 1.5*runif(2)) 
  }
  map
}

# - - - - - - - - - - - - - - - -
# Compile c code
compile.c<-function(){
  require("Rcpp")
  setwd("./c_code")
  #system("R CMD SHLIB c_model2.c")
  system("R CMD SHLIB c_model2_sr.c")
  #dyn.load("c_model2.so")
  dyn.load("c_model2_sr.so") # Note edit to remove ./ for cluster runs
  # sourceCpp("./cpp_steven.cpp")
  setwd("..")
}

# - - - - - - - - - - - - - - - -
# Define expected titre function

func1 <- function(x,titredat,dd,dd2,theta,testyear_index) {
  if (!is.numeric(x)){stop("argument x must be numeric")}
  out <- .C("c_model2_sr",
            n=as.integer(length(x)),
            itot=as.integer(sum(x)),
            nsample=as.integer(length(titredat)),
            x=as.double(x),
            x1=as.double(rep(0,length(x))),
            titre=as.double(titredat),
            titrepred=as.double(rep(0,length(titredat))),
            dd=as.double(dd),
            dd2=as.double(dd2),
            ntheta=as.integer(length(theta)),
            theta=as.double(theta),
            inputtestyr=as.integer(testyear_index)
  )
  # browser()
  # out$titrepred - out2$titrepred
  return(out$titrepred)
}

# - - - - - - - - - - - - - - - -
# Define likelihood function given expected titre and data - include uniform error term

likelihood.titre<-function(expect,titredat,theta){
  
  # DEBUG   expect=runif(20,min=0,max=10); titredat=floor(expect); theta[["error"]] = 0.1
  
  largett=(titredat > 8)  # Identify censored titres in data (>=8)
  smalltt=(titredat <= 0)  # Identify censored titres in data (>=8)
  
  
  # First sum up titres 0 < . < 8
  p_jkMID =  ( sum( log(  pnorm(as.numeric(titredat[!largett & !smalltt])+1, mean = expect[!largett & !smalltt], sd=theta[["error"]], log = FALSE) 
                     - pnorm(as.numeric(titredat[!largett & !smalltt]), mean = expect[!largett & !smalltt], sd=theta[["error"]], log = FALSE)  ) ) )
  
  # Add titres >=8
  p_jkSML = sum( (  pnorm(1, mean = expect[smalltt], sd=theta[["error"]], log = T) ) )
  
  # Add titres = 0
  p_jkLRG = sum( (  pnorm(8, mean = expect[largett], sd=theta[["error"]], log = T, lower.tail = F) ) ) 
  
  p_jk = p_jkSML + p_jkMID + p_jkLRG

  p_jk

}

# - - - - - - - - - - - - - - - -
# Calculate likelihood for given participant and test year

estimatelik<-function(ii,jj,historyii,dmatrix,dmatrix2,theta_star,test.list,testyearI){ # ii=participant | jj=test year
  
  # jj=jj_year[kk];historyii=as.numeric(history_star[ii,]);testyearI=testyear_index[kk]
  
  test.II=test.list[[ii]]
  test.jj=test.II[[jj]]
  
  # Check test data available
  if(length(test.jj[,1])==1){0}else{
    
    # Set up test strains
    test.part=as.numeric(test.jj[4,]) # index of sample strains data available for
    titredat=test.jj[2,] # Define titre data
    
    d.ij=dmatrix[test.part,] # Define cross-immunity matrix 1 for sample strain
    d_vector=melt(t(d.ij))$value #melt is by column
    
    d.ij2=dmatrix2[test.part,] # Define cross-immunity matrix 2 for sample strain
    d_vector2=melt(t(d.ij2))$value #melt is by column

    expect=func1(historyii,titredat,d_vector,d_vector2,theta_star,testyearI) # Output expectation

    lik = likelihood.titre(expect,titredat,theta_star)
    
    
    lik
    
  }
  
}


# - - - - - - - - - - - - - - - -
# Simulation infection history data

simulate_data<-function(test_years,
                        historytabPost=NULL, # This imposes a particular history
                        inf_years,strain_years,n_part=20,thetastar=theta0,p.inf=0.2,seedi=1,
                        roundv=F, # round expected titres to nearest integer?
                        linD=F, # use linear cross-reaction function?
                        antigenic.map.in=NULL,pmask=NULL,am.spline=NULL){ # ii=participant | jj=test year
  
  # DEBUG pickyr=1; test_years=test.yr[pickyr]; historytabPost=hist.sample; thetastar=theta.max; p.inf=0.1; antigenic.map.in=NULL ; linD=T; pmask=NULL
  
  # Variables needed: test_years,inf_years,strain_years,n_part
  #strain_years=seq(1968,2010,4)
  
  # Make adjustments depending on what is fitted and not
  if(sum(pmask=="muShort")>0){thetastar[["muShort"]]=1e-10} # Set short term boosting ~ 0 if waning not fitted
  if(sum(pmask=="map.fit")>0){ thetastar[["sigma"]]=1} # Set cross-reactivity = 1 and don't fit if antigenic map also fitted (to avoid overparameterisation)
  if(sum(pmask=="sigma2")>0){ thetastar[["sigma2"]]=thetastar[["sigma"]] } # Fix equal if sigma same for both 
  
  # Set year of birth
  age.yr=sample(1:80,n_part,replace = TRUE)
  
  test.n=length(test_years)
  inf.n=length(inf_years)
  nstrains=length(strain_years)
  sample.index=strain_years-min(inf_years)+1
  theta.sim.out=thetastar
  historytabSim2=historytabPost
  
  # Check inputs are correct
  if(sum(max(test_years)==inf_years)==0){
    stop("need infection years >= test years")
    return
  }
  
  # Define antigenic map
  xx=scalemap(inf_years,inf_years)
  yy=predict(am.spl,xx)$y
  antigenic.map.in=cbind(xx,yy)
  
  if(is.null(antigenic.map.in)){antigenic.map.in=inf_years} # If no specified antigenic map, use linear function by year
  
  # NOTE HARD CODED FOR LINEAR FUNCTION
  dmatrix = 1- thetastar[["sigma"]]*outputdmatrix.fromcoord(thetastar[["sigma"]],inf_years,antigenic.map.in,linearD=linD)
  dmatrix[dmatrix<0]=0
  dmatrix2 = 1- thetastar[["sigma2"]]*outputdmatrix.fromcoord(thetastar[["sigma2"]],inf_years,antigenic.map.in,linearD=linD)
  dmatrix2[dmatrix2<0]=0
  
  #Set per year incidence, to create correlation between participant infection histories
  log.sd=1
  if(length(p.inf)==1){
    attack.yr=rlnorm(inf.n,meanlog=log(p.inf)-log.sd^2/2,sdlog=log.sd)
  }else{
    attack.yr=p.inf
  }
  
  # Simulate random infection history for each infection year
  if(is.null(historytabPost)){
    historytabSim=matrix(0,ncol=inf.n,nrow=n_part)
    for(ii in 1:inf.n){
      #hist0=(runif(inf.n)<attack.yr)+0
      #alive=((max(test_years)-age.yr[ii])<=inf_years) - ignore age structure in simulations
      historytabSim[sample(n_part,round(n_part*attack.yr[ii])),ii]=1
    }
  }else{
    historytabSim=historytabPost
  }

  # Simulate titres for each participant
  
  test.list=list()
  
  for(ii in 1:n_part){
    
    subjectn=ii
    i.list=list()
    historyii=historytabSim[ii,]
    
    for(jj in 1:test.n){
      
      d.ij=dmatrix[sample.index,] # Define cross-immunity matrix for sample strain
      d_vector=melt(t(d.ij))$value
      
      d.ij2=dmatrix2[sample.index,] # Define cross-immunity matrix for sample strain
      d_vector2=melt(t(d.ij2))$value
      
      testyr=test_years[jj]
      testyearI=c(1:inf.n)[inf_years==testyr]
      
      expect=func1(historyii,sample.index,d_vector,d_vector2, thetastar,testyearI) # Output expectation
      
      if(roundv==T){
        titredat=sapply(expect,function(x){ floor( rnorm(1,mean=x,sd=thetastar[["error"]]) ) })
        titredat[titredat<0]=0
      }else{
        titredat=expect} # Generate test titre
      titredat=sapply(titredat,function(x){min(x,8)})

      i.list[[jj]]=rbind(test.year=rep(testyr,nstrains),
                         titredat,
                         strain_years,
                         sample.index
      )
    }
    #i.list[[1]][2,]

    test.list[[ii]]=i.list
  }
  test.listSim=test.list
  
  # Export data
  #browser()
  if(is.null(historytabPost)){
    save(test_years,inf_years,strain_years,n_part,test.listSim,theta.sim.out, age.yr,antigenic.map.in,historytabSim,file=paste("R_datasets/Simulated_data_",seedi,".RData",sep=""))
  }else{
    save(test_years,inf_years,strain_years,n_part,test.listSim,theta.sim.out, age.yr,antigenic.map.in,historytabSim2,file=paste("R_datasets/Simulated_dataPost_",seedi,".RData",sep=""))
  }
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Resample infection history - included ageA table in case needed later


SampleHistory<-function(historyA,pick,inf.n,ageA,inf_years,age.mask){

  for(ii in pick){

    #ls_pick=foreach(ii=(1:length(pick))) %dopar% {  # Parallel loop - slower to farm out
    rand1=runif(1)
    x=historyA[ii,age.mask[ii]:inf.n]
    
    infvector=c(1:length(x))
    infvector2=rev(infvector)
    
    # Remove infection
    if(rand1<1/3){
      infectID = infvector[(as.numeric(x)>0)]
      if(length(infectID)>0){
        x[sample(c(infectID),1)]=0 # 
      }
    }
    
    # Add new infection
    if(rand1>1/3 & rand1<2/3){
      ninfecID=infvector[(as.numeric(x)==0)]
      if(length(ninfecID)>0){
        x[sample(c(ninfecID),1)]=1
      }
    }
    
    # Move infection
    if(rand1>2/3){
      infectID=infvector[(as.numeric(x)>0)]
      ninfecID=infvector[(as.numeric(x)==0)]
      
      if(length(infectID)>0 & length(ninfecID)>0){
        x[sample(c(infectID),1)]=0
        x[sample(c(ninfecID),1)]=1
      }
    }

    
    historyA[ii,age.mask[ii]:inf.n]=x
    
  } # end loop over individuals
  historyA
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Resample age - add 1, 0, -1 with equal probability -- NOTE THIS IS NOT CURRENTLY ACTIVE

SampleAge<-function(pick,ageA){
  
  b1=sapply(ageA[pick],function(x){x+sample(c(-1:1),1)})
  ageA[pick]=b1
  ageA
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Resample antigenic location -- NOTE THIS IS NOT CURRENTLY ACTIVE

SampleAntigenicMap<-function(anti.map.star,epsilon.map,inf_years){
  Sigma0=(diag(1+0*inf_years))*epsilon.map
  #Sigma0[1,1]=0 # Anchor initial point
  if(length(anti.map.star)==length(inf_years)){
    sort(as.numeric(mvrnorm(1,anti.map.star, Sigma=Sigma0)))
  }else{
    apply(anti.map.star,2, function(x){as.numeric(mvrnorm(1,x, Sigma=Sigma0))})
  }
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Convert infection history to binary -- NOTE THIS IS NOT CURRENTLY ACTIVE

convert_binary <- function(x){sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

ComputeProbability<-function(marg_likelihood,marg_likelihood_star,theta_star, thetatabM,pmask,turn_off_likelihood){
  
  if(turn_off_likelihood==1){
    marg_likelihood = -1
    marg_likelihood_star = -1
  }
  
  # Identify fitted parameters
  match0 = match(pmask,names(theta_star)); match0 = match0[!is.na(match0)]
  
  # Flat priors on theta => symmetric update probability

  # Include comparison of values to account for lognormal resampling of theta
  calc.lik = exp(marg_likelihood_star-marg_likelihood + sum(log(theta_star[-match0])) - sum(log(thetatabM[-match0]))) # *theta_star[-match0][3] /thetatabM[-match0][3]  #+ sum(log(thetatabM[-match0])-log(theta_star[-match0])) )
  
  calc.lik[calc.lik>1]=1
  calc.lik
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample theta without reflective boundary conditions

SampleTheta2 <-function(theta_initial,m,covartheta,covarbasic,nparam,pmask){

  # sample from multivariate normal distribution - no adaptive sampling
  #theta_star = as.numeric(mvrnorm(1,theta_initial, Sigma=covarbasic))
  theta_star = as.numeric(exp(mvrnorm(1,log(theta_initial), Sigma=covarbasic)))
  
  names(theta_star)=names(theta_initial)

  likelihoodOK = 1
  
  # Check parameters biologically plausible
  if(theta_star[["wane"]]>1 & sum(pmask=="wane")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
  if(theta_star[["mu"]]>10 & sum(pmask=="mu")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
  if(theta_star[["muShort"]]>10 & sum(pmask=="muShort")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
  if(theta_star[["error"]]>10 & sum(pmask=="error")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
  if(theta_star[["sigma"]]>10 & sum(pmask=="sigma")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
  if(theta_star[["sigma2"]]>10 & sum(pmask=="sigma2")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
  if(theta_star[["tau2"]]>10 & sum(pmask=="tau2")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 

  return(list(thetaS=theta_star,likOK = likelihoodOK))
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Metropolis-Hastings algorithm

run_mcmc<-function(
  test.yr,
  test_years,
  inf_years,
  strain_years,
  n_part,
  test.list,
  theta,
  runs,
  varpart_prob,
  hist.true=NULL,
  switch1=2,
  seedi=1,
  pmask=NULL,
  linD=F, # toggles linear/exponential cross-reactivity function
  antigenic.map.in=NULL, # define specific map structure (or initial structure if fitting)
  am.spline=NULL, # fit antigenic map along defined spline function
  flu_type = NULL,
  turn_off_likelihood=0 # DEBUG turn off likelihood
  ){
  
  # DEBUG SIMULATION params <<<
   # test.yr=define.year; test.list=test.listSim; theta=theta0;runs=1e2; varpart_prob=vp1;hist.true=historytabSim;  switch1=10; pmask=pmask0;
   
   # seedi=loadseed; antigenic.map.in = NULL; am.spline = am.spl;  flu_type = flu.type; linD=F
  
  # DEBUG set params <<<
  # hist.true=NULL; test.yr=c(2009); runs=200; switch1=10; varpart_prob=0.05 ;   seedi=1; linD=T; antigenic.map.in=NULL; flu_type="H3HN"
  
  time.1 = Sys.time() # DEBUG TIME 1
  
  if(is.null(antigenic.map.in)){antigenic.map.in=inf_years} # if no input map, assume 1D
  test.n=length(test_years); inf.n=length(inf_years); nstrains=length(strain_years)
  sample.index=strain_years-min(inf_years)+1
  test.listPost=test.list

  # Predefine index variables to speed up code
  jj_year=match(test.yr,test_years); testyear_index = test.yr - min(inf_years) + 1
  sample.n=length(jj_year)
  
  # Extract ages and create mask for FluScape data
  if(flu_type=="H3FS" | flu_type=="H3FS_HI"){
    age.list <- array(unlist(test.list),dim=c(5,length(strain_years),n_part))[5,1,]
    age.mask <- sapply(age.list,function(x){if(is.na(x)){1}else{match(max(min(inf_years),test_years[1]-x),inf_years)  }  })
  }else{
    age.mask <- rep(1,n_part) # No mask if ages unknown
  }

  # Make adjustments depending on what is fitted and not
  if(sum(pmask=="muShort")>0){theta[["muShort"]]=1e-10} # Set short term boosting ~ 0 if waning not fitted
  if(sum(pmask=="sigma2")>0){ theta[["sigma2"]]=theta[["sigma"]] } # Fix equal if sigma same for both 
  if(sum(pmask=="error")>0){ theta[["error"]]=1e-10 } # Fix equal if sigma same for both 
  if(sum(pmask=="tau1")>0){ theta[["tau1"]]=1e-10 } # Fix equal if sigma same for both 
  if(sum(pmask=="tau2")>0){ theta[["tau2"]]=1e-10 } # Fix equal if sigma same for both 
  
  # Preallocate memory
  nparam=length(theta); npcov=rep(1,nparam); npcov[match(pmask,names(theta))]=0 # mask specified parameters
  cov_matrix_theta0 = diag(npcov)
  cov_matrix_thetaA=cov_matrix_theta0
  
  thetatab=matrix(NA,nrow=(runs+1),ncol=length(theta)); colnames(thetatab)=names(theta)
  thetatab[1,]=theta
  
  historytab=matrix(NA,nrow=n_part,ncol=inf.n)
  historytabCollect=historytab
  age.tab=matrix(NA,nrow=n_part,ncol=1)
  map.tab=antigenic.map.in
  map.tabCollect=list()
  
  dmatrix0 = outputdmatrix.fromcoord(theta[["sigma"]],inf_years,anti.map.in=map.tab,linearD=linD) # Arrange antigenic map into cross-reaction matrix
  dmatrix20 = outputdmatrix.fromcoord(theta[["sigma2"]],inf_years,anti.map.in=map.tab,linearD=linD) # Arrange antigenic map into cross-reaction matrix
  
  # Pick plausible initial conditions -- using all test years
  if(is.null(hist.true)){
    for(ii in 1:n_part){
      histIC=NULL
      for(kk in 1:length(jj_year)){
        histIC=rbind(histIC,setuphistIC_SA(ii,jj_year[kk],inf.n,test.list,testyear_index,test_years, inf_years))
      }
      histA=as.numeric(colSums(histIC)>0) # combine all histories
      histA0=histA*0
      histA0[c(age.mask[ii]:inf.n)]=histA[c(age.mask[ii]:inf.n)] # make sure no infections outside lifetime
      
      # Constrain to max of X infections in IC
      pick.X = 10
      if( sum(histA0) > pick.X ){ hist0_0 = histA0; pick_subset = sample(which(hist0_0==1),pick.X); histA0[-pick_subset]=0}
      
      historytab[ii,]=histA0
    }
  } else { historytab=hist.true }

  
  colnames(historytab)=as.character(inf_years)

  # Preallocate matrices
  likelihoodtab=matrix(-Inf,nrow=(runs+1),ncol=n_part)
  accepttabT=NULL
  accepttabH=NULL
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Run MCMC
  
  for (m in 1:runs){
    
    #time.3 = Sys.time() # DEBUG Set Time 3
    
    # Adaptive covariance matrix
    if(m==1){
      epsilon0=0.001 # DEBUG
      #cov_matrix_theta=epsilon0*cov_matrix_thetaA
      cov_matrix_basic=epsilon0*cov_matrix_theta0
      varpart_prob0=varpart_prob
    }else{
      epsilon0=max(0.00001,min(1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m)))
      #cov_matrix_theta=epsilon0*cov_matrix_thetaA
      cov_matrix_basic=epsilon0*cov_matrix_theta0
      #varpart_prob0=max(0.02,min(0.25,exp(log(varpart_prob0)+(accept_rateH-0.234)*0.999^m))) # resample max of 25%, min of 2%
    }
    
    # - - - - - - - - - - - - - - - -
    # Resample parameters

    if(m %% switch1==0 | m==1){ # m==1 condition as have to calculate all liks on first step
      if(m==1){
        theta_star = thetatab[m,]
        likelihood_sample_OK = 1 # Check parameters within boundaries
      }else{
        sample_theta = SampleTheta2(thetatab[m,], m,cov_matrix_theta,cov_matrix_basic,nparam=sum(cov_matrix_theta0),pmask) #resample theta
        theta_star = sample_theta$thetaS
        likelihood_sample_OK = sample_theta$likOK # Check parameters within boundaries
      }

      if(sum(pmask=="sigma2")>0){ theta_star[["sigma2"]]=theta_star[["sigma"]] } # Fix equal if sigma same for both 

      map_star=map.tab
      
      #age_star = age.tab
      history_star = historytab
      pickA=c(1:n_part)
      
    }else{
      pickA=NULL
      pickA=sample(n_part, ceiling(varpart_prob0*n_part)) # check that not length zero (i.e. at least one person sampled)
      #age_star = age.tab #SampleAge(pickA,age.tab) #resample age (not for now)
      history_star = SampleHistory(historytab,pickA,inf.n,age_star,inf_years,age.mask) #resample history
      theta_star = thetatab[m,]
      likelihood_sample_OK = 1 # Check parameters within boundaries
    }

    #time.4 = Sys.time() # DEBUG Set Time 3
    #print(paste(m,"/sample:",time.4 - time.3))

    #print(am.spline) # DEBUG
    dmatrix =  1-theta_star[["sigma"]] *dmatrix0;  dmatrix[dmatrix<0]=0 # Arrange antigenic map into cross-reaction matrix
    dmatrix2 = 1-theta_star[["sigma2"]]*dmatrix20; dmatrix2[dmatrix2<0]=0 # Arrange antigenic map into cross-reaction matrix

    # - - - - - - - - - - - - - - - -
    # LIKELIHOOD function - Only calculate for updated history
    
    lik_val=likelihoodtab[m,]
    
    if((turn_off_likelihood!=1 | m==1) & likelihood_sample_OK==1){ # Only calculate if likelihood running
      for(ii in pickA){
        # Set history to zero after test date
        lik.ii=rep(NA,sample.n)
        for(kk in 1:sample.n){
          #For DEBUG: set params <<<  ii=1;kk=2;historyii=as.numeric(history_star[ii,])
          lik.ii[kk]=estimatelik(ii,jj_year[kk],as.numeric(history_star[ii,]),dmatrix,dmatrix2,theta_star,test.list,testyear_index[kk])
          #if(lik.ii[kk]==-Inf){print(c(ii,kk))  } For DEBUG
        }
        lik_val[ii]=sum(lik.ii)
        #if(is.na(lik_val[ii])){ lik_val[ii]=-Inf} For DEBUG
      }
        
    }
    
    #time.6 = Sys.time() # DEBUG Set Time 3
    #print(paste("likelihood calc:",time.6-time.5))

    # - - - - - - - - - - - - - - - -
    # Metropolis Hastings step


    # History sample step
    if( (m %% switch1 != 0) & m>1){
      
      # Calculate piecewise likelihood
      output_prob = ComputeProbability(likelihoodtab[m,pickA],lik_val[pickA],theta_star, thetatab[m,],pmask,turn_off_likelihood) # Only calculate for selected
      pickCP = pickA[runif( length(pickA) ) < output_prob]
      
      historytab[pickCP,] = history_star[pickCP,]
      likelihoodtab[m+1,] = likelihoodtab[m,]
      likelihoodtab[m+1,pickCP] = lik_val[pickCP]
      thetatab[m+1,] = thetatab[m,]
      
    } # end history step
      
    # Theta sample step
    if( (m %% switch1==0) | m==1){
      
      # Estimate probability of update  
      output_prob = ComputeProbability(sum(likelihoodtab[m,]),sum(lik_val),theta_star, thetatab[m,],pmask,turn_off_likelihood) 
      if(is.na(output_prob) & m==1){stop(paste('check initial parameter values',theta_star[["error"]]))}
        
      if(likelihood_sample_OK==0){output_prob = 0}
      
      if(runif(1) < output_prob){
        
        thetatab[m+1,] = theta_star
        accepttabT=c(accepttabT,1)
        likelihoodtab[m+1,] = lik_val
        
      }else{
        thetatab[m+1,] = thetatab[m,]
        likelihoodtab[m+1,] = likelihoodtab[m,]
        if(m %% switch1==0){accepttabT=c(accepttabT,0)}
      }
    } # End theta sample step


    if(m<max(100)){
      accept_rateT=0.234 # target acceptance rate for theta
      #accept_rateH=0.234 # target acceptance rate for infection history
    }else{
      accept_rateT=sum(accepttabT)/length(accepttabT)
      #accept_rateH=sum(accepttabH)/length(accepttabH)
      #cov_matrix_thetaA=cov(thetatab[1:m,]) # Include adaptive covariance matrix for MCMC
    }

    
    # Store infection history
    if(m %% min(runs,20) ==0){
      historytabCollect=rbind(historytabCollect,historytab)
      #map.tabCollect[[round(m/20)]]=map.tab DEPRECATED
    }

    if(m %% min(runs,500) ==0){
      print(c(m,accept_rateT,varpart_prob0,round(sum(likelihoodtab[m,])))) # DEBUG HERE
      save(likelihoodtab,thetatab,inf_years,n_part,test.listPost,historytab,historytabCollect,map.tabCollect,age.tab,test.yr,switch1,file=paste("posterior_sero_runs/outputR_f",paste(test.yr,"_",collapse="",sep=""),"s",seedi,"_lin",linD,".RData",sep=""))
    }
    
    #time.8 = Sys.time() # DEBUG Set Time 3
    #print(paste("loop time:",time.8-time.3))
    
    
  } #End runs loop
  
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Inference using cross-sectional vs longitudinal data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

data.infer <- function(year_test,mcmc.iterations=1e3,loadseed=1,
                       flutype="H3HN",fix.param=NULL , fit.spline=NULL, 
                       switch0=2,linearFn=F,  vp1=0.2, # Probability resample history
                       turn_off_likelihood=0,
                       leave_out_10 = F
                       ) {
  
  #DEBUG  year_test=c(2007:2012); loadseed=1; switch0=2; vp1=0.2; mcmc.iterations=1e2; strain.fix=T; flutype="H3HN"; fix.param=NULL; linearFn=T; pmask=c("tau2"); turn_off_likelihood=0; leave_out_10 = T;
  
  # INFERENCE MODEL
  # Run MCMC for specific data set
  if(flutype=="H3HN" & leave_out_10==F){load("R_datasets/HaNam_data.RData")}
  if(flutype=="H3HN" & leave_out_10==T){load("R_datasets/HaNam_data_V_minus_10.RData")}
  if(flutype=="H3FS"){load("R_datasets/FluScapeH3_data.RData")}
  if(flutype=="H3FS_HI"){load("R_datasets/FluScapeH3_HI_data.RData")}
    #am.spl<-load.flu.map.data() # define spline from antigenic map data
  if(flutype=="B"){load("R_datasets/Fluscape_data.RData")}
  if(flutype=="H1"){load("R_datasets/HK_data.RData")}
  
  # Plot simulation data vs history
  if(flutype=="H3HN"){history_in = read.csv("R_datasets/hist_IC_H3.csv"); history_in = history_in[,-1]; colnames(history_in)=NULL; history_in= as.matrix(history_in) }else{history_in = NULL}
  
  # Set initial theta
  theta0=c(mu=NA,tau1=NA,tau2=NA,wane=NA,sigma=NA,muShort=NA,error=NA,disp_k=NA,sigma2=NA)
  theta0[["mu"]]=2 + if(sum(fix.param=="vary.init")>0){0.2*runif(1,c(-1,1))}else{0} # basic boosting
  theta0[["tau1"]]=0.05 # - NOT CURRENTLY USED
  theta0[["tau2"]]=0.05 + if(sum(fix.param=="vary.init")>0){0.01*runif(1,c(-1,1))}else{0} # suppression via AGS
  theta0[["wane"]]= 0.8 + if(sum(fix.param=="vary.init")>0){0.05*runif(1,c(-1,1))}else{0} # short term waning - half life of /X years -- add noise to IC if fitting
  theta0[["sigma"]]=0.1 + if(sum(fix.param=="vary.init")>0){0.02*runif(1,c(-1,1))}else{0} # cross-reaction
  theta0[["sigma2"]]=0.05 + if(sum(fix.param=="vary.init")>0){0.005*runif(1,c(-1,1))}else{0} # short-term cross-reaction
  theta0[["muShort"]]=2.5 + if(sum(fix.param=="vary.init")>0){0.2*runif(1,c(-1,1))}else{0} # short term boosting
  theta0[["error"]]=1.5 + if(sum(fix.param=="vary.init")>0){0.1*runif(1,c(-1,1))}else{0} # measurement error
  theta0[["disp_k"]]=1 # dispersion parameter - NOT CURRENTLY USED
  theta=theta0
  print(theta0)
  
  define.year=year_test # years to include in inference
  
  # browser()
  
  # Define antigenic map for fitting
  xx=scalemap(inf_years,inf_years)
  yy=predict(am.spl,xx)$y
  antigenic.map.in0=cbind(xx,yy)

  # RUN MCMC
  # Note: NEED TO RE-INITIALISE DATAFRAME IF REPEAT RUN (i.e. reload dataset above)
  run_mcmc(
    test.yr=define.year,
    test_years,
    inf_years,
    strain_years,
    n_part,
    test.list,
    theta=theta0,
    runs=mcmc.iterations, # number of MCMC runs
    varpart_prob=vp1,
    hist.true=history_in,
    switch1=switch0, # ratio of infection history resamples to theta resamples. This is fixed
    pmask=c(fix.param,"disp_k"), #c("wane"), #,"muShort"), # specify parameters to fix
    seedi=paste(loadseed,"_",flutype,sep=""), # record output
    antigenic.map.in = antigenic.map.in0, # define specific map structure (or initial structure if fitting)
    am.spline = NULL, # decide whether to fit antigenic map along "am.spl" spline function
    linD=linearFn,
    flu_type= flutype,
    turn_off_likelihood)
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define simulation model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


simulation.infer <- function(seed_i,mcmc.iterations=1e3, strain.fix=T,flu.type="H3HN", vp1=0.2,fit.spline = NULL, fix.param="vary.init",linearFn=F) {

  #DEBUG seed_i=1; mcmc.iterations=40; strain.fix=T; flu.type="H3HN"; fix.param ="vary.init"; linearFn= F
  
  # Load data
  loadseed=paste("SIM_",seed_i,sep="")
  if(flu.type=="H3HN"){load("R_datasets/HaNam_data.RData"); npartM=69; define.year=c(2007:2012); pmask0=c("tau1")}
  if(flu.type=="H3FS"){load("R_datasets/FluScapeH3_data.RData"); npartM=151; define.year=c(2009); pmask0=c("muShort","tau1")}
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # SIMULATION MODEL
  # Generate simulated data 
  #tau1=back-boost  / tau2=suppress / disp_k=dispersion (deprecated) 
  #sigma1=long-term cross-reactivity / sigma 2=short-term CR
  
  thetaSim = c(mu=2,tau1=0.02,tau2=0.05,wane=0.75,sigma=0.15,muShort=2,error=1.5,disp_k=1,sigma2=0.05)  # NOTE EDITED FOR SIMULATION RUNS
  
  if(strain.fix==T){
    strain_years0 = strain_years
  }else{
    strain_years0 = inf_years
  }
  inf_years.in = inf_years
  
  # Generate 2D map
  xx=scalemap(inf_years.in,inf_years.in)
  yy=predict(am.spl,xx)$y
  antigenic.map0=cbind(xx,yy)
  
  #attack.yr=read.csv("datasets/sim_attack.csv")[1:length(inf_years),1]
  # Simulate multiple
  sd.param=0.5
  attack.yr=rlnorm(inf_years.in,meanlog=log(0.15)-sd.param^2/2,sdlog=sd.param)
  attack.yr[1]=rlnorm(1,meanlog=log(0.5)-(sd.param/2)^2/2,sdlog=(sd.param/2))
  
  write.csv(attack.yr,paste("datasets/sim_attackS",seed_i,".csv",sep=""))
  
  simulate_data(test_years=define.year, # this needs to be vector
                inf_years=inf_years.in,strain_years=strain_years0,n_part=npartM, #leave strain years blank to use HaNam strains
                roundv=T, # Generate integer titre data
                thetastar=thetaSim,
                antigenic.map.in = antigenic.map0,
                #pmask=c("wane","sigma2"), # Specify what is included
                linD = linearFn,
                p.inf=attack.yr,seedi=loadseed)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # INFERENCE MODEL
  # Run MCMC for simulated data set
  
  load(paste("R_datasets/Simulated_data_",loadseed,".RData",sep="")) # Load simulation data for inference step that follows
  
  # Set initial theta
  theta0=c(mu=NA,tau1=NA,tau2=NA,wane=NA,sigma=NA,muShort=NA,error=NA,disp_k=1,sigma2=0.1)
  theta0[["mu"]]=2+ if(sum(fix.param=="vary.init")>0){runif(1,c(-1,1))}else{0} # basic boosting
  theta0[["tau1"]]=0.1+ if(sum(fix.param=="vary.init")>0){0.03*runif(1,c(-1,1))}else{0} # back-boost
  theta0[["tau2"]]=0.05+ if(sum(fix.param=="vary.init")>0){0.03*runif(1,c(-1,1))}else{0} # suppression via AGS
  theta0[["wane"]]= 0.8 + if(sum(fix.param=="vary.init")>0){0.1*runif(1,c(-1,1))}else{0} # -log(0.5)/1 # short term waning - half life of /X years
  theta0[["sigma"]]=0.15 + if(sum(fix.param=="vary.init")>0){0.1*runif(1,c(-1,1))}else{0} # long-term cross-reaction
  theta0[["sigma2"]]=0.05+ if(sum(fix.param=="vary.init")>0){0.05*runif(1,c(-1,1))}else{0} # short-term cross-reaction
  theta0[["muShort"]]=2 + if(sum(fix.param=="vary.init")>0){runif(1,c(-1,1))}else{0} # short term boosting
  theta0[["error"]]= 1.5 + if(sum(fix.param=="vary.init")>0){0.1*runif(1,c(-1,1))}else{0} # measurement error
  theta0[["disp_k"]]=0.1 # overdispersion (deprecated)
  theta=theta0

  print(theta0)

  # RUN MCMC
  # Note: NEED TO RE-INITIALISE DATAFRAME IF REPEAT RUN (i.e. reload dataset above)
  run_mcmc(
    test.yr=define.year,
    test_years,
    inf_years,
    strain_years,
    n_part,
    test.list=test.listSim, # use simulated data as input
    theta=theta0,
    runs=mcmc.iterations, # number of MCMC runs
    varpart_prob=vp1,
    hist.true= NULL, # True starting point?  # *** DEBUG ***  historytabSim
    switch1=2, # ratio of infection history resamples to theta resamples. This is fixed
    pmask=pmask0, # ,"map.fit" specify parameters to fix
    seedi=loadseed,
    antigenic.map.in = antigenic.map0, # Define random initial map to fit
    am.spline = NULL, # decide whether to fit antigenic map along "am.spl" spline function
    flu_type = flu.type,
    linD=linearFn)
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run code on network [SR code]
fn.network<-function(){
  system.time(
    for(ii in 1:1){  # Use multiple seeds for simulation code
      fnSeedLoop(ii)
    } # End loop over seeds
  )
  
  # Do some of these over the network
  library("didewin")
  didewin::didewin_config_global(credentials="~/.smbcredentials",
                                 home="~/dide/home",
                                 temp="~/dide/tmp")
  
  make_trees <- function(n, nspp) {
    lapply(seq_len(n), function(...) ape::rtree(nspp))
  }
}

# Just the likelihood function for DEBUGGING

likelihood_function <- function(theta_star,inf_years,test.yr,map_star,am.spline,test.list,history_star,n_part){
  
  testyear_index = test.yr - min(inf_years) + 1
  test_years = test.yr
  jj_year=match(test.yr,test_years); testyear_index = test.yr - min(inf_years) + 1
  sample.n=length(jj_year)
  
  dmatrix = outputdmatrix.fromcoord(theta_star[["sigma"]] ,inf_years,anti.map.in=map_star) # Arrange antigenic map into cross-reaction matrix
  dmatrix2 = outputdmatrix.fromcoord(theta_star[["sigma2"]],inf_years,anti.map.in=map_star) # Arrange antigenic map into cross-reaction matrix
  
  # - - - - - - - - - - - - - - - -
  # LIKELIHOOD function - Only calculate for updated history
  
  lik_val=NULL
  for(ii in 1:n_part){
    # Set history to zero after test date
    lik.ii=rep(NA,sample.n)
    for(kk in 1:sample.n){
      #For DEBUG: set params <<<  ii=1;kk=2;historyii=as.numeric(history_star[ii,])
      lik.ii[kk]=estimatelik(ii,jj_year[kk],as.numeric(history_star[ii,]),dmatrix,dmatrix2,theta_star,test.list,testyear_index[kk])
      #if(lik.ii[kk]==-Inf){print(c(ii,kk))  } For DEBUG
    }
    lik_val=c(lik_val,sum(lik.ii))
    #if(is.na(lik_val[ii])){ lik_val[ii]=-Inf} For DEBUG
  }
  
  sum(lik_val)
  
}


