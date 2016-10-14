# Stevens functions in R

# Declare some libraries that are needed here
library(Rcpp)

#' Load up the fluscape data and transform it into the same shape
#' as the Vietnam data.
#' make_fluscape_rdata()
make_fluscape_rdata <- function(
  pathfssvn="~/Dropbox/svn/fluscape/trunk/",
  sero="YAM") {
  
  # Source the main fluscape functions needed for loading the data
  source(paste(pathfssvn,"/source/R/GeneralUtility.r",sep=""))
  
  # Run the load and merge function for the first three visits
  fsd_tmp <- load.and.merge.part.V1.V2.V3(
    topdir=pathfssvn,
    convert.titers=TRUE)
  
  # Define a vectors of the assay results in the correct order
  if (sero=="YAM") {
    asTest <- c("HI.B.1988.V1","HI.B.2002.V1","HI.B.2006.V1",
                "HI.B.1988.V2","HI.B.2002.V2","HI.B.2006.V2")
  } else {
    stop("only yam implemented in load_fluscape")
  }
  
  # Subset the data for those for whom we have both flu B assays
  # and date of birth. Initially, just do the 
  subset_mask <- (  !is.na(fsd_tmp$HI.B.1987.V1) &
                    !is.na(fsd_tmp$HI.B.1987.V2) &
                    !is.na(fsd_tmp$HI.B.2004.V1) &
                    !is.na(fsd_tmp$HI.B.2004.V2) &
                    !is.na(fsd_tmp$HI.B.2008.V1) &
                    !is.na(fsd_tmp$HI.B.2008.V2) &
                    !is.na(fsd_tmp$HI.B.1988.V1) &
                    !is.na(fsd_tmp$HI.B.1988.V2) &
                    !is.na(fsd_tmp$HI.B.2002.V1) &
                    !is.na(fsd_tmp$HI.B.2002.V2) &
                    !is.na(fsd_tmp$HI.B.2006.V1) &
                    !is.na(fsd_tmp$HI.B.2006.V2) &
                    !is.na(fsd_tmp$PART_AGE.V1) &
                    !is.na(fsd_tmp$PART_AGE.V2)  
                    ) 
  fsd <- fsd_tmp[subset_mask,]
  
  test_years <- c(2011,2012)
  inf_years <- c(1987:2012)
  if (sero=="YAM") {
    strain_years=c(1987,1988,2002,2004,2006,2008)
  } else {
    stop("only yam implemented in load_fluscape")
  }

  # SR array function - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # setup return obects problem here with noP
  noP <- dim(fsd)[1]
  rtn_yam_td <- array(
    dim = c(noP,4,6),
    dimnames = list(1:noP,
      c("test.year","titredat","strain_years","sample.index"))
  )
  
  rtn_yam_td[,"test.year",] <- c(2011,2011,2011,2012,2012,2012)
  rtn_yam_td[,"strain_years",] <- c(1988,2002,2006,1988,2002,2006)
  rtn_yam_td[,"sample.index",] <- c(1,2,3,1,2,3) # AK: updated to make years from first strain
  
  # Cycle through people putting in the right tire values
  for (i in 1:noP) {
    rtn_yam_td[i,"titredat",] <- as.numeric(fsd[i,asTest])
  }
  
  npart <- noP
  test.list <- rtn_yam_td
  age.yr <- fsd$PART_AGE.V2
  
  # Return the list of required output
  save(test_years, strain_years, npart, test.list,
       file="R_datasets/Fluscape_data.RData")
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # AK: set up code to return as list by participant
  # setup return obects problem here with noP
  
  asTest1 <- c("HI.B.1987.V1","HI.B.1988.V1","HI.B.2002.V1","HI.B.2004.V1","HI.B.2006.V1","HI.B.2008.V1")
  asTest2 <- c("HI.B.1987.V2","HI.B.1988.V2","HI.B.2002.V2","HI.B.2004.V2","HI.B.2006.V2","HI.B.2008.V2")
  
  n_part <- noP
  test.n <- length(test_years)
  age.yr <- fsd$PART_AGE.V2
  nstrains <- 6
  
  fsd <- fsd_tmp[subset_mask,]
  noP <- dim(fsd)[1]
  rtn_yam_td <- array(
    dim = c(noP,4,6),
    dimnames = list(1:noP,
                    c("test.year","titredat","strain_years","sample.index"))
  )
  
  sample.index <- strain_years-min(inf_years)+1 # AK: updated to make years from first strain
  
  # Cycle through participants putting in the right tire values
  
  test.list <- list()
  for(ii in 1:n_part){
    
    subjectn=ii
    i.list=list()
    
    #   FORMAT
    #   test.year    2010 2010 2010 2010 2010 2010
    #   titredat        1    2    4    3    6    3
    #   strain_years 1990 1994 1998 2002 2006 2010 # Year of isolate
    #   sample.index    1    5    9   13   17   21 # Numerical index of strain isolate (start with first possible year of infection)
    #   age
    
    for(jj in 1:test.n){
      testyr=test_years[jj]
      testpick=if(jj==1){asTest1}else{asTest2}
      dataI=as.numeric(fsd[ii,testpick])
      i.list[[jj]]=rbind(test.year=rep(testyr,nstrains),
                         titredat=dataI,
                         strain_years,
                         sample.index,
                         age=rep(as.numeric(fsd[ii,"PART_AGE.V1"]))+jj-1 # record age
      )
    }
    test.list[[ii]]=i.list
  }
  
  # Return the list of required output
  save(test_years, inf_years,strain_years, n_part, test.list,
         file="R_datasets/Fluscape_data_List.RData")
  
}


# Load FluScape data and format to same shape

make_FluScape_rdata <- function(){

  # Load up the raw data
  data0 <- read.csv("datasets/Fluscape_SupplmentalDataS1.csv",stringsAsFactors=FALSE)

  
  # Define all the quads
  data1 <- data.frame(data0)
  aa <- unique(data1$neut.against)
  strain0 <- sapply(unique(data1$neut.against),function(x){x})
  
  n_part <- 151
  
  data1$titers=round(sapply(data1$titers,function(x){log2(exp(as.numeric(x))/10)+1}),6)  # Make titre log2 -- NOTE 0 to 8 scale
  
  
  
  # Use last 3 test years only to keep gaps uniform...
  test_years <- c(2009)
  test.n <- length(test_years)
  inf_years <- c(1968:2009)
  strain_years <- c(1968,1975,1979,1989,1995,2002,2003,2005,2008)
  #inf_years <- strain_years
  nstrains <- 1 #length(strain_years)
  
  test.list <- list()
  for(ii in 1:n_part){
    
    subjectn=ii; i.list=list()
    dataP <- data1[data1$id==ii,]
    
    #   FORMAT
    #   test.year    2010 2010 2010 2010 2010 2010
    #   titredat        1    2    4    3    6    3
    #   strain_years 1990 1994 1998 2002 2006 2010 # Year of isolate
    #   sample.index    1    5    9   13   17   21 # Numerical index of strain isolate (start with first possible year of infection)
    #   age
    
    for(jj in 1:test.n){
      
      testyr <- test_years[jj]
      dataI <- dataP$titers
      strain_year <- strain_years
      i.list[[jj]]=rbind(test.year=rep(testyr,nstrains),
                         titredat=dataI,
                         strain_year,
                         strain_years-min(strain_years)+1,
                         age=dataP[1,"age"] # record age
      )
    }
    test.list[[ii]]=i.list
    
  }
  
  # Return the list of required output
  save(test_years, inf_years,strain_years, n_part, test.list,
       file="R_datasets/FluScapeH3_data.RData")
  
}


# Load HK data and format to same shape

make_HK_rdata <- function(){

  #source("~/Documents/hkiss/R/hks.R")
  #dataHK <- hks.read.clean.study("~/Documents/hkiss/data/hkiss_main_v1.csv")
  
  # Load up the raw data
  rtn <- read.csv("~/Documents/hkiss/data/hkiss_main_v1.csv",stringsAsFactors=FALSE)
  
  # Define all the quads
  rtn$is.quad <- rtn$blood_1=="Yes" & rtn$blood_2=="Yes" & rtn$blood_3=="Yes" & (rtn$blood_4 %in% c("Yes","yes","YES"))
  rtn <- rtn[!is.na(rtn$is.quad),] # Remove NA
  rtn <- rtn[rtn$is.quad==T,] # 
  
  # Correct the dates of the bloods
  rtn$blood_1_date <- as.Date(as.character(rtn$blood_1_date),"%d/%m/%Y")
  rtn$blood_2_date <- as.Date(as.character(rtn$blood_2_date),"%d/%m/%Y")
  rtn$blood_3_date <- as.Date(as.character(rtn$blood_3_date),"%d/%m/%Y")
  rtn$blood_4_date <- as.Date(as.character(rtn$blood_4_date),"%d/%m/%Y")
  
  hist(c(rtn$blood_1_date,rtn$blood_2_date,rtn$blood_3_date,rtn$blood_4_date),breaks=as.Date(14000+(0:40)*35,origin="1970-01-01"))
  
  rtn$pT1 <- log((rtn$H1N1.T1 / 5),2)
  rtn$pT2 <- log((rtn$H1N1.T2 / 5),2)
  rtn$pT3 <- log((rtn$H1N1.T3 / 5),2)
  rtn$pT4 <- log((rtn$H1N1.T4 / 5),2)
  rtn$sT1 <- log((rtn$H3N2.T1 / 5),2)
  rtn$sT2 <- log((rtn$H3N2.T2 / 5),2)
  rtn$sT3 <- log((rtn$H3N2.T3 / 5),2)
  rtn$sT4 <- log((rtn$H3N2.T4 / 5),2)
  
  rtn <- rtn[!is.na(rtn$pT2),]
  
  HK.data <- rtn[,c("Age_1","pT2","pT3","pT4")]
  HK.data$Age_1 <- round(HK.data$Age_1)
  
  #rtn[,c("pT1","pT2","pT3","pT4")] %>% apply(2,function(x){x>2}) %>% colSums - Check titres DEBUG

  # Use last 3 test years only to keep gaps uniform...
  
  n_part <- length(rtn$sr.index)
  
  test_years <- c(2009:2011)
  test.n <- length(test_years)
  inf_years <- c(2009:2011)
  strain_years <- c(2009:2011)
  nstrains <- 1 #length(strain_years)
  
  test.list <- list()
  for(ii in 1:n_part){
    
    subjectn=ii; i.list=list()
    
    #   FORMAT
    #   test.year    2010 2010 2010 2010 2010 2010
    #   titredat        1    2    4    3    6    3
    #   strain_years 1990 1994 1998 2002 2006 2010 # Year of isolate
    #   sample.index    1    5    9   13   17   21 # Numerical index of strain isolate (start with first possible year of infection)
    #   age
    
    for(jj in 1:test.n){
      testyr <- test_years[jj]
      dataI <- as.numeric(HK.data[ii,jj+1]) # EDIT THIS TO MAKE REPEATS?
      strain_year <- strain_years[jj]
      i.list[[jj]]=rbind(test.year=rep(testyr,nstrains),
                         titredat=dataI,
                         strain_year,
                         jj,
                         age=rep(as.numeric(HK.data[ii,"Age_1"]))+jj-1 # record age
      )
    }
    test.list[[ii]]=i.list
    
  }
  
  # Return the list of required output
  save(test_years, inf_years,strain_years, n_part, test.list,
       file="R_datasets/HK_data.RData")

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# simLocations()
simLocations <- function(
  yrs=rep(2000:2018),
  mu_cs=4,
  c_cs=2,
  dr=c(2,6,2,6)) {
    noyears <- length(yrs)
    rtn <- matrix(nrow=noyears,ncol=2)
    next_cluster <- c_cs
    cur_x <- 0
    cur_y <- 0
    for (i in noyears:1) {
        rtn[i,] <- c(cur_x,cur_y)
        next_cluster <- next_cluster - 1
        if (next_cluster<1) {
          next_cluster <- 1 + rpois(1,mu_cs-1)
          cur_x <- cur_x + runif(1,min = dr[1],max=dr[2])
          cur_y <- cur_y + runif(1,min = dr[3],max=dr[4])
        }
    }
    rtn
}