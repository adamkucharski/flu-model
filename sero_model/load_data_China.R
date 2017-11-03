# Load FluScape data and format to same shape

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
  
