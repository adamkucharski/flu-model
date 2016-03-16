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
  
#   [,1] [,2] [,3] [,4] [,5] [,6]
#   test.year    2010 2010 2010 2010 2010 2010
#   titredat        1    2    4    3    6    3
#   strain_years 1990 1994 1998 2002 2006 2010
#   sample.index    1    5    9   13   17   21
  
  # setup return obects problem here with noP
  noP <- dim(fsd)[1]
  rtn_yam_td <- array(
    dim = c(noP,4,6),
    dimnames = list(1:noP,
      c("test.year","titredat","strain_years","sample.index"))
  )
  
  rtn_yam_td[,"test.year",] <- c(2011,2011,2011,2012,2012,2012)
  rtn_yam_td[,"strain_years",] <- c(1988,2002,2006,1988,2002,2006)
  rtn_yam_td[,"sample.index",] <- c(1,2,3,1,2,3)
  
  # Cycle through people putting in the right tire values
  for (i in 1:noP) {
    rtn_yam_td[i,"titredat",] <- as.numeric(fsd[i,asTest])
  }
  
  test_years <- c(2011,2012)
  inf_years <- 1987:2012
  if (sero=="YAM") {
    strain_years=c(1988,2002,2006)
  } else {
    stop("only yam implemented in load_fluscape")
  }
  npart <- noP
  test.list <- rtn_yam_td
  age.yr <- fsd$PART_AGE.V2
  
  # Return the list of required output
  save(test_years, strain_years, npart, test.list,
       file="R_datasets/Fluscape_data.RData")
  
}

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