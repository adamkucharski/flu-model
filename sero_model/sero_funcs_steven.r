# Stevens functions in R

#' Load up the fluscape data and transform it into the same shape
#' as the Vietnam data.
load_fluscape <- function(pathfssvn="~/Dropbox/svn/fluscape/trunk/") {
  
  # Source the main fluscape functions needed for loading the data
  source(paste(pathfssvn,"/source/R/GeneralUtility.r",sep=""))
  
  # Run the load and merge function for the first three visits
  fsd_tmp <- load.and.merge.part.V1.V2.V3(topdir=pathfssvn)
  
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
  rtn_yam_td[,"strain_years",] <- c(2011,2011,2011,2012,2012,2012)
  
  
}