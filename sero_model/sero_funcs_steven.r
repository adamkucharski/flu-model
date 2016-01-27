# Stevens functions in R

#' Load up the fluscape data and transform it into the same shape
#' as the Vietnam data.
load_fluscape <- function(pathfssvn="~/Dropbox/svn/fluscape/trunk/") {
  
  # source the main fluscape functions needed for loading the data
  source(paste(pathfssvn,"/source/R/GeneralUtility.r",sep=""))
  
  # Run the general utility amd load up the new data
  fsd <- load.and.merge.part.V1.V2.V3(topdir=pathfssvn)
  
  # 
  
}