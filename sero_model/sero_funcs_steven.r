# Stevens functions in R

#' Load up the fluscape data and transform it into the same shape
#' as the Vietnam data.
load_fluscape <- function(psthfssvn="~/Dropbox/svn/fluscape") {
  
  # source the main fluscape functions needed for loading the data
  source(paste(pathfssvn,"/trunk/source/R/GeneralUtility.r",sep=""))
  
}