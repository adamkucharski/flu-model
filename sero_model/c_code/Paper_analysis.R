# Set up parallel

library(foreach)
library(doMC)
#registerDoMC(4)  #change the 2 to your number of CPU cores
#getDoParWorkers()




setwd("~/Dropbox/Imperial/Fluscape_2/sero_model/c_code")



system("R CMD SHLIB c_model2.c")
dyn.load("./c_model2.so")


# Pick test year


dataplot1=data1[data1$Subject.number==subjectn,]


# Set up test strains

func1 <- function(x) {
  if (!is.numeric(x))
    stop("argument x must be numeric")
  out <- .C("c_model2",
            n=as.integer(length(x)),
            x=as.double(x))
  return(out$x)
}

func1(c(1:5))



#list<-foreach(iiG=1:10) %dopar% {
