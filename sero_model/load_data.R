# Model of serological dynamics - uses PLOS Biology model (Kucharski et al. 2015)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data (Fonville et al.)

options(StringsAsFactors=F)
data1=read.csv("datasets/HaNamCohort.csv", as.is=T)

# List test strains
nstrains=length(data1)-2 # remove subject and sample year
strain_names=names(data1)[3:(nstrains+2)]
test.index=c(1:nstrains)

# Convert to log titres and set missing data = NA
data1[data1=="*"]=NA
data1[,strain_names]=apply(data1[,strain_names],2,function(x){log2(as.numeric(x)/10)+1}) 

# Convert names into strain years
strain_years=as.numeric(sapply(strain_names,function(x){
a1=max(which(strsplit(x, "")[[1]]=="."))
lstr=nchar(x)
yr1=substr(x, a1+1, lstr)

if(nchar(yr1)>4){yr1=substr(yr1, 1, 4)}
year=yr1
if(nchar(yr1)==2 & as.numeric(yr1)>15){year=paste("19",yr1,sep="")}
if(nchar(yr1)==2 & as.numeric(yr1)<15){year=paste("20",yr1,sep="")}
year
}
))

strain_years_unique=sort(unique(strain_years)) # years of samples tested against

# Gather participants and infection years

n_part=max(data1$Subject.number) # number of participants

test_years=unique(data1$Sample.year) # year of testing
test.n=length(test_years) # number of test years

inf_years=seq(min(strain_years),max(c(test_years,strain_years))) #annual infection model
inf.n=length(inf_years) # number of possible infecting strains

# Set up list of test data for quick access

data.Test=data1[,strain_names]
test.list=list()

for(ii in 1:n_part){

subjectn=ii
i.list=list()

#   [,1] [,2] [,3] [,4] [,5] [,6]
#   test.year    2010 2010 2010 2010 2010 2010
#   titredat        1    2    4    3    6    3
#   strain_years 1990 1994 1998 2002 2006 2010 # Year of isolate
#   sample.index    1    5    9   13   17   21 # Numerical index of strain isolate (start with first possible year of infection)

for(jj in 1:test.n){

testyr=test_years[jj]
dataI=data.Test[data1$Subject.number==subjectn & data1$Sample.year==testyr,]
i.list[[jj]]=rbind(rep(testyr,nstrains),
      dataI[,!is.na(dataI)],
      strain_years[!is.na(dataI)],
      strain_years[test.index[!is.na(dataI)]]-min(strain_years)+1
      )

}


test.list[[ii]]=i.list

}

save(test_years,inf_years,strain_years,n_part,test.list,file=paste("R_datasets/HaNam_data_V.RData",sep=""))

