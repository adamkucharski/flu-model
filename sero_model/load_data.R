# Model of serological dynamics - uses PLOS Biology model (Kucharski et al. 2015)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data (Fonville et al.)

options(StringsAsFactors=F)
data1=read.csv("datasets/HaNamCohort.csv", as.is=T)

# List test strains
nstrains=length(data1)-2
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

strain_years_unique=sort(unique(strain_years))

# Gather participants and infection years

n_part=max(data1$Subject.number)
inf_years=seq(min(strain_years),max(strain_years)) #annual infection model
inf.n=length(inf_years)

test.years=unique(data1$Sample.year)
test.n=length(test.years)

# Set up list of test data for quick access

data.Test=data1[,strain_names]
test.list=list()

for(ii in 1:n_part){

subjectn=ii
i.list=list()

for(jj in 1:test.n){

testyr=test.years[jj]
dataI=data.Test[data1$Subject.number==subjectn & data1$Sample.year==testyr,]
i.list[[jj]]=rbind(rep(testyr,nstrains),
      dataI[,!is.na(dataI)],
      strain_years[!is.na(dataI)],
      strain_years[test.index[!is.na(dataI)]]-min(strain_years)+1
      )

}

test.list[[ii]]=i.list

}