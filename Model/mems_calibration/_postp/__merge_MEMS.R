rm(list=ls(all=TRUE))

library(foreign)
library(stringr)
library(data.table)

##################################################################################
#merge simulations arable
##################################################################################
MAOM<-fread("/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_out/MAOM_1.csv")
POM<-fread("/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_out/POM_1.csv")

for (j in 2:1000) {

MAOM_<-fread(paste('/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_out/MAOM_',j,'.csv', sep=""))
 MAOM<-rbind(MAOM, MAOM_)

 POM_<-fread(paste('/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_out/POM_',j,'.csv', sep=""))
 POM<-rbind(POM, POM_) 
 
 
  
################################################################################

}

 
  write.dbf(MAOM, file="/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_merge/MAOM_MEMS_LCS.dbf")
  write.dbf(POM, file="/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_merge/POM_MEMS_LCS.dbf")
  
