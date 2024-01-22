
library(foreign)
#library(hydroGOF)
library(simpleboot)
library(raster)
library(maptools)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(randomForest)

##### LUCAS point spatial join with over3
lucasSHP<-readShapePoints("D:\\SIT\\caprese\\soil\\LUCAS\\LUCAS-2009-data-as-distributed-by-ESDAC\\LAEA\\LCS_09_ALL_qc_nut2.shp")  ###LUCAS soil
#lucasSHP<-lucasSHP[which(!duplicated(lucasSHP$POINT_ID)), ]

 lucasDB<-as(lucasSHP, "data.frame")[,1:19]

lucasDB[lucasDB==-999]<-NA		####limit of detection
lucasDB$OC[lucasDB$OC<2]<-1
lucasDB$OC<-ifelse((substring(lucasDB$NUTS_2,1,2) == "RO" | substring(lucasDB$NUTS_2,1,2) == "BG") & lucasDB$OC == 6, 1, lucasDB$OC) ###min SOC =1

lucasDB$CaCO3[lucasDB$CaCO3<1]<-0.5
lucasDB$N[lucasDB$N<0.2]<-0.1
lucasDB$P[lucasDB$P<10]<-5


#LC1<-(ifelse((substring(lucasDB$NUTS_2,1,2) %in% c("RO","BG")), as.character(lucasDB$LC12), as.character(lucasDB$LC09)))
#lucasDB<-cbind(lucasDB, "LC1"=LC1)


##correction SOC for CACO3## NUTS definition 2013!!!

lucasDB[lucasDB$NUTS_2 == "ITD3" ,9]<- lucasDB[lucasDB$NUTS_2 == "ITD3" ,9] * (1+ (0.000006*lucasDB[lucasDB$NUTS_2 == "ITD3" ,10]^2 -0.0043*lucasDB[lucasDB$NUTS_2 == "ITD3" ,10]))
lucasDB[lucasDB$NUTS_2 == "ITD4" ,9]<- lucasDB[lucasDB$NUTS_2 == "ITD4" ,9] * (1+ (0.000006*lucasDB[lucasDB$NUTS_2 == "ITD4" ,10]^2 -0.0043*lucasDB[lucasDB$NUTS_2 == "ITD4" ,10]))


##correction BD##

 BDm<- 0.80806 + (0.823844*exp(-0.27993*lucasDB$OC/10))+(0.0014065*lucasDB$sand)-(0.0010299*lucasDB$clay) 
  BDm_c<-(100- lucasDB$coarse)/((lucasDB$coarse/2.65)+(100-lucasDB$coarse)/BDm)
  V_GR<- 1-BDm_c/BDm

 BDo<- 0.074 + 2.632* exp(-0.076*lucasDB$OC/10) 
# BDo<- 1.4903 - 0.33293 * log(lucasDB$OC/10)
  BDo_c<-(100- lucasDB$coarse)/((lucasDB$coarse/2.65)+(100-lucasDB$coarse)/BDo)
  Vo_GR<- 1-BDo_c/BDo

 BD<-ifelse(lucasDB$OC/10<20, BDm, BDo)
   VGR<-ifelse(lucasDB$OC/10<20, V_GR, Vo_GR)


lucasDB<-cbind(lucasDB,"NUT0"= substring(lucasDB$NUTS_2,1,2))
lucasDB<-cbind(lucasDB, "BD"=BD, "VGR"= VGR, "SOC_lcs"=0)


lucasDB[lucasDB$NUT0 == "RO" ,11]<-lucasDB[lucasDB$NUT0 == "RO" ,11]/100
lucasDB[lucasDB$NUT0 == "BG" ,11]<-lucasDB[lucasDB$NUT0 == "BG" ,11]/100


#####################################################MONTE CARLO for SOC
lucasDB$VGR[lucasDB$VGR<0]<-0


i=nrow(lucasDB)

MC<- array(data=NA, dim=c(i,5000))		


if(cm_prof==20) {

 lucasDB$SOC_lcs<- cm_prof*lucasDB$OC/10*BD*(1-VGR)			## cm depth!!
 
for (i in 1 : i) {
  MC[i,]<-  rnorm(5000, BD[i], BD[i]*0.0475) * rnorm(5000, cm_prof, cm_prof*0.0425) * rnorm(5000, lucasDB$OC[i]/10, lucasDB$OC[i]/10*0.115) * (1 - rnorm(5000, lucasDB$VGR[i], lucasDB$VGR[i]*0.4)) 
 }
 
}else{

#####SOC stratification grassland 20-30cm 1.15--1.27--1.34   @@@@@(1-exp(-0.04*30))/(1-exp(-0.04*20) sd +/- 50%		

lucasDB$SOC_lcs<-ifelse(lucasDB$LC09 %in% c("E10", "E20", "E30"), (cm_prof-10)*lucasDB$OC/10*BD*(1-VGR)*ifelse(lucasDB$OC<200,1.27,1.5), cm_prof*lucasDB$OC/10*BD*(1-VGR))

 for (i in 1 : i) {
  if(lucasDB$LC09[i] %in% c("E10", "E20", "E30")) {
   MC[i,]<-  rnorm(5000, BD[i], BD[i]*0.0475) * rnorm(5000, (cm_prof-10) , (cm_prof-10)*0.0425) * rnorm(5000, lucasDB$OC[i]/10, lucasDB$OC[i]/10*0.115) * (1 - rnorm(5000, lucasDB$VGR[i], lucasDB$VGR[i]*0.4)) * 
          ifelse(lucasDB$OC[i]<200,runif(5000, 1.15, 1.34), 1.5)#rnorm(5000, 1.23, 0.25)
  }else{
   MC[i,]<-  rnorm(5000, BD[i], BD[i]*0.0475) * rnorm(5000, cm_prof, cm_prof*0.0425) * rnorm(5000, lucasDB$OC[i]/10, lucasDB$OC[i]/10*0.115) * (1 - rnorm(5000, lucasDB$VGR[i], lucasDB$VGR[i]*0.4)) 
  }
 }
}




SOCsd_lcs<- -apply(MC, 1, sd)* 2 + lucasDB$SOC_lcs
SOCsdd_lcs<- apply(MC, 1, sd)* 2 + lucasDB$SOC_lcs

lucasDB<-cbind(lucasDB, "SOCsd_lcs"=SOCsd_lcs, "SOCsdd_lcs"=SOCsdd_lcs)

###############################################################ver
par(mfrow=c(1,3))
plot(lucasDB$sample_ID, lucasSHP$sample_ID)
plot(lucasDB$sample_ID - lucasSHP$sample_ID, lucasSHP$sample_ID)
plot(lucasDB$OC, BD)

