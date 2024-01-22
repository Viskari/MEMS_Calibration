call_MEMS <- function(ParVec,slcs) { 

#rm(list=ls(all=TRUE))

###set dir reading raster
dirSL<-("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/SptLYR/")
dirSC<-("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/statCRP/")
dirMEMS <- "/eos/jeodpp/data/projects/SOIL-NACA/MEMS/Calibration/Model/mems_calibration/MEMS/"
dirCal <- "/eos/jeodpp/data/projects/SOIL-NACA/MEMS/Calibration/"

setwd(dirMEMS)

###read LUCAS all LU 2009
LUsub<-c("C10", "C20", "C30", "E10", "E20", "E30", "CR")


####grid meteo############################################################################
  
  EOBS01<-read.dbf(paste0(dirSL, "EOBS_0.1reg_land.dbf"))
  EOBS11<-read.dbf(paste0(dirSL, "EOBS_11rot_land.dbf"))
 
  #  lat_long<-as.array(stack(paste0(dirSL,"lat_long.tif")))/100	##*100
  grid2D.rcm<- cbind(EOBS01$x, EOBS01$y) 		####0.1d regular
  grid2D.rcm_<- cbind(EOBS11$x, EOBS11$y) 		####11' rotated
  
  
  
  ###conversion LUCAS LU to MEMS plant type
  
  PLT<-read.table(paste0(dirCal,"Model/mems_calibration/DB/LU_SCH.txt"), header = TRUE, sep="\t")  

  ###read MEMS input/in_Control
  inCTR<-read.csv(paste0(dirMEMS,"input/in_Control_calibration.csv"), stringsAsFactors=FALSE)
  inST<-read.csv(paste0(dirMEMS,"input/in_Site.csv"), stringsAsFactors=FALSE)
  inSL<-read.csv(paste0(dirMEMS,"input/in_Soil.csv"), stringsAsFactors=FALSE)
  inWTH<-read.csv(paste0(dirMEMS,"input/in_Wth.csv"), stringsAsFactors=FALSE)  
  inPara<-read.csv(paste0(dirMEMS,"input/in_Basepars_init.csv"), stringsAsFactors=FALSE)  


###dim outputs
MAOM<- as.data.frame(array(data=NA, dim=c(length(slcs$NPP),1)))
POM<- as.data.frame(array(data=NA, dim=c(length(slcs$NPP),1)))

# Write the parameter file for the model runs

inPara[,2] <- c(1.01,ParVec)
write.csv(inPara,"input/in_Basepars.csv", quote = FALSE, row.names=FALSE)


###########################################################################################

 for (i in 1:length(slcs$NPP)) {
   #lcs<-subset(lcs, LU %in% LUsub)		###subset LU
   if (!(slcs$LU[i]  %in% c("C10", "C20", "C30", "E10", "E20", "E30", "CR"))) next	
   
### in_Control ###########################################################################
 SCH_N<-PLT[PLT$LU==as.character(slcs$LU[i]),2]
 inCTR[15,2]<-SCH_N
 inCTR[22,2]<-SCH_N
 
 write.csv(inCTR,"input/in_Control.csv", quote = FALSE, row.names=FALSE)
 
### in_Control ###########################################################################
 
 
### in_Soil ##############################################################################
 BDm<- 0.80806 + (0.823844 * exp(-0.27993 * slcs$OC[i] * 0.1))+(0.0014065 * slcs$sand[i])-(0.0010299 * slcs$clay[i])   #BD Fine Earth 
 BDm_c<- (100- slcs$coarse[i])/((slcs$coarse[i]/2.65)+(100-slcs$coarse[i])/BDm)							 			#BD FE without gravel
 PORm<- (1-BDm/2.65)
 sk_p<- (BDm - BDm_c)/BDm																 							#volume gravel

 WC33   <-   (0.2576 +(-0.002 * slcs$sand[i])+(0.0036 * slcs$clay[i])+ 0.0299 * min(slcs$OC[i]* 0.1 * 1.72, 6))
 WC1500 <-   (0.026 + (0.005 * slcs$clay[i])+ 0.0158 * min(slcs$OC[i]* 0.1 * 1.72, 6)) 
 
 Ks<- (1930*(PORm-WC33)^(3- log(WC1500/WC33)/log(33/1500)) )/10 													##Ks mm/h
 
 WC33s<- WC33 * (1-sk_p)																							##correction for gravel
 WC1500s <- WC1500 * (1-sk_p)
 Kss <- Ks *(1-sk_p)
 
 ### input In_soil 
 inSL[5, 4]<- BDm
 inSL[5, 5]<- slcs$sand[i]
 inSL[5, 6]<- slcs$clay[i]
 inSL[5, 7]<- slcs$coarse[i]
 inSL[5, 8]<- WC33s
 inSL[5, 9]<- WC1500s 
 inSL[5, 10]<- Ks
 inSL[5, 11]<- slcs$pH[i]
  
 write.csv(inSL,"input/in_Soil.csv", quote = FALSE, row.names=FALSE) 
  
 ### in_Soil ##############################################################################
 
 
 ### in_site ###########################################################################
 lat<- slcs$lat[i]
 long<-slcs$long[i]
 
 inST[2,2]<-lat
  inST[3,2]<-long
  
 write.csv(inST,"input/in_Site.csv", quote = FALSE, row.names=FALSE)
 
 ### in_site ###########################################################################
 
 
 
 ### in_WTH ###########################################################################
 f<-spDistsN1(pts=grid2D.rcm, pt=as.matrix(c(long,lat)), longlat=TRUE) 
 EOB01r<-EOBS01[which.min(f),1]
 
 f_<-spDistsN1(pts=grid2D.rcm_, pt=as.matrix(c(long,lat)), longlat=TRUE) 
 EOB11r<-EOBS11[which.min(f_),3]
 
 wth<- paste("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/meteo/_clima/EOBS0.1reg_day_v25/", EOB01r, ".wth", sep="") 
 if (!file.exists(wth)) next 
 
 
 obs<-read.table(wth)[10594:10958, 3:7]                 ###read 2009-2018
 
 NPPin<-PLT[PLT$LU==as.character(slcs$LU[i]),4]					###NPPin after harvest ratio
 
 NPPd<-dnorm(seq(1:365), 200, 50) * (slcs$NPP[i]*100)*NPPin		###NPP (g/m2)
 NPPdd<-dnorm(seq(1:366), 200, 50) * (slcs$NPP[i]*100)*NPPin
 obs[,6]<-c(NPPd)
 
 names(obs)<-names(inWTH)
 inWTH<-rbind(inWTH[1:3,], obs)
  write.csv( inWTH, "input/in_Wth.csv", quote = FALSE, row.names=FALSE)
 
  ### in_WTH ###########################################################################
 
 
 
###shell #################################################################################
bat<-(paste("java MEMS input/in_Control.csv > dirlog.log", sep=""))
system(bat, wait=TRUE)
###shell #################################################################################


###out ###################################################################################
if (!file.exists("input/yearlySOM.csv")) { next				###jump in case of missing sim
} else { 
out<-read.csv("input/yearlySOM.csv")
 file.remove("input/yearlySOM.csv")
 file.remove("input/yearlySurfOM.csv")
 } 

MAOM[i,]<- out[,8]/100                #MAOM (C t/ha)
 POM[i,]<-(out[,4] + out[,10])/100   #POM (C t/ha)

 
 # write.csv(out, "input/yearlySOM_.csv", quote = FALSE, row.names=FALSE)
###out ###################################################################################

}
 
###write ##################################################################################
# write.dbf(SOC, paste("/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas/out/SOC_", fold, ".dbf", sep=""))
# write.csv(MAOM, paste("/eos/jeodpp/data/projects/SOIL-NACA/MEMS//Calibration/Model/mems_calibration/_out/MAOM_", fold, ".csv", sep=""), quote = FALSE, row.names=FALSE)
# write.csv(POM, paste("/eos/jeodpp/data/projects/SOIL-NACA/MEMS//Calibration/Model/mems_calibration/_out/POM_", fold, ".csv", sep=""), quote = FALSE, row.names=FALSE)

setwd(dirCal)

return(MAOM+POM)
#return(MAOM)

warnings()

}