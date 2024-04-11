dirSL<-("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/SptLYR/")
dirCal <- "/eos/jeodpp/data/projects/SOIL-NACA/MEMS/Calibration/"

####grid meteo############################################################################

EOBS01<-read.dbf(paste0(dirSL, "EOBS_0.1reg_land.dbf"))
EOBS11<-read.dbf(paste0(dirSL, "EOBS_11rot_land.dbf"))

#  lat_long<-as.array(stack(paste0(dirSL,"lat_long.tif")))/100	##*100
grid2D.rcm<- cbind(EOBS01$x, EOBS01$y) 		####0.1d regular
grid2D.rcm_<- cbind(EOBS11$x, EOBS11$y) 		####11' rotated

PLT<-read.table(paste0(dirCal,"Model/mems_calibration/DB/LU_SCH.txt"), header = TRUE, sep="\t")

z <- length(slcs$NPP)

obs <- array(0.,c(365,6,z))

for (i in 1:length(slcs$NPP)) {
#for (i in 1:1) {

  if (!(slcs$LU[i]  %in% c("C10", "C20", "C30", "E10", "E20", "E30", "CR"))) next
   
  lat<- slcs$lat[i]
  long<-slcs$long[i]
  
  f<-spDistsN1(pts=grid2D.rcm, pt=as.matrix(c(long,lat)), longlat=TRUE) 
  EOB01r<-EOBS01[which.min(f),1]
  
  f_<-spDistsN1(pts=grid2D.rcm_, pt=as.matrix(c(long,lat)), longlat=TRUE) 
  EOB11r<-EOBS11[which.min(f_),3]
  
  wth<- paste("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/meteo/_clima/EOBS0.1reg_day_v25/", EOB01r, ".wth", sep="") 
  if (!file.exists(wth)) next
  
  clim<-as.matrix(read.table(wth)[10594:14245, 3:7])
  mc <- matrix(0.,nrow=365,ncol=5)
  
  for (k in 1:10) {
    mc <- mc + clim[(1+(k-1)*365):(365+(k-1)*365),1:5]
  }

  obs[,1:5,i] <- mc/10.
    
  NPPin<-PLT[PLT$LU==as.character(slcs$LU[i]),4]					###NPPin after harvest ratio
  
  NPPd<-dnorm(seq(1:365), 200, 50) * (slcs$NPP[i]*100)*NPPin		###NPP (g/m2)

  obs[,6,i]<-c(NPPd)    
}

saveRDS(obs,"driver.rds")

