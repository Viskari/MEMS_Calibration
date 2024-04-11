calib_MEMS <- function(Para,driver) {

  ii <- length(driver$clay)

  soil <- rep(1,7)
  LQ <- rep(1,4)
  env_in <- matrix(0.,nrow=365,ncol=4)

  POM <- rep(0.,ii)
  MOM <- rep(0.,ii)
        
  for(jj in 1:ii){
    C <- rep(0,11)
    
    soil[1] <- driver$sand[jj]
    soil[2] <- driver$clay[jj]
    soil[3] <- driver$coarse[jj]
    soil[4] <- driver$pH[jj]
    soil[5] <- driver$OC[jj]
    soil[6] <- driver$BDm[jj]
    soil[7] <- driver$depth[jj]
    
    LQ[1] <- driver$LitN[jj]
    LQ[2] <- driver$fsi[jj]
    LQ[3] <- driver$flig[jj]
    LQ[4] <- driver$fdoc[jj]
    
    env_in[,1] <- driver$maxTmp[,jj]
    env_in[,2] <- driver$minTmp[,jj]
    env_in[,3] <- driver$prec[,jj]
    env_in[,4] <- driver$aNPP[,jj]

    res <- .Fortran("mems",Para,env_in,soil,LQ,C)

    POM[jj] <- (res[[5]][5] +  res[[5]][10])/100.
    MOM[jj] <- res[[5]][9]/100.
  }

  soilC <- POM + MOM
  fracMOM <- MOM/soilC
  
  OutC <- cbind(soilC,fracMOM)
  return(OutC)
      
}