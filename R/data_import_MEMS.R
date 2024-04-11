# Import data for YASSO calibration. This function shapes the data as far as
# possible, so that likelihood does not have to do additional operations.

import_data <- function(path, datasets, datatype = "full") {
  # in: path - path to the calibration data folder
  # in: datasets - character vector with names of datasets to be imported
  # in: datatype - import full ("full"), training ("train") or testing ("test") data
  # out: list of imported datasets
  
  # Create a list element for each dataset
  list_data <- vector("list", length = length(datasets))
  names(list_data) <- datasets
  
  # Loop over datasets
  for (ds in datasets) {
        
    # Import data for standard calibration
    if (datatype == "full") {
      input_data <- read.csv(
        file = paste0(path, ds, "/full_data_", ds, ".csv")
      )
    # Import data for training
    } else if (datatype == "train") {
      input_data <- read.csv(
        file = paste0(path, ds, "/train_data_", ds, ".csv"),
        colClasses = "numeric"
      )
    # Import data for testing 
    } else if (datatype == "test") {
      input_data <- read.csv(
        file = paste0(path, ds, "/test_data_", ds, ".csv"),
        colClasses = "numeric"
      )
    }

    Env_dr <- readRDS(paste0(path, ds, "/driver.rds"))
    
    # Select which data points to use depending on if POM/MOM data 
    # or total carbon
    
    if (ds %in% c("LUCAS_PM")) {
      N_start <- 1
      N_end <- 348
    } else if (ds %in% c("LUCAS_tot")) {
      N_start <- 349
      N_end <- 19824
    }
            
              
    # Standardized inputs

    obs_soc <- input_data$SOC_lcs[N_start:N_end]
    obs_pom <- input_data$OC_pom_g_kg[N_start:N_end]
    obs_maom <- input_data$OC_sc_g_kg[N_start:N_end]
    obs_frac <- obs_maom/(obs_pom+obs_maom)
        
    unc_soc <- input_data$SOCsd_lcs[N_start:N_end]
    unc_pom <- 0.1*obs_pom
    unc_maom <- 0.1*obs_maom
    unc_frac <- 0.05*obs_maom_frac
 
    LU <- input_data$LU[N_start:N_end]
    clay <- input_data$clay[N_start:N_end]
    sand <- input_data$sand[N_start:N_end]
    coarse <- input_data$coarse[N_start:N_end]/100.
    OC <- input_data$OC[N_start:N_end]
    pH <- input_data$pH_in_H2O[N_start:N_end]

    lat <- input_data$GPS_LAT[N_start:N_end]
    long <- input_data$GPS_LONG[N_start:N_end]

    NPP <- input_data$NPP[N_start:N_end]
    
    maxTmp <- Env_dr[,3,N_start:N_end]
    minTmp <- Env_dr[,4,N_start:N_end]
    prec <- Env_dr[,5,N_start:N_end]    
    aNPP <- Env_dr[,6,N_start:N_end]

    cNPP <- colSums(aNPP)
    
    # Calculating inputs required by MEMS
    
    BDm<- 0.80806 + (0.823844 * exp(-0.27993 * OC * 0.1))+(0.0014065 * sand)-(0.0010299 * clay)   #BD Fine Earth 
    BDm_c<- (100- coarse)/((coarse/2.65)+(100-coarse)/BDm)							 			#BD FE without gravel
    PORm<- (1-BDm/2.65)
    sk_p<- (BDm - BDm_c)/BDm																 							#volume gravel
    
    WC33   <-   (0.2576 +(-0.002 * sand)+(0.0036 * clay)+ 0.0299 * min(OC* 0.1 * 1.72, 6))
    WC1500 <-   (0.026 + (0.005 * clay)+ 0.0158 * min(OC* 0.1 * 1.72, 6)) 
    
    Ks<- (1930*(PORm-WC33)^(3- log(WC1500/WC33)/log(33/1500)) )/10 													##Ks mm/h
    
    WC33s<- WC33 * (1-sk_p)																							##correction for gravel
    WC1500s <- WC1500 * (1-sk_p)

    LitN <- rep(1,length(NPP))
    fsi <- rep(1,length(NPP))
    flig <- rep(1,length(NPP))
    fdoc <- rep(1,length(NPP))
    
    depth <- rep(20,length(NPP))

    LitN[which(LU == "E10" | LU == "E20" | LU == "E30" | LU == "CR")] <- 1.1
    fsi[which(LU == "E10" | LU == "E20" | LU == "E30" | LU == "CR")] <- 0.35
    flig[which(LU == "E10" | LU == "E20" | LU == "E30" | LU == "CR")] <- 0.15
    fdoc[which(LU == "E10" | LU == "E20" | LU == "E30" | LU == "CR")] <- 0.35

    LitN[which(LU == "C10")] <- 1.32
    fsi[which(LU == "C10")] <- 0.4
    flig[which(LU == "C10")] <- 0.27
    fdoc[which(LU == "C10")] <- 0.15
    
    LitN[which(LU == "C20")] <- 0.41
    fsi[which(LU == "C20")] <- 0.35
    flig[which(LU == "C20")] <- 0.32
    fdoc[which(LU == "C20")] <- 0.15
    
    LitN[which(LU == "C30")] <- 0.87
    fsi[which(LU == "C30")] <- 0.375
    flig[which(LU == "C30")] <- 0.295
    fdoc[which(LU == "C30")] <- 0.15

    in_comp <- which(!(is.na(NPP) | LU == "SHR" |  cNPP == 0.))
            
    # Save imported variables to a list
    # Every dataset gets a single list with all the necessary data
    list_data[[ds]] <- list(
      obs_soc = obs_soc[in_comp],
      obs_pom = obs_pom[in_comp],
      obs_maom = obs_maom[in_comp],
      obs_frac = obs_frac[in_comp],
      unc_soc = unc_soc[in_comp],
      unc_pom = unc_pom[in_comp],
      unc_maom = unc_maom[in_comp],
      unc_frac = unc_frac[in_comp],
      LU = LU[in_comp],
      clay = clay[in_comp],
      sand = sand[in_comp],
      coarse = coarse[in_comp],
      OC = OC[in_comp],
      pH = pH[in_comp],
      depth = depth[in_comp],
      lat = lat[in_comp],
      long = long[in_comp],
      NPP = NPP[in_comp],
      BDm = BDm[in_comp],
      WC33s = WC33s[in_comp],
      WC1500s = WC1500s[in_comp],
      Ks = Ks[in_comp],
      maxTmp = maxTmp[,in_comp],
      minTmp = minTmp[,in_comp],
      prec = prec[,in_comp],
      aNPP = aNPP[,in_comp],
      LitN = LitN[in_comp],
      fsi = fsi[in_comp],
      flig = flig[in_comp],
      fdoc = fdoc[in_comp]
    )
    
  } # end for over datasets
  
  # Return the list with imported data for each dataset
  list_data
  
} # end function 
