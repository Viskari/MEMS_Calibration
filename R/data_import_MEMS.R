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
    
    # Standardized inputs

    N <- 2
    
    obs_soc <- input_data$SOC_lcs[1:N]
    obs_pom <- input_data$OC_pom_g_kg[1:N]
    obs_maom <- input_data$OC_sc_g_kg[1:N]
    
    unc <- input_data$SOCsd_lcs[1:N]
 
    LU <- input_data$LU[1:N]
    clay <- input_data$clay[1:N]
    sand <- input_data$sand[1:N]
    coarse <- input_data$coarse[1:N]
    OC <- input_data$OC[1:N]
    pH <- input_data$pH_in_H2O[1:N]

    lat <- input_data$GPS_LAT[1:N]
    long <- input_data$GPS_LONG[1:N]

    NPP <- input_data$NPP[1:N]

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
                    
    # Save imported variables to a list
    
    # Every dataset gets a single list with all the necessary data
    list_data[[ds]] <- list(
      obs_soc = obs_soc,
      obs_pom = obs_pom,
      obs_maom = obs_maom,
      unc = unc,
      LU = LU,
      clay = clay,
      sand = sand,
      coarse = coarse,
      OC = OC,
      pH = pH,
      lat = lat,
      long = long,
      NPP = NPP,
      BDm = BDm,
      WC33s = WC33s,
      WC1500s = WC1500s,
      Ks = Ks
    )
    
  } # end for over datasets
  
  # Return the list with imported data for each dataset
  list_data
  
} # end function 
