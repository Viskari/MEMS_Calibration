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
        file = paste0(path, ds, "/full_data_", ds, ".csv"),
        colClasses = "numeric"
      )
    # Import data for training YASSO
    } else if (datatype == "train") {
      input_data <- read.csv(
        file = paste0(path, ds, "/train_data_", ds, ".csv"),
        colClasses = "numeric"
      )
    # Import data for testing YASSO 
    } else if (datatype == "test") {
      input_data <- read.csv(
        file = paste0(path, ds, "/test_data_", ds, ".csv"),
        colClasses = "numeric"
      )
    }
    
    # Convert to matrix
    input_data <- as.matrix(input_data)
    
    # Define n_runs
    n_runs <- as.integer(nrow(input_data))
    
    # Define leaching
    leac <- as.double(0)
    
    
    # Standardized inputs
    
    # Initial values
    init <- input_data[, c("A_init", "W_init", "E_init", "N_init", "H_init")]
    
    # Times in years
    time <- input_data[, "time"]
    
    # Temperature
    temp <- input_data[, c(paste0("T_0", 1:9), paste0("T_", 10:12))]
    
    # Precipitation
    prec <- input_data[, "prec"]
    
    
    # Inputs that vary between data sets
    
    # ED1/2 observations and uncertainties are separate for AWEN (4 columns)
    if (ds %in% c("ed1", "ed2")) {
      
      # Observations
      obs <- input_data[, c("A_obs", "W_obs", "E_obs", "N_obs")]
      
      # Uncertainty
      unc <- input_data[, c("A_unc", "W_unc", "E_unc", "N_unc")]
      
    # Other observations and uncertainties are summed over AWEN (1 column)
    } else if (ds %in% c("cidet", "lidet", "gss", "tarasov", "makinen", "hob3")) {
      
      # Observations
      obs <- input_data[, "c_obs"]
      
      # Uncertainty
      unc <- input_data[, "c_unc"]
    }
    
    # GSS has multiple size and litter input values for each observation
    if (ds %in% c("gss")) {
      
      # Size
      size <- input_data[, c("size_nonw", "size_smallw", "size_largew")]
      
      # Litter input
      litter <- input_data[, c(
        "A_in_nonw", "W_in_nonw", "E_in_nonw","N_in_nonw", "H_in_nonw",
        "A_in_smallw", "W_in_smallw", "E_in_smallw", "N_in_smallw", "H_in_smallw",
        "A_in_largew", "W_in_largew", "E_in_largew", "N_in_largew", "H_in_largew"
        )]
      
    # Other datasets share a common format with a single value per observation
    } else if (ds %in% c("ed1", "ed2", "cidet", "lidet", "tarasov", "makinen", "hob3")) {
      
      # Size
      size <- input_data[, "size"]
      
      # Litter input
      litter <- input_data[, c("A_in", "W_in", "E_in", "N_in", "H_in")]
    }
 
    
    # Save imported variables to a list
    
    # Every dataset gets a single list with all the necessary data
    list_data[[ds]] <- list(
      obs = obs,
      unc = unc,
      n_runs = n_runs,
      time = time,
      temp = temp,
      prec = prec,
      init = init,
      litter = litter,
      size = size,
      leac = leac
    )
    
    # For GSS we use three sublists - one for each size group
    if (ds %in% c("gss")) {
      list_gss <- list(
        nonw = list_data[[ds]][!names(list_data[[ds]]) %in% c("obs", "unc")],
        smallw = list_data[[ds]][!names(list_data[[ds]]) %in% c("obs", "unc")],
        largew = list_data[[ds]][!names(list_data[[ds]]) %in% c("obs", "unc")],
        # Store obs and unc separately, since they are not YASSO inputs
        obs = obs,
        unc = unc
      )
      
      # Pick the correct woody size and litter input for each size group
      list_gss$nonw$size <- list_data[[ds]][["size"]][, 1]
      list_gss$smallw$size <- list_data[[ds]][["size"]][, 2]
      list_gss$largew$size <- list_data[[ds]][["size"]][, 3]
      
      list_gss$nonw$litter <- list_data[[ds]][["litter"]][, 1:5]
      list_gss$smallw$litter <- list_data[[ds]][["litter"]][, 6:10]
      list_gss$largew$litter <- list_data[[ds]][["litter"]][, 11:15]
      
      # Overwrite
      list_data[[ds]] <- list_gss
    }
    
  } # end for over datasets
  
  # Return the list with imported data for each dataset
  list_data
  
} # end function 
