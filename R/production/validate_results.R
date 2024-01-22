# This function validates the trained model. It runs YASSO with the calibrated
# MAPs and test datasets and compares the simulated results against
# observations.

# TODO: Currently this function only calculates RMSEs. Add other stats and
# plots, once it becomes clear what could be useful.


# -------------------------------------------------------------------------

validate_results <- function(path_main_dir, run_name, datasets, testsets, out, default_par, par_chosen, burn_in = 0) {

    
# -------------------------------------------------------------------------
# Setup output path and import data

# Path to save the validation results to
path_run_val <- paste0(path_main_dir, "results/", run_name, "/validation")

# Import testing data
list_data_test <- import_data(paste0(path_main_dir, "data/"), testsets, "test")


# -------------------------------------------------------------------------
# Setup trained parameters

# Define starting value for sampling based on sampler choice, burn-in
sampler <- out[[1]]$settings$sampler
if (sampler == "Metropolis"){
  start <- floor(burn_in) + 1
  # DE sampler subchains need division so sampling does not go out of bounds 
} else if (sampler %in% c("DEzs", "DREAMzs")) {
  start <- floor(burn_in / 3) + 1
}

# Obtain trained parameter values from MAPs, respect burn-in
trained_maps <- MAP(out, start = start)$parametersMAP

# Setup the entire parameter set
trained_par <- default_par$value
# Override the parameters which have trained values
trained_par[par_chosen] <- trained_maps

# Correct the humus flow: there is no flow to H pool in the litter bag
# measurements (that we are using for validation), hence pH is set to zero
trained_par[31] <- 0

# Correct the other p-parameters accordingly
# A-pool
trained_par[8] <- 1 # pAW = 1 - pH
# E-pool
trained_par[6] <- 1 - trained_par[9] # pEA = 1 - pEW - pH
# N-pool
trained_par[7] <- 1 # pNA = 1 - pH


# -------------------------------------------------------------------------
# Run YASSO with the calibrated MAPs and the testing data

# Create list for storing the simulation results and observations for each 
# dataset
list_validation <- vector(mode = "list", length = length(testsets))
names(list_validation) <- testsets

# Create dataframe for storing calculated RMSE values
df_diff <- data.frame(dataset = testsets, rmse = double(length(testsets)))

# Loop through the test datasets, similarly as in likelihood
for (ts in testsets) {

  # Setup the leaching parameters
  
  # If only one dataset was used in training, all test sets use the leaching for
  # that dataset
  if (all(datasets %in% "cidet")) {
    list_data_test[[ts]][["leac"]] <- as.double(trained_par[18])
  }
  if (all(datasets %in% "lidet")) {
    list_data_test[[ts]][["leac"]] <- as.double(trained_par[21])
  }
  if (all(datasets %in% c("ed1", "ed2"))) {
    list_data_test[[ts]][["leac"]] <- as.double(trained_par[17])
  }
  # If all datasets were used in training, all test sets use a trained leaching
  if (all(c("cidet", "lidet", "ed1", "ed2") %in% datasets)) {
    if (ts %in% c("ed1", "ed2")) {
      list_data_test[[ts]][["leac"]] <- as.double(trained_par[17])
    }
    if (ts == "cidet") {
      {list_data_test[[ts]][["leac"]] <- as.double(trained_par[18])}
    }
    if (ts == "lidet") {
      list_data_test[[ts]][["leac"]] <- as.double(trained_par[21])
    }
  }
  
  # TODO: Change above logic so that it works for calibrations with e.g. both
  # CIDET & LIDET. Could just check if the current ts dataset was calibrated.
  # But how to do the non-calibrated leachings in this case, e.g. how to choose
  # leaching for ED if calibrated with CIDET & LIDET?
  
  # ED1/2
  if (ts %in% c("ed1", "ed2")) {
    
    # Run trained model and calculate differences to observations
    
    # Run YASSO15 with the trained MAPs and the testset data
    simulated <- call_yasso(trained_par, list_data_test[[ts]])
    
    # The simulated N + H pool should be compared to the measured N pool
    simulated[, 4] <- simulated[, 4] + simulated[, 5]
    
    # Choose model results for the pools (AWEN) that have observations
    simulated <- simulated[, 1:4]
    
    # Difference between simulated and observed results
    diff <- simulated - list_data_test[[ts]][["obs"]]
    
    # Do statistics
    
    # Calculate RMSE
    rmse <- sqrt(mean(diff ^ 2))
    
    # Store the RMSE to dataframe
    df_diff[df_diff[, "dataset"] == ts, "rmse"] <- rmse
    
    # Store the difference between simulated and observed results to list
    list_validation[[ts]][["difference"]] <- diff
  }
  
  # CIDET, LIDET
  if (ts %in% c("cidet", "lidet","hob3")) {
    
    # Run trained model and calculate differences to observations
    
    # Run YASSO15 with the trained MAPs and the testset data
    simulated <- call_yasso(trained_par, list_data_test[[ts]])
    
    # Sum the results over AWENH, since observations are summed 
    simulated <- rowSums(simulated)
    
    # Difference between simulated and observed results
    diff <- simulated - list_data_test[[ts]][["obs"]]
    
    # Do statistics
    
    # Calculate RMSE
    rmse <- sqrt(mean(diff ^ 2))
    
    # Store the rmse to dataframe
    df_diff[df_diff[, "dataset"] == ts, "rmse"] <- rmse
    
    # Store the difference between simulated and observed results to list
    list_validation[[ts]][["difference"]] <- diff
  }
}


# -------------------------------------------------------------------------
# Calculate additional statistics

# Calculate the combined RMSE for ED1/2

# Combine the differences into a single matrix
ed_total_diff <- rbind(
  list_validation[["ed1"]][["difference"]],
  list_validation[["ed2"]][["difference"]]
)

# Calculate the combined RMSE
ed_rmse <- sqrt(mean(ed_total_diff ^ 2))

# Store to the RMSE dataframe
df_diff <- rbind(df_diff, data.frame(dataset = "eurodeco", rmse = ed_rmse))


# -------------------------------------------------------------------------
# Save the statistics

# Save calculated RMSEs
write.csv(
  df_diff, file = paste0(path_run_val, "/rmse_", run_name, ".csv"), 
  row.names = FALSE
)


# -------------------------------------------------------------------------
# TODO: Draw relevant plots based on the validation results

# At the time of writing, I could not think of any plots that would be useful in
# every situation. Once the inspiration and need comes up, they can be drawn
# here using the results obtained above. You can e.g. store necessary
# information to list_validation during the runs and use it here to plot.

}
