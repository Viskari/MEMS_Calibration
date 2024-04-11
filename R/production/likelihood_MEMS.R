# The MEMS likelihood are defined in this file.
library(raster)
library(foreign)
library(stringr)


# -------------------------------------------------------------------------
# Likelihood function

likelihood <- function(trial_set) {

  # Initialize likelihood sum
  ll_sum <- double(1)
  
  # Setup proposal: use the default parameter values as basis
  try_par <- default_par$value
  # Override the to-be-calibrated parameter values with the proposal (trial) set
  try_par[par_chosen] <- trial_set

  # Loop through the datasets
  for (ds in datasets) {
    
      # Run 
      # simulated <- as.matrix(call_MEMS(try_par, list_data[[ds]]))
      simulated <- as.matrix(calib_MEMS(try_par, list_data[[ds]]))
      
      # Difference between simulated and observed results
      diff_soc <- c(simulated[,1] - list_data[[ds]]$obs_soc)
      
      # Log-probability density for each element in diff, weighed by uncertainty
      ll_values_soc <- dnorm(diff_soc, sd = c(list_data[[ds]]$unc_soc), log = TRUE)

      # Difference between simulated and observed results
      diff_frac <- c(simulated[,2] - list_data[[ds]]$obs_frac)
      
      # Log-probability density for each element in diff, weighed by uncertainty
      ll_values_frac <- dnorm(diff_frac, sd = c(list_data[[ds]]$unc_frac), log = TRUE)
      
      # Increase the overall likelihood sum
      ll_sum <- ll_sum + sum(ll_values_soc) + sum(ll_values_frac)
    }
  
  # Return likelihood
  return(ll_sum)
}