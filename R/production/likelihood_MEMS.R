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
      simulated <- as.matrix(call_MEMS(try_par, list_data[[ds]]))
      
      # Difference between simulated and observed results
      diff <- c(simulated - list_data[[ds]]$obs_soc)
      
      # Log-probability density for each element in diff, weighed by uncertainty
      ll_values <- dnorm(diff, sd = c(list_data[[ds]]$unc), log = TRUE)
      
      # Increase the overall likelihood sum
      ll_sum <- ll_sum + sum(ll_values)
    }
  
  # Return likelihood
  return(ll_sum)
}