# The function for calling YASSO15 and the likelihood are defined in this file.


# -------------------------------------------------------------------------
# Function to call YASSO15 from likelihood. Can be freely changed, as long as
# likelihood is modified accordingly. By default the calibrate_yasso() function
# from Ryassofortran.

call_yasso <- function(try_par, ld, sspred = 0L) {
  
  # Run yasso (typeset parameters to be safe)
  soil_c <- calibrate_yasso(
    par = as.double(try_par),
    n_runs = ld$n_runs,
    time = ld$time,
    temp = ld$temp,
    prec = ld$prec,
    init = ld$init,
    litter = ld$litter,
    wsize = ld$size,
    leac = ld$leac,
    sspred = sspred
  )
  
  # Return the simulated soil carbon (AWENH)
  soil_c
}


# -------------------------------------------------------------------------
# Likelihood function

likelihood <- function(trial_set) {

  # Initialize likelihood sum
  ll_sum <- double(1)
  
  # Setup proposal: use the default parameter values as basis
  try_par <- default_par$value
  # Override the to-be-calibrated parameter values with the proposal (trial) set
  try_par[par_chosen] <- trial_set

  
  # Release constraint: YASSO is only allowed to release carbon from certain
  # pools. The p-parameters of the other pools are summed up to one.
  # A-pool
  try_par[8] <- 1 - try_par[31] # pAW = 1 - pH
  # E-pool
  try_par[6] <- max(1 - try_par[9] - try_par[31], 0) # pEA = 1 - pEW - pH (>= 0)
  # N-pool
  try_par[7] <- 1 - try_par[31] # pNA = 1 - pH

    
  # Sanity check: the p-parameters of a pool cannot sum up to more than one,
  # since this would indicate there is more than 100 % of carbon in the pool.
  # pA_WENH
  if (sum(try_par[8], try_par[11], try_par[14], try_par[31]) > 1) {
    return(-Inf) # reject proposal
  }
  # pW_AENH
  if (sum(try_par[5], try_par[12], try_par[15], try_par[31]) > 1) {
    return(-Inf)
  }
  # pE_AWNH
  if (sum(try_par[6], try_par[9], try_par[16], try_par[31]) > 1) {
    return(-Inf)
  }
  # pN_AWEH
  if (sum(try_par[7], try_par[10], try_par[13], try_par[31]) > 1) {
    return(-Inf)
  }
    
  # Humus flow correction: there is no flow to H pool in the litter bag
  # measurements, hence pH is set to zero
  litter_par <- try_par
  litter_par[31] <- 0
  
  # Correct the other p-parameters accordingly
  # A-pool
  litter_par[8] <- 1 # pAW = 1 - pH
  # E-pool
  litter_par[6] <- 1 - litter_par[9] # pEA = 1 - pEW - pH
  # N-pool
  litter_par[7] <- 1 # pNA = 1 - pH
  
  # TODO: Should the litter_par be sanity checked too? Overhaul this section.
  
  
  # Loop through the datasets
  for (ds in datasets) {
    
    # ED1/2
    if (ds %in% c("ed1", "ed2")) {
      
      # Write the value of the leaching parameter(s) to the input leac
      list_data[[ds]][["leac"]] <- as.double(try_par[17])
      
      # Run YASSO15 with the humus-corrected proposal set and the site data
      simulated <- call_yasso(litter_par, list_data[[ds]])
      
      # The simulated N + H pool should be compared to the measured N pool
      # TODO: can be removed, since there is no humus flow for litter bags
      simulated[, 4] <- simulated[, 4] + simulated[, 5]
      
      # Choose model results for the pools (AWEN) that have observations
      simulated <- simulated[, 1:4]
      
      # Difference between simulated and observed results
      diff <- c(simulated - list_data[[ds]][["obs"]])
      
      # Log-probability density for each element in diff, weighed by uncertainty
      ll_values <- dnorm(diff, sd = c(list_data[[ds]][["unc"]]), log = TRUE)
      
      # Increase the overall likelihood sum
      ll_sum <- ll_sum + sum(ll_values)
    }
    
    # CIDET, LIDET, Tarasov, Mäkinen
    if (ds %in% c("cidet", "lidet", "tarasov", "makinen", "hob3")) {

      # Write the value of the leaching parameter(s) to the input leac
      if (ds == "cidet") {list_data[[ds]][["leac"]] <- as.double(try_par[18])}
      if (ds == "lidet") {list_data[[ds]][["leac"]] <- as.double(try_par[21])}
      if (ds == "hob3") {list_data[[ds]][["leac"]] <- as.double(try_par[19])}
      
      # Run YASSO15 with the humus-corrected proposal set and the site data
      simulated <- call_yasso(litter_par, list_data[[ds]])
      
      # Sum the results over AWENH, since observations are summed 
      simulated <- rowSums(simulated)
      
      # Difference between simulated and observed results
      diff <- c(simulated - list_data[[ds]][["obs"]])
      
      # Log-probability density for each element in diff, weighed by uncertainty
      ll_values <- dnorm(diff, sd = c(list_data[[ds]][["unc"]]), log = TRUE)
      
      # Increase the overall likelihood sum
      ll_sum <- ll_sum + sum(ll_values)
    }
    
    # GSS
    if (ds %in% c("gss")) {
      
      # Run YASSO15 with the proposal set and the site data. For GSS three runs
      # have to be made and summed, since the observations consist of carbon of
      # three sizes.
      
      # Run YASSO15 for each size
      sim_nonw <- call_yasso(try_par, list_data[[ds]]$nonw, sspred = 1L)
      sim_smallw <- call_yasso(try_par, list_data[[ds]]$smallw, sspred = 1L)
      sim_largew <- call_yasso(try_par, list_data[[ds]]$largew, sspred = 1L)
      
      # Sum the results for the three sizes
      simulated <- (sim_nonw + sim_smallw + sim_largew)
      
      # Sum the results over AWENH, since observations are summed
      simulated <- rowSums(simulated)
      
      # Difference between simulated and observed results
      diff <- c(simulated - list_data[[ds]][["obs"]])
      
      # Log-probability density for each element in diff, weighed by uncertainty
      ll_values <- dnorm(diff, sd = c(list_data[[ds]][["unc"]]), log = TRUE)
      
      # Increase the overall likelihood sum
      ll_sum <- ll_sum + sum(ll_values)
    }
    
  }
  
  # Return likelihood
  return(ll_sum)
}