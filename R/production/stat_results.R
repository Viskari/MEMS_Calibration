# Calculate and save statistics for the results of multi-chain mcmc runs.

# Inputs:
# path_results_run - path to the results of the current run 
# run_name - name of the run
# out - outputs from the mcmc run
# n_samples - amount of samples to be taken from the posterior
# plot_chains - should individual chains be plotted


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Start of statistics logic

stat_mcmc_results <- 
  function(path_results_run, run_name, out, default_par, par_chosen, 
           inits, likelihood, burn_in = 0) {

    
# -------------------------------------------------------------------------
# Save starting values

# AM
if (out[[1]]$settings$sampler == "Metropolis") {
  for (j in 1:length(out)) {
    write.csv(
      t(inits[j, ]),
      file = paste0(
        path_results_run, "/ch", j,
        "/inits_", run_name, "_", "ch", j, ".txt"
      ),
      row.names = FALSE
    )
  }
  # DE samplers
} else if (out[[1]]$settings$sampler %in% c("DEzs", "DREAMzs")) {
  for (j in 1:length(out)) {
    write.csv(
      inits,
      file = paste0(
        path_results_run, "/ch", j,
        "/inits_", run_name, "_", "ch", j, ".txt"
      ),
      row.names = FALSE
    )
  }
}


# -------------------------------------------------------------------------
# Summaries

# Save overall light summary
print("Saving run summary.")
capture.output(
  light_summary(out, run_name, burn_in = burn_in),
  file = paste0(path_results_run, "/summary_", run_name, ".txt")
)

# Save summary for each chain (no respect for burn-in). TODO: light summary
for (ch in 1:length(out)) {
  print(paste0("Saving summary for chain ", ch, " out of ", n_chains, "."))
  ch_name <- paste0("ch", ch)
  capture.output(
    light_summary(out[[ch]], run_name, burn_in = burn_in),
    file = paste0(
      path_results_run, "/", ch_name,
      "/summary_", run_name, "_", ch_name, ".txt"
    )
  )
}


# -------------------------------------------------------------------------
# MAP likelihoods

# Define starting value for sampling based on burn-in
if (out[[1]]$settings$sampler == "Metropolis"){
  start <- floor(burn_in) + 1
  # DE sampler subchains need division so sampling does not go out of bounds 
} else if (out[[1]]$settings$sampler %in% c("DEzs", "DREAMzs")) {
  start <- floor(burn_in / 3) + 1
}

# Save likelihoods associated with produced MAP values
capture.output(
  print(sprintf(
    "Default parameters -- likelihood: %.2f",
    likelihood(default_par$value[par_chosen])
  )),
  print(sprintf(
    "All chains' MAP -- likelihood: %.2f",
    likelihood(MAP(out, start = start)$parametersMAP)
  )),
  invisible(sapply(1:length(out), function(ch) {
    print(sprintf(
      "Chain %d MAP -- likelihood: %.2f",
      ch, likelihood(MAP(out[[ch]], start = start)$parametersMAP)
    ))
  })),
  file = paste0(path_results_run, "/MAP_lhoods_", run_name, ".txt")
)

} # end function


# -------------------------------------------------------------------------

# This function provides a more focused and computationally light summary than
# the BayesianTools "summary()". The light summary also properly respects
# burn-in, unlike the BayesianTools version. Only works for multi-chain runs.
light_summary <- function(out, run_name, burn_in = 0, thin = 1) {
  
  if (length(out) < 2L) {stop("Light summary only works for multi-chain runs")}

  # Basic chain information
  if ("mcmcSamplerList" %in% class(out)) {
    sampler <- out[[1]]$settings$sampler
    runtime <- out[[1]]$settings$runtime[[3]]
    n_par <- out[[1]]$setup$numPars
    n_iter <- out[[1]]$settings$iterations
    n_chain <- length(out)
  } else if ("mcmcSampler" %in% class(out)) {
    sampler <- out$settings$sampler
    runtime <- out$settings$runtime[[3]]
    n_par <- out$setup$numPars
    n_iter <- out$settings$iterations
    n_chain <- 1L
  }
      
  # Define starting value for sampling based on burn-in
  if (sampler == "Metropolis"){
    start <- floor(burn_in) + 1
  # DE sampler subchains need division so sampling does not go out of bounds 
  } else if (sampler %in% c("DEzs", "DREAMzs")) {
    start <- floor(burn_in / 3) + 1
  }
  
  # Calculate MAPs for each parameter with BTools functions
  MAPs <- MAP(out, start = start, thin = thin)$parametersMAP
  # Calculate G-R convergence for cases with chains
  if ((sampler == "Metropolis") & ("mcmcSampler" %in% class(out))) {
    gel_rub_single <- matrix(rep(NA, n_par * 2), ncol = 2)
    gel_rub_multi <- NA
  } else {
    gel_rub <- gelmanDiagnostics(out, start = start, thin = thin)
    gel_rub_single <- gel_rub$psrf
    gel_rub_multi <- round(gel_rub$mpsrf, 3)
  }
  
  # Take sample from chain in coda format
  chain <- getSample(out, parametersOnly = T, coda = T, start = start, thin = thin)
    
  # Statistics from coda chain
  rej_rate <- round(mean(coda::rejectionRate(chain)), 3)
  eff_samp <- round(mean(coda::effectiveSize(chain)), 0)
  med <- numeric(n_par)
  for (j in 1:n_par) {
    med[j] <- median(unlist(chain[, j]))
  }
  
  # Summarize summary
  par_summary <- round(cbind(MAPs, med, gel_rub_single), 3)
  colnames(par_summary) <- c("MAP", "median", "G-R point", "G-R upper CI")
  
  # Print everything
  cat(rep("#", 25), "\n")
  cat("## MCMC light summary ##","\n")
  cat(rep("#", 25), "\n", "\n")
  cat("# Run name: ", run_name, "\n")
  cat("# MCMC sampler: ", sampler, "\n")
  if (sampler == "Metropolis") {
    cat("# Nr. Chains: ", n_chain, "\n")
  } else if (sampler %in% c("DEzs", "DREAMzs")) {
    cat("# Nr. Chains: ", n_chain, " (x 3)\n")
  }
  cat("# Chain runtime: ", runtime, " sec\n")
  cat("# Iterations per chain: ", n_iter, "\n")
  cat("# Burn-in per chain: ", floor(burn_in), "\n")
  cat("# Avg. rejection rate: ", rej_rate, "\n")
  cat("# Avg. effective sample size: ", eff_samp, "\n")
  cat("# Gelman-Rubin multivariate psrf: ", gel_rub_multi, "\n")
  cat("\n# Parameters and convergence: \n\n")
  print(par_summary)
  cat("\n")
}
