# Plot the summarized results of multi-chain mcmc runs with ggplot. The mcmc
# output should be saved as an object with the name "out".

# Inputs:
# path_results_run - path to the results of the current run 
# run_name - name of the run
# out - outputs from the mcmc run
# n_samples - amount of samples to be taken from the posterior
# plot_chains - should individual chains be plotted


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Start of plotting logic

plot_mcmc_results <-
  function(path_results_run, path_results_out, run_name, 
           out, n_samples = 10002, burn_in = 0, plot_chains = FALSE) {
    
# -------------------------------------------------------------------------
# Libraries

print("Loading additional packages for plotting.")
    
# Mapping
library(purrr)
# Basic shaping
library(dplyr)
# Pivotting
library(tidyr)
# Reading parameters, parsing factors
library(readr)
# Plotting
library(ggplot2)
# Correlation plots for individual chains
library(GGally)


# -------------------------------------------------------------------------
# Read basic chain settings

if (!("mcmcSamplerList" %in% class(out)) | (length(out) < 2)) {
  stop("Please provide an mcmcSamplerList with 2 or more chains for plotting.")
}
    
# Read sampler name
sampler <- out[[1]]$settings$sampler
    
# Read the amount of chains
n_chain <- length(out)

# Read the amount of iterations in a chain (should to be the same for each chain)
n_iter <- out[[1]]$settings$iterations

# Define starting value for sampling based on burn-in
if (sampler == "Metropolis"){
  start <- floor(burn_in) + 1
  # DE sampler subchains need division so sampling does not go out of bounds 
} else if (sampler %in% c("DEzs", "DREAMzs")) {
  start <- floor(burn_in / 3) + 1
}


# -------------------------------------------------------------------------
# Take a sample from each chain, combine results into a single tibble

df_samples <- out %>%
  # Take a sample from each chain
  map(
    getSample, 
    parametersOnly = F, start = start, thin = 1, numSamples = n_samples
  ) %>%
  # Convert to tibble
  map(as_tibble) %>%
  # Add an arbitrary iteration counter to each tibble (for trace plot x-axis)
  map(
    mutate,
    iteration = seq(from = floor(burn_in) + 1, to = n_iter, length.out = nrow(.[[1]]))
  ) %>%
  # Name each tibble with a chain number, rowbind tibbles
  set_names(nm = 1:n_chain) %>%
  map_dfr(as_tibble, .id = "chain")


# -------------------------------------------------------------------------
# Convert to long format for plotting

df_samples_long <- df_samples %>%
  # Select relevant columns
  select(
    -Lposterior,
    -Lprior,
    likelihood = Llikelihood
  ) %>%
  pivot_longer(
    cols = c(-chain, -iteration),
    names_to = "name",
    values_to = "value"
  ) %>%
  # Convert to factor to preserve order in ggplot, make x name more descriptive
  mutate(
    chain = parse_factor(chain),
    name = parse_factor(name)
  )


# -------------------------------------------------------------------------
# Import the priors (default parameters) used in the run to draw the parameter
# limits to plots.
# DEPRECATED: commented here and in the plots below. Can be restored, if 
# functionality is considered useful.

# Use readr, since it gets the factors levels correct by default.
# default_par <-
#   read_csv(
#     file = paste0(path_results_run, "/priors_", run_name, ".txt"),
#     col_types = "fddd"
#   ) %>% 
#   filter(name %in% unique(df_samples_long$name))


# -------------------------------------------------------------------------
# Plots over all chains

print("Drawing plots.")

# Individual chain densities for each parameter
facet_density <- df_samples_long %>% 
  ggplot() +
  geom_density(aes(x = value, col = chain)) +
  facet_wrap(vars(name), scales = "free") +
  theme(axis.text.x = element_text(angle = 17))
  # Draw parameter limits (DEPRECATED)
  # geom_vline(aes(xintercept = low), data = default_par, color = "blue") +
  # geom_vline(aes(xintercept = high), data = default_par, color = "blue")

# Overall density for each parameter
total_density <- df_samples_long %>% 
  ggplot() +
  geom_density(aes(x = value)) +
  facet_wrap(vars(name), scales = "free") +
  labs(title = "Overall density of all chains") +
  theme(axis.text.x = element_text(angle = 17))
  # geom_vline(aes(xintercept = low), data = default_par, color = "blue") +
  # geom_vline(aes(xintercept = high), data = default_par, color = "blue")

# Individual chain trace plots for each parameter
facet_trace <- df_samples_long %>%
  ggplot() +
  geom_line(aes(x = iteration, y = value, col = chain)) +
  facet_wrap(vars(name), scales = "free") +
  theme(axis.text.x = element_text(angle = 17))
  # geom_hline(aes(yintercept = low), data = default_par, color = "blue") +
  # geom_hline(aes(yintercept = high), data = default_par, color = "blue")


# -------------------------------------------------------------------------
# Save plots

# List plots
plots <- list(
  facet_density,
  total_density,
  facet_trace
)

# Name plots
plot_names <- paste0(
  c(
    "chain_density",
    "total_density",
    "chain_trace"
    # "chain_trace_all"
  ),
  "_", run_name, ".png"
)

# Save plots
walk2(
  plot_names, plots,
  ggsave,
  path = path_results_run,
  width = 12,
  height = 6
)

print("Summary plots have been saved!")


# -------------------------------------------------------------------------
# Plot for individual chains (if specified)

if (plot_chains) {
  
  # Loop over chains
  for (chain_id in 1:n_chain) {
    
    # Make sure that a subdirectory exists for each chain
    chain_dir <- paste0(path_results_run, "/ch", chain_id)
    if (!dir.exists(chain_dir)) { 
      dir.create(chain_dir) 
    }
    
    print(paste0("Plotting chain ", chain_id, " out of ", n_chain, "."))
    
    # Plot correlations
    corr_chain <- df_samples %>%
      filter(chain == chain_id) %>%
      select(-chain, -Lposterior:-iteration) %>%
      ggcorr(label = FALSE)
    
    
    # Plot density for each parameter in the chain
    density_chain <- df_samples_long %>% 
      filter(chain == chain_id) %>% 
      ggplot() +
      geom_density(aes(x = value)) +
      facet_wrap(vars(name), scales = "free")
    
    # Add parameter limits to the density plots
    # density_chain_limits <- density_chain +
    #   geom_vline(aes(xintercept = low), data = default_par, color = "blue") +
    #   geom_vline(aes(xintercept = high), data = default_par, color = "blue")
    
    
    # Plot trace for each parameter in the chain
    trace_chain <- df_samples_long %>% 
      filter(chain == chain_id) %>% 
      ggplot() +
      geom_line(aes(x = iteration, y = value)) +
      facet_wrap(vars(name), scales = "free")
    
    # Add parameter limits to the trace plots
    # trace_chain_limits <- trace_chain +
    #   geom_hline(aes(yintercept = low), data = default_par, color = "blue") +
    #   geom_hline(aes(yintercept = high), data = default_par, color = "blue")
    
    
    # Save the plots
    
    # List plots
    plots <- list(
      corr_chain,
      density_chain,
      # density_chain_limits,
      trace_chain
      # trace_chain_limits
    )
    
    # Name plots
    plot_names <- paste0(
      c(
        "correlation",
        "density", 
        # "density_limits",
        "trace" 
        # "trace_limits"
      ),
      "_", run_name, "_ch", chain_id, ".png"
    )
    
    # Save plots
    walk2(
      plot_names, plots,
      ggsave,
      path = chain_dir,
      width = 12,
      height = 6
    )
    
  } # end for
  print("Individual chain plots have been saved!")
} # end if

} # end function
