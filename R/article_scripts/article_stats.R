# This script creates the stats for the calibration article

# The code currently works for DEzs and DREAMzs based runs.

# -------------------------------------------------------------------------
# Libraries

# Standard tidyverse libraries
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
# Sampling from mcmc results
library(BayesianTools)


# -------------------------------------------------------------------------
# Define paths and names for the runs

# Define path to folder with the mcmc results (change)
path_results <- 
  "R/article_scripts/mcmc_out/"

# Define path for saving the statistics (change)
path_stats_out <- 
  "R/article_scripts/stats/"

# Define the run for each dataset
run_cidet <- "run-207"
run_lidet <- "run-208"
run_ed <- "run-206"
run_global <- "run-212"

# Name the runs and create paths to them
run_ids <- c("cidet", "lidet", "ed", "global")
run_names <- c(run_cidet, run_lidet, run_ed, run_global)
run_paths <- paste0(path_results, "out_", run_names, ".rds")


# -------------------------------------------------------------------------
# Recreate settings used for calibration

# Number of iterations used in each chain of the runs (do not change, unless
# changing runs)
n_iterations <- 1.5e6

# Burn-in used in each chain of the runs (can be changed)
burn_in <- 0.9 * n_iterations

# Change the burn-in format, so that the sampler understands DE subchains
# correctly
start <- floor(burn_in / 3) + 1

# Number of samples to be taken. I am not defining this for now, but instead
# using all chain values after burn-in.
# n_samples <- 10000


# -------------------------------------------------------------------------
# Import results for each run

# Create list for sampled data
sampled_list <- vector(mode = "list", length = 4L)
names(sampled_list) <- run_ids 

# We should load and sample one dataset at a time to avoid memory crashes. This
# can take a couple of minutes to complete, so run it once and then use the
# codes in sections below.
for (j in 1:length(run_names)) {
  
  # Load the run results
  run_results <- readRDS(run_paths[j])
  
  # Sample 
  sampled_results <- getSample(run_results, start = start)
  
  # Calculate MAPs with BTools function
  map_run <- MAP(run_results, start = start)$parametersMAP
  
  # Save to list
  sampled_list[[j]][["samples"]] <- sampled_results
  sampled_list[[j]][["maps"]] <- map_run
}


# -------------------------------------------------------------------------
# Calculate statistics

# Initialize list for stats of each dataset
stats_dataset <- vector(mode = "list", length = 4L)
names(stats_dataset) <- run_ids

# Loop over datasets
for (j in 1:length(run_names)) {

  # Create df with medians, CIs and MAPs
  stats_dataset[[j]] <- sampled_list[[j]][["samples"]] %>%
    as_tibble() %>%
    # To long format to calculate stats easily
    pivot_longer(cols = everything()) %>%
    # Preserves parameter order after summarising
    mutate(name = parse_factor(name)) %>% 
    # Group by parameters
    group_by(name) %>%
    # Calculate stats
    summarise(
      median = median(value),
      `97.5 %` = quantile(value, probs = 0.975),
      `2.5 %` = quantile(value, probs = 0.025)
    ) %>%
    # Add MAP and label datasets
    mutate(
      MAP = sampled_list[[j]][["maps"]],
      dataset = run_ids[j]
    ) %>% 
    # Rearrange columns
    relocate(dataset, name, MAP, `2.5 %`, median, `97.5 %`)
}

# This could all be done with map(), but for-loop is easier to read


# -------------------------------------------------------------------------
# Combine statistics

# Combine all the datasets into a single file that can be saved
stats_all <- stats_dataset %>% 
  map_dfr(as_tibble)
  

# -------------------------------------------------------------------------
# Save statistics

# Save the file with statistics for all datasets
write_csv(stats_all, path = paste0(path_stats_out, "stats_calibrations.csv"))
