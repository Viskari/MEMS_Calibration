# This script creates ensemble plots (posterior predictive simulation) for the
# calibration article

# The global calibration results are compared against the individual calibration
# results with each of the test datasets. A group of parameter sets is sampled
# from the results of each calibration and all of these sets are used to run
# YASSO. The model results are passed through the error model (normal
# distribution). Lastly, the simulated mean is calculated and plotted with
# predictive intervals.

# The code works for DEzs and DREAMzs based runs (AM would require adding an
# if-condition for burn-in). It is dependent on reading the test datasets from
# `YASSO15-calibration/data/` and the import and YASSO calling codes from
# `YASSO15-calibration/R/`. The data and likelihood need to be the exact same as
# in the calibrations whose ensembles are plotted.


# -------------------------------------------------------------------------
# Libraries

# Standard tidyverse libraries
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)

# Sampling from mcmc results
library(BayesianTools)
# # Yasso model (check that version is correct)
library(Ryassofortran)
if (packageVersion("Ryassofortran") != "0.4.0") {
  stop("Please install the latest version (v0.4.0) of Ryassofortran from github/YASSOmodel.")
}


# -------------------------------------------------------------------------
# Define paths and names for the runs

# Full path to the project directory (change unless opening with project file)
path_main_dir <- ""

# Define path to folder with the mcmc results
path_results <- paste0(path_main_dir, "R/article_scripts/mcmc_out/")

# Define path for saving the plots
path_plots_out <- paste0(path_main_dir, "R/article_scripts/plots/")

# Define path to the data folder
path_data <- paste0(path_main_dir, "data/")

# Define the run to be plotted for each dataset 
# NOTE: If you change the runs, you should also change the settings in the next
# section and check that the the default parameter file, the yasso-calling
# function and the data are match the changed run.
run_cidet <- "run-202"
run_lidet <- "run-203"
run_ed <- "run-201"
run_global <- "run-212"

# Name the runs and create paths to them
run_ids <- c("cidet", "lidet", "ed", "global")
run_names <- c(run_cidet, run_lidet, run_ed, run_global)
run_paths <- paste0(path_results, "out_", run_names, ".rds")


# -------------------------------------------------------------------------
# Recreate settings used in calibration

# Test datasets used to validate the model
testsets <- c("cidet", "lidet", "ed1")

# Calibrated parameters for each run
par_cidet <- c(1:5, 9, 15, 18, 22:32)
par_lidet <- c(1:5, 9, 15, 21:32)
par_ed <- c(1:5, 9, 15, 17, 22:32)
par_global <- c(1:5, 9, 15, 17:18, 21:35)
# Store the parameters into a list
par_list <- list(par_cidet, par_lidet, par_ed, par_global)
names(par_list) <- run_ids

# Number of iterations used in each chain of the run (do not change, unless
# changing runs)
n_iterations <- 1.5e6

# Burn-in used in each chain of the run
burn_in <- 0.9 * n_iterations

# Change the burn-in format, so that the sampler understands DE subchains
# correctly
start <- floor(burn_in / 3) + 1

# Number of samples to be taken (can be changed)
n_samples <- 1000

# Set seed for sampling reproducability
set.seed(666)


# -------------------------------------------------------------------------
# Import user-made functions

# Functions for importing data
source(paste0(path_main_dir, "R/data_import.R"))

# Function for calling yasso (ensure that it matches the one used in the run)
source(paste0(path_main_dir, "R/production/likelihood.R"))


# -------------------------------------------------------------------------
# Import test datasets and default parameters

# Import testing data
list_data_test <- import_data(path_data, testsets, "test")

# Import default parameters used for the run (if the run changes, a new
# parameter file might need to be added)
default_par <- read.csv(
  file = paste0(path_results, "defaults_and_priors_", run_global, ".txt"),
  colClasses = c("character", rep("numeric", 3))
)


# -------------------------------------------------------------------------
# Import results of each calibration

# Create list for sampled data
list_sampled <- vector(mode = "list", length = length(run_ids))
names(list_sampled) <- run_ids 

# We should load and sample one run at a time to avoid memory crashes. This
# can take a couple of minutes to complete, so do it once and then use the
# codes in sections below.
for (j in 1:length(run_paths)) {
  
  # Load the run results
  run_results <- readRDS(run_paths[j])
  
  # Sample
  sampled_results <- getSample(
    run_results, start = start, thin = 1, numSamples = n_samples
  )
  
  # Save samples to list
  list_sampled[[j]] <- sampled_results
}


# -------------------------------------------------------------------------
# Create full parameter sets and do humus corrections

# Create list for the corrected full parameter sets
list_corrected <- vector(mode = "list", length = length(run_ids))
names(list_corrected) <- run_ids

# Loop over the calibrations
for (run in run_ids) {
  
  # Create a matrix with the default parameter values as rows
  trained_par <- matrix(
    default_par$value, nrow = nrow(list_sampled[[run]]), ncol = nrow(default_par), 
    byrow = TRUE, dimnames = list(NULL, default_par$name)
  )
  
  # Override the parameters that have calibrated values
  trained_par[, par_list[[run]]] <- list_sampled[[run]]
  
  # Correct the humus flow: there is no flow to H pool in the litter bag
  # measurements (that we are using for validation), hence pH is set to zero
  trained_par[, 31] <- 0
  
  # Correct the other p-parameters accordingly
  trained_par[, 6] <- 1 - trained_par[, 9] # pEA = 1 - pEW - pH
  trained_par[, 7] <- 1 # pNA = 1 - pH
  trained_par[, 8] <- 1 # pAW = 1 - pH
  
  # Store the corrected full parameter sets to list
  list_corrected[[run]] <- trained_par
}


# -------------------------------------------------------------------------
# Do the ensemble runs, i.e. run YASSO with all the corrected parameter sets

# This might be the most horrible mess of nested for-loops I have
# ever created but it works. I tried a mapped solution but that was more
# challenging to implement and to read.

# The EURODECO N + H pool logic present in likelihood is not added here, since
# it is assumed that the simulated H-pool will be zero due to pH = 1. This makes
# the code more readable. However, changing the pH to something else will
# require adding an "if" for ed pools as in likelihood

# TODO: There is for sure a more efficient data structure than the existing one.
# I considered this a "use once" script so I did not spend much time optimizing
# it. If the code will be commonly used, overhaul the data structure.

# Setup a list for storing the ensemble results of each testset 
list_ensembles <- vector("list", length = length(testsets))
names(list_ensembles) <- testsets

# Loop over testsets (this can take 30 s, all the YASSO runs are done within)
for (ts in testsets) {
  
  # Loop over the calibrations
  for (run in run_ids) {
    
    # Only produce results for the given combinations
    if (
      (run == "global") | (run == "cidet" & ts == "cidet") | 
      (run == "lidet" & ts == "lidet") | (run == "ed" & ts %in% c("ed1", "ed2"))
      ) {
      
        # Setup a sublist for ensemble results of current calibration
        list_ensemble_res <- vector("list", length = nrow(list_corrected[[run]]))
        
        # Loop over sampled parameter sets in current calibration (ugly for-loop
        # over rows in a matrix, forced by leaching)
        for (j in 1:nrow(list_corrected[[run]])) {
          
          # Correct leaching before each run
          if (ts == "cidet") {
            list_data_test[[ts]][["leac"]] <- as.double(list_corrected[[run]][j, 18])
          }
          if (ts == "lidet") {
            list_data_test[[ts]][["leac"]] <- as.double(list_corrected[[run]][j, 21])
          }
          if (ts %in% c("ed1", "ed2")) {
            list_data_test[[ts]][["leac"]] <- as.double(list_corrected[[run]][j, 17])
          }
          
          # Run YASSO with a sampled parameter set and the testset data
          simulated <- call_yasso(list_corrected[[run]][j, ], list_data_test[[ts]])
  
          # Sum the results over AWENH, since observations are summed 
          simulated <- rowSums(simulated)
          
          # Shape observations and uncertainties depending on testset, sum over
          # ED pools
          if (ts %in% c("ed1", "ed2")) {
            observed <- rowSums(list_data_test[[ts]][["obs"]])
            uncertainty <- rowSums(list_data_test[[ts]][["unc"]])
          } else {
            observed <- list_data_test[[ts]][["obs"]]
            uncertainty <- list_data_test[[ts]][["unc"]]
          }
          
          # Simulate uncertainties for predictive intervals
          sim_predictive <- rnorm(
            n = length(simulated), 
            mean = simulated, 
            sd = uncertainty
          )
          
          # Store results for current sampled set
          list_ensemble_res[[j]] <- tibble(
            testset = ts,
            calibration = run,
            step = 1:length(simulated), 
            simulated,
            observed,
            sim_predictive
          )
        } # end for over sampled parameter sets
        
        # Bind results of all sampled sets of current run, store to list
        list_ensembles[[ts]][[run]] <- map_dfr(
          list_ensemble_res, as_tibble, .id = "sample"
        )
        
    } # end if (given combinations)
  } # end for over calibrations
} # end for over datasets


# -------------------------------------------------------------------------
# Shape and plot the ensemble results

# Combine the global & individual ensemble results in each testset
combined_ensemble <- map(list_ensembles, bind_rows)

# If ed1 and ed2 need to be combined, it can be done here
# TODO: write logic to combine results, observations
# x <- bind_rows(combined_ensemble[c("ed1", "ed2")])

# Loop over the ensemble results of each testset
# (Mapping this would give cleaner environment, but for-loop is quicker to
# implement and read)
for (setname in names(list_ensembles)) {

  # Choose ensemble results for testset
  ensemble <- combined_ensemble[[setname]]
  
  # Store the observations before summarizing
  obs <- ensemble %>% 
    select(step, observed) %>% 
    distinct()
  
  # Summarize to compute the quantiles
  df_plot <- ensemble %>% 
    # Factorize to preserve right order through summarize
    mutate(calibration = parse_factor(calibration)) %>% 
    group_by(step, calibration) %>% 
    # Calculate mean and predictive intervals
    summarise(
      mean = mean(simulated),
      `2.5%` = quantile(sim_predictive, 0.025),
      `97.5%` = quantile(sim_predictive, 0.975)
    ) %>% 
    pivot_longer(cols = c(mean, `2.5%`, `97.5%`), names_to = "measure") %>% 
    # Factorize to get right order in plots
    mutate(measure = parse_factor(measure)) %>% 
    # Join in the observations
    left_join(obs, by = "step")
  
  # Plot the results, save the plots
  
  # Basis for plots
  ensemble_base <- ggplot(df_plot, aes(x = step)) +
    geom_point(aes(y = observed), alpha = 0.1) +
    scale_linetype_manual(
      "measure",
      values = c("97.5%" = 1, "2.5%" = 1, "mean" = 3)
    ) +
    # Order the legends similarly in all plots
    guides(
      colour = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    labs(y = "carbon")
  
  # Plot with raw lines
  ensemble_lines <- ensemble_base +
    geom_line(
      aes(y = value, linetype = measure, col = calibration)
    )
  
  # Plot with smoothing
  ensemble_smooth <- ensemble_base +
    geom_smooth(
      aes(y = value, linetype = measure, col = calibration),
      se = FALSE,
      span = 0.18, # this control the "smoothness" of the approximation
      alpha = 0.1
    )
    
  # List plots
  plots <- list(
    ensemble_lines,
    ensemble_smooth
  )
  
  # Name plots
  plot_names <- paste0(
    c(
      paste0("ensemble_", setname, "_lines"),
      paste0("ensemble_", setname)
    ),
    ".png"
  )
  
  # Save plots
  walk2(
    plot_names, plots,
    ggsave,
    path = path_plots_out,
    width = 12,
    height = 6
  )
  
} # end for over ensemble results