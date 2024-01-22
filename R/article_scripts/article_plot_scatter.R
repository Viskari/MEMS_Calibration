# This script creates the scatter plots for the calibration article

# The global calibration MAPs are used to run YASSO with test datasets and the
# difference between simulations and observations is plotted versus weather
# data.

# The code currently works for DEzs and DREAMzs based runs. It is dependent on
# reading the test datasets from `YASSO15-calibration/data/`.

# -------------------------------------------------------------------------
# Libraries

# Standard tidyverse libraries
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
# Annotating plots with regression coefficients the easy way
library(ggpubr) 

# Sampling from mcmc results
library(BayesianTools)
# # Yasso model (check that version is correct)
library(Ryassofortran)
if (packageVersion("Ryassofortran") != "0.4.0") {
  stop("Please install the latest version (v0.4.0) of Ryassofortran from github/YASSOmodel.")
}


# -------------------------------------------------------------------------
# Define paths and name for the run

# Full path to the project directory (change unless opening with project file)
path_main_dir <- ""

# Define path to folder with the mcmc results
path_results <- paste0(path_main_dir, "R/article_scripts/mcmc_out/")

# Define path for saving the plots
path_plots_out <- paste0(path_main_dir, "R/article_scripts/plots/")

# Define path to the data folder
path_data <- paste0(path_main_dir, "data/")

# Define the run to be plotted
run_name <- "run-220"

# Path to run results
path_run <- paste0(path_results, "out_", run_name, ".rds")


# -------------------------------------------------------------------------
# Recreate settings used for calibration

# Datasets used for calibration
datasets <- c("ed1", "ed2", "cidet", "lidet", "gss", "makinen")

# Test datasets used to validate the model
testsets <- c("ed1", "ed2", "cidet", "lidet","hob3")

# Calibrated parameters (correspond to rows in default_par, e.g. 1 = aA, 2 = aW)
par_chosen <- c(1:5, 9, 15, 17:18, 21:35)

# Number of iterations used in each chain of the run (do not change, unless
# changing runs)
n_iterations <- 1.5e6

# Burn-in used in each chain of the run
burn_in <- 0.9 * n_iterations

# Change the burn-in format, so that the sampler understands DE subchains
# correctly
start <- floor(burn_in / 3) + 1


# -------------------------------------------------------------------------
# Import user-made functions

# Functions for importing data
source(paste0(path_main_dir, "R/data_import.R"))

# Functions for calling yasso and calculating likelihood
source(paste0(path_main_dir, "R/production/likelihood.R"))


# -------------------------------------------------------------------------
# Import test datasets, default parameters and results

# Import testing data
list_data_test <- import_data(path_data, testsets, "test")

# Import default parameters used for the run (if the run changes, a new parameter file
# needs to be added)
default_par <- read.csv(
  file = paste0(path_results, "defaults_and_priors_", run_name, ".txt"),
  colClasses = c("character", rep("numeric", 3))
)

# Import run results
run_results <- readRDS(path_run)


# -------------------------------------------------------------------------
# Calculate MAPs, correct parameters

# Obtain trained parameter values from MAPs, respect burn-in
trained_maps <- MAP(run_results, start = start)$parametersMAP

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

# Create list for storing the difference between simulation results and
# observations for each testset
list_validation <- vector(mode = "list", length = length(testsets))
names(list_validation) <- testsets

# Loop through the test datasets
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
    
    # Sum the difference over rows and store for plotting
    list_validation[[ts]][["difference"]] <- rowSums(diff)
  }
  
  # CIDET, LIDET
  if (ts %in% c("cidet", "lidet")) {
    
    # Run trained model and calculate differences to observations
    
    # Run YASSO15 with the trained MAPs and the testset data
    simulated <- call_yasso(trained_par, list_data_test[[ts]])
    
    # Sum the results over AWENH, since observations are summed 
    simulated <- rowSums(simulated)
    
    # Difference between simulated and observed results
    diff <- simulated - list_data_test[[ts]][["obs"]]
    
    # Store the difference for plotting
    list_validation[[ts]][["difference"]] <- diff
  }
  
  # Store weather data for plotting
#  list_validation[[ts]][["temp_mean"]] <- rowMeans(list_data_test[[ts]][["temp"]])
#  list_validation[[ts]][["temp_sd"]] <- apply(list_data_test[[ts]][["temp"]], 1, sd)
#  list_validation[[ts]][["precipitation"]] <- list_data_test[[ts]][["prec"]]
  
  list_validation[[ts]][["a"]] <- rowMeans(list_data_test[[ts]][["temp"]])
  list_validation[[ts]][["b"]] <- apply(list_data_test[[ts]][["temp"]], 1, sd)
  list_validation[[ts]][["c"]] <- list_data_test[[ts]][["prec"]]
}


# -------------------------------------------------------------------------
# Shape data for plotting, create plots

# Shape list into a tibble, combine ed1 and ed2
df_plot <- map_dfr(list_validation, as_tibble, .id = "dataset") %>% 
  mutate(
    dataset = if_else(dataset %in% c("ed1", "ed2"), "ed", dataset),
    dataset = parse_factor(dataset)
  ) %>% 
#  pivot_longer(cols = temp_mean:precipitation)
  pivot_longer(cols = a:c)

# Create basis for all plots
base_plot <- ggplot(df_plot, aes(x = value, y = difference)) +
  theme_grey(base_size = 12) +
  geom_point(aes(col = dataset, shape = dataset)) +
  geom_hline(yintercept = 0, size = 1.0) +
  geom_smooth(method = "lm", linetype = "dashed", col = "black", se = TRUE) +
  labs(y = "simulated - observed", x='')

# Create wide facet plot
scatter_facet_wide <- base_plot +
  stat_regline_equation(label.x.npc = 0.23, label.y.npc = 0.96) +
  facet_wrap(vars(name), scales = "free", dir = "v")

# Create wide facet plot
scatter_facet_long <- base_plot +
  stat_regline_equation(label.x.npc = 0.68) +
  facet_wrap(vars(name), scales = "free_x")


# -------------------------------------------------------------------------
# Save plots


# List plots (continue the list with any extra plots you make)
plots <- list(
  scatter_facet_wide,
  scatter_facet_long
)

# Name plots (continue the names with any extra plots you make)
plot_names <- paste0(
  c(
    "scatter_facet_wide",
    "scatter_facet_long"
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
