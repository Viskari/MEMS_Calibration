# Wrapper for calibrating YASSO with multiple parallel chains.

# Each chain is run at a different CPU core using the parallel package. The
# calibration results are plotted and summarized to "results/<run_name>/". The 
# raw results are saved to "results/out/". The shell script "copy_wrappers.sh" 
# can be used to make copies of this wrapper and "run_wrappers.sh" can be used 
# to run the wrappers simultaneously.

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# User-defined settings

## Reproducability and paths

# Seed and run name (matched for reproducability). The run name dictates the 
# folder where the results are saved, so make sure you do not overwrite existing 
# results. In other words, use a unique seed for each calibration.

seed <- 100L
run_name <- paste0("run-", seed)

# Full path to the project directory
# path_main_dir <- "" # local, with project file
#path_main_dir <- "~/YASSO-calibration-master/" # server
path_main_dir <- "/eos/jeodpp/data/projects/SOIL-NACA/MEMS/Calibration/" # server


## MCMC settings

# Datasets used for calibration
datasets <- c("LUCAS")

# Calibrated parameters (correspond to rows in default_par, e.g. 1 = aA, 2 = aW)
#par_chosen <- c(1:5, 9, 15, 17:18, 21:35)
par_chosen <- c(23:25)

# Mcmc sampler (AM, DEzs or DREAMzs)
sampler_name <- "DEzs"

# Iterations and chains for calibration, burn-in for postprocessing
# n_iter <- 1.5e6
n_iter <- 1.5e2
n_chains <- 3
burn_in <- 0.9 * n_iter

# Use the full ("full") or training ("train") data in the calibration. If 
# training data was chosen, define datasets used to validate the trained model.
# datatype <- "train"
datatype <- "full"
#testsets <- c("ed1", "ed2", "cidet", "lidet","hob3")
#testsets <- c("hob3")

# Run the chains in parallel (TRUE) or sequentially (FALSE). Parallelization is
# highly recommended, since it saves a considerable amount of time and works
# even on low-end systems. Sequential runs should only be done in the rare
# cases, where the system does not support parallelization.
# parallel <- TRUE
parallel <- FALSE


## Post-processing the results: statistics

# Calculate and save statistics for the results. If TRUE (recommended), stats
# will be saved to "results/<run-name>".
do_stats <- TRUE


## Post-processing the results: plotting
# Required packages: tidyverse, GGally

# Plot results. If TRUE, results will be plotted to "results/<run-name>".
plot_results <- TRUE

# Plot individual chains. If TRUE, creates correlation, trace and density plots
# for each chain to "results/<run-name>/<chain>".
plot_chains <- TRUE


# -------------------------------------------------------------------------
# Additional paths

# Full path to the folder to save mcmc output files in
path_results_out <- paste0(path_main_dir, "results/out")

# Full path to the folder to save plots and statistics in
path_results_run <- paste0(path_main_dir, "results/", run_name)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Start of calibration logic

# -------------------------------------------------------------------------
# Libraries

# MCMC algorithms
 library(BayesianTools)
# MEMS model wrapper
source(paste0(path_main_dir,"Model/mems_calibration/calibration_MEMS.R"))

# Parallelization
library(parallel)

# Functions for importing data
source(paste0(path_main_dir, "R/data_import_MEMS.R"))
# Functions for calling yasso and calculating likelihood
source(paste0(path_main_dir, "R/production/likelihood_MEMS.R"))


# -------------------------------------------------------------------------
# Create folder structure for saving results and plots

# Create folders for saving the calibration results if they do not exist
if (!dir.exists(path_results_out)) { dir.create(path_results_out) }
if (!dir.exists(path_results_run)) { dir.create(path_results_run) }
if (!dir.exists(paste0(path_results_out, "/tmp"))) {
  dir.create(paste0(path_results_out, "/tmp"))
}

# Create folders for saving plots and statistics, if they do not exist and are
# required
if (do_stats | plot_results) {
  # Create a subfolder for each chain, if one does not exist
  chain_dir <- paste0(path_results_run, "/ch", 1:n_chains)
  invisible(lapply(chain_dir, function(dir) {if (!dir.exists(dir)) {dir.create(dir)}}))
}

# Create folder for saving validation results
if ((datatype == "train") & !dir.exists(paste0(path_results_run, "/validation"))) {
  dir.create(paste0(path_results_run, "/validation"))
}

# Create tmp folder for saving results of each chain to protect against crashes 
if (!dir.exists(paste0(path_results_out, "/tmp"))) {
  dir.create(paste0(path_results_out, "/tmp"))
}


# -------------------------------------------------------------------------
# Import data and parameters

# Import the specified dataset(s)
 list_data <- import_data(paste0(path_main_dir, "data/"), datasets, datatype)


# Import default parameter values and prior ranges
default_par <- read.csv(
  file = paste0(path_main_dir, "parameters/defaults_and_priors_MEMS.csv"),
  colClasses = c("character", rep("numeric", 3))
)

# Import starting values for the chains
inits <- as.matrix(read.csv(
  file = paste0(path_main_dir, "parameters/starting_values_MEMS.csv"))
)

# Choose the starting values of the parameters in current calibration
inits <- inits[, default_par$name[par_chosen]]

# Loading data for the MEMS to use
dirSL<-("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/SptLYR/")
lat_long<-as.array(stack(paste0(dirSL,"lat_long.tif")))/100


# -------------------------------------------------------------------------
# Setup priors, define BayesianTools settings

# Create uniform priors for each parameter
prior <- createUniformPrior(
  lower = default_par$low[par_chosen],
  upper = default_par$high[par_chosen]
)

# Create Bayesian setup
b_setup <-
  createBayesianSetup(
    likelihood,
    prior,
    names = default_par$name[par_chosen],
    parallel = FALSE
  )


# -------------------------------------------------------------------------
# Setup sampler

# AM
if (sampler_name == "AM") {
  
  # Change the sampler name so BayesianTools understands it
  sampler_name <- "Metropolis"
  
  # Settings for the sampler
  sampler_settings <- list(
    iterations = n_iter,
    nrChains = 1,
    adapt = TRUE,
    # Turn off optimization so initial values can be set manually
    optimize = FALSE,
    # Setup inits for chain 1. The other chains' inits are set up later.
    startValue = inits[1, ],
    consoleUpdates = 5000
  )

# DEzs, DREAMzs
} else if (sampler_name %in% c("DEzs", "DREAMzs")) {
  
  # Settings for the sampler
  sampler_settings <- list(
    iterations = n_iter,
    nrChains = 1,
    startValue = inits,
    consoleUpdates = 5000
  )
} else {
  stop("Please provide a valid sampler name.")
}


# -------------------------------------------------------------------------
# Print settings to logs, check likelihood

# Basic settings for the calibration
print(paste0("Starting ", run_name, " with sampler " , sampler_name))
print("Used datasets:")
print(datasets)
print("Calibrated parameters:")
print(default_par[par_chosen, "name"])

# Calculate and print the likelihood with the default parameters
print(paste0(
  "Default likelihood for chosen datasets is ",
  round(likelihood(default_par$value[par_chosen]), 5)
  )
)


# -------------------------------------------------------------------------
# Perform the calibration

# If parallelization is used (recommended): create a cluster and run the chains
# in parallel
if (parallel == TRUE) {
  
  # Create a cluster
  cl <- makeCluster(n_chains, outfile = "") # windows + linux
  # cl <- makeCluster(n_chains, type = "FORK", outfile = "") # linux
  
  # Send libraries to the cluster
  invisible(clusterEvalQ(cl, library(BayesianTools)))
  invisible(clusterEvalQ(cl, library(Ryassofortran)))
  
  # Send variables and functions to the cluster for likelihood and saving
  clusterExport(
    cl = cl,
    varlist = c(
      "datasets", "default_par", "par_chosen", "list_data", "call_yasso", 
      "path_results_out", "run_name"),
    envir = environment()
  )
  
  # Set cluster seed
  clusterSetRNGStream(cl, seed)
  
  print("Starting parallelized MCMC.")
  
  # Run MCMC in parallel
  out_chains <- parLapply(
    cl = cl,
    X = 1:n_chains,
    fun = function(X, inits, b_setup, sampler_name, sampler_settings) {
      
      # Setup AM starting values, which are different for each chain
      if (sampler_name == "Metropolis") {
        # Currently, there are values for three chains in inits. In case of more
        # chains, use the first chain inits.
        if (X <= nrow(inits)) {
          sampler_settings[["startValue"]] <- inits[X, ]
        } else {
          sampler_settings[["startValue"]] <- inits[1, ]
        }
        print(paste0("Setting initials for chain ", X))
      }
      
      # Run the calibration
      ch_out <- runMCMC(
        bayesianSetup = b_setup,
        sampler = sampler_name,
        settings = sampler_settings
      )
      
      # Save the MCMC output temporarily within the loop for crash protection
      saveRDS(
        ch_out,
        file = paste0(path_results_out, "/tmp/out_", run_name, "_", X, ".rds")
      )
      
      # Return the MCMC output
      ch_out
      
    },
    inits,
    b_setup,
    sampler_name,
    sampler_settings
  )
  
  # Stop the cluster
  stopCluster(cl)
  
# If parallelization is not used (not recommended), run sequentially 
} else if (parallel == FALSE) {
  
  # Initialize a list for saving the calibration results
  out_chains <- vector(mode = "list", length = n_chains)
  
  # Set seed (might not work equally across OS's due to BayesianTools)
  set.seed(seed)
  
  # Perform sequential calibration
  for (X in 1:n_chains) {
    
    # Setup AM starting values, which are different for each chain
    if (sampler_name == "Metropolis") {
      # Currently, there are values for three chains in inits. In case of more
      # chains, use the first chain inits.
      if (X <= nrow(inits)) {
        sampler_settings[["startValue"]] <- inits[X, ]
      } else {
        sampler_settings[["startValue"]] <- inits[1, ]
      }
      print(paste0("Setting initials for chain ", X))
    }
    
    # Run the calibration
    ch_out <- runMCMC(
      bayesianSetup = b_setup,
      sampler = sampler_name,
      settings = sampler_settings
    )
    
    # Save the MCMC output temporarily within the loop for crash protection
    saveRDS(
      ch_out,
      file = paste0(path_results_out, "/tmp/out_", run_name, "_", X, ".rds")
    )
    
    # Save the MCMC output to list
    out_chains[[X]] <- ch_out
  }
}

# Combine results into a single mcmc object
out <- createMcmcSamplerList(out_chains)


# -------------------------------------------------------------------------
# Save results

# Save the calibration results
saveRDS(out, file = paste0(path_results_out, "/out_", run_name, ".rds"))

# Remove temporary files for the run, if the saving was successful
if (file.exists(paste0(path_results_out, "/out_", run_name, ".rds")) &
    dir.exists(paste0(path_results_out, "/tmp"))) {
  tmp_files <- list.files(
    paste0(path_results_out, "/tmp"),
    pattern = paste0("*_", run_name, "_*")
    )
  file.remove(paste0(path_results_out, "/tmp/", tmp_files))
}

# Save the defaults and prior ranges used for the run
write.csv(
  default_par,
  file = paste0(path_results_run, "/defaults_and_priors_", run_name, ".txt"),
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Optional: load old results

# Load old results, if stats and plots should be done done later/again
# out <- readRDS(file = paste0(path_results_out, "/out_", run_name, ".rds"))

# Load old results from out/tmp/, in case the combination phase crashed
# tmp_files <- list.files(
#   paste0(path_results_out, "/tmp"),
#   pattern = paste0("*_", run_name, "_*")
#   )
# out_chains <- lapply(paste0(path_results_out, "/tmp/", tmp_files), readRDS)
# out <- createMcmcSamplerList(out_chains)


# -------------------------------------------------------------------------
# Optional: Calculate and save statistics to "results/<run-name>".

# Calculate and save priors, initials and summaries
if (do_stats) {
  source(paste0(path_main_dir, "R/production/stat_results.R"))
  stat_mcmc_results(
    path_results_run, run_name, out, default_par, par_chosen,
    inits, likelihood, burn_in = burn_in
  )
}


# -------------------------------------------------------------------------
# Optional: Draw plots to "results/<run-name>"

if (plot_results) {
  # Load function for plotting results (requires additional libraries)
  source(paste0(path_main_dir, "R/production/plot_results.R"))
  # Draw and save plots
  plot_mcmc_results(
    path_results_run, path_results_out, run_name,
    out = out, plot_chains = plot_chains, burn_in = burn_in
  )
}


# -------------------------------------------------------------------------
# Optional: Test the trained model with various datasets.

if (datatype == "train") {
  source(paste0(path_main_dir, "R/production/validate_results.R"))
  validate_results(path_main_dir, run_name, datasets, testsets, out, default_par, par_chosen, burn_in = burn_in)
}


# -------------------------------------------------------------------------
# Save session info for reproducability (at the end, so all the attached
# packages are visible)

capture.output(
  sessionInfo(),
  file = paste0(path_results_run, "/session_info_", run_name, ".txt")
)

print("All done!")
