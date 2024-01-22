# Simple example for running a BayesianTools Adaptive Metropolis calibration for
# YASSO15. Also useful for debugging the likelihood logic.

# More information about BayesianTools can be found on their github page:
# https://github.com/florianhartig/BayesianTools/tree/master/BayesianTools 
# as well as the vignette:
# https://cran.r-project.org/web/packages/BayesianTools/vignettes/BayesianTools.html


# -------------------------------------------------------------------------
# Set paths

# Set generic path (empty if opened with project file)
path_main <- ""
# Set data path
path_data <- paste0(path_main, "data/")


# -------------------------------------------------------------------------
# Libraries

# MCMC algorithms
library(BayesianTools)
# Yasso model (check that version is correct)
library(Ryassofortran)
if (packageVersion("Ryassofortran") != "0.4.0") {
  stop("Please install the latest version (v0.4.0) of Ryassofortran from github/YASSOmodel.")
}

# Functions for importing data
source(paste0(path_main, "R/data_import.R"))

# Functions for calling yasso and calculating likelihood
source(paste0(path_main, "R/example/likelihood_example.R"))


# -------------------------------------------------------------------------
# Import data and parameters

# Choose datasets used for calibration
# datasets <- c("ed1", "ed2", "cidet", "lidet", "gss", "tarasov", "makinen", "hob3")
datasets <- "cidet"

# Choose which parameters to calibrate. Rest will be kept at default values.
par_chosen <- c(1:5, 9, 15, 18, 22:32)

# Import the specified dataset(s)
list_data <- import_data(path_data, datasets)

# Import default parameter values and prior ranges
default_par <- read.csv(
  file = paste0(path_main, "parameters/defaults_and_priors.csv"),
  colClasses = c("character", rep("numeric", 3))
)

# Import starting values for the chains
inits <- as.matrix(read.csv(
  file = paste0(path_main, "parameters/starting_values.csv"))
)

# Use only the inits that are needed in current run
inits <- inits[, default_par$name[par_chosen]]


# -------------------------------------------------------------------------
# Setup priors, define BayesianTools settings

# Create uniform priors for each parameter
prior <- createUniformPrior(
  lower = default_par$low[par_chosen],
  upper = default_par$high[par_chosen],
  best = default_par$value[par_chosen]
)

# Create a Bayesian setup (basically define likelihood and prior for mcmc)
b_setup <-
  createBayesianSetup(
    likelihood = likelihood,
    prior = prior,
    names = default_par$name[par_chosen],
    parallel = FALSE
  )


# -------------------------------------------------------------------------
# Setup sampler

# Iterations and chains
n_iter <- 1e4
n_chains <- 1

# Sampler name
sampler_name <- "Metropolis"

# Get initial values for the parameters. Setting the initial values manually is
# optional, but a good practice for reproducability. Here, we will pick one row
# from inits, since we are using AM. Other samplers require all three rows.
inits <- inits[1, ]

# Settings for the AM sampler (see BayesianTools vignette for explanation)
sampler_settings <- list(
  iterations = n_iter,
  nrChains = n_chains,
  adapt = TRUE,
  # Turn off optimization so initial values can be set manually
  optimize = FALSE,
  startValue = inits
)


# -------------------------------------------------------------------------
# Setup global variables for debugging likelihood

# If you are going through this file to learn about calibration, you need to run
# this section but you do not need to understand it - it is not a part of
# BayesianTools logic. Here, we set up various variables for debugging purposes.

# Counter for recording global variables
passed <- 0

# Save likelihood and cost function from each likelihood() access
lh_cf <- matrix(
  nrow = n_iter + 50, ncol = 2,
  dimnames = list(NULL, c("l_hood", "c_func"))
)

# Save all proposed parameter sets (calibration only saves accepted)
all_proposed <- matrix(
  nrow = n_iter + 50, ncol = length(par_chosen),
  dimnames = list(NULL, default_par$name[par_chosen])
)

# Track amount of rejections in p-checks
pA_fail <- 0
pW_fail <- 0
pE_fail <- 0
pN_fail <- 0


# -------------------------------------------------------------------------
# Likelihood debug opportunity

# Use this likelihood call to debug the likelihood logic
likelihood(default_par$value[par_chosen])


# -------------------------------------------------------------------------
# Run MCMC

# Set seed for reproducability
set.seed(150)

# Run the calibration
out <- runMCMC(
  bayesianSetup = b_setup,
  sampler = sampler_name,
  settings = sampler_settings
)


# -------------------------------------------------------------------------
# Examples for working with the results

# Calculate summary for the results
summary(out)

# Draw trace and density plot for each parameter
# tracePlot(out, numSamples = 5000)

# Get a sample from the chain(s)
# xx <- getSample(out, numSamples = 4000, parametersOnly = FALSE)



# -------------------------------------------------------------------------
# Likelihood sanity check

# If you are going through this file to learn about calibration, you do not need
# to read this section - it is not a part of BayesianTools logic. Here, we
# basically ensure that the contents of the datasets and the likelihood are as
# expected by computing a likelihood for every dataset with the default
# parameter values.

# This section is tricky to hide into a function, since BayesianTools demands
# certain variables in likelihood to be loaded from the global environment.

# A check comparing the MAPs of a known calibration was implemented here
# earlier, but it did not work across all operating systems. It seems like
# BayesianTools Metropolis algorithm does not respect set.seed() across OS's.


# Store datasets and data used in calibration (both variables will be used in
# sanity check)
datasets_tmp <- datasets
list_data_tmp <- list_data

# Define a known-good likelihood value for every dataset. If you make changes to
# the datasets, the default variables or the likelihood, you might have to
# refresh these values.
ds_lhoods <- c(
  "ed1" = -4397.2, 
  "ed2" = -8920.2, 
  "cidet" = -7016.0, 
  "lidet" = -15463.6,
  "gss" = -12858.5,
  "tarasov" = -628.4, 
  "makinen" = -8996.1, 
  "hob3" = -2388.0
)
# Setup a vector for storing comparison truth values
ds_ok <- vector(length = length(ds_lhoods))
names(ds_ok) <- names(ds_lhoods)

# Set tolerance for checking equality
tolerance <- 0.1

cat("\nLikelihood sanity check:\n")

# Perform sanity check: Loop over the datasets
for (datasets in names(ds_lhoods)) {
  
  # Import data for the dataset
  list_data <- import_data(path_data, datasets)
  # Compute likelihood with default parameters for the dataset
  l_hood <- likelihood(default_par$value[par_chosen])
  
  # Compare computed likelihood to known-good values
  # If comparison is good, produce a message
  if (abs(l_hood - ds_lhoods[datasets]) < tolerance) {
    print(
      paste0("PASS ", datasets, " check with l.hood ", round(l_hood, 2))
    )
    # If comparison is ok, mark it as TRUE
    ds_ok[datasets] <- TRUE
  # If the computed likelihood is off, produce a message and a warning
  } else {
    print(
      paste0(
        "FAIL ", datasets, " check with l.hood ", round(l_hood, 2),
        " vs expected ", ds_lhoods[datasets]
      )
    )
    warning(paste0(datasets, " failed the likelihood sanity check!"))
  }
}

# Restore datasets and data used in calibration
datasets <- datasets_tmp
list_data <- list_data_tmp

# Let the user know if everything is good
if (exists("out") & all(ds_ok)) {
  cat("\nCalibration is working as intended!")
} else {
  warning("Calibration is not working as intended! Check for any warnings.")
}
