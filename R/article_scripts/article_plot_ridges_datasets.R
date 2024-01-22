# This script creates the dataset-divided ridge plots for the calibration article

# The code currently works for DEzs and DREAMzs based runs.

# -------------------------------------------------------------------------
# Libraries

# Standard tidyverse libraries
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)
# Sampling from mcmc results
library(BayesianTools)
# Ridge plots
library(ggridges)


# -------------------------------------------------------------------------
# Define paths and names for the runs

# Define path to folder with the mcmc results (change)
path_results <- 
  "R/article_scripts/mcmc_out/"

# Define path for saving the plots (change)
path_plots_out <- 
  "R/article_scripts/plots/"

# Define the run to be plotted for each dataset
run_cidet <- "run-202"
run_lidet <- "run-203"
run_ed <- "run-201"
run_global <- "run-220"

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

# Number of samples to be taken (can be changed)
n_samples <- 10000

# Set seed for sampling reproducability
set.seed(666)


# -------------------------------------------------------------------------
# Import results for each run

# Create list for sampled data
sampled_list <- vector(mode = "list", length = 4L)
names(sampled_list) <- run_ids 

# We should load and sample one dataset at a time to avoid memory crashes. This
# can take a couple of minutes to complete, so run it once and then use the
# codes in sections below.
for (j in 1:length(run_paths)) {

  # Load the run results
  run_results <- readRDS(run_paths[j])

  # Sample
  sampled_results <- getSample(
    run_results, start = start, thin = 1, numSamples = n_samples
  )

  # Save samples to list
  sampled_list[[j]] <- sampled_results
}


# -------------------------------------------------------------------------
# Shape the data into ggplot format

df_plot <- sampled_list %>%
  # Shape into a single dataframe
  map_dfr(as_tibble, .id = "dataset") %>%
  # Shape into long format for ggplot
  pivot_longer(cols = -dataset) %>%
  # Convert chars to factors to preserve correct order in ggplot
  mutate(
    # Reversing levels here to get them right in ridge plot
    dataset = parse_factor(
      dataset,
      levels = rev(c("global", "ed", "lidet", "cidet"))
    ),
    # Set parameter levels manually to ensure the facet plots have correct order
    name = parse_factor(
      name,
      levels = c(
        "aA", "aW", "aE", "aN", "pWA", "pEW", "pWN", "w1", "w2", "w5", "b1",
        "b2", "bN1", "bN2", "bH1", "bH2", "g", "gN", "gH", "pH", "aH", "th1",
        "th2", "r"
      )
    )
  )

  df_plot$name <- as.character(df_plot$name)

  df_plot$name[df_plot$name == 'aA'] <- 'alpha[A]'
  df_plot$name[df_plot$name == 'aW'] <- 'alpha[W]'
  df_plot$name[df_plot$name == 'aE'] <- 'alpha[E]'
  df_plot$name[df_plot$name == 'aN'] <- 'alpha[N]'
  df_plot$name[df_plot$name == 'aH'] <- 'alpha[H]'
  df_plot$name[df_plot$name == 'b1'] <- 'beta[1]'  
  df_plot$name[df_plot$name == 'b2'] <- 'beta[2]'  
  df_plot$name[df_plot$name == 'bN1'] <- 'beta[N1]'
  df_plot$name[df_plot$name == 'bN2'] <- 'beta[N2]'
  df_plot$name[df_plot$name == 'bH1'] <- 'beta[H1]'
  df_plot$name[df_plot$name == 'bH2'] <- 'beta[H2]'
  df_plot$name[df_plot$name == 'g'] <- 'gamma'
  df_plot$name[df_plot$name == 'gN'] <- 'gamma[N]'
  df_plot$name[df_plot$name == 'gH'] <- 'gamma[H]'
  df_plot$name[df_plot$name == 'th1'] <- 'phi[1]'
  df_plot$name[df_plot$name == 'th2'] <- 'phi[2]'
  df_plot$name[df_plot$name == 'w1'] <- 'w[ED]'
  df_plot$name[df_plot$name == 'w2'] <- 'w[CIDET]'
  df_plot$name[df_plot$name == 'w5'] <- 'w[LIDET]'
  df_plot$name[df_plot$name == 'pH'] <- 'p[H]'
  df_plot$name[df_plot$name == 'pWA'] <- 'p[WA]'
  df_plot$name[df_plot$name == 'pWN'] <- 'p[WN]'
  df_plot$name[df_plot$name == 'pEW'] <- 'p[EW]'

  mylevels <- c('alpha[A]', 'alpha[W]', 'alpha[E]', 'alpha[N]', 'p[WA]', 'p[EW]', 'p[WN]', 'w[ED]', 'w[CIDET]', 'w[LIDET]', 'beta[1]',
                'beta[2]', 'beta[N1]', 'beta[N2]', 'beta[H1]', 'beta[H2]', 'gamma', 'gamma[N]', 'gamma[H]', 'p[H]', 'alpha[H]', 'phi[1]',
                'phi[2]', 'r')
  df_plot$name <- factor(df_plot$name,levels=mylevels)
    
# Now the data is in the correct shape for ggplot and can be used for most types
# of plots. I have written the code for ridge plots below.


# -------------------------------------------------------------------------
# Plot

# Basic tweaks for all ridge plots
  ridge_base <- ggplot(df_plot, aes(x = value, y = dataset, fill = dataset)) +
  # Tick labels closer to x-axis
  scale_y_discrete(expand = c(0, 0)) +
  # Thicken lines overlapping x-axis
  coord_cartesian(clip = "off") +
  # Getting the legend to the same order as y-labels takes some hacking
  scale_fill_cyclical(
    # Define colors
    values = c("#009E73", "#CC79A7", "#0072B2", "#D55E00"),
    name = "dataset",
    labels = c("cidet", "lidet", "ed", "global"),
    guide = guide_legend(reverse = TRUE)) +
  # Font size
  theme_grey(base_size = 12) +
  # Remove legend title, set text angle, make facet strips narrower
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 19),
    strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"))
  ) +
  # Setup facets according to parameters
  facet_wrap(vars(name), scales = "free", labeller=label_parsed)

# "Standard" ridge density plot
ridge_densities <- ridge_base +
  geom_density_ridges(alpha = 0.5)

# Ridge plot with stat_density, cutoffs according to data
ridge_stat_densities <- ridge_base +
  geom_density_ridges(
    aes(height = stat(density)),
    alpha = 0.5, stat = "density", trim = TRUE
  )
  
# Ridge histogram
ridge_hist <- ridge_base +
  geom_density_ridges(
    alpha = 0.5, stat = "binline", bins = 35, scale = 0.95, 
    draw_baseline = FALSE
  )


# -------------------------------------------------------------------------
# Save plots

# Sometimes saving/showing the plots gives random errors related to
# UseMethod("depth"). Just try again and it works properly.

# List plots (continue the list with any extra plots you make)
plots <- list(
  ridge_densities,
  ridge_stat_densities,
  ridge_hist
)

# Name plots (continue the names with any extra plots you make)
plot_names <- paste0(
  c(
    "ridge_density",
    "ridge_density_stat",
    "ridge_histogram"
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
