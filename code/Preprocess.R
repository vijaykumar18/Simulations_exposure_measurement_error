# Load necessary libraries
library(stringr)

# Define the base directory
base_dir <- "~/Downloads/DiffMeasError/exp"

# Function to find the latest experiment number
find_latest_experiment <- function(base_dir) {
  experiment_dirs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
  experiment_numbers <- as.numeric(str_extract(experiment_dirs, "\\d{3}"))
  latest_experiment <- max(experiment_numbers, na.rm = TRUE)
  return(sprintf("%03d", latest_experiment))
}

# Get the latest experiment number and increment by 1
latest_exp <- find_latest_experiment(base_dir)
new_exp <- sprintf("%03d", as.numeric(latest_exp) + 1)

# Validate the new experiment number
if (as.numeric(new_exp) > 50) {
  stop("Error: Maximum number of experiments (50) reached.")
}

# Create new experiment directory and subdirectory
new_exp_dir <- file.path(base_dir, paste0("exp_", new_exp))
cat("Creating new experiment directory:", new_exp_dir, "\n")
dir.create(new_exp_dir)

# Create Results subdirectory
results_dir <- file.path(new_exp_dir, "Results")
cat("Creating Results directory:", results_dir, "\n")
dir.create(results_dir)

# Create a readme.txt file explaining the experiment number
readme_content <- paste("This is Experiment", new_exp, "of the project.")
writeLines(readme_content, file.path(new_exp_dir, "readme.txt"))

# Copy R scripts from the latest experiment directory to the new one
latest_exp_dir <- file.path(base_dir, paste0("exp_", latest_exp))
cat("Copying scripts from:", latest_exp_dir, "\n")
for (i in 1:10) {
  file.copy(file.path(latest_exp_dir, paste0("sims", i, ".R")), new_exp_dir)
}

# Copy VisualizeResults.R from the previous experiment's Results folder
previous_exp <- sprintf("%03d", as.numeric(latest_exp) - 1)
previous_results_dir <- file.path(base_dir, paste0("exp_", previous_exp), "Results")
file.copy(file.path(previous_results_dir, "VisualizeResults.R"), results_dir)

# Update the folder name in the copied scripts
for (i in 1:10) {
  script_path <- file.path(new_exp_dir, paste0("sims", i, ".R"))
  script_content <- readLines(script_path)
  script_content <- gsub(paste0("exp_", latest_exp), paste0("exp_", new_exp), script_content)
  writeLines(script_content, script_path)
}

# Update the first line of VisualizeResults.R to reflect the new experiment directory
visualize_script_content <- readLines(file.path(results_dir, "VisualizeResults.R"))
visualize_script_content[1] <- paste0('setwd("', results_dir, '")')
writeLines(visualize_script_content, file.path(results_dir, "VisualizeResults.R"))

cat("Experiment", new_exp, "preprocessed successfully. Results directory:", results_dir, "\n")
