#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

# Define paths
# Inputs
chem_features_path <- file.path("data", "processed", "02_prepare_features", "chemical_features.csv")
aop_features_path  <- file.path("data", "processed", "03_import_epi_features", "aopwin_features.csv")
mpbp_features_path <- file.path("data", "processed", "03_import_epi_features", "mpbpwin_features.csv")

# Output
# We'll save it in 03_import_epi_features directory for consistency
output_dir <- file.path("data", "processed", "03_import_epi_features")
output_path <- file.path(output_dir, "chemical_merged_features.csv")

ensure_directories <- function() {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

main <- function() {
  ensure_directories()
  
  # Check inputs
  if (!file.exists(chem_features_path)) {
    stop(sprintf("Base file not found: %s", chem_features_path))
  }
  
  message("Reading chemical features...")
  chem <- read_csv(chem_features_path, show_col_types = FALSE)
  
  # Read AOPWIN if exists
  if (file.exists(aop_features_path)) {
    message("Reading AOPWIN features...")
    aop <- read_csv(aop_features_path, show_col_types = FALSE)
  } else {
    warning(sprintf("AOPWIN file not found: %s. Creating empty placeholder.", aop_features_path))
    aop <- tibble(molecule_id = numeric())
  }
  
  # Read MPBPWIN if exists
  if (file.exists(mpbp_features_path)) {
    message("Reading MPBPWIN features...")
    mpbp <- read_csv(mpbp_features_path, show_col_types = FALSE)
  } else {
    warning(sprintf("MPBPWIN file not found: %s. Creating empty placeholder.", mpbp_features_path))
    mpbp <- tibble(molecule_id = numeric())
  }
  
  # Merge
  # We use left_join onto the 'chem' dataset because it contains the master list of molecules
  # derived from the raw data.
  message("Merging features by molecule_id...")
  
  merged_data <- chem %>%
    left_join(aop, by = "molecule_id") %>%
    left_join(mpbp, by = "molecule_id")
  
  # Save
  write_csv(merged_data, output_path)
  message(sprintf("Successfully merged data. Rows: %d, Cols: %d", nrow(merged_data), ncol(merged_data)))
  message(sprintf("Saved to: %s", output_path))
}

if (identical(environment(), globalenv())) {
  main()
}
