#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

# Define paths
raw_dir <- file.path("data", "raw")
processed_dir <- file.path("data", "processed", "02_prepare_features")
molecules_path <- file.path(raw_dir, "molecules.csv")
properties_path <- file.path(raw_dir, "properties.csv")
output_csv_path <- file.path(processed_dir, "chemical_features.csv")
output_smiles_path <- file.path(processed_dir, "canonical_smiles.txt")

ensure_directories <- function() {
  dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
}

main <- function() {
  ensure_directories()
  
  # Read data
  properties <- readr::read_csv(properties_path, show_col_types = FALSE)
  molecules <- readr::read_csv(molecules_path, show_col_types = FALSE)
  
  # Rename IDs as requested
  # properties.csv: id -> property_id
  properties <- properties %>%
    rename(property_id = id)
  
  # molecules.csv: id -> molecule_id
  molecules <- molecules %>%
    rename(molecule_id = id)
  
  # Merge tables
  # Using full_join to keep all data. 
  # Common columns other than key (molecule_id) will be suffixed with .x and .y
  merged_data <- properties %>%
    full_join(molecules, by = "molecule_id") %>%
    select(-any_of(c("created_at.x", "updated_at.x", "created_at.y", "updated_at.y")))
  
  # Filter out problematic molecules causing errors in EPISuite
  problematic_ids <- c(546623, 685189, 71730)
  merged_data <- merged_data %>%
    filter(!molecule_id %in% problematic_ids)
  
  # Save merged file
  readr::write_csv(merged_data, output_csv_path)
  message(sprintf("Merged data saved to %s", output_csv_path))
  
  # Generate canonical_smiles.txt
  # Format: "canonical_smiles molecule_id"
  # We need to filter out rows where canonical_smiles or molecule_id might be missing to ensure valid output
  smiles_data <- merged_data %>%
    filter(!is.na(canonical_smiles) & !is.na(molecule_id)) %>%
    select(canonical_smiles, molecule_id) %>%
    distinct() # Avoid duplicates if multiple properties map to same molecule
    
  # Create the formatted string
  # "canonical_smiles string comes first and ends with the first blank space, molecule_id follow the bank space"
  # Assuming this means just pasting them together with a space
  lines <- paste(smiles_data$canonical_smiles, smiles_data$molecule_id)
  
  # Split into chunks of 10,000 lines
  chunk_size <- 100000
  total_lines <- length(lines)
  num_chunks <- ceiling(total_lines / chunk_size)
  
  for (i in seq_len(num_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, total_lines)
    chunk <- lines[start_idx:end_idx]
    
    chunk_filename <- sprintf("canonical_smiles_part_%d.txt", i)
    chunk_path <- file.path(processed_dir, chunk_filename)
    
    readr::write_lines(chunk, chunk_path)
    message(sprintf("Saved chunk %d to %s", i, chunk_path))
  }
}

if (identical(environment(), globalenv())) {
  main()
}
