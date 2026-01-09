#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))

# Define paths
processed_dir <- file.path("data", "processed", "03_import_epi_features")
input_path <- file.path(processed_dir, "apowin_out_part_1.txt")
output_path <- file.path(processed_dir, "aopwin_features.csv")

ensure_directories <- function() {
  dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
}

# Helper to parse scientific notation from text (e.g., "1.23 E-12" -> 1.23e-12)
parse_sci <- function(text) {
  # Remove units and everything after
  val_str <- str_extract(text, "[0-9.]+\\s+E[-+]?[0-9]+")
  if (is.na(val_str)) return(NA_real_)
  as.numeric(str_replace(val_str, "\\s", ""))
}

# Helper to extract value based on regex patterns with capture group
extract_val <- function(lines, pattern) {
  # pattern should have one capture group for the number
  # Find lines matching the pattern
  matching_lines <- grep(pattern, lines, value = TRUE)
  if (length(matching_lines) == 0) return(NA_real_)
  
  # Extract the first match
  # str_match returns a matrix where column 2 is the first capture group
  val_str <- str_match(matching_lines[1], pattern)[1, 2]
  as.numeric(val_str)
}

# Helper to parse the first half-life line and return hours
parse_first_half_life_hours <- function(lines) {
  # Find all HALF-LIFE lines
  hl_lines <- grep("HALF-LIFE =", lines, value = TRUE)
  if (length(hl_lines) == 0) return(NA_real_)
  
  # Take the first one
  line <- hl_lines[1]
  
  # Extract value and unit
  val_part <- str_extract(line, "HALF-LIFE =\\s+([0-9.]+)\\s+(Hrs|Days|Min)")
  if (is.na(val_part)) return(NA_real_)
  
  val <- as.numeric(str_extract(val_part, "[0-9.]+"))
  unit <- str_extract(val_part, "Hrs|Days|Min")
  
  if (unit == "Hrs") {
    return(val)
  } else if (unit == "Min") {
    return(val / 60)
  } else if (unit == "Days") {
    # Check for specific 12-hr day assumption string
    # String: "(12-hr day; 1.5E6 OH/cm3)"
    # We use fixed string matching or regex
    if (str_detect(line, fixed("(12-hr day; 1.5E6 OH/cm3)"))) {
      return(val * 12)
    } else {
      return(val * 24)
    }
  }
  
  NA_real_
}

process_block <- function(block_text) {
  lines <- str_split(block_text, "\n")[[1]]
  
  # Extract Molecule ID (CHEM :)
  chem_line <- grep("^CHEM\\s+:", lines, value = TRUE)
  if (length(chem_line) == 0) return(NULL)
  molecule_id <- str_extract(chem_line, "[0-9]+") %>% as.numeric()
  
  # Split sections
  oh_start <- grep("HYDROXYL RADICALS", lines)
  ozone_start <- grep("OZONE REACTION", lines)
  
  # Ensure we have at least an OH section to proceed (or handled gracefully)
  if (length(oh_start) == 0) {
    # If no OH section, return NAs or skip? 
    # Usually AOP output has it. If not, we return basic NAs.
    return(tibble(
        molecule_id = molecule_id,
        aop_oh_rate_constant = NA_real_,
        aop_oh_half_life_hours = NA_real_,
        aop_ozone_rate_constant = NA_real_,
        aop_ozone_half_life_hours = NA_real_
    ))
  }
  
  oh_end <- if (length(ozone_start) > 0) ozone_start[1] - 1 else length(lines)
  oh_lines <- lines[oh_start[1]:oh_end]
  
  ozone_lines <- if (length(ozone_start) > 0) lines[ozone_start[1]:length(lines)] else character(0)

  # --- OH Section ---
  # OH Rate
  oh_rate_line <- grep("OVERALL OH Rate Constant =", oh_lines, value = TRUE)
  oh_rate <- if(length(oh_rate_line) > 0) parse_sci(oh_rate_line[1]) else NA_real_
  
  # OH Half Life (Hours)
  oh_hl_hours <- parse_first_half_life_hours(oh_lines)
  
  # --- Ozone Section ---
  ozone_rate <- NA_real_
  ozone_hl_hours <- NA_real_
  
  if (length(ozone_lines) > 0) {
    if (any(str_detect(ozone_lines, "NO OZONE REACTION ESTIMATION"))) {
       ozone_rate <- 0
    } else {
       oz_rate_line <- grep("OVERALL OZONE Rate Constant =", ozone_lines, value = TRUE)
       if (length(oz_rate_line) > 0) {
         ozone_rate <- parse_sci(oz_rate_line[1])
         ozone_hl_hours <- parse_first_half_life_hours(ozone_lines)
       }
    }
  }

  tibble(
    molecule_id = molecule_id,
    aop_oh_rate_constant = oh_rate,
    aop_oh_half_life_hours = oh_hl_hours,
    aop_ozone_rate_constant = ozone_rate,
    aop_ozone_half_life_hours = ozone_hl_hours
  )
}

main <- function() {
  suppressPackageStartupMessages(library(furrr))
  suppressPackageStartupMessages(library(future))

  ensure_directories()
  
  if (!file.exists(input_path)) {
    stop(sprintf("Input file not found: %s", input_path))
  }
  
  # Read full text
  content <- read_file(input_path)
  
  # Split by "SMILES :" which marks start of new entry
  # But the first entry starts with SMILES :, subsequent also.
  # We can split by regex "(?=\nSMILES :)" or just "SMILES :"
  # Use Lookahead to keep the delimiter? Or just prepend it back.
  # simple split by "SMILES :" gives empty first element.
  
  blocks <- str_split(content, "SMILES :")[[1]]
  blocks <- blocks[nchar(str_trim(blocks)) > 0] # remove empty
  
  # Set up parallel processing
  # Use available cores minus 1 to keep system responsive
  workers <- max(1, parallel::detectCores(logical = FALSE) - 1)
  plan(multisession, workers = workers)
  message(sprintf("Using %d workers for processing...", workers))

  # Process each block in parallel
  features <- future_map_dfr(blocks, process_block, .progress = TRUE)
  
  # Save
  readr::write_csv(features, output_path)
  message(sprintf("Saved %d features to %s", nrow(features), output_path))
}

if (identical(environment(), globalenv())) {
  main()
}
