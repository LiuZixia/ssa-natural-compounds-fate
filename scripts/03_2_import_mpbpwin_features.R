#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))

# Define paths
processed_dir <- file.path("data", "processed", "03_import_epi_features")
input_path <- file.path(processed_dir, "mpbpwin_out_part_1.txt")
output_path <- file.path(processed_dir, "mpbpwin_features.csv")

main <- function() {
  if (!file.exists(input_path)) {
    stop(sprintf("Input file not found: %s", input_path))
  }
  
  message("Reading MPBPWIN output...")
  lines <- read_lines(input_path)
  
  # Filter out empty lines
  lines <- lines[str_trim(lines) != ""]
  
  message(sprintf("Processing %d lines...", length(lines)))
  
  # --- Step 1: Merge split error lines ---
  # Pattern to identify ANY error line start
  patt_err_start <- "^\\s*\\*+>"
  # Pattern to remove header from continuation lines (matches up to first colon and spaces)
  patt_header_remove <- "^\\s*\\*+>.*?:\\s*"
  
  merged_lines <- character(length(lines))
  write_idx <- 1
  i <- 1
  n <- length(lines)
  
  while (i <= n) {
    current_line <- lines[i]
    
    if (str_detect(current_line, patt_err_start)) {
      # It is an error line. Check for continuations.
      j <- i + 1
      while (j <= n && str_detect(lines[j], patt_err_start)) {
        
        # Prepare to check if we should merge lines[j] into current_line
        
        # 1. Check if lines[j] content starts with a digit (Split ID case)
        # Remove the header prefix to get the content
        content_j <- str_remove(lines[j], patt_header_remove)
        starts_with_digit <- str_detect(content_j, "^[0-9]")
        
        # 2. Check if current_line looks like a complete entry (Ends with " ID")
        current_ends_with_id <- str_detect(current_line, "\\s+[0-9]+\\s*$")
        
        # DECISION: Merge if it's a split ID OR if the previous line was incomplete
        if (starts_with_digit || !current_ends_with_id) {
          # Merge!
          current_line <- paste0(current_line, content_j) # content_j is already trimmed of header
          j <- j + 1
        } else {
          # It's a new error entry (consecutive errors)
          break
        }
      }
      merged_lines[write_idx] <- current_line
      write_idx <- write_idx + 1
      i <- j # Skip consumed lines
    } else {
      # Normal line
      merged_lines[write_idx] <- current_line
      write_idx <- write_idx + 1
      i <- i + 1
    }
  }
  
  lines <- merged_lines[1:(write_idx - 1)]
  message(sprintf("Collapsed to %d lines after merging split errors", length(lines)))
  
  # --- Step 2: Parse Features ---
  
  # Regex patterns
  # Valid: Capture groups 1=BP, 2=MP, 3=VP, 4=ID
  regex_valid <- "^\\s*([0-9.E+-]+)\\(BP est\\)\\s+([0-9.E+-]+)\\(MP est\\)\\s+([0-9.E+-]+)\\(VP est\\)\\s+\\S+\\s+([0-9]+)\\s*$"
  
  # Error: Generic capture for ANY line starting with *>, ending in ID
  regex_error_generic <- "^\\s*\\*+>.*\\s([0-9]+)\\s*$"
  
  # 1. Match Valid Lines
  valid_matches <- str_match(lines, regex_valid)
  
  # 2. Match Error Lines (Generic)
  error_matches <- str_match(lines, regex_error_generic)
  
  # Construct data frames
  
  # Valid entries
  valid_indices <- !is.na(valid_matches[, 1])
  valid_df <- tibble(
    molecule_id = as.numeric(valid_matches[valid_indices, 5]),
    bp_est = as.numeric(valid_matches[valid_indices, 2]),
    mp_est = as.numeric(valid_matches[valid_indices, 3]),
    vp_est = as.numeric(valid_matches[valid_indices, 4])
  )
  
  # Error entries
  error_indices <- !is.na(error_matches[, 1])
  if (sum(error_indices) > 0) {
     error_df <- tibble(
      molecule_id = as.numeric(error_matches[error_indices, 2]),
      bp_est = NA_real_,
      mp_est = NA_real_,
      vp_est = NA_real_
    )
  } else {
    error_df <- tibble(molecule_id = numeric(), bp_est = numeric(), mp_est = numeric(), vp_est = numeric())
  }

  
  # Combine
  features <- bind_rows(valid_df, error_df) %>%
    arrange(molecule_id)
  
  # Check if we missed any lines (optional debugging)
  # total_parsed <- nrow(features)
  # total_lines <- length(lines)
  # if (total_parsed < total_lines) {
  #   warning(sprintf("Failed to parse %d lines", total_lines - total_parsed))
  # }
  
  # Write output
  write_csv(features, output_path)
  message(sprintf("Saved %d features to %s", nrow(features), output_path))
}

if (identical(environment(), globalenv())) {
  main()
}
