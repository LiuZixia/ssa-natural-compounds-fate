#!/usr/bin/env Rscript

# 05_2_clustering_transform.R
# Transform + standardize variables for clustering

if (!require("tidyverse", quietly = TRUE)) install.packages("tidyverse")
library(tidyverse)

# --- 1. Configuration ---
input_dir  <- "data/processed/05_clustering"
output_dir <- "data/processed/05_clustering"

if (!dir.exists(input_dir)) {
  input_dir  <- "../data/processed/05_clustering"
  output_dir <- "../data/processed/05_clustering"
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files_to_process <- list(
  enrichment     = "clustering_properties_enrichment.csv",
  aerosolization = "clustering_properties_aerosolization.csv"
)

SKEW_THRESHOLD <- 1.0
VP_EPS_FLOOR   <- 1e-30  # safe floor for log10

# --- Helper Functions ---
calc_skewness <- function(x) {
  x <- na.omit(x)
  if (length(x) < 3) return(0)
  s <- sd(x)
  if (s == 0) return(0)
  n <- length(x)
  m <- mean(x)
  sum((x - m)^3) / ((n - 1) * s^3)
}

is_binary_01 <- function(x) {
  ux <- unique(na.omit(x))
  length(ux) <= 2 && all(ux %in% c(0, 1))
}

is_vp_col <- function(colname) {
  # catches vp_est, vapor_pressure, vp, vpmmhg, etc.
  grepl("^vp(_|$)", colname, ignore.case = TRUE) ||
    grepl("vapor.*press", colname, ignore.case = TRUE)
}

process_file <- function(filename, label) {
  in_path <- file.path(input_dir, filename)
  if (!file.exists(in_path)) {
    message("Skipping ", label, " (File not found: ", in_path, ")")
    return(invisible(NULL))
  }
  
  message("\nProcessing: ", label, " (", in_path, ")")
  df <- read_csv(in_path, show_col_types = FALSE)
  
  # Drop rows with any NA across feature columns (keeps your current behavior)
  feature_cols <- setdiff(names(df), "molecule_id")
  before_n <- nrow(df)
  df <- df %>% drop_na(all_of(feature_cols))
  after_n <- nrow(df)
  if (before_n != after_n) {
    message("Dropped ", before_n - after_n, " rows with missing values.")
  }
  if (nrow(df) == 0) {
    message("No data left after dropping NAs.")
    return(invisible(NULL))
  }
  
  # --- Transform ---
  for (col in feature_cols) {
    x <- df[[col]]
    if (!is.numeric(x)) next
    
    # binary 0/1 -> keep as is (no log)
    if (is_binary_01(x)) next
    
    # Special handling for vapor pressure (vp_est etc.)
    if (is_vp_col(col)) {
      x_pos <- x[is.finite(x) & x > 0]
      eps <- if (length(x_pos) > 0) min(x_pos) / 10 else VP_EPS_FLOOR
      eps <- max(eps, VP_EPS_FLOOR)
      
      n_zero <- sum(x == 0, na.rm = TRUE)
      df[[col]] <- log10(pmax(x, eps))
      
      message(sprintf(
        "  VP: log10(pmax(%s, eps)) applied (eps=%.2e; zeros replaced=%d)",
        col, eps, n_zero
      ))
      next
    }
    
    # Skewness-based transform for other numeric columns
    sk <- calc_skewness(x)
    if (abs(sk) > SKEW_THRESHOLD) {
      min_val <- min(x, na.rm = TRUE)
      
      if (min_val > 0) {
        df[[col]] <- log10(x)
        message(sprintf("  Log10(x) applied to %s (Skew=%.2f, Min=%.2e)", col, sk, min_val))
        
      } else if (min_val == 0) {
        # OK for count-like variables, but NOT for vapor pressure (handled above)
        df[[col]] <- log10(x + 1)
        message(sprintf("  Log10(x+1) applied to %s (Skew=%.2f, Min=0)", col, sk))
        
      } else {
        message(sprintf("  Skipping Log for %s (Skew=%.2f) due to negatives (Min=%.2f)", col, sk, min_val))
      }
    }
  }
  
  # --- Standardize (Z-score): scale only non-binary numeric columns ---
  num_cols <- feature_cols[sapply(df[feature_cols], is.numeric)]
  scale_cols <- num_cols[!sapply(df[num_cols], is_binary_01)]
  if (length(scale_cols) > 0) {
    df[scale_cols] <- scale(df[scale_cols])
  }
  
  # Save
  out_name <- paste0("transformed_", filename)
  out_path <- file.path(output_dir, out_name)
  write_csv(df, out_path)
  message("Saved transformed data to: ", out_path)
  
  invisible(df)
}

# --- Execution ---
for (lbl in names(files_to_process)) {
  process_file(files_to_process[[lbl]], lbl)
}

message("\nDone.")
