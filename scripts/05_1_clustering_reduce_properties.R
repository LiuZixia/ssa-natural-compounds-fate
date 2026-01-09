#!/usr/bin/env Rscript

# 05_1_clustering_reduce_properties.R
# This script performs statistical filtering and feature selection for clustering.
# It processes two views: Enrichment and Aerosolization.
# Steps:
# 1. Basic statistical filtering
#    1.1 Near-zero variance filter (>95-99% same value)
#    1.2 Missingness filter (>40% missing dropped, unless critical)
#    1.3 Correlation/collinearity pruning (|rho| > 0.9)

suppressPackageStartupMessages({
    library(tidyverse)
    library(caret)
    library(reshape2)
})

# --- Configuration ---
input_path <- "data/processed/03_import_epi_features/chemical_merged_features.csv"
output_dir <- "data/processed/05_clustering"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Thresholds
nzv_freq_cut <- 95 / 5 # 95% cutoff
missing_threshold <- 0.40 # 40% missingness
correlation_cutoff <- 0.95

# Critical variables to keep even if missingness is high
critical_vars <- c("vp_est")

# Variable Definitions (from theory labels)
vars_both <- c(
    "total_atom_count", "heavy_atom_count", "molecular_weight",
    "exact_molecular_weight", "alogp", "formal_charge",
    "number_of_minimal_rings", "contains_sugar",
    "contains_ring_sugars", "contains_linear_sugars"
)

vars_enrichment_only <- c(
    "topological_polar_surface_area", "rotatable_bond_count",
    "hydrogen_bond_acceptors", "hydrogen_bond_donors",
    "aromatic_rings_count", "fractioncsp3", "van_der_walls_volume"
)

vars_aerosolization_only <- c("bp_est", "mp_est", "vp_est")

# Define Views
views <- list(
    enrichment = c(vars_enrichment_only, vars_both),
    aerosolization = c(vars_aerosolization_only, vars_both)
)

# Priority rules for collinearity (Higher score = keep)
# We prioritize "Both" variables often as they are core physchem,
# but mechanism-specific ones might be preferred if the conflict is within the set.
# However, the user says "Prefer the one with clearer mechanism... theory labels".
# We'll define a simple hierarchy or score.
# Here we just implement a heuristic: prioritize by lower missingness first, then manual preference.
# Manual preference list (higher is better):
manual_priority <- c(
    "vp_est" = 10,
    "molecular_weight" = 5,
    "exact_molecular_weight" = 4, # Prefer MW over Exact MW (simpler)
    "heavy_atom_count" = 3,
    "total_atom_count" = 3
)

# --- Helper Functions ---

get_priority <- function(var_name, missing_pct) {
    # Base priority from manual list (default 0)
    p <- ifelse(var_name %in% names(manual_priority), manual_priority[var_name], 0)
    # Penalize heavy missingness heavily
    p <- p - (missing_pct * 100)
    return(p)
}

process_view <- function(df, view_vars, view_name) {
    message(sprintf("\n--- Processing View: %s ---", view_name))

    # Select variables + ID
    # Ensure variables exist in df
    valid_vars <- intersect(view_vars, names(df))
    missing_vars <- setdiff(view_vars, names(df))
    if (length(missing_vars) > 0) message("Warning: Missing variables in data: ", paste(missing_vars, collapse = ", "))

    # Start with relevant subset
    # We assume molecule_id is key
    df_sub <- df %>% select(molecule_id, all_of(valid_vars))

    # --- 1.1 Near-Zero Variance ---
    # Only apply to numeric columns (vars_both etc are usually numeric or binary/counts)
    # Convert logical to numeric if any
    df_numeric <- df_sub %>%
        select(-molecule_id) %>%
        mutate(across(everything(), as.numeric))

    nzv_info <- nearZeroVar(df_numeric, saveMetrics = TRUE, freqCut = nzv_freq_cut)
    drop_nzv <- rownames(nzv_info)[nzv_info$nzv]

    if (length(drop_nzv) > 0) {
        message("Dropping NZV variables: ", paste(drop_nzv, collapse = ", "))
        df_sub <- df_sub %>% select(-all_of(drop_nzv))
        valid_vars <- setdiff(valid_vars, drop_nzv)
    }

    # --- 1.2 Missingness Filter ---
    missing_pct <- colMeans(is.na(df_sub %>% select(all_of(valid_vars))))
    drop_missing <- names(missing_pct)[missing_pct > missing_threshold]

    # Check critical variables
    kept_critical <- intersect(drop_missing, critical_vars)
    drop_missing <- setdiff(drop_missing, critical_vars)

    if (length(drop_missing) > 0) {
        message("Dropping High Missingness variables: ", paste(drop_missing, collapse = ", "))
        df_sub <- df_sub %>% select(-all_of(drop_missing))
        valid_vars <- setdiff(valid_vars, drop_missing)
    }

    # Handle Critical Variables (Add Missing Flag)
    for (cv in kept_critical) {
        if (cv %in% names(df_sub)) {
            flag_col <- paste0(cv, "_missing")
            message("Adding missing flag for critical variable: ", cv)
            df_sub[[flag_col]] <- as.integer(is.na(df_sub[[cv]]))
            # We do NOT impute the value itself unless needed for correlation.
            # Correlation calc below uses pairwise.complete.obs, so NAs are ok there.
        }
    }

    # --- 1.3 Correlation Pruning ---
    # Compute correlation on remaining numeric vars
    # (Exclude missing flags from correlation check? Usually we check feature-feature corr)

    # Re-evaluate valid vars (excluding flags for now, but flags might be useful?
    # Usually we prune the main features)
    current_vars <- valid_vars

    if (length(current_vars) > 1) {
        # Calculate cor matrix (Spearman)
        # Use pairwise complete obs
        cor_mat <- cor(df_sub %>% select(all_of(current_vars)), use = "pairwise.complete.obs", method = "spearman")
        cor_mat[is.na(cor_mat)] <- 0 # Handle cases where cor cannot be computed (e.g. 0 variance but passed nzv?)

        # Find high correlations
        high_corr <- which(abs(cor_mat) > correlation_cutoff & lower.tri(cor_mat), arr.ind = TRUE)

        vars_to_remove <- c()

        if (nrow(high_corr) > 0) {
            message("Found ", nrow(high_corr), " highly correlated pairs (|rho| > ", correlation_cutoff, "). Pruning...")

            # Sort pairs by correlation magnitude (highest first)
            pairs_df <- data.frame(
                v1 = rownames(cor_mat)[high_corr[, 1]],
                v2 = colnames(cor_mat)[high_corr[, 2]],
                rho = cor_mat[high_corr]
            ) %>% arrange(desc(abs(rho)))

            for (i in 1:nrow(pairs_df)) {
                v1 <- pairs_df$v1[i]
                v2 <- pairs_df$v2[i]

                # If one is already removed, skip
                if (v1 %in% vars_to_remove || v2 %in% vars_to_remove) next

                # Decide which to keep
                # 1. Check missingness in current df_sub
                m1 <- mean(is.na(df_sub[[v1]]))
                m2 <- mean(is.na(df_sub[[v2]]))

                # Calculate Priority Score
                p1 <- get_priority(v1, m1)
                p2 <- get_priority(v2, m2)

                remove_var <- if (p1 >= p2) v2 else v1
                keep_var <- if (p1 >= p2) v1 else v2

                message(sprintf("  Pair: %s vs %s (rho=%.2f). Keeping %s.", v1, v2, pairs_df$rho[i], keep_var))
                vars_to_remove <- c(vars_to_remove, remove_var)
            }
        }

        if (length(vars_to_remove) > 0) {
            df_sub <- df_sub %>% select(-all_of(vars_to_remove))
            message("Removed ", length(vars_to_remove), " features due to collinearity.")
        }
    }

    return(df_sub)
}

# --- Main Execution ---

if (!file.exists(input_path)) {
    stop("Input file not found: ", input_path)
}

message("Reading data from ", input_path)
df_full <- read_csv(input_path, show_col_types = FALSE)

# Process Enrichment View
df_enrich <- process_view(df_full, views$enrichment, "Enrichment")
enrich_out <- file.path(output_dir, "clustering_properties_enrichment.csv")
write_csv(df_enrich, enrich_out)
message("Saved Enrichment View to: ", enrich_out)

# Process Aerosolization View
df_aerosol <- process_view(df_full, views$aerosolization, "Aerosolization")
aerosol_out <- file.path(output_dir, "clustering_properties_aerosolization.csv")
write_csv(df_aerosol, aerosol_out)
message("Saved Aerosolization View to: ", aerosol_out)

message("\nDone.")
