#!/usr/bin/env Rscript

# 05_3_clustering_method_selection.R
# Compare clustering algorithms.
# Optimization Level: High (64GB RAM Tuning) + Dynamic K Range
# - Dynamic K: Starts at N_features + 1, grows geometrically (factor ~1.5).

packages <- c("tidyverse", "cluster", "factoextra", "fastcluster", "foreach", "doParallel", "ClusterR")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(cluster)
library(factoextra)
library(fastcluster)
library(foreach)
library(doParallel)
library(ClusterR)

# --- Configuration ---
input_dir <- "data/processed/05_clustering"
output_dir <- "reports/05_clustering_comparison"

if (!dir.exists(input_dir)) input_dir <- "../data/processed/05_clustering"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files_to_process <- list(
    enrichment = "transformed_clustering_properties_enrichment.csv",
    aerosolization = "transformed_clustering_properties_aerosolization.csv"
)

# Tuning for 64GB RAM
MAX_SIL_SAMPLE <- 10000 # Reduced slightly for multiple runs
MAX_HC_SAMPLE <- 30000
N_BOOT <- 25 # Number of bootstrap iterations

# K Range Settings
K_GROWTH_FACTOR <- 1.25
K_MAX_LIMIT <- 50 # Optimization cutoff

# Parallel Setup
n_cores <- 28
if (n_cores < 1) n_cores <- 1
message("Using ", n_cores, " cores for parallel processing.")
registerDoParallel(cores = n_cores)

# --- Dynamic K Generator ---
generate_k_sequence <- function(n_features) {
    k_start <- round(n_features / 2)
    k_seq <- c(k_start)

    curr <- k_start
    while (TRUE) {
        next_k <- round(curr * K_GROWTH_FACTOR)
        if (next_k == curr) next_k <- curr + 1 # Ensure growth

        if (next_k > K_MAX_LIMIT) break

        k_seq <- c(k_seq, next_k)
        curr <- next_k
    }
    return(unique(k_seq))
}

compare_methods <- function(df, label) {
    message("\nComparing methods for: ", label)

    df_model <- df %>%
        select(-matches("molecule_id")) %>%
        drop_na()

    if (nrow(df_model) < 50) {
        message("  Not enough data. Skipping.")
        return(NULL)
    }

    mat <- as.matrix(df_model)
    n_samples <- nrow(mat)
    n_features <- ncol(mat)
    message("  Data dimensions: ", n_samples, " x ", n_features)

    # Generate K Range
    k_range <- generate_k_sequence(n_features)
    message("  Testing K grid: ", paste(k_range, collapse = ", "))

    # Parameter grid: k * bootstrap_iterations
    grid <- expand.grid(k = k_range, boot = 1:N_BOOT)

    # --- Parallel Loop ---
    # Iterate over both K and Bootstrap samples
    results_list <- foreach(i = 1:nrow(grid), .packages = c("cluster", "fastcluster", "ClusterR"), .export = c("MAX_HC_SAMPLE", "MAX_SIL_SAMPLE")) %dopar% {
        k <- grid$k[i]

        # --- Sampling for this iteration ---
        # We resample evaluation set AND clustering set (for HC) to capture full variance
        if (n_samples > MAX_SIL_SAMPLE) {
            sil_idx <- sample(seq_len(n_samples), MAX_SIL_SAMPLE)
        } else {
            sil_idx <- seq_len(n_samples)
        }
        mat_sil <- mat[sil_idx, , drop = FALSE]

        # Local Dist for Metric (re-calculated per bootstrap)
        d_sil <- dist(mat_sil, method = "euclidean")

        if (n_samples > MAX_HC_SAMPLE) {
            hc_idx <- sample(seq_len(n_samples), MAX_HC_SAMPLE)
            mat_hc <- mat[hc_idx, , drop = FALSE]
        } else {
            mat_hc <- mat
            hc_idx <- seq_len(n_samples)
        }


        # 1. K-Means
        # Run on full data (fast), evaluate on subsample
        km_res <- ClusterR::KMeans_rcpp(mat, clusters = k, num_init = 5, max_iters = 50)
        labels_km <- km_res$clusters[sil_idx]

        ss_km <- silhouette(labels_km, d_sil)
        avg_sil_km <- mean(ss_km[, 3])

        # 2. HC (Ward)
        hc_ward <- fastcluster::hclust.vector(mat_hc, method = "ward")
        grp_ward <- cutree(hc_ward, k = k)

        # Map back to sil_idx.
        # If sil_idx is not subset of hc_idx, we have a problem.
        # Strategy: Cluster on mat_hc. Evaluate on intersection(hc_idx, sil_idx) OR just use hc_idx as sil_idx for HC?
        # Better: For HC, we must evaluate on the data we clustered if we can't project.
        # But we want comparable metrics.
        # Solution: Use sil_idx for BOTH HC and Evaluation.
        # i.e. subsample data -> cluster that subsample -> evaluate that subsample.
        # This is strictly "clustering stability on subset".

        hc_sub_ward <- fastcluster::hclust.vector(mat_sil, method = "ward")
        grp_sub_ward <- cutree(hc_sub_ward, k = k)
        ss_ward <- silhouette(grp_sub_ward, d_sil)
        avg_sil_ward <- mean(ss_ward[, 3])

        # 3. CLARA
        # CLARA handles sampling internally, but here we fix the sample to sil_idx for direct comparison?
        # No, let CLARA do its thing on full data, check labels on sil_idx.
        clara_res <- clara(mat, k = k, metric = "euclidean", samples = 50, pamLike = TRUE)
        labels_clara <- clara_res$clustering[sil_idx]

        if (any(is.na(labels_clara))) {
            avg_sil_clara <- NA
        } else {
            ss_clara <- silhouette(labels_clara, d_sil)
            avg_sil_clara <- mean(ss_clara[, 3])
        }

        # Return
        data.frame(
            k = k,
            boot = grid$boot[i],
            Method = c("K-Means", "HC-Ward", "CLARA"),
            Silhouette = c(avg_sil_km, avg_sil_ward, avg_sil_clara)
        )
    }

    results <- do.call(rbind, results_list) %>% drop_na()

    # --- Aggregation ---
    results_agg <- results %>%
        group_by(k, Method) %>%
        summarise(
            mean_sil = mean(Silhouette),
            sd_sil = sd(Silhouette),
            se_sil = sd(Silhouette) / sqrt(n()),
            .groups = "drop"
        )

    # 95% CI
    results_agg <- results_agg %>%
        mutate(
            ci_lower = mean_sil - 1.96 * se_sil,
            ci_upper = mean_sil + 1.96 * se_sil
        )

    # --- Plot ---
    p <- ggplot(results_agg, aes(x = k, y = mean_sil, color = Method, group = Method)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, alpha = 0.6) +
        theme_minimal() +
        labs(
            title = paste0("Clustering Performance: ", label),
            subtitle = paste0("Features = ", n_features, " | Bootstraps = ", N_BOOT, " | Error Bars = 95% CI"),
            y = "Avg Silhouette Width",
            x = "Number of Clusters (k)"
        ) +
        scale_x_continuous(breaks = k_range)

    plot_file <- file.path(output_dir, paste0("comparison_", label, ".png"))
    ggsave(plot_file, p, width = 8, height = 6)
    message("  Saved plot to: ", plot_file)

    csv_file <- file.path(output_dir, paste0("comparison_metrics_", label, ".csv"))
    write_csv(results, csv_file) # Save full results

    agg_file <- file.path(output_dir, paste0("comparison_summary_", label, ".csv"))
    write_csv(results_agg, agg_file) # Save summary

    best_res <- results_agg %>%
        group_by(Method) %>%
        filter(mean_sil == max(mean_sil)) %>%
        slice(1) %>%
        ungroup()
    print(best_res)
}

tryCatch({
    for (lbl in names(files_to_process)) {
        f <- files_to_process[[lbl]]
        f_path <- file.path(input_dir, f)
        if (file.exists(f_path)) compare_methods(read_csv(f_path, show_col_types = FALSE), lbl)
    }
}, finally = {
    stopImplicitCluster()
})

message("\nDone.")
