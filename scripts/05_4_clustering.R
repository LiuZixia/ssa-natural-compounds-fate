#!/usr/bin/env Rscript

# 05_4_clustering.R
# Perform clustering using K-means (k=8/10) + Sub-clustering of wide clusters.

packages <- c("tidyverse", "cluster", "factoextra", "ClusterR")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(cluster)
library(factoextra)
library(ClusterR)

# --- Configuration ---
input_dir <- "data/processed/05_clustering"
output_dir <- "data/processed/05_clustering"

if (!dir.exists(input_dir)) input_dir <- "data/processed/05_clustering"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files_to_process <- list(
    enrichment = "transformed_clustering_properties_enrichment.csv",
    aerosolization = "transformed_clustering_properties_aerosolization.csv"
)

# --- Constants ---
INITIAL_K_MAP <- list(
    enrichment = 10,
    aerosolization = 10
)
SUB_CLUSTER_K_RANGE <- 2:5
VARIANCE_QUANTILE_THRESHOLD <- 0.75 # Top 25% variance clusters are candidates for splitting
MIN_CLUSTER_SIZE <- 50 # Sub-cluster only if cluster is large enough

# --- Functions ---

calculate_cluster_variance <- function(data, labels) {
    # Calculate geometric mean of variances of all features within each cluster
    # Or simply the mean standard deviation (trace of covariance / p)

    unique_labels <- unique(labels)
    cluster_stats <- data.frame(cluster = unique_labels, avg_var = NA)

    for (i in 1:nrow(cluster_stats)) {
        lbl <- cluster_stats$cluster[i]
        subset_data <- data[labels == lbl, ]
        if (nrow(subset_data) < 2) {
            cluster_stats$avg_var[i] <- 0
            next
        }

        # Calculate variance for each column and take the mean
        # Since data is standardized (Z-score), we expect var ~ 1 for whole dataset.
        cluster_vars <- apply(subset_data, 2, var)
        cluster_stats$avg_var[i] <- mean(cluster_vars)
    }
    return(cluster_stats)
}

perform_clustering <- function(df, label) {
    message("\nProcessing: ", label)

    # 1. Prepare Data
    df_model <- df %>%
        select(-matches("molecule_id")) %>%
        drop_na()

    ids <- df$molecule_id

    mat <- as.matrix(df_model)
    n_samples <- nrow(mat)

    message("  Data shape: ", n_samples, " rows, ", ncol(mat), " cols")

    # 2. Initial K-Means
    current_k <- INITIAL_K_MAP[[label]]
    if (is.null(current_k)) current_k <- 10

    set.seed(123)
    km_res <- KMeans_rcpp(mat, clusters = current_k, num_init = 20, max_iters = 100)
    initial_labels <- km_res$clusters

    message("  Initial K-means (k=", current_k, ") complete.")

    # 3. Identify Wide Clusters
    # Calculate variance per cluster
    cluster_vars <- calculate_cluster_variance(mat, initial_labels)

    # Use a threshold to decide which to split.
    # Logic: If a cluster's internal variance is significantly higher than the median or specific quantile.
    threshold <- quantile(cluster_vars$avg_var, VARIANCE_QUANTILE_THRESHOLD)
    message("  Variance Threshold (", VARIANCE_QUANTILE_THRESHOLD * 100, "%): ", round(threshold, 3))

    final_labels <- as.character(initial_labels)

    clusters_to_check <- cluster_vars %>%
        filter(avg_var > threshold) %>%
        pull(cluster)

    message("  Checking ", length(clusters_to_check), " clusters for potential sub-clustering...")

    for (clust_id in clusters_to_check) {
        idx <- which(initial_labels == clust_id)
        sub_mat <- mat[idx, , drop = FALSE]

        if (nrow(sub_mat) < MIN_CLUSTER_SIZE) {
            message("    Skipping cluster ", clust_id, " (Too small: ", nrow(sub_mat), ")")
            next
        }

        # Attempt sub-clustering
        best_k <- NA
        best_sil <- -1
        best_sub_labels <- NULL

        # Calculate dist once for silhouette
        # Warning: dist on >30k points is slow. If cluster is large (>5k), maybe subsample for silhouette?
        # But sub-clusters usually aren't that huge if we start with k=10.
        # Let's assume sub-clusters are manageable.
        if (nrow(sub_mat) > 5000) {
            # fast approximate silhouette or just use simple metric?
            # For now, just perform sub-clustering without robust sil check if too huge,
            # OR just take a sample for evaluation.
            sil_idx <- sample(1:nrow(sub_mat), 5000)
            d_sub <- dist(sub_mat[sil_idx, ], method = "euclidean")
        } else {
            sil_idx <- 1:nrow(sub_mat)
            d_sub <- dist(sub_mat, method = "euclidean")
        }

        found_better <- FALSE

        for (k_sub in SUB_CLUSTER_K_RANGE) {
            set.seed(123)
            # Use kmeans (base R) or KMeans_rcpp - use base for simplicity on small sets or consistency
            sub_km <- kmeans(sub_mat, centers = k_sub, nstart = 10)

            # Evaluate Silhouette
            # Note: silhouette returns object.
            ss <- silhouette(sub_km$cluster[sil_idx], d_sub)
            avg_sil <- mean(ss[, 3])

            if (!is.na(avg_sil) && avg_sil > best_sil) {
                best_sil <- avg_sil
                best_k <- k_sub
                best_sub_labels <- sub_km$cluster
                found_better <- TRUE
            }
        }

        # Validation: Is the split actually "good"?
        # A simple Silhouette > 0.1 or 0.2 might be a safety check.
        # But if the user wants to split "wide" clusters, we probably should prioritize reducing variance.
        # Let's apply if we found a valid k.

        if (found_better && best_sil > 0.1) { # 0.05 is a very loose safety check
            message("    Splitting Cluster ", clust_id, " -> ", best_k, " sub-clusters (Sil: ", round(best_sil, 3), ", Size: ", nrow(sub_mat), ")")

            # Update labels
            # New label format: "1_1", "1_2" etc.
            new_labs <- paste0(clust_id, "_", best_sub_labels)
            final_labels[idx] <- new_labs
        } else {
            message("    Cluster ", clust_id, ": No better structure found (Best Sil: ", round(best_sil, 3), "). Keeping as is.")
        }
    }

    # 4. Save Results
    results_df <- data.frame(molecule_id = ids, cluster = final_labels, stringsAsFactors = FALSE)

    # Calculate final stats
    n_final_clusters <- length(unique(results_df$cluster))
    message("  Final cluster count: ", n_final_clusters)

    out_file <- file.path(output_dir, paste0("clustering_results_", label, ".csv"))
    write_csv(results_df, out_file)
    message("  Saved results to: ", out_file)

    return(results_df)
}

# --- Execution ---

for (lbl in names(files_to_process)) {
    f_name <- files_to_process[[lbl]]
    f_path <- file.path(input_dir, f_name)

    if (file.exists(f_path)) {
        tryCatch(
            {
                df <- read_csv(f_path, show_col_types = FALSE)
                perform_clustering(df, lbl)
            },
            error = function(e) {
                message("Error processing ", lbl, ": ", e$message)
            }
        )
    } else {
        message("File not found: ", f_path)
    }
}

message("\nDone.")
