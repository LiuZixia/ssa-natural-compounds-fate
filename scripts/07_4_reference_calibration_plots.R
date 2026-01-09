#!/usr/bin/env Rscript

# scripts/07_4_reference_calibration_Figure5.R
# Figure 5 (3.3): Literature-anchored calibration
#   A) PCA overlays (Aerosolization + Enrichment)
#   B) Reference factor distributions by cluster (points + median + bootstrap 95% CI)
#   C) Enrichment LOO rank stability (optional)
# Output: figures/Figure_5_ABC.png + CSVs in reports/07_module_results/

suppressPackageStartupMessages({
  pkgs <- c("tidyverse", "ggrepel", "ggbeeswarm", "scales", "cowplot")
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
  library(tidyverse)
  library(ggrepel)
  library(ggbeeswarm)
  library(scales)
  library(cowplot)
})

set.seed(123)

# -------------------------
# Config
# -------------------------
input_dir_factors  <- "data/raw/known_enrichment_factors"
input_dir_clusters <- "data/processed/05_clustering"
output_dir_figs    <- "figures"
output_dir_reports <- "reports/07_module_results"

ref_file     <- file.path(input_dir_factors,  "latest.csv")
clust_e_file <- file.path(input_dir_clusters, "clustering_results_enrichment.csv")
clust_a_file <- file.path(input_dir_clusters, "clustering_results_aerosolization.csv")
props_e_file <- file.path(input_dir_clusters, "clustering_properties_enrichment.csv")
props_a_file <- file.path(input_dir_clusters, "clustering_properties_aerosolization.csv")

figure_out <- file.path(output_dir_figs, "Figure_5_ABC.png")

# Visual controls
pca_bg_n          <- 25000  # background sample (still used for bin2d)
pca_bins          <- 140    # bin resolution (higher = finer density)
pca_zoom_quantile <- c(0.01, 0.99)  # zoom limits for PC axes
boot_B            <- 2000

include_panel_C <- TRUE  # set FALSE to remove panel C

dir.create(output_dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_reports, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(ref_file), file.exists(clust_e_file), file.exists(clust_a_file),
          file.exists(props_e_file), file.exists(props_a_file))

# Nature-ish accents (match PCA vibe)
COL_ENRICH <- "#3C5488FF"
COL_AERO   <- muted("red")  # softer than #E64B35FF, closer to PCA palette feel

# -------------------------
# Helpers
# -------------------------
safe_log10 <- function(x) ifelse(is.na(x) | x <= 0, NA_real_, log10(x))

boot_median_ci <- function(x, B = 2000, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(tibble(med = NA_real_, lo = NA_real_, hi = NA_real_, n = 0))
  if (length(x) == 1) return(tibble(med = x, lo = x, hi = x, n = 1))
  meds <- replicate(B, median(sample(x, replace = TRUE)))
  a <- (1 - conf) / 2
  tibble(
    med = median(x),
    lo  = unname(quantile(meds, probs = a, na.rm = TRUE)),
    hi  = unname(quantile(meds, probs = 1 - a, na.rm = TRUE)),
    n   = length(x)
  )
}

# -------------------------
# Load data
# -------------------------
refs <- read_csv(ref_file, show_col_types = FALSE) %>%
  mutate(molecule_id = as.character(molecule_id))

clusters_e <- read_csv(clust_e_file, show_col_types = FALSE) %>%
  transmute(molecule_id = as.character(molecule_id), cluster_e = as.character(cluster))

clusters_a <- read_csv(clust_a_file, show_col_types = FALSE) %>%
  transmute(molecule_id = as.character(molecule_id), cluster_a = as.character(cluster))

props_e <- read_csv(props_e_file, show_col_types = FALSE) %>%
  mutate(molecule_id = as.character(molecule_id))

props_a <- read_csv(props_a_file, show_col_types = FALSE) %>%
  mutate(molecule_id = as.character(molecule_id))

# -------------------------
# Process reference factors
# -------------------------
refs_proc <- refs %>%
  mutate(
    EF = case_when(
      input_type == "Concentrations" & !is.na(conc_sml) & !is.na(conc_bsw) & conc_bsw > 0 ~ conc_sml / conc_bsw,
      str_detect(input_type, regex("Direct", TRUE)) & str_detect(coalesce(direct_factor_type, ""), regex("Enrichment", TRUE)) ~ direct_factor_value,
      TRUE ~ NA_real_
    ),
    AF = case_when(
      input_type == "Concentrations" & !is.na(conc_ssa) & !is.na(conc_bsw) & conc_bsw > 0 ~ conc_ssa / conc_bsw,
      input_type == "Concentrations" & !is.na(conc_ssa) & !is.na(conc_sml) & conc_sml > 0 ~ conc_ssa / conc_sml,
      str_detect(input_type, regex("Direct", TRUE)) & str_detect(coalesce(direct_factor_type, ""), regex("Aerosol", TRUE)) ~ direct_factor_value,
      TRUE ~ NA_real_
    ),
    logEF = safe_log10(EF),
    logAF = safe_log10(AF),
    has_EF = is.finite(logEF),
    has_AF = is.finite(logAF),
    evidence_source = case_when(
      str_detect(coalesce(source_type, ""), regex("field", TRUE)) ~ "Field",
      str_detect(coalesce(source_type, ""), regex("wet|lab", TRUE)) ~ "Wetlab",
      TRUE ~ "Other"
    ),
    evidence_input = case_when(
      str_detect(coalesce(input_type, ""), regex("Concen", TRUE)) ~ "Concentrations",
      str_detect(coalesce(input_type, ""), regex("Direct", TRUE)) ~ "Direct factor",
      TRUE ~ "Other"
    )
  )

write_csv(
  tibble(
    total_records = nrow(refs_proc),
    n_with_EF = sum(refs_proc$has_EF, na.rm = TRUE),
    n_with_AF = sum(refs_proc$has_AF, na.rm = TRUE)
  ),
  file.path(output_dir_reports, "reference_set_summary.csv")
)

refs_e <- refs_proc %>% filter(has_EF) %>% inner_join(clusters_e, by = "molecule_id")
refs_a <- refs_proc %>% filter(has_AF) %>% inner_join(clusters_a, by = "molecule_id")

write_csv(refs_e, file.path(output_dir_reports, "references_mapped_enrichment.csv"))
write_csv(refs_a, file.path(output_dir_reports, "references_mapped_aerosolization.csv"))

message("Reference Set Summary:")
message("  Total Records: ", nrow(refs_proc))
message("  With Enrichment Evidence (finite log10 EF): ", sum(refs_proc$has_EF, na.rm = TRUE))
message("  With Aerosolization Evidence (finite log10 AF): ", sum(refs_proc$has_AF, na.rm = TRUE))

# -------------------------
# Stats for factor plots
# -------------------------
cluster_factor_stats <- function(ref_df, cluster_col, log_col, module) {
  ref_df %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(ci = list(boot_median_ci(.data[[log_col]], B = boot_B)), .groups = "drop") %>%
    unnest(ci) %>%
    rename(cluster = !!cluster_col) %>%
    mutate(module = module) %>%
    arrange(desc(med))
}

stats_e <- cluster_factor_stats(refs_e, "cluster_e", "logEF", "Enrichment")
stats_a <- cluster_factor_stats(refs_a, "cluster_a", "logAF", "Aerosolization")

write_csv(bind_rows(stats_a, stats_e), file.path(output_dir_reports, "cluster_factor_bootstrap_summary.csv"))

# Export ranked clusters for 07_5 usage
stats_e %>%
  rename(cluster_e = cluster, median_factor = med, n_refs = n) %>%
  write_csv(file.path(output_dir_reports, "ranked_enrichment_clusters_by_ref.csv"))

stats_a %>%
  rename(cluster_a = cluster, median_factor = med, n_refs = n) %>%
  write_csv(file.path(output_dir_reports, "ranked_aerosolization_clusters_by_ref.csv"))

# -------------------------
# Figure 5A: PCA overlay (density background)
# -------------------------
pca_overlay <- function(props_df, clusters_df, ref_df, cluster_col, log_col, module_label) {
  
  d <- props_df %>%
    inner_join(clusters_df, by = "molecule_id") %>%
    mutate(cluster = as.character(.data[[cluster_col]]))
  
  num_cols <- d %>% select(where(is.numeric)) %>% names()
  stopifnot(length(num_cols) >= 2)
  
  X <- d %>% select(all_of(num_cols)) %>% as.matrix()
  ok <- complete.cases(X)
  d_ok <- d[ok, , drop = FALSE]
  X_ok <- X[ok, , drop = FALSE]
  
  pca <- prcomp(X_ok, center = TRUE, scale. = TRUE)
  scores <- as_tibble(pca$x[, 1:2]) %>% setNames(c("PC1", "PC2"))
  d_plot <- bind_cols(d_ok %>% select(molecule_id, cluster), scores)
  
  # sample background for speed
  bg <- if (nrow(d_plot) > pca_bg_n) d_plot %>% slice_sample(n = pca_bg_n) else d_plot
  
  # zoom to avoid outlier compression
  xlim <- quantile(bg$PC1, probs = pca_zoom_quantile, na.rm = TRUE)
  ylim <- quantile(bg$PC2, probs = pca_zoom_quantile, na.rm = TRUE)
  
  ref_overlay <- ref_df %>%
    mutate(cluster = as.character(.data[[cluster_col]]), logF = .data[[log_col]]) %>%
    inner_join(d_plot %>% select(molecule_id, PC1, PC2), by = "molecule_id")
  
  # centroid labels (clean + stable)
  lab_df <- ref_overlay %>%
    group_by(cluster) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  
  # export overlay table
  out_csv <- file.path(output_dir_reports, paste0("pca_ref_overlay_", tolower(module_label), ".csv"))
  write_csv(ref_overlay, out_csv)
  
  ggplot() +
    geom_bin2d(
      data = bg,
      aes(PC1, PC2, fill = after_stat(count)),
      bins = pca_bins
    ) +
    scale_fill_gradient(low = "white", high = "grey65", guide = "none") +
    geom_point(
      data = ref_overlay,
      aes(PC1, PC2, fill = logF, shape = evidence_source),
      size = 3.2, alpha = 0.95, color = "black", stroke = 0.35
    ) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(PC1, PC2, label = cluster),
      size = 3.0, max.overlaps = 50,
      box.padding = 0.25, point.padding = 0.15
    ) +
    scale_shape_manual(values = c(Field = 21, Wetlab = 24, Other = 22)) +
    scale_fill_gradient2(
      low = muted("blue"), mid = "grey90", high = muted("red"),
      midpoint = 0, name = "log10(factor)",
      oob = squish
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = "PC1", y = "PC2", shape = "Evidence source") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.box = "horizontal")
}

pA_a <- pca_overlay(props_a, clusters_a, refs_a, "cluster_a", "logAF", "Aerosolization") + ggtitle("Aerosolization Module")
pA_e <- pca_overlay(props_e, clusters_e, refs_e, "cluster_e", "logEF", "Enrichment") + ggtitle("Enrichment Module")

# Extract shared legend from Enrichment (Enrichment has Field + Wetlab, Aerosolization might only have Field)
legend_A <- get_legend(pA_e + theme(legend.box.margin = margin(0, 0, 0, 0)))

pA_row <- plot_grid(
  pA_e + theme(legend.position = "none"), # Enrichment left
  pA_a + theme(legend.position = "none"), # Aerosolization right
  nrow = 1,
  labels = c("A", ""),
  label_size = 14
)

panel_A <- plot_grid(pA_row, legend_A, ncol = 1, rel_heights = c(1, 0.1))

# -------------------------
# Figure 5B: Factor distributions (bootstrap CI)
# -------------------------
factor_plot <- function(ref_df, stats_df, cluster_col, log_col, module_label, dot_col) {
  pts <- ref_df %>%
    mutate(cluster = as.character(.data[[cluster_col]]), logF = .data[[log_col]]) %>%
    filter(is.finite(logF))
  
  st <- stats_df %>% mutate(cluster = as.character(cluster))
  
  ord <- st %>% arrange(med) %>% pull(cluster)
  pts$cluster <- factor(pts$cluster, levels = ord)
  st$cluster  <- factor(st$cluster,  levels = ord)
  
  ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
    ggbeeswarm::geom_quasirandom(
      data = pts,
      aes(x = logF, y = cluster, shape = evidence_input),
      width = 0.25, alpha = 0.75, size = 2.1, color = "grey25"
    ) +
    geom_errorbarh(
      data = st,
      aes(y = cluster, xmin = lo, xmax = hi),
      height = 0.18, linewidth = 1.1, alpha = 0.9
    ) +
    geom_point(
      data = st,
      aes(y = cluster, x = med, size = n),
      shape = 21, fill = dot_col, color = "black", stroke = 0.6
    ) +
    geom_text(
      data = st,
      aes(y = cluster, x = hi + 0.05, label = paste0("n=", n)),
      hjust = 0, size = 3.0
    ) +
    scale_size(range = c(2.2, 7.5)) +
    scale_shape_manual(values = c("Concentrations" = 16, "Direct factor" = 17, "Other" = 15)) +
    labs(
      x = "log10(reference factor)",
      y = "Cluster ID",
      shape = "Evidence input",
      size = "n refs"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.margin = margin(5.5, 60, 5.5, 5.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    coord_cartesian(clip = "off")
}

pB_a <- factor_plot(refs_a, stats_a, "cluster_a", "logAF", "Aerosolization", COL_ENRICH) + ggtitle("Aerosolization Module") # Use same blue color
pB_e <- factor_plot(refs_e, stats_e, "cluster_e", "logEF", "Enrichment",     COL_ENRICH) + ggtitle("Enrichment Module")

legend_B <- get_legend(pB_a + theme(legend.box.margin = margin(0, 0, 0, 0)))

pB_row <- plot_grid(
  pB_e + theme(legend.position = "none"), # Enrichment left
  pB_a + theme(legend.position = "none"), # Aerosolization right
  nrow = 1,
  labels = c("B", ""),
  label_size = 14
)

panel_B <- plot_grid(pB_row, legend_B, ncol = 1, rel_heights = c(1, 0.1))

# -------------------------
# Figure 5C: Enrichment LOO rank stability (prettified + centered)
# -------------------------
loo_rank_intervals <- function(ref_df, cluster_col, log_col) {
  df <- ref_df %>%
    mutate(cluster = as.character(.data[[cluster_col]]), logF = .data[[log_col]]) %>%
    filter(is.finite(logF))
  
  if (nrow(df) < 4 || n_distinct(df$cluster) < 2) return(tibble())
  
  base <- df %>%
    group_by(cluster) %>%
    summarise(med = median(logF), n = n(), .groups = "drop")
  
  loo <- map_dfr(seq_len(nrow(df)), function(i) {
    d2 <- df[-i, , drop = FALSE]
    d2 %>%
      group_by(cluster) %>%
      summarise(med = median(logF), n = n(), .groups = "drop") %>%
      mutate(rank = dense_rank(desc(med)))
  })
  
  loo %>%
    group_by(cluster) %>%
    summarise(rank_med = median(rank), rank_min = min(rank), rank_max = max(rank), .groups = "drop") %>%
    left_join(base, by = "cluster") %>%
    arrange(rank_med)
}

loo_e <- loo_rank_intervals(refs_e, "cluster_e", "logEF")
loo_a <- loo_rank_intervals(refs_a, "cluster_a", "logAF") # Aerosolization LOO

write_csv(loo_e, file.path(output_dir_reports, "loo_rank_intervals_enrichment.csv"))
write_csv(loo_a, file.path(output_dir_reports, "loo_rank_intervals_aerosolization.csv"))

rank_stability_plot <- function(loo_df, title_label, color_val) {
  if (nrow(loo_df) == 0) {
    ggplot() + theme_void() + labs(title = paste0(title_label, " (insufficient data)"))
  } else {
    d <- loo_df %>% mutate(cluster = fct_reorder(cluster, -rank_med))
    max_rank <- max(d$rank_max, na.rm = TRUE)
    
    ggplot(d, aes(y = cluster)) +
      geom_vline(xintercept = seq(1, max_rank), linewidth = 0.3, alpha = 0.15) +
      geom_segment(aes(x = rank_min, xend = rank_max, yend = cluster),
                   linewidth = 2.4, alpha = 0.35, lineend = "round", color = "grey30") +
      geom_point(aes(x = rank_med, size = n),
                 shape = 21, fill = color_val, color = "black", stroke = 0.6) +
      scale_x_reverse(breaks = seq(1, max_rank)) +
      scale_size(range = c(2.6, 7.5)) +
      labs(
        title = title_label,
        x = "Rank (1 = highest median factor)",
        y = "Cluster ID",
        size = "n refs"
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0)
      )
  }
}

pC_e <- rank_stability_plot(loo_e, "Enrichment Module", COL_AERO)
pC_a <- rank_stability_plot(loo_a, "Aerosolization Module", COL_AERO)

panel_C <- if (!include_panel_C) NULL else {
  # Shared legend for size
  legend_C <- get_legend(pC_e + theme(legend.box.margin = margin(0,0,0,0)))
  
  pC_row <- plot_grid(
    pC_e + theme(legend.position = "none"),
    pC_a + theme(legend.position = "none"),
    nrow = 1,
    labels = c("C", ""),
    label_size = 14
  )
  
  plot_grid(pC_row, legend_C, ncol = 1, rel_heights = c(1, 0.1))
}

# -------------------------
# Assemble + Save
# -------------------------

# Visual separator
separator_line <- ggdraw() + draw_line(x = c(0.05, 0.95), y = c(0.5, 0.5), color = "grey80", linewidth = 1)

# Descriptions for Console
cat("\n--- Figure 5 Descriptions ---\n")
cat("Panel A: References overlaid on descriptor space (PCA).\n")
cat(sprintf("  Background = binned candidate density (n approx %d sampled); labels = reference centroids.\n", pca_bg_n))
cat("Panel B: Reference factor distributions by cluster.\n")
cat("  Points = references; dot = cluster median; whiskers = bootstrap 95% CI (median).\n")
cat("Panel C: Rank stability under leave-one-out (LOO).\n")
cat("  Dot = median rank; line = min-max rank across LOO (1 = highest median factor). Red = Enrichment & Aerosolization.\n")
cat("-----------------------------\n\n")

# Combine with separators
# We want: A, sep, B, sep, C
# Heights need adjustment.
figure5 <- if (include_panel_C) {
  plot_grid(
    panel_A,
    separator_line,
    panel_B,
    separator_line,
    panel_C,
    ncol = 1,
    rel_heights = c(1.15, 0.05, 1.15, 0.05, 0.8),
    align = "v"
  )
} else {
  plot_grid(
    panel_A,
    separator_line,
    panel_B,
    ncol = 1,
    rel_heights = c(1.15, 0.05, 1.15),
    align = "v"
  )
}

ggsave(figure_out, figure5, width = 12, height = if (include_panel_C) 10 else 9, dpi = 300, bg = "white")

message("Done. Figure saved to: ", normalizePath(figure_out))
message("CSVs written to: ", normalizePath(output_dir_reports))
