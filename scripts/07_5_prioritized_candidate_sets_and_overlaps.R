#!/usr/bin/env Rscript

# scripts/07_5_prioritized_candidate_sets_and_overlaps.R
# Identify prioritized candidate sets based on top-ranked clusters and analyze their overlaps
# with Atmospheric Stability.
# Generates Figure 6.

packages <- c("tidyverse", "ggVennDiagram", "patchwork", "ggsci")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

library(tidyverse)
library(ggVennDiagram)
library(patchwork)
library(ggsci)

# --- Configuration ---
input_dir_clusters <- "data/processed/05_clustering"
input_dir_ranks <- "reports/07_module_results"
input_dir_refs <- "data/raw/known_enrichment_factors"
input_features <- "data/processed/03_import_epi_features/chemical_merged_features.csv"

output_dir_figs <- "figures"
output_dir_reports <- "reports/07_module_results"

if (!dir.exists(output_dir_figs)) dir.create(output_dir_figs, recursive = TRUE)
if (!dir.exists(output_dir_reports)) dir.create(output_dir_reports, recursive = TRUE)

STABILITY_CUTOFF_HOURS <- 0.5 # 30 minutes

# --- 1. Load Data ---

# Ranked Cluster Stats (from 07_4)
rank_enrich_file <- file.path(input_dir_ranks, "ranked_enrichment_clusters_by_ref.csv")
rank_aero_file <- file.path(input_dir_ranks, "ranked_aerosolization_clusters_by_ref.csv")

if (!file.exists(rank_enrich_file) || !file.exists(rank_aero_file)) {
    stop("Ranked cluster files not found. Please run scripts/07_4_reference_coverage.R first.")
}

rank_e <- read_csv(rank_enrich_file, show_col_types = FALSE)
rank_a <- read_csv(rank_aero_file, show_col_types = FALSE)

# Full Cluster Members (from 05_4/07_2)
clust_enrich_file <- file.path(input_dir_clusters, "clustering_results_enrichment.csv")
clust_aero_file <- file.path(input_dir_clusters, "clustering_results_aerosolization.csv")

clusters_e <- read_csv(clust_enrich_file, show_col_types = FALSE) %>%
    mutate(cluster = as.character(cluster)) %>%
    rename(cluster_e = cluster)

clusters_a <- read_csv(clust_aero_file, show_col_types = FALSE) %>%
    mutate(cluster = as.character(cluster)) %>%
    rename(cluster_a = cluster)

# Reference List (for checking exemplars)
ref_file <- file.path(input_dir_refs, "latest.csv")
refs <- if (file.exists(ref_file)) read_csv(ref_file, show_col_types = FALSE) else NULL

# Features for Stability
if (!file.exists(input_features)) stop("Chemical features file not found.")
df_features <- read_csv(input_features, show_col_types = FALSE) %>%
    select(molecule_id, aop_oh_half_life_hours)

# --- 2. Identify Candidate Sets ---

# Set 1: Enrichment (Top 20 clusters by Median Factor)
top_e_clusters <- rank_e %>%
    arrange(desc(median_factor)) %>%
    head(20) %>%
    pull(cluster_e)

cand_e <- clusters_e %>%
    filter(cluster_e %in% top_e_clusters) %>%
    pull(molecule_id)
set_e <- unique(cand_e)

# Set 2: Aerosolization (Top 20 clusters by Median Factor)
top_a_clusters <- rank_a %>%
    arrange(desc(median_factor)) %>%
    head(20) %>%
    pull(cluster_a)

cand_a <- clusters_a %>%
    filter(cluster_a %in% top_a_clusters) %>%
    pull(molecule_id)
set_a <- unique(cand_a)

# Set 3: Stability (OH Half-life > 0.5h)
cand_s <- df_features %>%
    filter(aop_oh_half_life_hours > STABILITY_CUTOFF_HOURS) %>%
    pull(molecule_id)
set_s <- unique(cand_s)

# --- 3. Intersection Analysis ---

sets_list <- list(
    Enrichment = set_e,
    Aerosolization = set_a,
    Stable = set_s
)

# Intersections
intersect_ea <- intersect(set_e, set_a)
intersect_final <- intersect(intersect_ea, set_s) # E & A & S

message("Candidate Set Summary:")
message("  Enrichment (Top 20): ", length(set_e))
message("  Aerosolization (Top 20): ", length(set_a))
message("  Stable (>30m): ", length(set_s))
message("  E & A Overlap: ", length(intersect_ea))
message("  E & A & Stable (Final Candidates): ", length(intersect_final))

# --- 4. Cluster Contribution Analysis ---
# Analyzing which clusters drive the E & A Intersection (prior to stability filter, to see origins)
# We can also flag which ones result in Stable candidates.

intersection_df <- data.frame(molecule_id = intersect_ea) %>%
    inner_join(clusters_e, by = "molecule_id") %>%
    inner_join(clusters_a, by = "molecule_id") %>%
    left_join(df_features, by = "molecule_id") %>%
    mutate(is_stable = molecule_id %in% set_s)

# Top contributing Enrichment clusters (to E&A intersection)
top_contrib_e <- intersection_df %>%
    count(cluster_e, name = "n_intersection") %>%
    arrange(desc(n_intersection)) %>%
    head(10)

# Top contributing Aerosolization clusters (to E&A intersection)
top_contrib_a <- intersection_df %>%
    count(cluster_a, name = "n_intersection") %>%
    arrange(desc(n_intersection)) %>%
    head(10)

# Save Stats
write_csv(intersection_df %>% filter(is_stable), file.path(output_dir_reports, "prioritized_molecules_triple_intersection.csv"))
write_csv(top_contrib_e, file.path(output_dir_reports, "top_intersection_contributors_enrichment.csv"))
write_csv(top_contrib_a, file.path(output_dir_reports, "top_intersection_contributors_aerosolization.csv"))

# --- 5. Identify Exemplars ---
# Final candidates (E & A & S) that are known references
if (!is.null(refs)) {
    exemplars_known <- intersection_df %>%
        filter(is_stable) %>%
        inner_join(refs, by = "molecule_id")

    write_csv(exemplars_known, file.path(output_dir_reports, "prioritized_exemplars_known_refs_stable.csv"))
    message("  Known references in Triple Intersection: ", nrow(exemplars_known))
}

# --- 6. Generate Figure 6 (improved aesthetics + more informative bars) ---

# Colors (match your scheme)
col_enrich <- "#3C5488FF"
col_aero <- "#E64B35FF"
col_stable <- "#00A087FF"

# A cleaner paper theme (with generous margins to avoid clipping)
theme_paper <- theme_classic(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", size = 15, margin = margin(b = 6)),
        plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
        axis.text = element_text(color = "black"),
        plot.margin = margin(10, 18, 10, 18)
    )

# --------------------------
# Panel A: Venn (fix clipping)
# --------------------------
p_venn <- ggVennDiagram(
    sets_list,
    label = "both", # count + percent
    label_sep = "\n",
    label_alpha = 0, # no label background
    label_size = 4,
    set_size = 5,
    edge_size = 0.9,
    set_color = c(col_enrich, col_aero, col_stable)
) +
    # neutral fill (avoid "Stable looks blue" issue)
    scale_fill_gradient(low = scales::alpha("#F7F7F7", 0.25), high = scales::alpha("#D6D6D6", 0.25)) +
    labs(
        title = "Prioritized candidate overlap",
        subtitle = paste0("Stable defined as OH half-life > ", STABILITY_CUTOFF_HOURS, " h")
    ) +
    coord_equal(clip = "off") +
    theme_void(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", size = 12, margin = margin(t = 25, b = 5)),
        plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
        plot.margin = margin(18, 40, 10, 40), # <- key: prevents truncation of set labels
        legend.position = "none"
    )

# -------------------------------------------------------
# Panel B/C data: cluster contributions (total + stable)
# -------------------------------------------------------
contrib_e <- intersection_df %>%
    count(cluster_e, is_stable, name = "n") %>%
    group_by(cluster_e) %>%
    summarise(
        total = sum(n),
        stable = sum(n[is_stable]),
        .groups = "drop"
    ) %>%
    mutate(stable = ifelse(is.na(stable), 0, stable)) %>%
    arrange(desc(total)) %>%
    slice_head(n = 10) %>%
    mutate(cluster_e = factor(cluster_e, levels = rev(cluster_e)))

contrib_a <- intersection_df %>%
    count(cluster_a, is_stable, name = "n") %>%
    group_by(cluster_a) %>%
    summarise(
        total = sum(n),
        stable = sum(n[is_stable]),
        .groups = "drop"
    ) %>%
    mutate(stable = ifelse(is.na(stable), 0, stable)) %>%
    arrange(desc(total)) %>%
    slice_head(n = 10) %>%
    mutate(cluster_a = factor(cluster_a, levels = rev(cluster_a)))

# Optional: write the improved contributor tables (more useful than only totals)
write_csv(contrib_e, file.path(output_dir_reports, "top_intersection_contributors_enrichment_total_and_stable.csv"))
write_csv(contrib_a, file.path(output_dir_reports, "top_intersection_contributors_aerosolization_total_and_stable.csv"))

# -------------------------------------------------------
# Merged Panel B/C: Faceted Bar Chart
# -------------------------------------------------------
# -------------------------------------------------------
# Merged Panel B/C: Faceted Bar Chart
# -------------------------------------------------------
# Prepare data for stacked bar chart (Stability Class mapping)
contrib_data <- bind_rows(
    contrib_e %>% rename(cluster = cluster_e) %>% mutate(Module = "Enrichment Module"),
    contrib_a %>% rename(cluster = cluster_a) %>% mutate(Module = "Aerosolization Module")
) %>%
    mutate(
        unstable = total - stable,
        Module = factor(Module, levels = c("Enrichment Module", "Aerosolization Module"))
    ) %>%
    select(cluster, Module, stable, unstable, total) %>%
    pivot_longer(cols = c(stable, unstable), names_to = "type", values_to = "count") %>%
    mutate(
        # Define factor levels for stacking order and legend
        fill_label = case_when(
            type == "unstable" ~ "E ∩ A only (not stable)",
            type == "stable" ~ "Stable subset (E ∩ A ∩ S)"
        ),
        fill_label = factor(fill_label, levels = c("E ∩ A only (not stable)", "Stable subset (E ∩ A ∩ S)"))
    )

p_bars <- ggplot(contrib_data, aes(y = cluster, x = count, fill = fill_label)) +
    geom_col(position = "stack", width = 0.72) +
    # Add text label for total count (sum of stacked bars)
    geom_text(
        data = contrib_data %>% group_by(cluster, Module) %>% summarise(total = sum(count), .groups = "drop"),
        aes(x = total, y = cluster, label = scales::comma(total)),
        inherit.aes = FALSE,
        hjust = -0.08, size = 3.6, color = "black"
    ) +
    facet_wrap(~Module, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(
        "E ∩ A only (not stable)" = "#B0B0B099", # Grey, semi-transparent
        "Stable subset (E ∩ A ∩ S)" = scales::alpha(col_stable, 0.85)
    )) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    theme_paper +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    labs(
        x = "Count of molecules",
        y = "Cluster ID"
    )

# --------------------------
# Combine layout
# --------------------------
# --------------------------
# Combine layout
# --------------------------
layout <- (p_venn | p_bars) +
    plot_layout(widths = c(1, 2)) + # Left/Right arrangement, 1:3 ratio
    plot_annotation(tag_levels = "A")

ggsave(file.path(output_dir_figs, "Figure_6.png"),
    layout,
    width = 12, height = 6, dpi = 300, bg = "white"
)

ggsave(file.path(output_dir_figs, "Figure_6.png"),
    layout,
    width = 12, height = 6, bg = "white"
)


cat("\n--- Figure 6 Descriptions ---\n")
cat("Panel A: Prioritized candidate overlap (Enrichment, Aerosolization, Stable).\n")
cat(sprintf("  Stable defined as OH half-life > %s h.\n", STABILITY_CUTOFF_HOURS))
cat("Panel B: Top cluster drivers of the Enrichment-Aerosolization overlap.\n")
cat("  Light bar = all molecules in E-A intersection from this cluster; Green bar = subset that is also stable.\n")
cat("  Faceted by module: Enrichment Module (top) vs Aerosolization Module (bottom).\n")
cat("-----------------------------\n\n")

message("Analysis Complete.")
message("Figure saved to ", file.path(output_dir_figs, "Figure_6.png"))
