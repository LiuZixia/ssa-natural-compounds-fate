#!/usr/bin/env Rscript

# scripts/07_3_atmospheric_stability.R
# Analyze atmospheric stability (OH and O3 half-lives)
# - Missingness stats
# - Figure 4: Density plots with cutoff
# - Stability stats per cluster

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggsci)
    library(cowplot)
    library(grid)
    library(gridExtra)
})

# --- Configuration ---
STABILITY_CUTOFF_HOURS <- 0.5

input_features <- "data/processed/03_import_epi_features/chemical_merged_features.csv"
input_clusters_enrich <- "data/processed/05_clustering/clustering_results_enrichment.csv"
input_clusters_aerosol <- "data/processed/05_clustering/clustering_results_aerosolization.csv"

output_report_dir <- "reports/07_module_results"
output_plot_dir <- "figures"


# --- Theme ---
theme_nature <- function(base_size = 12) {
    theme_classic(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            panel.grid = element_blank(),
            legend.title = element_text(face = "bold"),
            legend.text = element_text(color = "black")
        )
}

# --- 1. Load Data ---
message("Loading data...")
df_features <- read_csv(input_features, show_col_types = FALSE)

# Check cols
target_cols <- c("aop_oh_half_life_hours", "aop_ozone_half_life_hours")
if (!all(target_cols %in% names(df_features))) {
    stop("Missing half-life columns in features file.")
}

df_stab <- df_features %>%
    select(molecule_id, all_of(target_cols))

# --- 2. Missingness Analysis ---
n_total <- nrow(df_stab)
usage_stats <- df_stab %>%
    summarise(
        n_missing_oh = sum(is.na(aop_oh_half_life_hours)),
        pct_missing_oh = mean(is.na(aop_oh_half_life_hours)) * 100,
        n_missing_oz = sum(is.na(aop_ozone_half_life_hours)),
        pct_missing_oz = mean(is.na(aop_ozone_half_life_hours)) * 100
    )

message("\nMissingness Stats:")
print(usage_stats)

write_csv(usage_stats, file.path(output_report_dir, "stability_missingness.csv"))

# --- 3. Figure 4: Density Plots ---
# Prepare data for plotting (long format)
# Filter unreasonable values if necessary (e.g. negative?) Half-lives should be > 0.
df_plot <- df_stab %>%
    pivot_longer(cols = all_of(target_cols), names_to = "Type", values_to = "HalfLife") %>%
    drop_na(HalfLife) %>%
    mutate(
        Type = ifelse(Type == "aop_oh_half_life_hours", "OH Half-Life", "Ozone Half-Life"),
        HalfLife = HalfLife + 1e-4 # Small offset prevents log(0) issues
    )

# Plot
# Using ggsci colors and log10 scale with readable labels
cutoff_label <- if (STABILITY_CUTOFF_HOURS < 1) {
    paste0(STABILITY_CUTOFF_HOURS * 60, " min Cutoff")
} else {
    paste0(STABILITY_CUTOFF_HOURS, "h Cutoff")
}

p_dens <- ggplot(df_plot, aes(x = HalfLife, fill = Type)) +
    geom_density(alpha = 0.6, color = NA) +
    geom_vline(xintercept = STABILITY_CUTOFF_HOURS, linetype = "dashed", color = "black", linewidth = 0.8) +
    annotate("text",
        x = STABILITY_CUTOFF_HOURS * 1.5, y = 0.05, label = cutoff_label,
        angle = 90, vjust = 1, hjust = 0
    ) +
    scale_x_continuous(
        trans = "log2",
        breaks = c(1 / 60, 10 / 60, 1, 6, 24, 168),
        labels = c("1 min", "10 min", "1 hr", "6 hr", "1 d", ">1 wk"),
        limits = c(1 / 60, 168),
        oob = scales::squish
    ) +
    scale_fill_manual(values = c("#E64B35FF", "#3C5488FF")) +
    theme_nature(base_size = 14) +
    labs(
        title = "Atmospheric Stability Distributions",
        x = "Half-Life",
        y = "Density"
    ) +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(file.path(output_plot_dir, "Figure_4.png"), p_dens, width = 8, height = 6, dpi = 300)
message("Saved Figure 4 to figures/Figure_4.png")


# --- 4. Cluster Integration ---

calculate_stability_stats <- function(cluster_file, label) {
    if (!file.exists(cluster_file)) {
        warning("Cluster file not found: ", cluster_file)
        return(NULL)
    }

    df_clust <- read_csv(cluster_file, show_col_types = FALSE) %>%
        mutate(cluster = as.character(cluster))

    # Merge
    merged <- inner_join(df_clust, df_stab, by = "molecule_id")

    # Calculate stats per cluster
    # Stable if OH > cutoff AND Ozone > cutoff?
    # Or just based on OH? Usually global transport relies more on OH.
    # User asked for "stability outputs", usually implying persistence.
    # Let's calculate for both individually and combined.
    # "stable" generally means resistant to degradation.
    # We will compute:
    # - % OH Stable (> cutoff)
    # - % Ozone Stable (> cutoff)
    # - % Both Stable

    stats <- merged %>%
        group_by(cluster) %>%
        summarise(
            n_total = n(),
            n_oh_stable = sum(aop_oh_half_life_hours > STABILITY_CUTOFF_HOURS, na.rm = TRUE),
            pct_oh_stable = mean(aop_oh_half_life_hours > STABILITY_CUTOFF_HOURS, na.rm = TRUE) * 100,
            n_oz_stable = sum(aop_ozone_half_life_hours > STABILITY_CUTOFF_HOURS, na.rm = TRUE),
            pct_oz_stable = mean(aop_ozone_half_life_hours > STABILITY_CUTOFF_HOURS, na.rm = TRUE) * 100,
            # For combined, handle NAs carefully: if either is missing, sum is NA?
            # Assuming if one is missing, we can't be sure, OR ignore?
            # Let's count available data only for percentages?
            # Simple approach: count valid comparisons.

            median_oh_hl = median(aop_oh_half_life_hours, na.rm = TRUE),
            median_oz_hl = median(aop_ozone_half_life_hours, na.rm = TRUE)
        ) %>%
        arrange(cluster) # Or numeric sort

    # Sort clusters numerically if possible
    stats$cluster <- factor(stats$cluster, levels = stringr::str_sort(unique(stats$cluster), numeric = TRUE))
    stats <- stats %>% arrange(cluster)

    outfile <- file.path(output_report_dir, paste0("stability_by_cluster_", label, ".csv"))
    write_csv(stats, outfile)
    message("Saved stability stats for ", label, " to ", outfile)

    return(stats)
}

# Enrichment
s_enrich <- calculate_stability_stats(input_clusters_enrich, "enrichment")

# Aerosolization
s_aerosol <- calculate_stability_stats(input_clusters_aerosol, "aerosolization")


# --- 5. Generate Panels B & C (Bar Plots) ---

plot_cluster_stability <- function(stats, title) {
    if (is.null(stats)) {
        return(NULL)
    }

    # Reshape for plotting
    df_long <- stats %>%
        select(cluster, pct_oh_stable, pct_oz_stable) %>%
        pivot_longer(cols = c(pct_oh_stable, pct_oz_stable), names_to = "Type", values_to = "PctStable") %>%
        mutate(
            Type = ifelse(Type == "pct_oh_stable", "OH Stable", "Ozone Stable")
        )

    p <- ggplot(df_long, aes(x = cluster, y = PctStable, fill = Type)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
        theme_nature(base_size = 12) +
        scale_fill_manual(values = c("#E64B35FF", "#3C5488FF")) +
        scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
        labs(
            title = title,
            x = "Cluster ID",
            y = "Fraction Stable (%)",
            fill = NULL
        ) +
        theme(
            legend.position = "top",
            legend.box.margin = margin(0, 0, 0, 0),
            legend.margin = margin(0, 0, 0, 0),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    return(p)
}

p_b <- plot_cluster_stability(s_enrich, "B. Enrichment Clusters")
p_c <- plot_cluster_stability(s_aerosol, "C. Aerosolization Clusters")

# --- 6. Combine and Save ---

if (!is.null(p_dens) && !is.null(p_b) && !is.null(p_c)) {
    # Layout: A on top (full width), B and C on bottom (half width each)
    # Use cowplot or gridExtra. gridExtra is loaded.

    # Remove legends from B and C if they are redundant?
    # They use the same colors as A, but labels are "OH Stable" vs "OH Half-Life".
    # Colors logic: Red = OH, Blue = Ozone in A?
    # Let's check A:
    # Type = ifelse(Type == "aop_oh_half_life_hours", "OH Half-Life", "Ozone Half-Life")
    # Plot A fill: Red/Blue.
    # Plot B/C fill: Red/Blue.
    # It seems consistent.

    # Combining
    lay <- rbind(
        c(1, 1),
        c(2, 3)
    )

    g <- arrangeGrob(p_dens, p_b, p_c, layout_matrix = lay, heights = c(1, 1))

    outfile_fig4 <- file.path(output_plot_dir, "Figure_4.png")
    ggsave(outfile_fig4, g, width = 12, height = 10, dpi = 300)
    message("Saved combined Figure 4 to ", outfile_fig4)
} else {
    warning("Could not generate all panels for Figure 4.")
}

message("\nDone.")
