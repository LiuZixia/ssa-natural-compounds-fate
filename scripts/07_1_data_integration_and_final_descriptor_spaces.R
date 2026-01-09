#!/usr/bin/env Rscript

# scripts/07_1_data_integration_and_final_descriptor_spaces.R
# Generate summaries (N0-N2, PE/PA, NE/NA, etc.) and Refactored Figures 1 and S1

suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(grid)
    library(gridExtra)
    library(reshape2)
    library(ggsci)
    library(scales)
})


# --- Configuration ---
# Paths
paths <- list(
    raw_molecules = "data/raw/molecules.csv",
    filtered_features = "data/processed/02_prepare_features/chemical_features.csv",
    merged_features = "data/processed/03_import_epi_features/chemical_merged_features.csv",
    enrichment_props = "data/processed/05_clustering/clustering_properties_enrichment.csv",
    aerosol_props = "data/processed/05_clustering/clustering_properties_aerosolization.csv",
    transformed_enrich = "data/processed/05_clustering/transformed_clustering_properties_enrichment.csv",
    transformed_aerosol = "data/processed/05_clustering/transformed_clustering_properties_aerosolization.csv"
)

output_plot_dir <- "figures"
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)

# --- Helper Functions ---
getPEMA <- function(path) {
    if (!file.exists(path)) {
        return(list(N = NA, P = NA, m = NA, rMax = NA, df = NULL, miss = NULL))
    }

    df <- read_csv(path, show_col_types = FALSE)
    N <- nrow(df)
    cols <- setdiff(names(df), "molecule_id")
    # Identify non-flag cols for correlation
    # We assume flags end in "_missing".
    feat_cols <- cols[!grepl("_missing$", cols)]
    P <- length(feat_cols)

    # Missingness
    miss_pct <- colMeans(is.na(df[cols])) * 100
    m_median <- median(miss_pct)

    # Correlation (rMax)
    num_cols <- feat_cols[sapply(df[feat_cols], is.numeric)]

    rMax <- 0
    if (length(num_cols) > 1) {
        cor_mat <- cor(df[num_cols], use = "pairwise.complete.obs", method = "spearman")
        diag(cor_mat) <- 0
        rMax <- max(abs(cor_mat), na.rm = TRUE)
    }

    list(N = N, P = P, m = m_median, rMax = rMax, df = df, miss = miss_pct, cor_mat = cor_mat, num_cols = num_cols)
}

plot_heatmap <- function(cor_mat, title_text, highlight_thr = 0.90, show_legend = TRUE) {
    if (is.null(cor_mat)) {
        return(ggplot() +
            labs(title = "No Data") +
            theme_void())
    }

    # order by hierarchical clustering
    dd <- as.dist((1 - cor_mat) / 2)
    hc <- hclust(dd)
    cmat <- cor_mat[hc$order, hc$order]

    # upper triangle only
    cmat[lower.tri(cmat, diag = TRUE)] <- NA
    melted <- melt(cmat, na.rm = TRUE) %>%
        mutate(
            highlight = abs(value) >= highlight_thr,
            Var1 = factor(Var1, levels = colnames(cmat)),
            Var2 = factor(Var2, levels = rownames(cmat))
        )

    p <- ggplot(melted, aes(Var1, Var2, fill = value)) +
        geom_tile(color = "grey95", linewidth = 0.35) +
        geom_text(
            data = subset(melted, highlight),
            aes(label = "×"),
            color = "black", size = 4, alpha = 0.65
        ) +
        scale_fill_gradient2(
            low = COL_BLUE, mid = "white", high = COL_RED,
            midpoint = 0, limits = c(-1, 1),
            name = expression(Spearman ~ rho)
        ) +
        coord_fixed() +
        scale_x_discrete(labels = nice_label) +
        scale_y_discrete(labels = nice_label) +
        labs(title = title_text, x = NULL, y = NULL) +
        theme_nature(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.key.height = unit(18, "pt"),
            legend.key.width = unit(8, "pt")
        )

    if (!show_legend) p <- p + theme(legend.position = "none")
    p
}

plot_distributions <- function(df, cols, title_text, fill_color = COL_BLUE) {
    if (is.null(df) || length(cols) == 0) {
        return(ggplot() +
            labs(title = "No Data") +
            theme_void())
    }

    df_long <- df %>%
        select(all_of(cols)) %>%
        pivot_longer(everything(), names_to = "Descriptor", values_to = "Value") %>%
        mutate(Descriptor = factor(Descriptor, levels = rev(cols)))

    # trim extremes for display only (keeps readable axes)
    xr <- quantile(df_long$Value, probs = c(0.01, 0.99), na.rm = TRUE)

    ggplot(df_long, aes(x = Value, y = Descriptor)) +
        geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
        geom_violin(fill = fill_color, color = NA, alpha = 0.55, scale = "width", trim = TRUE) +
        geom_boxplot(
            width = 0.13, outlier.shape = NA, fill = "white",
            color = COL_GREY, linewidth = 0.35
        ) +
        coord_cartesian(xlim = xr) +
        scale_y_discrete(labels = nice_label) +
        labs(title = title_text, x = "Transformed Value (z-score)", y = NULL) +
        theme_nature(base_size = 12) +
        theme(axis.text.y = element_text(size = 10))
}

# Nature-ish (NPG-like) palette
COL_BLUE <- "#3C5488FF"
COL_RED <- "#E64B35FF"
COL_GREY <- "grey30"

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

nice_label <- function(x) {
    map <- c(
        "fractioncsp3" = "Fsp3",
        "alogp" = "ALogP",
        "hydrogen_bond_donors" = "HBD",
        "hydrogen_bond_acceptors" = "HBA",
        "topological_polar_surface_area" = "TPSA",
        "rotatable_bond_count" = "RotB",
        "total_atom_count" = "Atoms",
        "molecular_weight" = "MW",
        "aromatic_rings_count" = "AromRings",
        "number_of_minimal_rings" = "MinRings",
        "vp_est" = "VP (est.)"
    )
    ifelse(x %in% names(map), unname(map[x]), x)
}

# --- 1. Data Counts (N0 - N2) ---
if (file.exists(paths$raw_molecules)) {
    N0 <- nrow(read_csv(paths$raw_molecules, show_col_types = FALSE))
} else {
    N0 <- NA
}

if (file.exists(paths$filtered_features)) {
    N1 <- nrow(read_csv(paths$filtered_features, show_col_types = FALSE))
} else {
    N1 <- NA
}

if (file.exists(paths$merged_features)) {
    d2 <- read_csv(paths$merged_features, show_col_types = FALSE)
    N2 <- nrow(d2)
    P0 <- ncol(d2) - 1
} else {
    N2 <- NA
    P0 <- NA
}


# --- 2. Module Analysis ---
enrich_stats <- getPEMA(paths$enrichment_props)
NE <- enrich_stats$N
PE <- enrich_stats$P
mE <- enrich_stats$m
rMax <- enrich_stats$rMax

aerosol_stats <- getPEMA(paths$aerosol_props)
NA_count <- aerosol_stats$N
PA <- aerosol_stats$P
mA <- aerosol_stats$m
rMaxA <- aerosol_stats$rMax

message("Generating combined Figure (A–D)...")

# --- Figure panels A/B: correlation heatmaps ---
pA <- NULL
pB <- NULL
if (!is.null(enrich_stats$cor_mat)) {
    pA <- plot_heatmap(enrich_stats$cor_mat, "A. Enrichment Module", show_legend = FALSE)
}
if (!is.null(aerosol_stats$cor_mat)) {
    pB <- plot_heatmap(aerosol_stats$cor_mat, "B. Aerosolization Module", show_legend = TRUE)
}

# --- Figure panels C/D: descriptor distributions ---
df_trans_e <- if (file.exists(paths$transformed_enrich)) read_csv(paths$transformed_enrich, show_col_types = FALSE) else NULL
df_trans_a <- if (file.exists(paths$transformed_aerosol)) read_csv(paths$transformed_aerosol, show_col_types = FALSE) else NULL

pC <- NULL
pD <- NULL

if (!is.null(df_trans_e)) {
    is_binary <- sapply(df_trans_e, function(x) {
        if (is.logical(x)) {
            return(TRUE)
        }
        ux <- unique(x)
        ux <- ux[!is.na(ux)]
        length(ux) <= 2 && all(ux %in% c(0, 1))
    })
    cont_cols <- setdiff(names(df_trans_e), c("molecule_id", names(which(is_binary))))
    pC <- plot_distributions(df_trans_e, cont_cols, "C. Enrichment descriptors (continuous)", fill_color = COL_BLUE)
}

if (!is.null(df_trans_a)) {
    is_binary_a <- sapply(df_trans_a, function(x) {
        if (is.logical(x)) {
            return(TRUE)
        }
        ux <- unique(x)
        ux <- ux[!is.na(ux)]
        length(ux) <= 2 && all(ux %in% c(0, 1))
    })
    cont_cols_a <- setdiff(names(df_trans_a), c("molecule_id", names(which(is_binary_a))))
    pD <- plot_distributions(df_trans_a, cont_cols_a, "D. Aerosolization descriptors (continuous)", fill_color = COL_RED)
}

# Safety: if any panel missing, replace with blank
blank_panel <- ggplot() +
    theme_void()
if (is.null(pA)) pA <- blank_panel
if (is.null(pB)) pB <- blank_panel
if (is.null(pC)) pC <- blank_panel
if (is.null(pD)) pD <- blank_panel

# Arrange 2x2
fig_ABCD <- cowplot::plot_grid(
    pA, pB,
    pC, pD,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    rel_heights = c(1.05, 0.95)
)

ggsave(file.path(output_plot_dir, "Figure_1_ABCD.png"),
    fig_ABCD,
    width = 14, height = 10, dpi = 300, bg = "white"
)



message("\nFigures saved in ", output_plot_dir)

# --- Summary Output ---
# --- Final Summary Text for Manuscript ---
message("\n--- Variable Values ---")

vals <- list(
    "N0 (Initial MarineNP candidate compounds)" = N0,
    "N1 (After valid SMILES + exclude problematic IDs)" = N1,
    "N2 (Integrated dataset after EPI Suite merge)" = N2,
    "P0 (Total descriptors in integrated dataset)" = P0,
    "PE (Descriptors retained: enrichment module)" = PE,
    "NE (Compounds retained: enrichment module)" = NE,
    "PA (Descriptors retained: aerosolization module)" = PA,
    "NA_count (Compounds retained: aerosolization module)" = NA_count,
    "mE (Median missingness %: enrichment module)" = mE,
    "mA (Median missingness %: aerosolization module)" = mA,
    "rMax (Max |Spearman rho| after pruning: enrichment)" = rMax,
    "rMaxA (Max |Spearman rho| after pruning: aerosolization)" = rMaxA
)

for (nm in names(vals)) {
    message(sprintf("%s: %s", nm, as.character(vals[[nm]])))
}
