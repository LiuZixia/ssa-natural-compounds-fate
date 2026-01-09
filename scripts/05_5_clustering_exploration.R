#!/usr/bin/env Rscript

# 05_5_clustering_exploration.R
# Shiny app to explore clustering results, search molecules, and analyze clusters.

# Check for required packages
required_packages <- c("shiny", "tidyverse", "plotly", "shinyjs", "DT")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(shiny)
library(tidyverse)
library(plotly)
library(shinyjs)
library(DT)

# --- Data Loading ---

# Define paths
data_dir_base <- "data/processed"
if (!dir.exists(data_dir_base)) data_dir_base <- "../data/processed"

features_path <- file.path(data_dir_base, "03_import_epi_features/chemical_merged_features.csv")
clustering_dir <- file.path(data_dir_base, "05_clustering")

clustering_enrichment_path <- file.path(clustering_dir, "clustering_results_enrichment.csv")
clustering_aerosolization_path <- file.path(clustering_dir, "clustering_results_aerosolization.csv")

# Validation
if (!file.exists(features_path)) stop("Features file not found: ", features_path)
if (!file.exists(clustering_enrichment_path)) warning("Enrichment clustering results not found.")
if (!file.exists(clustering_aerosolization_path)) warning("Aerosolization clustering results not found.")

message("Loading data...")
df_features <- read_csv(features_path, show_col_types = FALSE)

# Load clustering results
load_clustering <- function(path) {
    if (file.exists(path)) {
        read_csv(path, show_col_types = FALSE) %>%
            mutate(cluster = as.character(cluster)) # Ensure cluster is character
    } else {
        NULL
    }
}

df_clust_enrich <- load_clustering(clustering_enrichment_path)
df_clust_aero <- load_clustering(clustering_aerosolization_path)

# Join clustering to features (Left join to keep all features, but typically we analyze clustered items)
# slightly different approach: We create two main frames, or one frame with two cluster columns?
# One frame with two cluster columns is cleaner.

df_full <- df_features

if (!is.null(df_clust_enrich)) {
    df_full <- df_full %>%
        left_join(df_clust_enrich %>% rename(cluster_enrichment = cluster), by = "molecule_id")
} else {
    df_full$cluster_enrichment <- NA
}

if (!is.null(df_clust_aero)) {
    df_full <- df_full %>%
        left_join(df_clust_aero %>% rename(cluster_aerosolization = cluster), by = "molecule_id")
} else {
    df_full$cluster_aerosolization <- NA
}

# --- Variables for Analysis ---
vars_numeric <- names(df_full)[sapply(df_full, is.numeric)]
vars_enrichment_view <- c(
    "topological_polar_surface_area", "rotatable_bond_count",
    "hydrogen_bond_acceptors", "hydrogen_bond_donors",
    "aromatic_rings_count", "fractioncsp3", "van_der_walls_volume", "alogp", "molecular_weight"
)
vars_aerosolization_view <- c("bp_est", "mp_est", "vp_est", "alogp", "molecular_weight")

# Ensure they exist
vars_enrichment_view <- intersect(vars_enrichment_view, names(df_full))
vars_aerosolization_view <- intersect(vars_aerosolization_view, names(df_full))


# --- UI ---
ui <- fluidPage(
    useShinyjs(),
    titlePanel("Clustering Exploration & Analysis"),
    sidebarLayout(
        sidebarPanel(
            width = 3,
            h4("Search & Navigation"),
            selectizeInput("molecule_search", "Search Molecule (ID or Key):",
                choices = NULL, # Server-side load
                multiple = FALSE,
                options = list(placeholder = "Type to search...")
            ),
            hr(),
            h4("Clustering Context"),
            selectInput("cluster_type", "Select Clustering Type:",
                choices = c("Enrichment", "Aerosolization")
            ),

            # Dynamic UI for Cluster Selection (updates when type changes or molecule found)
            uiOutput("cluster_selector_ui"),
            hr(),
            h4("Visualization Settings"),
            selectInput("plot_x_var", "X Variable:", choices = vars_numeric, selected = "alogp"),
            selectInput("plot_y_var", "Y Variable:", choices = vars_numeric, selected = "molecular_weight"),
            checkboxInput("log_scale", "Log Scale (X & Y)", value = FALSE)
        ),
        mainPanel(
            tabsetPanel(
                id = "main_tabs",
                tabPanel(
                    "Cluster Overview",
                    br(),
                    h4(textOutput("cluster_title")),
                    fluidRow(
                        column(6, h5("Feature Distribution (Boxplots)"), plotlyOutput("cluster_boxplot", height = "500px")),
                        column(6, h5("Cluster Membership"), DTOutput("cluster_members_table"))
                    )
                ),
                tabPanel(
                    "Global Visualization",
                    br(),
                    h4("Scatter Plot (Colored by Cluster)"),
                    plotlyOutput("global_scatter", height = "600px"),
                    p("Click points to inspect molecules. Current cluster is highlighted.")
                ),
                tabPanel(
                    "Molecule Details",
                    br(),
                    verbatimTextOutput("molecule_info")
                )
            )
        )
    )
)

# --- Server ---
server <- function(input, output, session) {
    # 1. Update Molecule Search (Server-side)
    # Prepare search data
    # Check which columns are available
    search_cols <- c("molecule_id", "standard_inchi_key", "name", "iupac_name")
    existing_cols <- intersect(search_cols, names(df_full))

    # Create a display label: "Name (ID)" or just "ID"
    # We construct a dataframe for selectize
    search_df <- df_full %>%
        select(all_of(existing_cols)) %>%
        distinct(molecule_id, .keep_all = TRUE)

    # Create a clean label
    if ("name" %in% existing_cols) {
        search_df$display_label <- ifelse(!is.na(search_df$name) & search_df$name != "",
            paste0(search_df$name, " (", search_df$molecule_id, ")"),
            search_df$molecule_id
        )
    } else {
        search_df$display_label <- search_df$molecule_id
    }

    # Configure Selectize
    updateSelectizeInput(session, "molecule_search",
        choices = search_df,
        server = TRUE,
        options = list(
            placeholder = "Type name, InChiKey, or ID...",
            searchField = existing_cols,
            labelField = "display_label",
            valueField = "molecule_id",
            render = I("{
                option: function(item, escape) {
                    var label = escape(item.display_label);
                    var extra = '';
                    if (item.name && item.name !== item.display_label) extra += '<br><small>Name: ' + escape(item.name) + '</small>';
                    if (item.standard_inchi_key) extra += '<br><small>Key: ' + escape(item.standard_inchi_key) + '</small>';
                    return '<div><strong>' + label + '</strong>' + extra + '</div>';
                }
            }")
        )
    )

    # Reactive state for currently selected molecule (from search or click)
    selected_molecule <- reactiveVal(NULL)

    # Reactive state for current cluster (can be changed by dropdown OR by search)
    # using reactiveValues to allow bi-directional update pattern
    state <- reactiveValues(
        cluster_id = NULL
    )

    # Helper to get current cluster column depending on type
    get_cluster_col <- reactive({
        if (input$cluster_type == "Enrichment") "cluster_enrichment" else "cluster_aerosolization"
    })

    # 2. Render Cluster Selector
    output$cluster_selector_ui <- renderUI({
        col <- get_cluster_col()
        clusters <- sort(unique(df_full[[col]]))
        clusters <- clusters[!is.na(clusters)]

        selectInput("selected_cluster", "Select Cluster:",
            choices = clusters,
            selected = state$cluster_id
        ) # Uses state to set initial/updated value
    })

    # Observer: When user manually changes input$selected_cluster, update state
    observeEvent(input$selected_cluster,
        {
            state$cluster_id <- input$selected_cluster
        },
        ignoreInit = TRUE
    )

    # Observer: When Clustering Context changes, reset or try to maintain?
    observeEvent(input$cluster_type, {
        # If a molecule is selected, jump to its cluster in the NEW type
        mol_id <- selected_molecule()
        if (!is.null(mol_id)) {
            row <- df_full[df_full$molecule_id == mol_id, ]
            col <- get_cluster_col()
            new_clust <- row[[col]]
            if (!is.na(new_clust) && length(new_clust) > 0) {
                state$cluster_id <- new_clust
            }
        } else {
            # Default to first cluster if none selected
            col <- get_cluster_col()
            clusters <- sort(unique(df_full[[col]]))
            if (length(clusters) > 0) state$cluster_id <- clusters[1]
        }
    })

    # 3. Handle Molecule Search / Selection
    observeEvent(input$molecule_search, {
        req(input$molecule_search)
        mol_id <- input$molecule_search
        selected_molecule(mol_id)

        # Switch to Molecule Tab? Or just update state?
        # Requirement: "jump to the cluster it belongs to"

        row <- df_full[df_full$molecule_id == mol_id, ]
        if (nrow(row) == 0) {
            return()
        }

        col <- get_cluster_col()
        clust <- row[[col]]

        if (!is.na(clust)) {
            state$cluster_id <- clust
            showNotification(paste("Jumped to Cluster", clust), type = "message")
        } else {
            showNotification("Molecule not assigned to any cluster in this view.", type = "warning")
        }

        updateTabsetPanel(session, "main_tabs", selected = "Molecule Details")
    })

    # Handle Plot Click
    observeEvent(event_data("plotly_click", source = "global_scatter"), {
        d <- event_data("plotly_click", source = "global_scatter")
        if (is.null(d)) {
            return()
        }

        mol_id <- d$customdata
        if (!is.null(mol_id)) {
            selected_molecule(mol_id)
            # Also update search box to reflect this
            updateSelectizeInput(session, "molecule_search", selected = mol_id)

            # Update cluster if needed
            row <- df_full[df_full$molecule_id == mol_id, ]
            col <- get_cluster_col()
            clust <- row[[col]]
            if (!is.na(clust)) state$cluster_id <- clust
        }
    })

    # 4. Outputs

    output$cluster_title <- renderText({
        paste(input$cluster_type, "Cluster:", state$cluster_id)
    })

    output$molecule_info <- renderPrint({
        req(selected_molecule())
        row <- df_full %>%
            filter(molecule_id == selected_molecule()) %>%
            as.list()
        print(row)
    })

    # Boxplots: Show distribution of key features for THIS cluster vs Background
    output$cluster_boxplot <- renderPlotly({
        req(state$cluster_id)
        col <- get_cluster_col()

        # Identify features to plot
        features_to_plot <- if (input$cluster_type == "Enrichment") vars_enrichment_view else vars_aerosolization_view

        # Prepare data: "Current Cluster" vs "Others"
        plot_df <- df_full %>%
            filter(!is.na(.data[[col]])) %>% # Remove unclustered for this view
            mutate(Group = ifelse(.data[[col]] == state$cluster_id, "Selected Cluster", "Background")) %>%
            select(Group, all_of(features_to_plot)) %>%
            pivot_longer(cols = -Group, names_to = "Feature", values_to = "Value")

        # Normalize features? Maybe just raw values.
        # Since scales differ wildly (MW vs fractions), standard boxplot might be hard to read.
        # Faceting is better.

        p <- ggplot(plot_df, aes(x = Group, y = Value, fill = Group)) +
            geom_boxplot() +
            facet_wrap(~Feature, scales = "free_y") +
            theme_minimal() +
            theme(legend.position = "none") +
            labs(x = NULL, y = NULL)

        ggplotly(p)
    })

    output$cluster_members_table <- renderDT({
        req(state$cluster_id)
        col <- get_cluster_col()

        df_sub <- df_full %>%
            filter(.data[[col]] == state$cluster_id)

        # Select columns to display - Check if name exists
        disp_cols <- c("molecule_id", "standard_inchi_key", "molecular_weight", "alogp")
        if ("name" %in% names(df_full)) disp_cols <- c("molecule_id", "name", "standard_inchi_key", "molecular_weight", "alogp")

        df_sub <- df_sub %>% select(any_of(disp_cols))

        datatable(df_sub, selection = "single", options = list(pageLength = 10))
    })

    # Table click -> Select molecule
    observeEvent(input$cluster_members_table_rows_selected, {
        idx <- input$cluster_members_table_rows_selected
        req(idx)
        # Re-fetch data to get ID
        col <- get_cluster_col()
        df_sub <- df_full %>%
            filter(.data[[col]] == state$cluster_id)

        mol_id <- df_sub$molecule_id[idx]

        selected_molecule(mol_id)
        updateSelectizeInput(session, "molecule_search", selected = mol_id)
        updateTabsetPanel(session, "main_tabs", selected = "Molecule Details")
    })

    output$global_scatter <- renderPlotly({
        req(input$plot_x_var, input$plot_y_var)
        col <- get_cluster_col()

        # Colors: Highlight selected cluster, gray out others? Or just discrete colors?
        # Discrete colors for all clusters is nice, but if many clusters (e.g. 20), might be messy.
        # Let's highlight specific cluster.

        plot_data <- df_full %>%
            filter(!is.na(.data[[col]])) %>%
            mutate(
                IsSelected = ifelse(.data[[col]] == state$cluster_id, "Selected", "Other"),
                ClusterLabel = as.character(.data[[col]])
            )

        # If log scale
        if (input$log_scale) {
            plot_data[[input$plot_x_var]] <- log10(plot_data[[input$plot_x_var]] + 1e-6)
            plot_data[[input$plot_y_var]] <- log10(plot_data[[input$plot_y_var]] + 1e-6)
        }

        # Create colors
        # We want Selected to be Red/Distinct, Others to be Gray or lighter colors

        p <- plot_ly(plot_data,
            x = ~ .data[[input$plot_x_var]],
            y = ~ .data[[input$plot_y_var]],
            color = ~ClusterLabel,
            colors = "Paired", # Or manual
            type = "scatter",
            mode = "markers",
            text = ~ paste("ID:", molecule_id, "<br>Cluster:", ClusterLabel),
            customdata = ~molecule_id,
            marker = list(opacity = 0.7, size = 6),
            source = "global_scatter"
        ) %>%
            layout(
                title = paste(input$plot_y_var, "vs", input$plot_x_var),
                xaxis = list(title = input$plot_x_var),
                yaxis = list(title = input$plot_y_var),
                dragmode = "select"
            )

        p
    })
}

shinyApp(ui, server)
