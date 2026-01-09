#!/usr/bin/env Rscript

# 06_estimation_labeling.R
# Shiny app to input known enrichment and aerosolization factors.

# Check for required packages
required_packages <- c("shiny", "tidyverse", "DT", "shinyjs", "uuid")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(shiny)
library(tidyverse)
library(DT)
library(shinyjs)
library(uuid)

# --- Config & Data Loading ---

# Paths
features_path <- "../data/processed/03_import_epi_features/chemical_merged_features.csv"
clustering_enrichment_path <- "../data/processed/05_clustering/clustering_results_enrichment.csv"
clustering_aerosolization_path <- "../data/processed/05_clustering/clustering_results_aerosolization.csv"
known_factors_dir <- "../data/raw/known_enrichment_factors"
latest_factors_path <- file.path(known_factors_dir, "latest.csv")

# Ensure directory exists
if (!dir.exists(known_factors_dir)) dir.create(known_factors_dir, recursive = TRUE)

# Initialize latest.csv if missing
if (!file.exists(latest_factors_path)) {
    write_csv(tibble(
        molecule_id = integer(),
        source_type = character(),
        reference_doi = character(),
        remarks = character(),
        timestamp = character(),
        input_type = character(),
        conc_bsw = numeric(),
        conc_sml = numeric(),
        conc_ssa = numeric(),
        direct_factor_value = numeric(),
        direct_factor_type = character()
    ), latest_factors_path)
}

# Load static data
message("Loading static data...")
if (file.exists(features_path)) {
    df_features <- read_csv(features_path, show_col_types = FALSE) %>%
        select(molecule_id, standard_inchi, standard_inchi_key, canonical_smiles, name, iupac_name)
} else {
    stop("Features file not found: ", features_path)
}

load_clustering <- function(path) {
    if (file.exists(path)) {
        read_csv(path, show_col_types = FALSE) %>%
            mutate(cluster = as.character(cluster))
    } else {
        NULL
    }
}

df_clim_enrich <- load_clustering(clustering_enrichment_path)
df_clim_aero <- load_clustering(clustering_aerosolization_path)

# Normalize molecule_id in features for joining
df_features <- df_features %>% mutate(molecule_id = as.integer(molecule_id))

# --- UI ---

ui <- fluidPage(
    useShinyjs(),
    titlePanel("Estimation Labeling: Enrichment & Aerosolization"),
    sidebarLayout(
        sidebarPanel(
            width = 4,
            h4("1. Select Molecule"),
            selectizeInput("molecule_search", "Search Molecule:",
                choices = NULL,
                options = list(
                    placeholder = "Type name, InChIKey, or ID...",
                    maxOptions = 10
                )
            ),
            uiOutput("molecule_info_ui"),
            hr(),
            h4("2. Input Data"),
            selectInput("source_type", "Source Type:",
                choices = c("Field Sample", "Wetlab")
            ),
            textInput("reference_doi", "Reference DOI/Citation:"),
            textAreaInput("remarks", "Remarks (Optional):", rows = 2),
            radioButtons("input_type", "Input Mode:",
                choices = c("Concentrations" = "Concentrations", "Direct Factor" = "Direct Factor"),
                inline = TRUE
            ),

            # Contextual inputs based on mode
            conditionalPanel(
                condition = "input.input_type == 'Concentrations'",
                numericInput("conc_bsw", "Bulk Seawater (BSW):", value = NA),
                numericInput("conc_sml", "Sea Surface Microlayer (SML):", value = NA),
                numericInput("conc_ssa", "Sea Spray Aerosol (SSA):", value = NA)
            ),
            conditionalPanel(
                condition = "input.input_type == 'Direct Factor'",
                selectInput("direct_factor_type", "Factor Type:",
                    choices = c(
                        "Enrichment Factor (SML/BSW)",
                        "Aerosolization Factor (SSA/BSW)",
                        "Aerosolization Factor (SSA/SML)"
                    )
                ),
                numericInput("direct_factor_value", "Factor Value:", value = NA)
            ),
            br(),
            actionButton("save_btn", "Save Record", class = "btn-primary", width = "100%"),
            br(), br(),
            h4("3. Manage Backups"),
            actionButton("show_rollback_btn", "Show Rollback Options", class = "btn-warning")
        ),
        mainPanel(
            width = 8,
            tabsetPanel(
                tabPanel(
                    "Current Data",
                    br(),
                    h4("Registered Factors (latest.csv)"),
                    DTOutput("factors_table")
                ),
                tabPanel(
                    "Cluster Context",
                    br(),
                    h4("Cluster Information for Selected Molecule"),
                    verbatimTextOutput("cluster_details_basic"),
                    hr(),
                    fluidRow(
                        column(
                            6,
                            h4("Enrichment Cluster"),
                            verbatimTextOutput("cluster_id_enrich"),
                            h5("Labeled Data in this Cluster"),
                            DTOutput("cluster_data_enrich")
                        ),
                        column(
                            6,
                            h4("Aerosolization Cluster"),
                            verbatimTextOutput("cluster_id_aero"),
                            h5("Labeled Data in this Cluster"),
                            DTOutput("cluster_data_aero")
                        )
                    )
                ),
                tabPanel(
                    "Labeling Coverage",
                    br(),
                    h4("Cluster Labeling Status"),
                    fluidRow(
                        column(6, h5("Enrichment Clusters"), DTOutput("coverage_enrichment")),
                        column(6, h5("Aerosolization Clusters"), DTOutput("coverage_aerosolization"))
                    )
                )
            )
        )
    )
)

# --- Server ---

server <- function(input, output, session) {
    # Reactive Values
    rv <- reactiveValues(
        data = tibble(), # The current content of latest.csv
        selected_mol_id = NULL
    )

    # 1. Load Data
    load_data <- function() {
        if (file.exists(latest_factors_path)) {
            rv$data <- read_csv(latest_factors_path,
                show_col_types = FALSE,
                col_types = cols(
                    molecule_id = col_integer(),
                    source_type = col_character(),
                    reference_doi = col_character(),
                    remarks = col_character(),
                    timestamp = col_character(),
                    input_type = col_character(),
                    conc_bsw = col_double(),
                    conc_sml = col_double(),
                    conc_ssa = col_double(),
                    direct_factor_value = col_double(),
                    direct_factor_type = col_character()
                )
            )
        }
    }
    load_data()

    # 2. Molecule Search Setup
    # Create a compact search dictionary
    # Combining name, inchikey, id into a search string?
    # For performance with large lists, server-side selectize is best.

    # We'll use the 'name' column as label and 'molecule_id' as value
    # But filtering should happen on multiple columns.

    updateSelectizeInput(session, "molecule_search",
        choices = df_features,
        server = TRUE,
        options = list(
            labelField = "name",
            valueField = "molecule_id",
            searchField = c("name", "standard_inchi_key", "molecule_id", "canonical_smiles"),
            render = I("{
                option: function(item, escape) {
                    return '<div><strong>' + escape(item.name) + '</strong><br>' +
                           '<small>' + escape(item.standard_inchi_key) + ' (ID: ' + escape(item.molecule_id) + ')</small></div>';
                }
            }")
        )
    )

    # 3. Observe Selection
    observeEvent(input$molecule_search, {
        req(input$molecule_search)
        rv$selected_mol_id <- as.integer(input$molecule_search)
    })

    # Render Molecule Info
    output$molecule_info_ui <- renderUI({
        req(rv$selected_mol_id)
        mol <- df_features %>%
            filter(molecule_id == rv$selected_mol_id) %>%
            slice(1)

        # Find clusters
        clust_enrich <- df_clim_enrich %>%
            filter(molecule_id == rv$selected_mol_id) %>%
            pull(cluster)
        clust_aero <- df_clim_aero %>%
            filter(molecule_id == rv$selected_mol_id) %>%
            pull(cluster)

        tagList(
            strong("Name: "), mol$name, br(),
            strong("InChiKey: "), mol$standard_inchi_key, br(),
            strong("Enrichment Cluster: "), if (length(clust_enrich) > 0) clust_enrich else "N/A", br(),
            strong("Aerosolization Cluster: "), if (length(clust_aero) > 0) clust_aero else "N/A"
        )
    })

    # Render Cluster Basic Info
    output$cluster_details_basic <- renderPrint({
        req(rv$selected_mol_id)
        cat("Selected Molecule ID:", rv$selected_mol_id, "\n")
    })

    # Render Enrichment Context
    output$cluster_id_enrich <- renderPrint({
        req(rv$selected_mol_id)
        cid <- df_clim_enrich %>%
            filter(molecule_id == rv$selected_mol_id) %>%
            pull(cluster)
        if (length(cid) > 0) cat("Cluster ID:", cid) else cat("Cluster ID: N/A")
    })

    output$cluster_data_enrich <- renderDT(
        {
            req(rv$selected_mol_id)

            cid <- df_clim_enrich %>%
                filter(molecule_id == rv$selected_mol_id) %>%
                pull(cluster)

            if (length(cid) > 0) {
                # Find all molecules in this cluster
                cluster_mols <- df_clim_enrich %>%
                    filter(cluster == cid) %>%
                    pull(molecule_id)

                # Find data for these molecules
                rv$data %>%
                    filter(molecule_id %in% cluster_mols) %>%
                    arrange(desc(timestamp)) %>%
                    select(molecule_id, everything()) # Ensure ID is first
            } else {
                tibble()
            }
        },
        options = list(pageLength = 5, scrollX = TRUE)
    )

    # Render Aerosolization Context
    output$cluster_id_aero <- renderPrint({
        req(rv$selected_mol_id)
        cid <- df_clim_aero %>%
            filter(molecule_id == rv$selected_mol_id) %>%
            pull(cluster)
        if (length(cid) > 0) cat("Cluster ID:", cid) else cat("Cluster ID: N/A")
    })

    output$cluster_data_aero <- renderDT(
        {
            req(rv$selected_mol_id)

            cid <- df_clim_aero %>%
                filter(molecule_id == rv$selected_mol_id) %>%
                pull(cluster)

            if (length(cid) > 0) {
                # Find all molecules in this cluster
                cluster_mols <- df_clim_aero %>%
                    filter(cluster == cid) %>%
                    pull(molecule_id)

                # Find data for these molecules
                rv$data %>%
                    filter(molecule_id %in% cluster_mols) %>%
                    arrange(desc(timestamp)) %>%
                    select(molecule_id, everything())
            } else {
                tibble()
            }
        },
        options = list(pageLength = 5, scrollX = TRUE)
    )

    # 4. Save Logic
    observeEvent(input$save_btn, {
        req(rv$selected_mol_id)

        # Validation
        if (input$input_type == "Concentrations") {
            if (is.na(input$conc_bsw) && is.na(input$conc_sml) && is.na(input$conc_ssa)) {
                showNotification("Please enter at least one concentration.", type = "error")
                return()
            }
        } else {
            if (is.na(input$direct_factor_value)) {
                showNotification("Please enter a factor value.", type = "error")
                return()
            }
        }

        # Create record
        new_row <- tibble(
            molecule_id = rv$selected_mol_id,
            source_type = input$source_type,
            reference_doi = input$reference_doi,
            remarks = input$remarks,
            timestamp = as.character(Sys.time()),
            input_type = input$input_type,
            conc_bsw = if (input$input_type == "Concentrations") input$conc_bsw else NA,
            conc_sml = if (input$input_type == "Concentrations") input$conc_sml else NA,
            conc_ssa = if (input$input_type == "Concentrations") input$conc_ssa else NA,
            direct_factor_value = if (input$input_type == "Direct Factor") input$direct_factor_value else NA,
            direct_factor_type = if (input$input_type == "Direct Factor") input$direct_factor_type else NA
        )

        # Backup
        backup_file <- file.path(known_factors_dir, paste0("backup_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv"))
        file.copy(latest_factors_path, backup_file)

        # Append and Save
        updated_data <- bind_rows(rv$data, new_row)
        write_csv(updated_data, latest_factors_path)

        # Update State
        rv$data <- updated_data
        showNotification("Record saved successfully!", type = "message")

        # Reset inputs (optional, keeping some might be useful)
        updateTextAreaInput(session, "remarks", value = "")
        updateNumericInput(session, "conc_bsw", value = NA)
        updateNumericInput(session, "conc_sml", value = NA)
        updateNumericInput(session, "conc_ssa", value = NA)
        updateNumericInput(session, "direct_factor_value", value = NA)
    })

    # 5. Rollback Logic
    observeEvent(input$show_rollback_btn, {
        backups <- list.files(known_factors_dir, pattern = "backup_.*\\.csv")
        showModal(modalDialog(
            title = "Rollback Data",
            selectInput("backup_to_restore", "Select Backup to Restore:", choices = sort(backups, decreasing = TRUE)),
            footer = tagList(
                modalButton("Cancel"),
                actionButton("confirm_rollback", "Restore Selected Backup", class = "btn-danger")
            )
        ))
    })

    observeEvent(input$confirm_rollback, {
        req(input$backup_to_restore)
        backup_path <- file.path(known_factors_dir, input$backup_to_restore)

        if (file.exists(backup_path)) {
            # Make a backup of CURRENT before rolling back? Yes.
            pre_rollback_backup <- file.path(known_factors_dir, paste0("backup_PREROLLBACK_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv"))
            file.copy(latest_factors_path, pre_rollback_backup)

            # Restore
            file.copy(backup_path, latest_factors_path, overwrite = TRUE)
            load_data() # Reload to reactive val

            showNotification(paste("Restored from", input$backup_to_restore), type = "message")
            removeModal()
        }
    })

    # Render Main Table
    output$factors_table <- renderDT({
        rv$data %>% arrange(desc(timestamp))
    }, )

    # 6. Coverage Logic
    output$coverage_enrichment <- renderDT(
        {
            req(rv$data)

            # Count total molecules per cluster
            total_counts <- df_clim_enrich %>%
                count(cluster, name = "total_molecules")

            # Count labeled molecules per cluster
            # Filter rv$data for unique molecules that have at least one record
            labeled_mols <- unique(rv$data$molecule_id)

            labeled_counts <- df_clim_enrich %>%
                filter(molecule_id %in% labeled_mols) %>%
                count(cluster, name = "labeled_molecules")

            # Join
            coverage <- total_counts %>%
                left_join(labeled_counts, by = "cluster") %>%
                replace_na(list(labeled_molecules = 0)) %>%
                mutate(
                    pct_covered = round(100 * labeled_molecules / total_molecules, 1)
                ) %>%
                arrange(desc(pct_covered), desc(total_molecules))

            coverage
        },
        options = list(pageLength = 10, scrollX = TRUE)
    )

    output$coverage_aerosolization <- renderDT(
        {
            req(rv$data)

            total_counts <- df_clim_aero %>%
                count(cluster, name = "total_molecules")

            labeled_mols <- unique(rv$data$molecule_id)

            labeled_counts <- df_clim_aero %>%
                filter(molecule_id %in% labeled_mols) %>%
                count(cluster, name = "labeled_molecules")

            coverage <- total_counts %>%
                left_join(labeled_counts, by = "cluster") %>%
                replace_na(list(labeled_molecules = 0)) %>%
                mutate(
                    pct_covered = round(100 * labeled_molecules / total_molecules, 1)
                ) %>%
                arrange(desc(pct_covered), desc(total_molecules))

            coverage
        },
        options = list(pageLength = 10, scrollX = TRUE)
    )
}

shinyApp(ui = ui, server = server)
