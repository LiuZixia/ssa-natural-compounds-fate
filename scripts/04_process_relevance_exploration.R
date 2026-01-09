#!/usr/bin/env Rscript

# Check for required packages and install if missing
required_packages <- c("shiny", "tidyverse", "plotly", "shinyjs", "DT")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(shiny)
library(tidyverse)
library(plotly)
library(shinyjs)
library(DT)

# --- Data Loading ---
# Define path relative to the script location or project root
# Assuming script is run from project root, or we try to locate it
data_path <- "data/processed/03_import_epi_features/chemical_merged_features.csv"
if (!file.exists(data_path)) {
  # Fallback if running from scripts directory
  data_path <- "../data/processed/03_import_epi_features/chemical_merged_features.csv"
}

if (!file.exists(data_path)) {
  stop("Could not find 'chemical_merged_features.csv'. Please run this script from the project root or 'scripts' directory.")
}

message("Loading data from: ", data_path)
df_full <- read_csv(data_path, show_col_types = FALSE)

# Pre-calculate variable types for efficiency
# Identify numeric vs categorical columns
col_types <- sapply(df_full, function(x) {
  if (is.numeric(x)) "numeric" else "categorical"
})

# --- UI ---
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Chemical Relevance Exploration"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,

      
      hr(),
      h4("Dynamic Filters"),
      p("Select variables below to create filters for them."),
      selectizeInput("filter_vars", "Add Filters", 
                     choices = names(df_full), 
                     multiple = TRUE,
                     options = list(placeholder = 'Type to search variables...')),
      
      # This area will hold the generated filter UI elements
      uiOutput("dynamic_filters_ui"),
      
      hr(),
      actionButton("reset_filters", "Reset All Filters", icon = icon("refresh")),
      
      # Plot Tab Controls
      conditionalPanel(
        condition = "input.main_tabs == 'Plot'",
        hr(),
        h4("Plot Settings"),
        selectInput("x_var", "Select Variable (X-Axis)", 
                    choices = names(df_full), 
                    selected = "molecular_weight"),
        selectInput("x_transform", "Transformation",
                    choices = c("None", "Log10", "Sqrt"),
                    selected = "None"),
        checkboxInput("remove_outliers", "Remove Extreme Outliers (Top 0.1%)", value = FALSE)
      ),
      
      # PCA Tab Controls
      conditionalPanel(
        condition = "input.main_tabs == 'PCA'",
        hr(),
        h4("PCA Settings"),
        selectInput("pca_group", "Select Feature Set for PCA:", 
                    choices = c("Enrichment", "Aerosolization")),
        uiOutput("pca_vars_ui")
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "main_tabs",
        tabPanel("Plot", 
                 br(),
                 h4("Visual Exploration"),

                 plotlyOutput("dist_plot", height = "600px"),
                 br(),
                 h5("Instructions:"),
                 tags$ul(
                   tags$li("Click on a point to open the Coconut database entry for that chemical."),
                   tags$li("Hover to see details.")
                 )
        ),
        tabPanel("PCA",
                 br(),
                 # PCA Controls

                 hr(),
                 plotlyOutput("pca_plot", height = "600px"),
                 br(),
                 p("Note: PCA is performed on the selected subset of variables (including 'Both'). Missing values are removed.")
        ),
        tabPanel("Data Table",
                 br(),
                 DTOutput("data_table")
        )
      )
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  # Reactive values to track filter states if needed (or rely on input directly)
  # using a named list to store filter inputs dynamically might be tricky with just 'input'
  # so we rely on the specific naming convention: "filter_{col_name}"
  
  # Variable sets
  vars_both <- c("total_atom_count", "heavy_atom_count", "molecular_weight",
                 "exact_molecular_weight", "alogp", "formal_charge",
                 "number_of_minimal_rings", "contains_sugar",
                 "contains_ring_sugars", "contains_linear_sugars")

  vars_enrichment <- c("topological_polar_surface_area", "rotatable_bond_count",
                       "hydrogen_bond_acceptors", "hydrogen_bond_donors",
                       "aromatic_rings_count", "fractioncsp3", "van_der_walls_volume")

  vars_aerosolization <- c("bp_est", "mp_est", "vp_est")

  # 0. PCA Vars UI
  output$pca_vars_ui <- renderUI({
    req(input$pca_group)
    
    # Define available variables based on group
    available_vars <- if (input$pca_group == "Enrichment") {
      unique(c(vars_enrichment, vars_both))
    } else {
      unique(c(vars_aerosolization, vars_both))
    }
    
    # Filter to those actually in df_full
    available_vars <- intersect(available_vars, names(df_full))
    
    selectizeInput("pca_vars", "Select Variables to Include:",
                   choices = available_vars,
                   selected = available_vars,
                   multiple = TRUE,
                   options = list(plugins = list('remove_button')))
  })

  # 1. Generate Filter UI
  output$dynamic_filters_ui <- renderUI({
    req(input$filter_vars)
    
    # Create a list of UI elements
    lapply(input$filter_vars, function(var) {
      curr_type <- col_types[[var]]
      data_vec <- df_full[[var]]
      ns_id <- paste0("filter_", var) # Name space for input
      
      div(
        class = "panel panel-default",
        style = "padding: 10px; margin-bottom: 10px; background: #f9f9f9;",
        if (curr_type == "numeric") {
          rng <- range(data_vec, na.rm = TRUE)
          sliderInput(ns_id, label = var, min = floor(rng[1]), max = ceiling(rng[2]), 
                      value = rng)
        } else {
          # Categorical: limit choices if too many
          choices <- unique(data_vec)
          choices <- choices[!is.na(choices)]
          if (length(choices) > 50) {
             # If too many, maybe use selectize or just show top 50?
             # For exploration, standard selectize is ok.
             selectizeInput(ns_id, label = var, choices = choices, multiple = TRUE)
          } else {
             checkboxGroupInput(ns_id, label = var, choices = choices, selected = choices)
          } 
        }
      )
    })
  })
  
  # 2. Filter Data
  filtered_data <- reactive({
    res <- df_full
    
    # Iterate through active filter variables
    # We only filter if the input for that filter actually exists
    req_vars <- input$filter_vars
    
    for (var in req_vars) {
      inp_id <- paste0("filter_", var)
      val <- input[[inp_id]]
      
      # If filter hasn't rendered yet, skip
      if (is.null(val)) next 
      
      curr_type <- col_types[[var]]
      
      if (curr_type == "numeric") {
        # Slider returns length 2 vector (min, max)
        res <- res %>% filter(.data[[var]] >= val[1] & .data[[var]] <= val[2])
      } else {
        # Checkbox/Select returns vector of selected values
        # If nothing selected, typically means 'none', but for filters usually 'show these'
        # If user unchecks all, result is empty.
        if (length(val) > 0) {
          res <- res %>% filter(.data[[var]] %in% val)
        } else {
          # If categorical filter active but nothing selected -> show nothing? 
          # Or treat as "all"? Usually in shiny checkboxGroup, empty = empty.
          res <- res[0, ] 
        }
      }
    }
    
    res
  })
  
  # 3. Render Plot
  output$dist_plot <- renderPlotly({
    req(input$x_var)
    data <- filtered_data()
    
    if (nrow(data) == 0) {
      return(plot_ly() %>% layout(title = "No data selected"))
    }
    
    x_val <- data[[input$x_var]]
    is_num <- is.numeric(x_val)
    
    if (is_num) {
      # --- Numeric: Scatter/Jitter Plot ---
      
      # Handle Outliers
      if (input$remove_outliers) {
        q999 <- quantile(x_val, 0.999, na.rm = TRUE)
        keep_mask <- x_val <= q999
        keep_mask[is.na(keep_mask)] <- FALSE
        data <- data[keep_mask, ]
        x_val <- data[[input$x_var]]
      }
      
      # Handle Transformation
      plot_x <- x_val
      x_label <- input$x_var
      
      if (nrow(data) > 0) {
        if (input$x_transform == "Log10") {
          plot_x <- log10(x_val)
          x_label <- paste0("Log10(", input$x_var, ")")
        } else if (input$x_transform == "Sqrt") {
          plot_x <- sqrt(x_val)
          x_label <- paste0("Sqrt(", input$x_var, ")")
        }
      }
      
      has_link <- "standard_inchi_key" %in% names(data)
      
      p <- plot_ly(data, x = plot_x, y = ~rnorm(nrow(data)), 
                   type = 'scatter', mode = 'markers',
                   text = ~paste0("ID: ", molecule_id, 
                                  if("name" %in% names(data)) paste0("<br>Name: ", name) else "",
                                  if(has_link) paste0("<br>Key: ", standard_inchi_key) else "",
                                  "<br>Value: ", signif(x_val, 4)),
                   customdata = if(has_link) ~standard_inchi_key else ~molecule_id, 
                   marker = list(opacity = 0.6, size = 8)) %>%
        layout(
          title = paste("Distribution of", input$x_var),
          xaxis = list(title = x_label),
          yaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
          hovermode = "closest"
        )
      return(p)
      
    } else {
      # --- Categorical: Pie Chart ---
      # Calculate counts
      counts <- as.data.frame(table(x_val, useNA = "ifany"))
      names(counts) <- c("Category", "Count")
      
      p <- plot_ly(counts, labels = ~Category, values = ~Count, type = 'pie',
                   textinfo = 'label+percent',
                   hoverinfo = 'label+value+percent') %>%
        layout(title = paste("Distribution of", input$x_var))
      return(p)
    }
  })
  
  # 4. Handle Click
  observeEvent(event_data("plotly_click"), {
    click_data <- event_data("plotly_click")
    if (is.null(click_data)) return()
    
    # key is in customdata
    key <- click_data$customdata
    
    if (!is.null(key) && !is.na(key) && key != "") {
      # Construct URL
      url <- paste0("https://coconut.naturalproducts.net/search?q=", key)
      
      # Use shinyjs to open window
      runjs(sprintf("window.open('%s', '_blank');", url))
    }
  })
  
  # 5. Data Table
  output$data_table <- renderDT({
    filtered_data()
  }, options = list(
    pageLength = 8, 
    scrollX = TRUE,
    columnDefs = list(list(
      targets = "_all",
      render = JS(
        "function(data, type, row, meta) {",
        "  if (type === 'display' && data != null && typeof data === 'string' && data.length > 50) {",
        "    return '<span class=\"truncated\">' + data.substr(0, 50) + '...</span>' +",
        "           '<span class=\"full-text\" style=\"display:none;\">' + data + '</span>' +",
        "           ' <span class=\"toggle-btn\" style=\"cursor:pointer; color:blue; font-size:0.8em;\">[+]</span>';",
        "  }",
        "  return data;",
        "}"
      )
    ))
  ), callback = JS(
    "table.on('click', '.toggle-btn', function() {",
    "  var btn = $(this);",
    "  var truncated = btn.siblings('.truncated');",
    "  var full = btn.siblings('.full-text');",
    "  if (full.is(':visible')) {",
    "    full.hide();",
    "    truncated.show();",
    "    btn.text('[+]');",
    "  } else {",
    "    truncated.hide();",
    "    full.show();",
    "    btn.text('[-]');",
    "  }",
    "});"
  ))
  

  # Helper to prepare data for PCA
  get_pca_data <- function(df, group) {
    # This helper is deprecated as we moved logic to renderPlotly with input$pca_vars
    # But we can keep it for reference or remove it. 
    # For cleanliness, we should remove it if not used, but let's just ignore.
    return(NULL)
  }

  output$pca_plot <- renderPlotly({
    req(input$pca_vars) # Only proceed if variables are selected
    
    target_vars <- input$pca_vars
    
    if (length(target_vars) < 2) return(plot_ly() %>% layout(title = "Please select at least 2 variables"))

    # Filter mainly to keep track of IDs
    # We need to ensure we valid_vars check again or handle missing cols gracefully
    valid_vars <- intersect(target_vars, names(df_full))
    
    if (length(valid_vars) < 2) return(plot_ly() %>% layout(title = "Variables not found in data"))

    # Create a temp dataframe with ID and variables
    df_work <- df_full %>% 
      select(molecule_id, standard_inchi_key, all_of(valid_vars))
    
    # Convert known booleans in place
    for (v in valid_vars) {
       if (is.logical(df_work[[v]])) df_work[[v]] <- as.numeric(df_work[[v]])
    }
    
    # Keep only complete cases for the variables used in PCA
    df_work <- df_work %>% drop_na(all_of(valid_vars))
    
    if (nrow(df_work) < 5) return(plot_ly() %>% layout(title = "Not enough data for PCA"))
    
    # Compute PCA
    pca_res <- prcomp(df_work[, valid_vars], scale. = TRUE, center = TRUE)
    
    # Extract Scores
    df_pca <- data.frame(
      PC1 = pca_res$x[, 1],
      PC2 = pca_res$x[, 2],
      molecule_id = df_work$molecule_id,
      standard_inchi_key = df_work$standard_inchi_key
    )
    
    # Variance explained
    var_exp <- summary(pca_res)$importance[2, 1:2]
    
    # --- Biplot Arrows (Loadings) ---
    loadings <- as.data.frame(pca_res$rotation)
    
    # Scaling factor to make arrows visible over the data points
    # Heuristic: Scale arrows to 80% of the max range of scores
    max_score <- max(abs(c(df_pca$PC1, df_pca$PC2)))
    max_load <- max(abs(c(loadings$PC1, loadings$PC2)))
    scale_factor <- (max_score / max_load) * 0.8
    
    loadings$PC1_scaled <- loadings$PC1 * scale_factor
    loadings$PC2_scaled <- loadings$PC2 * scale_factor
    loadings$VarName <- rownames(loadings)
    
    # Create segment lines for arrows (from 0,0 to x,y)
    # Interleave with NAs to draw unconnected lines in one trace
    arrow_x <- c()
    arrow_y <- c()
    for (i in 1:nrow(loadings)) {
      arrow_x <- c(arrow_x, 0, loadings$PC1_scaled[i], NA)
      arrow_y <- c(arrow_y, 0, loadings$PC2_scaled[i], NA)
    }
    
    # --- Plotting ---
    p <- plot_ly() %>%
      # 1. Chemicals (Points)
      add_trace(
        data = df_pca, 
        x = ~PC1, y = ~PC2, 
        type = 'scatter', mode = 'markers',
        name = 'Chemicals',
        text = ~paste("ID:", molecule_id),
        customdata = ~standard_inchi_key,
        marker = list(opacity = 0.6)
      ) %>%
      # 2. Variable Arrows (Lines)
      add_trace(
        x = arrow_x, y = arrow_y,
        type = 'scatter', mode = 'lines',
        name = 'Variables',
        legendgroup = 'vars',
        line = list(color = 'red', width = 1)
      ) %>%
      # 3. Variable Labels (Text)
      add_trace(
        data = loadings,
        x = ~PC1_scaled, y = ~PC2_scaled,
        type = 'scatter', mode = 'text',
        name = 'Variables',
        legendgroup = 'vars',
        showlegend = FALSE,
        text = ~VarName,
        textposition = 'top center',
        textfont = list(color = 'red')
      ) %>%
      layout(
        title = paste0("PCA Biplot: ", input$pca_group),
        xaxis = list(title = paste0("PC1 (", scales::percent(var_exp[1]), ")")),
        yaxis = list(title = paste0("PC2 (", scales::percent(var_exp[2]), ")")),
        hovermode = "closest",
        legend = list(title = list(text = "Click to Hide:"))
      )
      
    p
  })
  
  # Link click for PCA plot as well
  observeEvent(event_data("plotly_click", source = "pca_plot"), {
    click_data <- event_data("plotly_click", source = "pca_plot")
    if (is.null(click_data)) return()
    key <- click_data$customdata
    if (!is.null(key) && !is.na(key) && key != "") {
      url <- paste0("https://coconut.naturalproducts.net/search?q=", key)
      runjs(sprintf("window.open('%s', '_blank');", url))
    }
  })

  # Reset filters
  observeEvent(input$reset_filters, {
    updateSelectizeInput(session, "filter_vars", selected = character(0))
  })
}

# Run the app
shinyApp(ui = ui, server = server)
