#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

raw_dir <- file.path("data", "raw")
zenodo_zip_url <- "https://zenodo.org/records/15641716/files/marinenp_v1_r20250608_csv.zip?download=1"
zenodo_zip_path <- file.path(raw_dir, "marinenp_v1_r20250608_csv.zip")
molecules_path <- file.path(raw_dir, "molecules.csv")
properties_path <- file.path(raw_dir, "properties.csv")

ensure_directories <- function() {
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
}

read_current_raw <- function() {
  if (file.exists(molecules_path) && file.exists(properties_path)) {
    list(
      molecules = readr::read_csv(molecules_path, show_col_types = FALSE),
      properties = readr::read_csv(properties_path, show_col_types = FALSE),
      found = TRUE
    )
  } else {
    list(molecules = tibble(), properties = tibble(), found = FALSE)
  }
}

download_archive <- function(url, dest) {
  message("[download] Fetching dataset from ", url)
  tryCatch(
    download.file(url, destfile = dest, mode = "wb", quiet = TRUE),
    error = function(e) {
      stop("Download failed: ", e$message, "\nIf downloads are blocked, manually place molecules.csv and properties.csv in ", raw_dir)
    }
  )
  message("[download] Saved archive to ", dest)
}

extract_required_files <- function(zip_path) {
  contents <- tryCatch(unzip(zip_path, list = TRUE), error = function(e) {
    stop("Unable to inspect downloaded archive: ", e$message)
  })

  required <- c("molecules.csv", "properties.csv")
  missing <- setdiff(required, contents$Name)
  if (length(missing) > 0) {
    stop("Archive is missing required files: ", paste(missing, collapse = ", "))
  }

  unzip(zip_path, files = required, exdir = raw_dir, overwrite = TRUE)
  invisible(required)
}

main <- function() {
  ensure_directories()

  existing <- read_current_raw()

  if (!existing$found) {
    download_archive(zenodo_zip_url, zenodo_zip_path)
    extract_required_files(zenodo_zip_path)
    file.remove(zenodo_zip_path)
    existing <- read_current_raw()
  } else {
    message("[download] Using already downloaded raw files in ", raw_dir)
  }

  readr::write_csv(existing$molecules, molecules_path)
  readr::write_csv(existing$properties, properties_path)
  message("[download] Molecules and properties saved to ", raw_dir)
}

if (identical(environment(), globalenv()) && !length(sys.frames())) {
  main()
}
